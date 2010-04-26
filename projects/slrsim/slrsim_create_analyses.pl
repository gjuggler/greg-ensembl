#!/usr/bin/env perl

use warnings;
use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::Greg::EslrUtils;
use File::Path;
use File::Basename;

my ( $clean, $url ) = undef;
GetOptions(
  'clean' => \$clean,
  'url=s' => \$url
);
Bio::EnsEMBL::Registry->no_version_check(1);

$url = 'mysql://slrsim:slrsim@mysql-greg.ebi.ac.uk:4134/slrsim_anisimova' if ( !$url );
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new( -url => $url );
my $hive_dba = Bio::EnsEMBL::Hive::DBSQL::DBAdaptor->new( -url => $url );
my $dbc = $dba->dbc;

clean_tables();

load();
simulate();
align();
scores();
omegas();
collect_stats();
plots();

connect_analysis( "LoadTrees",   "PhyloSim" );
connect_analysis( "PhyloSim",    "Align" );
connect_analysis( "Align",       "AlignScores" );
connect_analysis( "AlignScores", "Omegas" );
connect_analysis( "Omegas",      "CollectStats" );

wait_for("Plots",["Omegas","CollectStats"]);

sub clean_tables {
  if ($clean) {

    my @truncate_tables = qw^
      protein_tree_member protein_tree_node protein_tree_tag member sequence
      sequence sequence_cds
      analysis analysis_job dataflow_rule hive
      sitewise_omega
      stats_sites stats_genes
      parameter_set
      node_set_member node_set
      ^;
    map {
      print "$_\n";
      eval { $dba->dbc->do("truncate table $_"); }
    } @truncate_tables;
    eval { $dba->dbc->do("drop table stats_sites"); };
    eval { $dba->dbc->do("drop table stats_genes"); };
  }
}

sub load {
  my $logic_name  = "LoadTrees";
  my $module      = "Bio::Greg::Slrsim::LoadTrees";
  my $params      = { simulation_set => "filter_sweeps" };
  my $analysis_id = _create_analysis( $logic_name, $module, $params );

  _add_job_to_analysis( "LoadTrees", {} );    # Add a dummy job to run and load the trees.
}

sub simulate {
  my $logic_name = "PhyloSim";
  my $module     = "Bio::Greg::Hive::PhyloSim";
  my $params     = {};
  _create_analysis( $logic_name, $module, $params, 30, 1 );
}

sub align {
  my $logic_name = "Align";
  my $module     = "Bio::Greg::Hive::Align";
  my $params     = {
    # These params will be filled in by the LoadTree simulation definitions.
    t_coffee_executable => '/homes/greg/src/T-COFFEE_distribution_Version_8.06/bin/binaries/linux/t_coffee'
  };
  _create_analysis( $logic_name, $module, $params, 100, 1 );
}

sub scores {
  my $logic_name = "AlignScores";
  my $module     = "Bio::Greg::Hive::AlignmentScores";
  my $params     = {

    # These params will be filled in by the LoadTree simulation definitions.
  };
  _create_analysis( $logic_name, $module, $params, 100, 1 );
}

sub omegas {
  my $logic_name = "Omegas";
  my $module     = "Bio::Greg::Hive::PhyloAnalysis";
  my $params     = {

    # These params will be filled in by the LoadTree simulation definitions.
  };
  _create_analysis( $logic_name, $module, $params, 500, 1 );
}

sub collect_stats {
  my $logic_name = "CollectStats";
  my $module     = "Bio::Greg::Slrsim::CollectSlrsimStats";
  my $params     = {};
  _create_analysis( $logic_name, $module, $params, 50, 1 );
}

sub plots {
  my $logic_name = "Plots";
  my $module     = "Bio::Greg::Slrsim::Plots";
  my $params     = {};
  _create_analysis( $logic_name, $module, $params, 50, 1 );

  _add_job_to_analysis( "Plots", {} );    # Add a dummy job to plot at the end.
}

########*********########
#-------~~~~~~~~~-------#
########*********########

sub _combine_hashes {
  my @hashes = @_;

  my $new_hash = {};
  foreach my $hash (@hashes) {
    foreach my $key ( keys %{$hash} ) {
      $new_hash->{$key} = $hash->{$key};
    }
  }
  return $new_hash;
}

our $param_set_counter;

sub _add_parameter_set {
  my $params = shift;

  $param_set_counter = 1 if ( !$param_set_counter );
  my $parameter_set_id = $params->{'parameter_set_id'} || $param_set_counter++;
  $params->{'parameter_set_id'} = $parameter_set_id;

  if ( exists $params->{'parameter_set_name'} ) {
    my $parameter_set_name = $params->{'parameter_set_name'};
    my $name_cmd =
      "REPLACE INTO parameter_set VALUES ('$parameter_set_id','name',\"$parameter_set_name\");";
    $dbc->do($name_cmd);
  }

  if ( exists $params->{'parameter_set_shortname'} ) {
    my $parameter_set_shortname = $params->{'parameter_set_shortname'} || '';
    my $shortname_cmd =
      "REPLACE INTO parameter_set VALUES ('$parameter_set_id','parameter_set_shortname',\"$parameter_set_shortname\");";
    $dbc->do($shortname_cmd);
  }

  my $param_string = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
  my $cmd = "REPLACE INTO parameter_set VALUES ('$parameter_set_id','params',\"$param_string\");";
  $dbc->do($cmd);
}

our $analysis_counter = 0;

sub _create_analysis {
  my $logic_name    = shift;
  my $module        = shift;
  my $params        = shift;
  my $hive_capacity = shift || 500;
  my $batch_size    = shift || 1;

  my $analysis_id = ++$analysis_counter;

  my $param_string = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
  my $cmd          = qq{REPLACE INTO analysis SET
		 created=now(),
		 analysis_id=$analysis_id,
		 logic_name="$logic_name",
		 module="$module",
		 parameters="$param_string"
		 ;};
  $dbc->do($cmd);
  $cmd = qq{REPLACE INTO analysis_stats SET
	      analysis_id=$analysis_id,
	      hive_capacity=$hive_capacity,
	      batch_size=$batch_size,
	      failed_job_tolerance=20000
	      ;};
  $dbc->do($cmd);
  return $analysis_id;
}

sub connect_analysis {
  my $from_name   = shift;
  my $to_name     = shift;
  my $branch_code = shift;
  $branch_code = 1 unless ( defined $branch_code );

  my $dataflow_rule_adaptor = $hive_dba->get_DataflowRuleAdaptor;
  my $analysis_adaptor      = $hive_dba->get_AnalysisAdaptor;

  my $from_analysis = $analysis_adaptor->fetch_by_logic_name($from_name);
  my $to_analysis   = $analysis_adaptor->fetch_by_logic_name($to_name);

  if ( $from_analysis and $to_analysis ) {
    $dataflow_rule_adaptor->create_rule( $from_analysis, $to_analysis, $branch_code );
    warn "Created DataFlow rule: [$branch_code] $from_name -> $to_name\n";
  } else {
    die "Could not fetch analyses $from_analysis -> $to_analysis to create a dataflow rule";
  }
}

sub wait_for {
  my $waiting_name  = shift;
  my $wait_for_list = shift;

  my $ctrl_rule_adaptor = $hive_dba->get_AnalysisCtrlRuleAdaptor;
  my $analysis_adaptor  = $hive_dba->get_AnalysisAdaptor;

  my $waiting_analysis = $analysis_adaptor->fetch_by_logic_name($waiting_name);

  foreach my $wait_for_name (@$wait_for_list) {
    my $wait_for_analysis = $analysis_adaptor->fetch_by_logic_name($wait_for_name);

    if ( $waiting_analysis and $wait_for_analysis ) {
      $ctrl_rule_adaptor->create_rule( $wait_for_analysis, $waiting_analysis );
      warn "Created Control rule: $waiting_name will wait for $wait_for_name\n";
    } else {
      die "Could not fetch $waiting_name -> $wait_for_name to create a control rule";
    }
  }
}

sub _add_job_to_analysis {
  my $analysis_name = shift;
  my $input_id_hash = shift;

  my $analysis_adaptor = $hive_dba->get_AnalysisAdaptor;
  my $analysis         = $analysis_adaptor->fetch_by_logic_name($analysis_name);

  my $job_id = Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor->CreateNewJob(
    -input_id     => $input_id_hash,
    -analysis     => $analysis,
    -input_job_id => 0,
  );
  return $job_id;
}

sub _add_nodes_to_analysis {
  my $analysis_id   = shift;
  my $params        = shift || {};
  my $node_arrayref = shift;

  my @node_ids = @{$node_arrayref};
  my $sth      = $dbc->prepare(
    qq{REPLACE INTO analysis_job SET
			       analysis_id=?,
			       input_id=?;}
  );
  foreach my $node_id (@node_ids) {
    $params->{'node_id'} = $node_id;
    my $input_id = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
    $sth->execute( $analysis_id, $input_id );
    my $analysis_job_id = $sth->{'mysql_insertid'};
    print "New AnalysisJob: $analysis_id  $analysis_job_id  $input_id\n";
  }
  $sth->finish;
}

sub _select_node_ids {
  my $cmd = shift;
  if ( !defined $cmd ) {
    $cmd = "SELECT node_id FROM protein_tree_node WHERE parent_id=1 AND root_id=1";
  }
  my $sth = $dbc->prepare($cmd);
  $sth->execute();

  my $array_ref = $sth->fetchall_arrayref( [0] );
  my @node_ids = @{$array_ref};
  @node_ids =
    map { @{$_}[0] } @node_ids;    # Some weird mappings to unpack the numbers from the arrayrefs.
  $sth->finish;
  return @node_ids;
}
