#!/usr/bin/env perl

use warnings;
use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
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
my $dbc = $dba->dbc;

my $LIMIT            = '';
my $SKIP_ADDING_JOBS = 0;

clean_tables();

simulate_alignments();
align_sequences();
alignment_scores();
calculate_omegas();
collect_stats();

connect_analysis( "PhyloSim",    "Align",        1 );
connect_analysis( "Align",       "AlignScores",  1 );
connect_analysis( "AlignScores", "Omegas",       1 );
connect_analysis( "Omegas",      "CollectStats", 1 );

sub clean_tables {
  if ($clean) {
    my @truncate_tables = qw^
      sequence
      analysis analysis_job dataflow_rule hive
      sitewise_omega
      stats_slrsim
      parameter_set
      node_set_member node_set
      ^;
    map {
      print "$_\n";
      eval { $dba->dbc->do("truncate table $_"); }
    } @truncate_tables;
  }
}

sub simulate_alignments {
  my $analysis_id = 1;

  my $logic_name = "PhyloSim";
  my $module     = "Bio::Greg::PhyloSim";
  my $params     = {

  };
  _create_analysis( $analysis_id, $logic_name, $module, $params, 30, 1 );

  my $cmd   = "SELECT node_id FROM protein_tree_node WHERE parent_id=0 AND root_id=0 $LIMIT;";
  my @nodes = _select_node_ids($cmd);
  _add_nodes_to_analysis( $analysis_id, $params, \@nodes );
}

sub align_sequences {
  my $analysis_id = 2;

  my $logic_name = "Align";
  my $module     = "Bio::EnsEMBL::Compara::RunnableDB::MCoffee";
  my $params     = {
    alignment_table => 'aln_mcoffee',
    executable => '/homes/greg/src/T-COFFEE_distribution_Version_8.06/bin/binaries/linux/t_coffee'
  };
  _create_analysis( $analysis_id, $logic_name, $module, $params, 100, 1 );
}

sub alignment_scores {
  my $analysis_id = 3;

  my $logic_name = "AlignScores";
  my $module     = "Bio::Greg::AlignmentScores";
  my $params     = {
    alignment_table         => 'aln_mcoffee',
    alignment_scores_action => 'gblocks prank trimal'
  };
  _create_analysis( $analysis_id, $logic_name, $module, $params, 100, 1 );
}

sub calculate_omegas {
  my $analysis_id = 4;
  my $logic_name  = "Omegas";
  my $module      = "Bio::EnsEMBL::Compara::RunnableDB::Sitewise_dNdS";
  my $params      = { sitewise_parameter_sets => 'all' };
  _create_analysis( $analysis_id, $logic_name, $module, $params, 500, 1 );

  $params = {
    parameter_set_id => 1,
    name             => "Everything",
  };
  _add_parameter_set($params);
}

sub collect_stats {
  my $analysis_id = 5;

  my $logic_name = "CollectStats";
  my $module     = "Bio::Greg::Slrsim::CollectStats";
  my $params     = {
    alignment_table              => 'aln_mcoffee',
    collect_stats_parameter_sets => 'all'
  };
  _create_analysis( $analysis_id, $logic_name, $module, $params, 50, 1 );

#my @nodes = _select_node_ids("SELECT distinct(node_id) FROM protein_tree_tag where tag='sim_name'");
#_add_nodes_to_analysis($analysis_id,{},\@nodes);
}

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

sub _add_parameter_set {
  my $params           = shift;
  my $parameter_set_id = $params->{'parameter_set_id'};
  if ( exists $params->{'name'} ) {
    my $parameter_set_name = $params->{'name'};
    delete $params->{'name'};
    my $name_cmd =
      "REPLACE INTO parameter_set VALUES ('$parameter_set_id','name',\"$parameter_set_name\");";
    $dbc->do($name_cmd);
  }

  my $param_string = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
  my $cmd = "REPLACE INTO parameter_set VALUES ('$parameter_set_id','params',\"$param_string\");";
  $dbc->do($cmd);
}

my $analysis_hash;

sub _create_analysis {
  my $analysis_id   = shift;
  my $logic_name    = shift;
  my $module        = shift;
  my $params        = shift;
  my $hive_capacity = shift || 500;
  my $batch_size    = shift || 1;

  $analysis_hash->{$logic_name} = $analysis_id;

  my $param_string = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);

  my $cmd = qq{REPLACE INTO analysis SET
		 created=now(),
		 analysis_id=$analysis_id,
		 logic_name="$logic_name",
		 module="$module",
		 parameters="$param_string"
		 ;};
  $dbc->do($cmd);

  $cmd = qq{REPLACE INTO analysis_stats SET
	      hive_capacity=$hive_capacity,
	      batch_size=$batch_size,
	      analysis_id=$analysis_id
	      ;};
  $dbc->do($cmd);
}

sub connect_analysis {
  my $from_name   = shift;
  my $to_name     = shift;
  my $branch_code = shift;

  my $from_id = $analysis_hash->{$from_name};

  $branch_code = 1 unless ( defined $branch_code );

  my $cmd = qq{REPLACE INTO dataflow_rule SET
		 from_analysis_id=$from_id,
		 to_analysis_url="$to_name",
		 branch_code=$branch_code
		 ;};
  $dbc->do($cmd);
}

sub _add_nodes_to_analysis {
  return if ($SKIP_ADDING_JOBS);

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
