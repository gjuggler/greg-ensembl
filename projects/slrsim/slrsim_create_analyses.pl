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

my ($clean,$url) = undef;
GetOptions('clean' => \$clean,
           'url=s' => \$url
	   );
Bio::EnsEMBL::Registry->no_version_check(1);

$url = 'mysql://greg:TMOqp3now@mysql-greg.ebi.ac.uk:4134/gj1_slrsim_test' if (!$url);
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-url => $url);
my $dbc = $dba->dbc;

my $LIMIT='';
my $SKIP_ADDING_JOBS=0;

clean_tables();

simulate_alignments();
align_sequences();
#alignment_scores();
calculate_omegas();
collect_stats();


connect_analysis("PhyloSim","Align",1);
connect_analysis("Align","Omegas",1);
#connect_analysis("AlignScores","Omegas",1);
connect_analysis("Omegas","CollectStats",1);

sub clean_tables {
  if ($clean) {

    $dbc->do("create table if not exists omega_tr LIKE sitewise_omega");
    $dbc->do("create table if not exists omega_mc LIKE sitewise_omega");
    $dbc->do("create table if not exists aln_mcoffee LIKE protein_tree_member;");
    $dbc->do("create table if not exists aln_mcoffee_score LIKE protein_tree_member_score;");
    $dbc->do("create table if not exists aln_mcoffee_prank LIKE protein_tree_member_score;");
    $dbc->do("create table if not exists aln_mcoffee_trimal LIKE protein_tree_member_score;");
    $dbc->do("create table if not exists aln_mcoffee_gblocks LIKE protein_tree_member_score;");

    my @truncate_tables = qw^
sequence
aln_mcoffee aln_mcoffee_score aln_mcoffee_prank aln_mcoffee_trimal aln_mcoffee_gblocks
analysis analysis_job dataflow_rule hive
sitewise_omega sitewise_stats
omega_mc omega_tr
parameter_set
node_set_member node_set
      ^;
    map {print "$_\n";eval {$dba->dbc->do("truncate table $_");}} @truncate_tables;
  }
}

sub simulate_alignments {
  my $analysis_id=1;

  my $logic_name = "PhyloSim";
  my $module = "Bio::Greg::PhyloSim";
  my $params = {
    
  };
  _create_analysis($analysis_id,$logic_name,$module,$params,20,1);

  my $cmd = "SELECT node_id FROM protein_tree_node WHERE parent_id=0 AND root_id=0 $LIMIT;";
  my @nodes = _select_node_ids($cmd);
  _add_nodes_to_analysis($analysis_id,$params,\@nodes);  
}

sub align_sequences {
  my $analysis_id=2;

  my $logic_name = "Align";
  my $module = "Bio::EnsEMBL::Compara::RunnableDB::MCoffee";
  my $params = {
    input_table_base => 'protein_tree',
    output_table => 'aln_mcoffee',
    executable => '/homes/greg/src/T-COFFEE_distribution_Version_8.06/bin/binaries/linux/t_coffee'
  };
  _create_analysis($analysis_id,$logic_name,$module,$params,50,1);
}

sub alignment_scores {
  my $analysis_id=3;
  my $scores_table = "aln_mcoffee";

  my $logic_name = "AlignScores";
  my $module = "Bio::EnsEMBL::Compara::RunnableDB::AlignmentScores";
  my $params = {
    input_table => 'aln_mcoffee',
    action => 'gblocks prank trimal'
  };
  _create_analysis($analysis_id,$logic_name,$module,$params,30,1);
}

sub calculate_omegas {
  my $analysis_id=4;
  my $logic_name = "Omegas";
  my $module = "Bio::EnsEMBL::Compara::RunnableDB::Sitewise_dNdS";
  my $params = {
    parameter_sets => "1,2,3"
  };
  _create_analysis($analysis_id,$logic_name,$module,$params,400,1);

  $params = {
    parameter_set_id => 1,
    name => "SLR (true alignment)",
    input_table => 'protein_tree_member',
    output_table => 'omega_tr'
  };
  _add_parameter_set($params);
  
  my $bp = {
    input_table => 'aln_mcoffee',
    output_table => 'omega_mc'
  };
  
  $params = {
    parameter_set_id => 2,
    name => "SLR"
  };
  _add_parameter_set(_combine_hashes($bp,$params));

  $params = {
    parameter_set_id => 3,
    name => "PAML M3",
    action => 'paml',
    model => 'M3'
  };
  _add_parameter_set(_combine_hashes($bp,$params));

#  $params = {
#    parameter_set_id => 4,
#    name => "PAML M2",
#    action => 'paml',
#    model => 'M2'
#  };
#  _add_parameter_set(_combine_hashes($bp,$params));

#  $params = {
#    parameter_set_id => 5,
#    name => "SLR+trimAl",
#    alignment_score_filtering => 1,
#    alignment_score_table => 'aln_mcoffee_trimal',
#    alignment_score_threshold => 5,
#  };
#  _add_parameter_set(_combine_hashes($bp,$params));
#
#  $params = {
#    parameter_set_id => 6,
#    name => "SLR+Gblocks",
#    alignment_score_filtering => 1,
#    alignment_score_table => 'aln_mcoffee_gblocks',
#    alignment_score_threshold => 5,
#  };
#  _add_parameter_set(_combine_hashes($bp,$params));
#
#  $params = {
#    parameter_set_id => 7,
#    name => "SLR+Prank",
#    alignment_score_filtering => 1,
#    alignment_score_table => 'aln_mcoffee_prank',
#    alignment_score_threshold => 5,
#  };
#  _add_parameter_set(_combine_hashes($bp,$params));
#
#  $params = {
#    parameter_set_id => 8,
#    name => "SLR+MCoffee",
#    alignment_score_filtering => 1,
#    alignment_score_table => 'aln_mcoffee_score',
#    alignment_score_threshold => 5,
#  };
#  _add_parameter_set(_combine_hashes($bp,$params));

#  $params = {
#    parameter_set_id => 9,
#    name => "PAML LRT",
#    action => 'paml_lrt',
#    model => 'M7',
#    model_b => 'M8'
#  };
#  _add_parameter_set(_combine_hashes($bp,$params));
  
}

sub collect_stats {
  my $analysis_id=5;

  my $logic_name = "CollectStats";
  my $module = "Bio::Greg::Slrsim::CollectStats";
  my $params = {
    alignment_table => 'aln_mcoffee',
    parameter_sets => '2,3'
  };
  _create_analysis($analysis_id,$logic_name,$module,$params,30,1);

  #my @nodes = _select_node_ids("SELECT distinct(node_id) FROM protein_tree_tag where tag='sim_name'");
  #_add_nodes_to_analysis($analysis_id,{},\@nodes);
}

sub _combine_hashes {
  my @hashes = @_;

  my $new_hash = {};
  foreach my $hash (@hashes) {
    foreach my $key (keys %{$hash}) {
      $new_hash->{$key} = $hash->{$key};
    }
  }
  return $new_hash;
}


sub _add_parameter_set {
  my $params = shift;
  my $parameter_set_id = $params->{'parameter_set_id'};
  if (exists $params->{'name'} ) {
    my $parameter_set_name = $params->{'name'};
    delete $params->{'name'};
    my $name_cmd = "REPLACE INTO parameter_set VALUES ('$parameter_set_id','name',\"$parameter_set_name\");";
    $dbc->do($name_cmd);
  }
  
  my $param_string = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
  my $cmd = "REPLACE INTO parameter_set VALUES ('$parameter_set_id','params',\"$param_string\");";
  $dbc->do($cmd);
}

my $analysis_hash;
sub _create_analysis {
  my $analysis_id = shift;
  my $logic_name = shift;
  my $module = shift;
  my $params = shift;
  my $hive_capacity = shift || 500;
  my $batch_size = shift || 1;
  
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
  my $from_name = shift;
  my $to_name = shift;
  my $branch_code = shift;

  my $from_id = $analysis_hash->{$from_name};

  $branch_code = 1 unless (defined $branch_code);

  my $cmd = qq{REPLACE INTO dataflow_rule SET
		 from_analysis_id=$from_id,
		 to_analysis_url="$to_name",
		 branch_code=$branch_code
		 ;};
  $dbc->do($cmd);
}

sub _add_nodes_to_analysis {
  return if ($SKIP_ADDING_JOBS);

  my $analysis_id = shift;
  my $params = shift || {};
  my $node_arrayref = shift;

  my @node_ids = @{$node_arrayref};  
  my $sth = $dbc->prepare( qq{REPLACE INTO analysis_job SET
			       analysis_id=?,
			       input_id=?;}
			   );
  foreach my $node_id (@node_ids) {
    $params->{'node_id'} = $node_id;
    my $input_id = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
    $sth->execute($analysis_id,$input_id);
    my $analysis_job_id = $sth->{'mysql_insertid'};
    print "New AnalysisJob: $analysis_id  $analysis_job_id  $input_id\n";
  }
  $sth->finish;
}

sub _select_node_ids {
  my $cmd = shift;

  if (!defined $cmd) {
    $cmd = "SELECT node_id FROM protein_tree_node WHERE parent_id=1 AND root_id=1";
  }

  my $sth = $dbc->prepare($cmd);
  $sth->execute();

  my $array_ref = $sth->fetchall_arrayref([0]);
  my @node_ids = @{$array_ref};
  @node_ids = map {@{$_}[0]} @node_ids;  # Some weird mappings to unpack the numbers from the arrayrefs.
  $sth->finish;
  return @node_ids;
}
