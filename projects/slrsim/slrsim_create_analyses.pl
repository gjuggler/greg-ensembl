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

my ($clean) = undef;
GetOptions('clean' => \$clean
	   );

Bio::EnsEMBL::Registry->no_version_check(1);

my $url = 'mysql://greg:TMOqp3now@mysql-greg.ebi.ac.uk:4134/gj1_slrsim_1';
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-url => $url);
my $dbc = $dba->dbc;
my $dbh = $dbc->db_handle;
my $pta = $dba->get_ProteinTreeAdaptor();
my $mba = $dba->get_MemberAdaptor();

my $LIMIT='';
my $SKIP_ADDING_JOBS=0;

create_tables();
simulate_alignments();
align_sequences();
calculate_omegas();

connect_analysis("PhyloSim","Align",1);
connect_analysis("Align","Omegas",1);

sub create_tables {
  $dbc->do("CREATE TABLE IF NOT EXISTS aln_mcoffee LIKE protein_tree_member");
  $dbc->do("CREATE TABLE IF NOT EXISTS aln_mcoffee_score LIKE protein_tree_member_score");
  $dbc->do("CREATE TABLE IF NOT EXISTS omega_tr LIKE sitewise_aln");
  $dbc->do("CREATE TABLE IF NOT EXISTS omega_mc LIKE sitewise_aln");
  $dbc->do("TRUNCATE TABLE analysis_job") if ($clean);
}

sub simulate_alignments {
  my $analysis_id=1;
  my $logic_name = "PhyloSim";
  my $module = "Bio::Greg::PhyloSim";
  my $params = {
    
  };
  _create_analysis($analysis_id,$logic_name,$module,$params,30,1);

  my $cmd = "SELECT node_id FROM protein_tree_node WHERE parent_id=0 AND root_id=0 $LIMIT;";
  my @nodes = _select_node_ids($cmd);
  $dbc->do("DELETE from analysis_job WHERE analysis_id=$analysis_id;");
  _add_nodes_to_analysis($analysis_id,$params,\@nodes);  
}

sub align_sequences {
  my $tree_table = "aln_mcoffee";
  if ($clean) {
    # Delete all trees and reset the counter.
    $dba->dbc->do("truncate table  ${tree_table};");
    $dba->dbc->do("truncate table ${tree_table}_score;");
  }

  my $analysis_id=2;
  my $logic_name = "Align";
  my $module = "Bio::EnsEMBL::Compara::RunnableDB::MCoffee";
  my $params = {
    input_table_base => 'protein_tree',
    output_table => 'aln_mcoffee',
    executable => '/homes/greg/src/T-COFFEE_distribution_Version_8.06/bin/binaries/linux/t_coffee'
  };
  _create_analysis($analysis_id,$logic_name,$module,$params,30,1);
}

sub calculate_omegas {
  my $analysis_id=3;
  my $logic_name = "Omegas";
  my $module = "Bio::EnsEMBL::Compara::RunnableDB::Sitewise_dNdS";
  my $params = {
    parameter_sets => "1,2,3",
  };
  _create_analysis($analysis_id,$logic_name,$module,$params,30,1);

  $params = {
    parameter_set_id => 1,
    input_table => 'protein_tree_member',
    output_table => 'omega_tr'
  };
  _add_parameter_set($params);

   $params = {
    parameter_set_id => 2,
    input_table => 'aln_mcoffee',
    output_table => 'omega_mc'
  };
  _add_parameter_set($params);

  $params = {
    parameter_set_id => 3,
    trimal_filtering => 1,
    input_table => 'aln_mcoffee',
    output_table => 'omega_mc'
  };
  _add_parameter_set($params);
  
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
  
  if (exists $params->{'parameter_set_name'} ) {
    my $parameter_set_name = $params->{'parameter_set_name'};
    delete $params->{'parameter_set_name'};
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

  $cmd = qq{UPDATE analysis_stats SET
	      hive_capacity=$hive_capacity,
	      batch_size=$batch_size
	      WHERE analysis_id=$analysis_id
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
