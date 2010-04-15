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

my ($url) = undef;
GetOptions('url=s' => \$url);
my $clean = 1;
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-url => $url);
my $dbc = $dba->dbc;
my $analysis_hash; # Hash to store mapping between analysis names and ID numbers.

# Clean up our mess.
clean_tables();

# Define parameters (species sets, filtering options, etc).
parameter_sets();

# Create analyses.
node_sets();
align();
omegas();
mapping();
collect_stats();

# Connect the dots.
connect_analysis("NodeSets","Align");
connect_analysis("Align","SplitByParameterSet");
connect_analysis("SplitByParameterSet","SplitByProteinDomain");
connect_analysis("SplitByProteinDomain","GeneOmegas");
connect_analysis("GeneOmegas","SitewiseOmegas");
connect_analysis("SitewiseOmegas","CollectStats");

connect_analysis("Align","Mapping");
wait_for("CollectStats",["Mapping"]);
wait_for("OutputTabularData",["CollectStats"]);

sub parameter_sets {
  my $params;
  my $base_params={};

  # Subroutines to return a list of taxon IDs with specific features.
  sub clade_taxon_ids {
    my $clade = shift || 1;
    my @genomes = Bio::EnsEMBL::Compara::ComparaUtils->get_genomes_within_clade($dba,$clade);
    my @taxon_ids = map {$_->taxon_id} @genomes;
    return @taxon_ids;
  }
  sub coverage_taxon_ids {
    my $coverage = shift;

    my @output;
    my @all_gdb = Bio::EnsEMBL::Compara::ComparaUtils->get_genomes_within_clade($dba,1);
    foreach my $gdb (@all_gdb) {
      # This is finicky: we need to call the "db_adaptor" method to get the Bio::EnsEMBL::DBSQL::DBAdaptor object, and then the meta container.
      my $meta = $gdb->db_adaptor->get_MetaContainer;
      my $str = @{$meta->list_value_by_key('assembly.coverage_depth')}[0];
      push @output, $gdb->taxon_id if ($str eq $coverage);
    }
    return @output;
  }
  sub subtract {
    my $list_a = shift;
    my @remove_us = @_;
    my $hash;
    map {$hash->{$_}=1} @$list_a;
    foreach my $list_b (@remove_us) {
      map {delete $hash->{$_}} @$list_b;
    }
    return keys %$hash;
  }

  my @all_arr = clade_taxon_ids();
  my @mammals_arr = clade_taxon_ids("Eutheria");
  my @primates_arr = clade_taxon_ids("Primates");
  my @homininae_arr = clade_taxon_ids("Homininae");
  my @hominidae_arr = clade_taxon_ids("Hominidae");

  my $everything = join(",",clade_taxon_ids());

  my $primates = join(",",clade_taxon_ids("Primates"));
  my $homininae = join(",",clade_taxon_ids("Homininae"));
  my $hominidae = join(",",clade_taxon_ids("Hominidae"));
  my $non_homininae = join(",",subtract(\@primates_arr,\@homininae_arr));
  my $non_hominidae = join(",",subtract(\@primates_arr,\@hominidae_arr));
  my $non_gorilla = join(",",subtract(\@hominidae_arr,[9593]));
  my $simiiformes = join(",",clade_taxon_ids("Simiiformes"));
  my $haplorrhini = join(",",clade_taxon_ids("Haplorrhini"));

  my $glires = join(",",clade_taxon_ids("Glires"));
  my $laurasiatheria = join(",",clade_taxon_ids("Laurasiatheria"));
  my $afrotheria = join(",",clade_taxon_ids("Afrotheria"));

  my $mammals = join(",",clade_taxon_ids("Eutheria"));

  my $non_primates = join(",",subtract(\@mammals_arr,\@primates_arr));
 
  $params = {
    parameter_set_name => "Homininae",
    parameter_set_shortname => 'hmn',
    keep_species => $homininae,
    gorilla_count_species => '9593,9598,9606'
    };
  _add_parameter_set($params);
  
  $params = {
    parameter_set_name => "Hominidae",
    parameter_set_shortname => 'hmd',
    keep_species => $hominidae
    };
  _add_parameter_set($params);
  
  $params = {
    parameter_set_name => "NonHominidPrimates",
    parameter_set_shortname => 'nonhom',
    keep_species => $non_hominidae
    };
  _add_parameter_set($params);
  
  $params = {
    parameter_set_name => "Haplorrhini",
    parameter_set_shortname => 'hpl',
    keep_species => $haplorrhini
    };
  _add_parameter_set($params);
  
  $params = {
    parameter_set_name => "Primates",
    parameter_set_shortname => 'p',
    keep_species => $primates
    };
  _add_parameter_set($params);
  
  $params = {
    parameter_set_name => "Mammals",
    parameter_set_shortname => 'm',
    keep_species => $mammals
    };
  _add_parameter_set($params);
  
  $params = {
    parameter_set_name => "NonPrimateMammals",
    parameter_set_shortname => 'nonprm',
    keep_species => $non_primates
    };
  _add_parameter_set($params);
}

sub node_sets {
  my $logic_name = "NodeSets";
  my $module = "Bio::Greg::Hive::NodeSets";
  my $params = {
    flow_node_set => 'Primates'
  };
  my $analysis_id = _create_analysis($logic_name,$module,$params,50,1);

  # Add all root nodes to this analysis.
  $params = {};
  my $cmd = "SELECT node_id FROM protein_tree_node WHERE parent_id=1;";
  my @nodes = _select_node_ids($cmd);
  _add_nodes_to_analysis($analysis_id,$params,\@nodes);  
}

sub align {
  my $logic_name = "Align";
  my $module = "Bio::EnsEMBL::Compara::RunnableDB::MCoffee";
  my $params = {
    alignment_method => 'prank_f',
  };

  _create_analysis($logic_name,$module,$params,400,1);
}

sub sequence_quality {
  my $logic_name = "SequenceQuality";
  my $module = "Bio::Greg::SequenceQualityLoader";
  my $params = {};
  
  _create_analysis($logic_name,$module,$params,100,1);
}

sub omegas {
  my $logic_name = "Omegas";
  my $module = "Bio::EnsEMBL::Compara::RunnableDB::Sitewise_dNdS";
  my $base_params = {
    parameter_sets => "all",
    sequence_quality_filtering => 0,
    alignment_score_filtering => 1,
    };
  _create_analysis($logic_name,$module,$base_params,500,1);
}

sub mapping {
  my $logic_name = "Mapping";
  my $module = "Bio::Greg::SitewiseMapper";
  my $params = {
  };
  _create_analysis($logic_name,$module,$params,200,1);
}


sub collect_stats {
  my $logic_name = "CollectStats";
  my $module = "Bio::Greg::Eslr::CollectEslrStats";
  my $params = {
  };
  _create_analysis($logic_name,$module,$params,80,1);
}


sub clean_tables {
  if ($clean) {    
    my @truncate_tables = qw^
      analysis analysis_job analysis_stats dataflow_rule hive
      parameter_set
      node_set_member node_set
      sitewise_omega sitewise_tag sitewise_genome
      go_terms
      stats_sites stats_genes
      ^;
    map {
      print "$_\n";
      eval {$dba->dbc->do("truncate table $_");}} @truncate_tables;
  }
}


########*********########
#-------~~~~~~~~~-------#
########*********########

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

our $param_set_counter;
sub _add_parameter_set {
  my $params = shift;

  $param_set_counter = 1 if (!$param_set_counter);
  my $parameter_set_id = $params->{'parameter_set_id'} || $param_set_counter++;
  $params->{'parameter_set_id'} = $parameter_set_id;

  if (exists $params->{'parameter_set_name'} ) {
    my $parameter_set_name = $params->{'parameter_set_name'};
    my $name_cmd = "REPLACE INTO parameter_set VALUES ('$parameter_set_id','name',\"$parameter_set_name\");";
    $dbc->do($name_cmd);
  }

  if (exists $params->{'parameter_set_shortname'} ) {
    my $parameter_set_shortname = $params->{'parameter_set_shortname'} || '';
    my $shortname_cmd = "REPLACE INTO parameter_set VALUES ('$parameter_set_id','parameter_set_shortname',\"$parameter_set_shortname\");";
    $dbc->do($shortname_cmd);
  }
  
  my $param_string = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
  my $cmd = "REPLACE INTO parameter_set VALUES ('$parameter_set_id','params',\"$param_string\");";
  $dbc->do($cmd);
}

our $analysis_counter = 0;
sub _create_analysis {
  my $logic_name = shift;
  my $module = shift;
  my $params = shift;
  my $hive_capacity = shift || 500;
  my $batch_size = shift || 1;

  my $analysis_id = ++$analysis_counter;
  
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
	      analysis_id=$analysis_id,
	      hive_capacity=$hive_capacity,
	      batch_size=$batch_size,
	      failed_job_tolerance=20000
	      ;};
  $dbc->do($cmd);
  return $analysis_id;
}

sub connect_analysis {
  my $from_name = shift;
  my $to_name = shift;
  my $branch_code = shift;
  $branch_code = 1 unless (defined $branch_code);

  my $dataflow_rule_adaptor = $hive_dba->get_DataflowRuleAdaptor;
  my $analysis_adaptor = $hive_dba->get_AnalysisAdaptor;

  my $from_analysis = $analysis_adaptor->fetch_by_logic_name($from_name);
  my $to_analysis = $analysis_adaptor->fetch_by_logic_name($to_name);
  
  if($from_analysis and $to_analysis) {
    $dataflow_rule_adaptor->create_rule( $from_analysis, $to_analysis, $branch_code);
    warn "Created DataFlow rule: [$branch_code] $from_name -> $to_name\n";
  } else {
    die "Could not fetch analyses $from_analysis -> $to_analysis to create a dataflow rule";
  }
}

sub wait_for {
  my $waiting_name = shift;
  my $wait_for_list = shift;
  
  my $ctrl_rule_adaptor = $hive_dba->get_AnalysisCtrlRuleAdaptor;
  my $analysis_adaptor = $hive_dba->get_AnalysisAdaptor;

  my $waiting_analysis = $analysis_adaptor->fetch_by_logic_name($waiting_name);

  foreach my $wait_for_name (@$wait_for_list) {
    my $wait_for_analysis = $analysis_adaptor->fetch_by_logic_name($wait_for_name);

    if($waiting_analysis and $wait_for_analysis) {
      $ctrl_rule_adaptor->create_rule( $wait_for_analysis, $waiting_analysis);
      warn "Created Control rule: $waiting_name will wait for $wait_for_name\n";
    } else {
      die "Could not fetch $waiting_name -> $wait_for_name to create a control rule";
    }
  }
}

sub _add_nodes_to_analysis {
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
