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

my ($url) = 'mysql://ensadmin:ensembl@ens-research/gj1_eslr_nodesets';
#GetOptions('url=s' => \$url);
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
#omegas();
mapping();
collect_stats();

# Connect the dots.
connect_analysis("NodeSets","Mapping",1);
#connect_analysis("Omegas","Mapping",1);
connect_analysis("Mapping","CollectStats",1);

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

  my @all = clade_taxon_ids();
  my @mamms = clade_taxon_ids("Eutheria");
  my $not_mammals = join(",",subtract(\@all,\@mamms));

  my $everything = join(",",clade_taxon_ids());
  my $mammals = join(",",clade_taxon_ids("Eutheria"));
  my $primates = join(",",clade_taxon_ids("Primates"));
  my $glires = join(",",clade_taxon_ids("Glires"));
  my $laurasiatheria = join(",",clade_taxon_ids("Laurasiatheria"));
  my $afrotheria = join(",",clade_taxon_ids("Afrotheria"));

  my $sauria = join(",",clade_taxon_ids("Sauria"));
  my $fish = join(",",clade_taxon_ids("Clupeocephala"));
  
  # Get only hi-coverage genomes.
  my $hi_coverage = join(",",coverage_taxon_ids("high"));
  my $lo_coverage = join(",",coverage_taxon_ids("low"));

  $params = {
    parameter_set_name => "Mammals",
    shortname => 'm',
    keep_species => $mammals,
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "Primates",
    shortname => 'p',
    keep_species => $primates
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "Glires",
    shortname => 'g',
    keep_species => $glires
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "Laurasiatheria",
    shortname => 'l',
    keep_species => $laurasiatheria
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "No 2x",
    shortname => 'n2',
    keep_species => $hi_coverage,
    remove_species => $not_mammals
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "2x Only",
    shortname => '2o',
    keep_species => $lo_coverage,
    remove_species => $not_mammals
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "No Primates",
    shortname => 'np',
    keep_species => $mammals,
    remove_species => $primates
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "No Glires",
    shortname => 'ng',
    keep_species => $mammals,
    remove_species => $glires
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

}

sub node_sets {
  my $analysis_id=101;
  my $logic_name = "NodeSets";
  my $module = "Bio::Greg::NodeSets";
  my $params = {
    flow_node_set => 10
  };
  _create_analysis($analysis_id,$logic_name,$module,$params,30,1);

  # Add all root nodes to this analysis.
  $params = {};
  my $cmd = "SELECT node_id FROM protein_tree_node WHERE parent_id=1;";
  my @nodes = _select_node_ids($cmd);
  _add_nodes_to_analysis($analysis_id,$params,\@nodes);  
}

sub align {
  my $analysis_id = 102;
  my $logic_name = "Align";
  my $module = "Bio::EnsEMBL::Compara::RunnableDB::MCoffee";
  my $params = {
    method => 'cmcoffee'
  };

  _create_analysis($analysis_id,$logic_name,$module,$params,400,1);
}

sub sequence_quality {
  my $analysis_id=103;
  my $logic_name = "SequenceQuality";
  my $module = "Bio::Greg::SequenceQualityLoader";
  my $params = {};
  
  _create_analysis($analysis_id,$logic_name,$module,$params,100,1);
}

sub omegas {
  my $analysis_id = 105;
  my $logic_name = "Omegas";
  my $module = "Bio::EnsEMBL::Compara::RunnableDB::Sitewise_dNdS";
  my $base_params = {
    parameter_sets => "all",
    sequence_quality_filtering => 0,
    alignment_quality_filtering => 0,
    };
  _create_analysis($analysis_id,$logic_name,$module,$base_params,500,1);

}

sub mapping {
  my $analysis_id=106;
  my $logic_name = "Mapping";
  my $module = "Bio::Greg::Eslr::SitewiseMapper";
  my $params = {
  };
  _create_analysis($analysis_id,$logic_name,$module,$params,200,1);
}


sub collect_stats {
  my $analysis_id=107;
  my $logic_name = "CollectStats";
  my $module = "Bio::Greg::Eslr::CollectEslrStats";
  my $params = {
  };
  _create_analysis($analysis_id,$logic_name,$module,$params,50,1);
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

our $param_set_counter;
sub _add_parameter_set {
  my $params = shift;

  $param_set_counter = 1 if (!$param_set_counter);
  my $parameter_set_id = $params->{'parameter_set_id'} || $param_set_counter++;
  $params->{'parameter_set_id'} = $parameter_set_id;

  if (exists $params->{'parameter_set_name'} ) {
    my $parameter_set_name = $params->{'parameter_set_name'};
    my $name_cmd = "REPLACE INTO parameter_set VALUES (\"$parameter_set_id\",'name',\"$parameter_set_name\");";
    $dbc->do($name_cmd);
  }

  if (exists $params->{'shortname'} ) {
    my $shortname = $params->{'shortname'} || '';
    my $shortname_cmd = "REPLACE INTO parameter_set VALUES (\"$parameter_set_id\",'shortname',\"$shortname\");";
    $dbc->do($shortname_cmd);
  }

  
  my $param_string = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
  my $cmd = "REPLACE INTO parameter_set VALUES ('$parameter_set_id','params',\"$param_string\");";
  $dbc->do($cmd);
}

sub clean_tables {
  if ($clean) {    
    my @truncate_tables = qw^
      analysis analysis_job analysis_stats dataflow_rule hive
      parameter_set
      node_set_member node_set
      sitewise_omega sitewise_tag sitewise_genome sitewise_pfam
      go_terms      
      ^;
    map {
      #print "$_\n";
      eval {$dba->dbc->do("truncate table $_");}} @truncate_tables;
  }
}



########*********########
#-------~~~~~~~~~-------#
########*********########

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
	      analysis_id=$analysis_id,
	      hive_capacity=$hive_capacity,
	      batch_size=$batch_size
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
