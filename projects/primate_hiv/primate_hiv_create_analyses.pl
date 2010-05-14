#!/usr/bin/env perl

use warnings;
use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::Greg::EslrUtils;
use File::Path;
use File::Basename;

use Bio::Greg::Hive::HiveLoaderUtils;
use Bio::Greg::Hive::ComparaHiveLoaderUtils;

my ($url) = undef;
GetOptions('url=s' => \$url);
my $clean = 1;

my $h = new Bio::Greg::Hive::ComparaHiveLoaderUtils;
$h->init($url);

# Clean up our mess.
if ($clean) {
  $h->clean_compara_analysis_tables;
  $h->clean_hive_tables;
}

# Define parameters (species sets, filtering options, etc).
parameter_sets();

# Create analyses.
node_sets();
align();
sequence_quality();
mapping();
split_by_windows();
split_by_parameter_set();
gene_omegas();
sitewise_omegas();
collect_stats();
output_data();

# Connect the dots.
$h->connect_analysis("NodeSets","Align");
$h->connect_analysis("Align","SequenceQuality");
$h->connect_analysis("SequenceQuality","SplitAlignments");
$h->connect_analysis("SplitAlignments","SplitParams");
$h->connect_analysis("SplitParams","GeneOmegas");
$h->connect_analysis("GeneOmegas","SitewiseOmegas");
$h->connect_analysis("SitewiseOmegas","CollectStats");

$h->connect_analysis("NodeSets","Mapping");
$h->wait_for("CollectStats",["Mapping"]);
$h->wait_for("OutputTabularData",["CollectStats"]);

my @genes = gene_list();
$h->add_genes_to_analysis("NodeSets",\@genes);

sub gene_list {
  my @list =  qw(CYBB);
  return @list;
}

sub parameter_sets {
  my $params;
  my $base_params={};

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
 
#  $params = {
#    parameter_set_name => "Homininae",
#    parameter_set_shortname => 'hmn',
#    keep_species => $homininae,
#    gorilla_count_species => '9593,9598,9606'
#    };
#  $h->add_parameter_set($params);
  
#  $params = {
#    parameter_set_name => "Hominidae",
#    parameter_set_shortname => 'hmd',
#    keep_species => $hominidae
#    };
#  $h->add_parameter_set($params);
  
#  $params = {
#    parameter_set_name => "NonHominidPrimates",
#    parameter_set_shortname => 'nonhom',
#    keep_species => $non_hominidae
#    };
#  _add_parameter_set($params);
  
  $params = {
    parameter_set_name => "Primates",
    parameter_set_shortname => 'p',
    keep_species => $primates
    };
  $h->add_parameter_set($params);

#  $params = {
#    parameter_set_name => "Laurasiatheria",
#    parameter_set_shortname => 'l',
#    keep_species => $laurasiatheria
#    };
#  _add_parameter_set($params);
  
#  $params = {
#    parameter_set_name => "Mammals",
#    parameter_set_shortname => 'm',
#    keep_species => $mammals
#    };
#  _add_parameter_set($params);
  
#  $params = {
#    parameter_set_name => "NonPrimateMammals",
#    parameter_set_shortname => 'nonprm',
#    keep_species => $non_primates
#    };
#  _add_parameter_set($params);
}

sub node_sets {
  my $logic_name = "NodeSets";
  my $module = "Bio::Greg::Hive::NodeSets";
  my $params = {
    flow_node_set => 'Primates'
  };
  $h->create_analysis($logic_name,$module,$params,50,1);
}

sub align {
  my $logic_name = "Align";
  my $module = "Bio::Greg::Hive::Align";
  my $params = {
    alignment_method => 'none',
  };

  $h->create_analysis($logic_name,$module,$params,400,1);
}

sub sequence_quality {
  my $logic_name = "SequenceQuality";
  my $module = "Bio::Greg::Hive::SequenceQualityLoader";
  my $params = {
    sequence_quality_filtering => 1
  };
  
  $h->create_analysis($logic_name,$module,$params,50,1);
}

sub split_by_domains {
  my $logic_name = "SplitAlignments";
  my $module = "Bio::Greg::Hive::SplitByProteinDomain";
  my $params = {
  };
  $h->create_analysis($logic_name,$module,$params,50,1);
}

sub split_by_windows {
  my $logic_name = "SplitAlignments";
  my $module = "Bio::Greg::Hive::SplitBySlidingWindow";
  my $params = {
    window_width => 200,
    window_step => 50
  };
  $h->create_analysis($logic_name,$module,$params,50,1);
}

sub split_by_parameter_set {
  my $logic_name = "SplitParams";
  my $module = "Bio::Greg::Hive::SplitByParameterSet";
  my $params = {
    flow_parameter_sets => 'all'
  };
  $h->create_analysis($logic_name,$module,$params,50,1);
}

sub gene_omegas {
  my $logic_name = "GeneOmegas";
  my $module = "Bio::Greg::Hive::PhyloAnalysis";
  my $base_params = {
    analysis_action => 'hyphy_dnds',
  };
  $h->create_analysis($logic_name,$module,$base_params,500,1);
}

sub sitewise_omegas {
  my $logic_name = "SitewiseOmegas";
  my $module = "Bio::Greg::Hive::PhyloAnalysis";
  my $base_params = {
    analysis_action => 'slr',
  };
  $h->create_analysis($logic_name,$module,$base_params,500,1);
}

sub mapping {
  my $logic_name = "Mapping";
  my $module = "Bio::Greg::Hive::SitewiseMapper";
  my $params = {
  };
  $h->create_analysis($logic_name,$module,$params,200,1);
}

sub collect_stats {
  my $logic_name = "CollectStats";
  my $module = "Bio::Greg::PrimateHIV::CollectPrimateHIVStats";
  my $params = {
  };
  $h->create_analysis($logic_name,$module,$params,80,1);
}

sub output_data {
  my $logic_name = "OutputTabularData";
  my $module = "Bio::Greg::PrimateHIV::OutputPrimateHIVData";
  my $params = {
  };
  $h->create_analysis($logic_name,$module,$params,50,1);
}

########*********########
#-------~~~~~~~~~-------#
########*********########

# Subroutines to return a list of taxon IDs with specific features.
sub clade_taxon_ids {
  my $clade = shift || 1;

  my $dba = $h->dba;
  my @genomes = Bio::EnsEMBL::Compara::ComparaUtils->get_genomes_within_clade($dba,$clade);
  my @taxon_ids = map {$_->taxon_id} @genomes;
  return @taxon_ids;
}
sub coverage_taxon_ids {
  my $coverage = shift;
  
  my $dba = $h->dba;
  my @output;
  my @all_gdb = Bio::EnsEMBL::Compara::ComparaUtils->get_genomes_within_clade($dba,1);
  foreach my $gdb (@all_gdb) {
    # This is finicky: we need to call the "db_adaptor" method to get the Bio::EnsEMBL::DBSQL::DBAdaptor object, and then the meta container.
    my $meta = $gdb->db_adaptor->get_MetaContainer;
    my $str = @{$meta->list_value_by_key('assembly.coverage_depth')}[0];
    push @output, $gdb->taxon_id if ($str eq $coverage);
  }  return @output;
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

sub _select_node_ids {
  my $cmd = shift;
  if ( !defined $cmd ) {
    $cmd = "SELECT node_id FROM protein_tree_node WHERE parent_id=1 AND root_id=1";
  }
  my $dbc = $h->dbc;
  my $sth = $dbc->prepare($cmd);
  $sth->execute();

  my $array_ref = $sth->fetchall_arrayref( [0] );
  my @node_ids = @{$array_ref};
  @node_ids =
    map { @{$_}[0] } @node_ids;    # Some weird mappings to unpack the numbers from the arrayrefs.
  $sth->finish;
  return @node_ids;
}

sub _add_nodes_to_analysis {
  my $logic_name = shift;
  my $node_id_arrayref = shift;

  my @node_ids = @{$node_id_arrayref};
  foreach my $node_id (@node_ids) {
    my $input_id = {
      node_id => $node_id
    };
    $h->add_job_to_analysis($logic_name,$input_id);
  }
}
