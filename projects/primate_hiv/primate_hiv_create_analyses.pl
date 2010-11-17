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

my $url = 'mysql://ensadmin:ensembl@ens-research:3306/gj1_hiv_58';
my $clean = 1;

my $h = new Bio::Greg::Hive::ComparaHiveLoaderUtils;
$h->init($url);

# Clean up our mess.
if ($clean) {
  $h->clean_compara_analysis_tables;
  $h->clean_hive_tables;
  $h->drop_tables(['stats_windows']);
}

# Define parameters (species sets, filtering options, etc).
parameter_sets();

# Create analyses.
node_sets();
split_params();
window_analysis();
output_data();

# Connect the dots.
$h->connect_analysis("NodeSets","SplitParams");
$h->connect_analysis("SplitParams","WindowAnalysis");

$h->wait_for("OutputData",["NodeSets","SplitParams","WindowAnalysis"]);

add_all_genes();

#my @genes = gene_list();
#@genes = @genes[1..30];
#@genes = grep {$_ =~ m/(ENSG00000197123|ENSG00000119977)/gi} @genes;
#$h->add_genes_to_analysis("NodeSets",\@genes);
#add_all_nodes();

sub add_all_genes {
  my $cmd = "SELECT stable_id FROM member m join protein_tree_member ptm USING(member_id) WHERE taxon_id=9606 and source_name='ENSEMBLPEP'";

  my $dbc = $h->dbc;
  my $sth = $dbc->prepare($cmd);
  $sth->execute();
  my $array_ref = $sth->fetchall_arrayref( [0] );
  my @protein_ids = @{$array_ref};
  @protein_ids =
    map { @{$_}[0] } @protein_ids;    # Some weird mappings to unpack the numbers from the arrayrefs.
  $sth->finish;

  $h->add_genes_to_analysis("NodeSets",\@protein_ids);
}

sub add_all_nodes {
  my $cmd = "SELECT node_id FROM protein_tree_node WHERE parent_id=1;";
  my @nodes = _select_node_ids($cmd);
  my $logic_name = "NodeSets";
  _add_nodes_to_analysis($logic_name,\@nodes);
}


sub gene_list {
  open INPUT, "<scale_up_genes.txt";
  my @lines = <INPUT>;
  close INPUT;
  return @lines;
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

  foreach my $aln_type ('genomic_mammals') {
    my $aln_short = 'gm';
    $aln_short = 'c' if ($aln_type eq 'compara');
    $aln_short = 'ga' if ($aln_type eq 'genomic_all');

      $params = {
        parameter_set_name => "Primates".' '.$aln_type,
        parameter_set_shortname => 'p_'.$aln_short,
        quality_threshold => 30,
        keep_species => $primates,
        aln_type => $aln_type
      };
      $h->add_parameter_set($params);

#      $params = {
#        parameter_set_name => "Mammals".' '.$aln_type,
#        parameter_set_shortname => 'm_'.$aln_short,
#        quality_threshold => 30,
#        keep_species => $mammals,
#        aln_type => $aln_type
#      };
#      $h->add_parameter_set($params);
      
  }
}

sub node_sets {
  my $logic_name = "NodeSets";
  my $module = "Bio::Greg::Hive::NodeSets";
  my $params = {
    flow_node_set => 'Primates'
  };
  $h->create_analysis($logic_name,$module,$params,20,1);
}

sub sequence_quality {
  my $logic_name = "SequenceQuality";
  my $module = "Bio::Greg::Hive::SequenceQualityLoader";
  my $params = {
    sequence_quality_filtering => 1
  };
  
  $h->create_analysis($logic_name,$module,$params,50,1);
}

sub split_params {
  my $logic_name = "SplitParams";
  my $module = "Bio::Greg::Hive::SplitByParameterSet";
  my $params = {
    flow_parameter_sets => 'all'
  };
  $h->create_analysis($logic_name,$module,$params,20,1);
}

sub window_analysis {
 my $logic_name = "WindowAnalysis";
 my $module = "Bio::Greg::PrimateHIV::WindowAnalysis";
 my $params = {
 };
 $h->create_analysis($logic_name,$module,$params,50,1);

}

sub output_data {
  my $logic_name = "OutputData";
  my $module = "Bio::Greg::PrimateHIV::OutputPrimateHIVData";
  my $params = {
  };
  $h->create_analysis($logic_name,$module,$params,50,1);
  $h->add_job_to_analysis($logic_name,{});
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
