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

my $url = 'mysql://ensadmin:ensembl@ens-research/gj1_hiv_62';
my $clean = 1;

my $h = new Bio::Greg::Hive::ComparaHiveLoaderUtils;
$h->init($url);

# Clean up our mess.
if ($clean) {
  $h->clean_hive_tables;
}

# Create analyses.
window_analysis();
output_data();

add_all_genes();

sub add_all_genes {
  my $cmd = "SELECT distinct(mg.stable_id) FROM member mg, member m WHERE mg.member_id=m.gene_member_id AND m.taxon_id=9606 and m.source_name='ENSEMBLPEP'";

  my $compara_dba = $h->compara_dba;
  my $sth = $compara_dba->dbc->prepare($cmd);
  $sth->execute();
  my $array_ref = $sth->fetchall_arrayref( [0] );
  my @protein_ids = @{$array_ref};
  @protein_ids =
    map { @{$_}[0] } @protein_ids;    # Some weird mappings to unpack the numbers from the arrayrefs.
  $sth->finish;

  foreach my $id (@protein_ids) {
    $h->add_job_to_analysis("WindowAnalysis", {
      gene_id => $id
                            });
  }
}

sub window_analysis {
 my $logic_name = "WindowAnalysis";
 my $module = "Bio::Greg::PrimateHIV::WindowAnalysis";
 my $params = {
 };
 $h->create_analysis($logic_name,$module,$params,600,1);
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

  my $dba = $h->compara_dba;
  my @genomes = Bio::EnsEMBL::Compara::ComparaUtils->get_genomes_within_clade($dba,$clade);
  my @taxon_ids = map {$_->taxon_id} @genomes;
  return @taxon_ids;
}
sub coverage_taxon_ids {
  my $coverage = shift;
  
  my $dba = $h->compara_dba;
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
  my $dbc = $h->compara_dba;
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
