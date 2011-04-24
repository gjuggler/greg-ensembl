#!/usr/bin/env perl

use warnings;
use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::ComparaUtils;
use File::Path;
use File::Basename;

use Bio::Greg::Hive::HiveLoaderUtils;
use Bio::Greg::Hive::ComparaHiveLoaderUtils;

# First create a database for the experiment.
my $url = 'mysql://ensadmin:ensembl@ens-research/gj1_gorilla';
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new( -url => $url );

$url = 'mysql://ensadmin:ensembl@ens-research/' . "gj1_gorilla";
my $h = new Bio::Greg::Hive::ComparaHiveLoaderUtils;
$h->init($url);

$h->clean_hive_tables;
my @truncate_tables = qw^
      sites genes
  ^;
$h->truncate_tables(\@truncate_tables);

my $gene_id_hash;

full_tests();
branch_tests();

add_top_genes();
add_all_genes();

sub add_top_genes {
  my @genes = ();

  open INPUT, "< gorilla_genes.txt";
  my @lines = <INPUT>;
  close INPUT;
  @genes =  @lines;

  foreach my $gene (@genes) {
    chomp $gene;
    next if ($gene eq '');
    $h->add_job_to_analysis( 'BranchSiteTests', {gene_id => $gene});
    $gene_id_hash->{$gene} = 1;
  }
}

sub add_all_genes {
  my @genes = ();

  open INPUT, "< gorilla_genes_all.txt";
  my @lines = <INPUT>;
  close INPUT;
  @genes =  @lines;

  foreach my $gene (@genes) {
    chomp $gene;
    next if ($gene eq '');
    next if (defined $gene_id_hash->{$gene});
    $h->add_job_to_analysis( 'BranchTests', {gene_id => $gene});
  }
}

sub full_tests {
  my $logic_name = "BranchSiteTests";
  my $module     = "Bio::Greg::Gorilla::BranchTests";
  my $params     = {
    run_branchsite => 1
  };
  $h->create_analysis( $logic_name, $module, $params, 700, 1 );
}

sub branch_tests {
  my $logic_name = "BranchTests";
  my $module     = "Bio::Greg::Gorilla::BranchTests";
  my $params     = {
    run_branchsite => 0
  };
  $h->create_analysis( $logic_name, $module, $params, 700, 1 );
}

sub plots {
  my $logic_name = "Plots";
  my $module     = "Bio::Greg::Gorilla::Plots";
  my $params     = {};
  $h->create_analysis( $logic_name, $module, $params, 1, 1 );

  $h->add_job_to_analysis( "Plots", $params );    # Add a dummy job to plot at the end.
}

########*********########
#-------~~~~~~~~~-------#
########*********########
