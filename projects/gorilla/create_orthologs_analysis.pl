#!/usr/bin/env perl

use warnings;
use strict;
use DBI;
use Bio::EnsEMBL::Compara::ComparaUtils;

use Bio::Greg::Hive::HiveLoaderUtils;
use Bio::Greg::Hive::ComparaHiveLoaderUtils;

# First create a database for the experiment.
my $url = 'mysql://ensadmin:ensembl@ens-research/gj1_orthologs';
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new( -url => $url );

$url = 'mysql://ensadmin:ensembl@ens-research/' . "gj1_orthologs";
my $h = new Bio::Greg::Hive::ComparaHiveLoaderUtils;
$h->init($url);

$h->clean_hive_tables;
my @truncate_tables = qw^
      orthologs
  ^;
$h->truncate_tables(\@truncate_tables);

count_orthologs();
add_all_genes();

sub add_all_genes {
  my $cmd = "SELECT distinct(m2.stable_id) FROM member m, member m2 WHERE m.gene_member_id=m2.member_id and m.source_name='ENSEMBLPEP' and m2.source_name='ENSEMBLGENE' and m.taxon_id=9606 and m2.taxon_id=9606";

  my $dbc = $h->compara_dba;
  my $sth = $dbc->prepare($cmd);
  $sth->execute();
  my $array_ref = $sth->fetchall_arrayref( [0] );
  my @protein_ids = @{$array_ref};
  @protein_ids =
    map { @{$_}[0] } @protein_ids;    # Some weird mappings to unpack the numbers from the arrayrefs.
  $sth->finish;

  foreach my $id (@protein_ids) {
    $h->add_job_to_analysis( 'CountOrthologs', {gene_id => $id});
  }
}

sub count_orthologs {
  my $logic_name = "CountOrthologs";
  my $module     = "Bio::Greg::Hive::CountOrthologs";
  my $params     = {
  };
  $h->create_analysis( $logic_name, $module, $params, 300, 1 );
}

########*********########
#-------~~~~~~~~~-------#
########*********########
