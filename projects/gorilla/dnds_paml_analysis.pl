#!/usr/bin/env perl
use warnings;
use strict;
use DBI;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::ComparaUtils;

use Bio::Greg::Hive::HiveLoaderUtils;
use Bio::Greg::Hive::ComparaHiveLoaderUtils;

# First create a database for the experiment.
my $url = 'mysql://ensadmin:ensembl@ens-research/gj1_dnds';
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new( -url => $url );

$url = 'mysql://ensadmin:ensembl@ens-research/' . "gj1_dnds";
my $h = new Bio::Greg::Hive::ComparaHiveLoaderUtils;
$h->init($url);

$h->clean_hive_tables;
my @truncate_tables = qw^
      
  ^;
$h->truncate_tables(\@truncate_tables);

paml_tests();
add_lrt_quantiles();

sub add_lrt_quantiles {
  my $slice_a = Bio::EnsEMBL::Registry->get_adaptor( 'Human', 'core', 'slice' );

  my $intervals = [
    [-209.746, -43.601],
    [-43.601, -23.754],
    [-23.754, -11.566],
    [-11.566, 0.8505],
    [0.8505, 118.981]
    ];
  foreach my $interval (@$intervals) {
    my ($lo, $hi) = @$interval;
    $h->add_job_to_analysis('PAMLTests', 
                            {
                              lrt_lo => $lo,
                              lrt_hi => $hi
                            });
  }
}

sub paml_tests {
  my $logic_name = "PAMLTests";
  my $module     = "Bio::Greg::Gorilla::PAMLTests";
  my $params     = {
  };
  $h->create_analysis( $logic_name, $module, $params, 50, 1 );
}
