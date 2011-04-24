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
my $url = 'mysql://ensadmin:ensembl@ens-research/gj1_ils';
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new( -url => $url );

$url = 'mysql://ensadmin:ensembl@ens-research/' . "gj1_ils";
my $h = new Bio::Greg::Hive::ComparaHiveLoaderUtils;
$h->init($url);

$h->clean_hive_tables;
my @truncate_tables = qw^
      sites
  ^;
$h->truncate_tables(\@truncate_tables);

ils_counts();

add_mb_regions();

sub add_mb_regions {
  my $slice_a = Bio::EnsEMBL::Registry->get_adaptor( 'Human', 'core', 'slice' );

  my $win_width = 1_000_000;
  my @chrs = @{$slice_a->fetch_all('chromosome')};  
  foreach my $chr (@chrs) {
    my $win_lo = 0;
    while ($win_lo <= $chr->length) {
      my $win_hi = $win_lo + $win_width;
      $win_hi = $chr->length if ($win_hi > $chr->length);
      $h->add_job_to_analysis('ILSCounts', 
                              {
                                chr_name => $chr->seq_region_name,
                                chr_start => $win_lo,
                                chr_end => $win_hi
                              });
      $win_lo += $win_width;
    }
  }
}

sub ils_counts {
  my $logic_name = "ILSCounts";
  my $module     = "Bio::Greg::Gorilla::ILSCounts";
  my $params     = {
  };
  $h->create_analysis( $logic_name, $module, $params, 10, 1 );
}
