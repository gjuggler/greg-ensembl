#!/usr/bin/env perl

use warnings;
use strict;

use Bio::Greg::Hive::ComparaHiveLoaderUtils;

my $url = 'mysql://ensadmin:ensembl@ens-research:3306/gj1_hum_chimp';
my $clean = 1;

my $h = new Bio::Greg::Hive::ComparaHiveLoaderUtils;
$h->init($url);

# Clean up our mess.
if ($clean) {
  $h->clean_hive_tables;
}

# Create analyses.
split_by_genes();
collect_substitutions();

# Connect the dots.
$h->connect_analysis("SplitByGenes","CollectSubstitutions");

my $feeder_id =
  {
   genome_taxon_id => 9606,
   genome_window_width => 1000000
};
$h->add_job_to_analysis("SplitByGenes",$feeder_id);

sub split_by_genes {
  my $logic_name = "SplitByGenes";
  my $module = "Bio::Greg::Hive::SplitByGenes";
  my $params = {
  };
  $h->create_analysis($logic_name,$module,$params,1,1);
}

sub collect_substitutions {
  my $logic_name = "CollectSubstitutions";
  my $module = "Bio::Greg::Hive::CollectSubstitutions";
  my $params = {
                keep_species => '9606,9598'
  };
  $h->create_analysis($logic_name,$module,$params,30,1);
}
