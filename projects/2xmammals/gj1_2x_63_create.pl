#!/usr/bin/env perl

use warnings;
use strict;
use Bio::Greg::Hive::ComparaHiveLoaderUtils;

my $url = 'mysql://ensadmin:ensembl@ens-research/gj1_2x_63';
my $h = new Bio::Greg::Hive::ComparaHiveLoaderUtils;
$h->init($url);

my $clean = 1;

# Clean up our mess.
if ($clean) {
  $h->clean_compara_analysis_tables;
  $h->clean_hive_tables;
#  $h->truncate_tables([qw(meta genes sites)]);
}

# Create analyses.
sitewise_mammals();

# Connect the dots.

# Add some trees.
add_nodes();

sub add_nodes {
  my $cmd = "SELECT tree_id as node_id FROM gj1_orthologs_63.trees where method_id='Eutheria';";
  my $dba = $h->dba;
  my $sth = $dba->dbc->prepare($cmd);
  $sth->execute();
  my $array_ref = $sth->fetchall_arrayref( [0] );
  my @node_ids = @{$array_ref};
  my @node_ids =
    map { @{$_}[0] } @node_ids;    # Some weird mappings to unpack the numbers from the arrayrefs.
  $sth->finish;

  foreach my $node_id (@node_ids) {
    $h->add_job_to_analysis('SitewiseMammals',{
      node_id => $node_id
    });
  }
}

### Process definitions.
###

sub sitewise_mammals {
  my $logic_name = "SitewiseMammals";
  my $module = "Bio::Greg::Mammals::SitewiseMammals";
  my $base_params = {
  };
  $h->create_analysis($logic_name,$module,$base_params,50,1);
}
