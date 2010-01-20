#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;

my ($url) = undef;
GetOptions('url=s' => \$url);
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-url => $url);

  my $nhx = Bio::EnsEMBL::Compara::ComparaUtils->get_genome_tree_nhx($dba,{labels => 'mnemonics',
									   images => 1});
  print $nhx."\n";
