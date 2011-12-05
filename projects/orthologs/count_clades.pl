#!/usr/bin/env perl

use warnings;
use strict;
use DBI;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::Greg::EslrUtils;

use Bio::Greg::Hive::HiveLoaderUtils;
use Bio::Greg::Hive::ComparaHiveLoaderUtils;

my $url = 'mysql://ensadmin:ensembl@ens-research/gj1_orthologs_63';
my $clean = 1;

my $h = new Bio::Greg::Hive::ComparaHiveLoaderUtils;
$h->init($url);

count_clades();

sub count_clades {
  my $compara_dba = $h->compara_dba;

  my @clades = ('Fungi/Metazoa group',
                'Chordata',
                'Vertebrata',
                'Amniota',
                'Mammalia',
                'Eutheria',
                'Primates',
                'Glires',
                'Laurasiatheria',
                'Clupeocephala',
                'Sauria'
    );

  foreach my $clade (@clades) {
    my $tree = Bio::EnsEMBL::Compara::ComparaUtils->get_genome_taxonomy_below_level($compara_dba,
                                                                                    $clade);
    my $species = scalar($tree->leaves);
    my @leaves = $tree->leaves;
    print "$clade  ".$species."\n";
    my @aliases = map {$_->ensembl_alias} @leaves;
    my @taxids = map {$_->taxon_id} @leaves;
    print join(', ', sort {$a <=> $b} @taxids)."\n";
    print join(', ', sort {$a cmp $b} @aliases)."\n";
    #print $tree->ascii;
  }
}
