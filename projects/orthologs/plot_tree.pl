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

use Bio::Greg::Hive::HiveLoaderUtils;
use Bio::Greg::Hive::ComparaHiveLoaderUtils;

Bio::EnsEMBL::Compara::ComparaUtils->load_registry;

my $compara_dba = Bio::EnsEMBL::Registry->get_DBAdaptor( 'multi', 'compara' );

my $tree = Bio::EnsEMBL::Compara::ComparaUtils->get_genome_tree_with_extras($compara_dba,
                                                                            {});
my $taxid_map;
map {$taxid_map->{$_->id} = $_->taxon_id if ($_->id ne '')} $tree->nodes;

my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
print $treeI->ascii(1,1,1);

my @internal_includes = (
  'Primates', 'Glires', 'Laurasiatheria',
  'Eutheria', 'Theria', 'Mammalia', 'Amniota', 'Tetrapoda', 'Euteleostomi', 'Vertebrata', 'Chordata', 
  'Clupeocephala', 'Sauria', 'Atlantogenata',
  'Eukaryota'
  );

foreach my $node ($treeI->nodes) {
  my $cvrg = $node->get_tag_value('low_coverage');
  $node->remove_all_tags;
  if (defined $cvrg) {
    $node->set_tag_value('low_coverage', $cvrg);
    #print $node->id."  ". $node->get_tag_value('low_coverage')."\n";
  }

  if (!$node->is_leaf) {
#    $node->id('') unless (grep {$node->id eq $_} @internal_includes);
  }
}

my $nhx_str = $treeI->root->as_text('nhx');
#my $nhx_str = $treeI->root->as_text('newick');
open(OUT, ">compara_63_tree.nhx");
print OUT $nhx_str."\n";
close(OUT);

# Output a taxon_id tree.
my $taxids = $treeI;
$taxids->translate_ids($taxid_map);
#print $taxids->ascii;
my $tax_str = $taxids->root->as_text('newick');
open(OUT, ">compara_63_taxids.nh");
print OUT $tax_str."\n";
close(OUT);

my $plot_tree_file = "compara_63.pdf";

my $cmd = qq^
library(devtools)
load_all("~/lib/greg-ensembl/projects/ggphylo", reset=TRUE)
  phylo <- tree.read.nhx("compara_63_tree.nhx")

  phylo <- tree.scale.to(phylo, 1)
  phylo <- tree.normalize.branchlengths(phylo, push.to.tips=F)
  phylo <- ladderize(phylo)

  pdf(file="${plot_tree_file}", width=8, height=12)
  p <- ggphylo(phylo,
    do.plot=F,
    x.expand = c(0.5, 0.5),
    internal.label.angle = 10,
    line.alpha = 0.5,
    line.color.by='low_coverage',
    line.color.scale = scale_colour_gradient("Sequencing Coverage", low='black', high='red', breaks=c(0, 1), labels=c('High', 'Low'))
  )
  p <- p + theme_bw()
  print(p)
  dev.off()
^;

Bio::Greg::EslrUtils->run_r($cmd);
