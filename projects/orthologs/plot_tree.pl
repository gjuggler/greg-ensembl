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

my $url = 'mysql://ensadmin:ensembl@ens-research/gj1_orthologs_63';
my $clean = 1;

my $h = new Bio::Greg::Hive::ComparaHiveLoaderUtils;
$h->init($url);

my $compara_dba = $h->compara_dba;

my $tree = Bio::EnsEMBL::Compara::ComparaUtils->get_genome_tree_with_extras($compara_dba,
                                                                            {});
my $taxid_map;
map {$taxid_map->{$_->id} = $_->taxon_id if ($_->id ne '')} $tree->nodes;

my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
#print $treeI->ascii;

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
print $taxids->ascii;
my $tax_str = $taxids->root->as_text('newick');
open(OUT, ">compara_63_taxids.nh");
print OUT $tax_str."\n";
close(OUT);

my $plot_tree_file = "compara_63.pdf";

my $cmd = qq^
library(phylosim);
source("~/src/greg-ensembl/projects/phylosim/PhyloSimPlots.R")

  sim <- PhyloSim()
  str <- readLines(con=file("compara_63_tree.nhx"))
  phylo <- read.nhx.tree(str)

  phylo <- remove.branchlengths(phylo)
  phylo <- ladderize(phylo)

  sim\$.phylo <- phylo
  pdf(file="${plot_tree_file}", width=10, height=12)
  xlab <- ""
  obj <- plotTree(sim, tree.do.plot=F, color.by='low_coverage', tree.xlab=xlab, tree.xlim.expand=c(1, 2), line.width=1, include.internal.labels=T, internal.angle=15, internal.size=1, internal.color=rgb(.1, .1, .7))
  p <- obj[['grob']]
  p <- p + theme_bw()
  p <- p + scale_colour_gradient("Sequencing Coverage",
    low="black", high="red",
    breaks=c(0, 1),
    labels=c('High', 'Low')
  )
  p <- p + scale_y_continuous("")
  p <- p + opts(
    title="Ensembl Compara Species v63",
    axis.text.x = theme_blank(),
    axis.text.y = theme_blank()
  )
  print(p)
  dev.off()
^;

Bio::Greg::EslrUtils->run_r($cmd);
