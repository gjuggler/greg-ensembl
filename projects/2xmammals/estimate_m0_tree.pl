#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use Bio::EnsEMBL::Compara::TreeUtils;
use Bio::EnsEMBL::Compara::AlignUtils;
use Bio::Greg::Hive::PhyloAnalysis;

my $tree_in;
my $model;
my $tree_out;
my $nhx_out;
my $file_out;
my $aln_in;
my $clean;
GetOptions('tree=s' => \$tree_in,
           'model=s' => \$model,
           'tree_out=s' => \$tree_out,
           'nhx_out=s' => \$nhx_out,
           'file_out=s' => \$file_out,
           'aln=s' => \$aln_in,
           'clean' => \$clean);

print "$tree_in\n$model\n$tree_out\n$nhx_out\n$file_out\n$aln_in\n";

my $phylo = new Bio::Greg::Hive::PhyloAnalysis;
my $tree = Bio::EnsEMBL::Compara::TreeUtils->from_file($tree_in);
my $lines;

if (!-e $file_out || $clean) {
  my $aln = Bio::EnsEMBL::Compara::AlignUtils->from_file($aln_in);
  my $pep_aln = $phylo->_tx_aln($aln);
  print $tree->ascii;
  $phylo->pretty_print($pep_aln);

  my $params;
  if ($model eq 'm0') {
    $params = {
      model => 0,
      fix_blength => 0,
      method => 1,
      cleandata => 0,
      getSE => 0,
      RateAncestor => 0,
      Small_Diff => 1e-6
    };
  } else {
    $params = {
      model => 1,
      fix_blength => 2,
      method => 1,
      cleandata => 0,
      RateAncestor => 0,
      getSE => 1,
      Small_Diff => 1e-6
    };      
  }
  print "  running PAML $model...\n";
  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
  foreach my $node ($treeI->nodes) {
    if (!$node->is_leaf) {
      $node->name('');
    }
  }
  my $results = Bio::Greg::Codeml->branch_model_likelihood( $treeI, $aln, $phylo->worker_temp_directory, $params );
  $lines = $results->{lines};

  open(my $fh, ">", $file_out) or die $!;
  print $fh join("", @{$lines})."\n";
  close $fh;
} else {

  open(my $fh, $file_out);
  my @in_lines = <$fh>;
  $lines = \@in_lines;
  close($fh);
}

my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
my $codeml_tree = Bio::Greg::Codeml->parse_codeml_results($lines);
Bio::EnsEMBL::Compara::TreeUtils->transfer_branchlengths($codeml_tree, $treeI);
$treeI->root->branch_length(0.01); # Set small b.l. on the root node.

my $newick_str = $treeI->as_text('newick');
print $newick_str."\n";

open(my $fh, '>', $tree_out) or die $!;
print $fh $newick_str."\n";
close($fh);

my @nodes = $codeml_tree->nodes;
foreach my $node (@nodes) {
  $node->remove_tag('substitutions');
  $node->remove_tag('parent');
  $node->remove_tag('node_a');
  $node->remove_tag('node_b');
  $node->remove_tag('id');
  $node->remove_tag('branch_label');
}
Bio::EnsEMBL::Compara::TreeUtils->transfer_annotations($codeml_tree, $treeI);
my $nhx_str = $treeI->as_text('nhx');
print $nhx_str."\n";

open(OUT, '>', $nhx_out);
print OUT $nhx_str."\n";
close(OUT);

exit(0);
