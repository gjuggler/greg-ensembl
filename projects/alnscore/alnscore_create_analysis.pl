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

use Bio::TreeIO;

# First create a database for the experiment.
my $mysql_base;
if (Bio::Greg::EslrUtils->is_ebi) {
  $mysql_base = 'mysql://slrsim:slrsim@mysql-greg.ebi.ac.uk:4134/';
} else {
  $mysql_base  = 'mysql://ensadmin:ensembl@ens-research/';
}
print "  MySQL base [$mysql_base]\n";

my $url = $mysql_base . "gj1_alnscore";
my $h = new Bio::Greg::Hive::ComparaHiveLoaderUtils;
$h->init($url);

print "  cleaning tables...\n";
$h->clean_hive_tables;
my @truncate_tables = qw^
      aln
  ^;
$h->truncate_tables(\@truncate_tables);

make_ladder();
make_balanced();
aln_score();

sub aln_score {
  my $logic_name  = "AlnScore";
  my $module      = "Bio::Greg::AlnScore::AlnScore";
  my $params      = {};
  $h->create_analysis( $logic_name, $module, $params, 500, 1 );

#  my @aligners = ('clustalw', 'mafft', 'muscle', 'probcons', 'fsa', 'prank', 'prank_codon', 'pagan');
  my @aligners = ('mafft', 'clustalw', 'probcons', 'fsa', 'prank_codon', 'pagan');
  #my @aligners = ('mafft');
  my @mpls = (1.0, 2.0, 0.2, 0.5, 3.0, 6.0);
  #my @mpls = (0.5);
  #my @trees = ('balanced.nh', 'ladder.nh');
  my @trees = ('balanced.nh', 'ladder.nh');

  my @n_domains = (1, 2, 4, 8);
  my @linker_ratios = (0.2);
  my @domain_ins_rate_mults = (0.1);

  #my $n_reps = 1;
  my $n_reps = 100;

  my $indel_rate = 0.1;

  my $i=0;
  foreach my $tree (@trees) {
    foreach my $n_domain (@n_domains) {
      foreach my $linker_ratio (@linker_ratios) {
        foreach my $domain_ins_rate_mult (@domain_ins_rate_mults) {
          foreach my $mpl (@mpls) {
            foreach my $aligner (@aligners) {
              foreach my $rep (1 .. $n_reps) {
                $i++;
                my $cur_indel = $indel_rate;
                $cur_indel = $indel_rate / $mpl / 2 if ($mpl >= 3);
                $cur_indel = sprintf("%.4f", $cur_indel);
                my $p = {
                  data_id => $i,
                  replicate => $rep,
                  random_seed => int(rand(999999)),

                  aligner => $aligner,
                  mpl => $mpl,
                  tree => $tree,
                  indel => $cur_indel,

                  n_domains => $n_domain,
                  linker_ratio => $linker_ratio,
                  domain_ins_rate_mult => $domain_ins_rate_mult,
                };
                $h->add_job_to_analysis("AlnScore", $p);
              }
            }
          }
        }
      }
    }
  }
}

sub make_ladder {
  # Create a ladder-like tree.
  my $node = new Bio::Tree::Node();
  my $tree = new Bio::Tree::Tree(-root => $node);
  my $n = 32;
  foreach my $i (1 .. $n) {
    if ($i < $n) {
      my $child = new $node;
      $child->name($i);
      $child->branch_length($n - $i);
      $node->add_child($child);
      
      my $internal_child = new $node;
      $node->add_child($internal_child);
      $node->branch_length(1);
      $node = $internal_child;
    } else {
      $node->name($i);
    }
  }
  print $tree->ascii(1, 1, 0)."\n";
  my $newick_str = $tree->to_newick;
  open(OUT, ">ladder.nh");
  print OUT $newick_str."\n";
  close(OUT);
  $tree->root->scale_mean_path_to(1);
  print $tree->root->total_branch_length."\n";
}

sub make_balanced {

  # Create a balanced tree.
  my $node = new Bio::Tree::Node();
  my $tree = new Bio::Tree::Tree(-root => $node);
  my $n = 32;
  my $cur_n = 1;
  my @outer_nodes = ($node);
  while (scalar(@outer_nodes) < $n) {
    foreach my $outer_node (@outer_nodes) {
      my $child_a = new $outer_node;
      my $child_b = new $outer_node;
      $child_a->branch_length(1);
      $child_b->branch_length(1);
      $outer_node->add_child($child_a);
      $outer_node->add_child($child_b);
    }
    @outer_nodes = $tree->leaves;
  }
  my $i=1;
  foreach my $leaf ($tree->leaves) {
    $leaf->name($i++);
  }
  print $tree->ascii(1, 1, 0)."\n";
  my $newick_str = $tree->to_newick;
  open(OUT, ">balanced.nh");
  print OUT $newick_str."\n";
  close(OUT);
  $tree->root->scale_mean_path_to(1);
  print $tree->root->total_branch_length."\n";
}


########*********########
#-------~~~~~~~~~-------#
########*********########
