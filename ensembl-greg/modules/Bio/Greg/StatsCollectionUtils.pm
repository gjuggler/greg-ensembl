package Bio::Greg::StatsCollectionUtils;

use strict;

use Bio::EnsEMBL::Compara::NestedSet;

our @ISA = qw();

sub max_path {
  my $class = shift;
  my $tree = shift;
  my $max = $tree->max_distance;
  return sprintf "%.3f", $max;
}

sub tree_length {
  my $class = shift;
  my $tree = shift;
  my $dist = 0;
  map {$dist += $_->distance_to_parent} $tree->nodes;
  return sprintf "%.3f", $dist;
}

sub mean_path {
  my $class = shift;
  my $tree = shift;
  my $dist = 0;
  map {$dist += $class->dist_to_root($_);} $tree->leaves;
  $dist = $dist / scalar($tree->leaves);
  return sprintf "%.3f", $dist;
}

sub dist_to_root {
  my $class = shift;
  my $leaf = shift;

  my $d = $leaf->distance_to_parent;
  my $p = $leaf->parent;
  while ($p) {
    $d += $p->distance_to_parent;
    $p = $p->parent;
  }
  return $d;
}

sub max_branch {
  my $class = shift;
  my $tree = shift;

  my $max = 0;
  foreach my $node ($tree->nodes) {
    $max = $node->distance_to_parent if ($node->distance_to_parent > $max);
  }
  return $max;
}

sub mean_branch {
  my $class = shift;
  my $tree = shift;

  my $total = 0;
  foreach my $node ($tree->nodes) {
    $total += $node->distance_to_parent;
  }
  return $total / scalar($tree->nodes);
}

1;
