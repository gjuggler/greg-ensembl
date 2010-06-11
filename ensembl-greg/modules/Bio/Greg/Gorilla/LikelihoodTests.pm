package Bio::Greg::Gorilla::LikelihoodTests;

use strict;
use Bio::Greg::Codeml;

use base ('Bio::Greg::Hive::Process');

my $TREE = 'Bio::EnsEMBL::Compara::TreeUtils';

sub fetch_input {
  my ($self) = @_;

  ### DEFAULT PARAMETERS ###
  my $params = {
  };
  ##########################

  # Fetch parameters from all possible locations.
  $self->load_all_params($params);

}


sub run {
  my $self = shift;

  my $full_tree = $self->get_tree;
  my $tree = $self->good_primate_tree($full_tree);

  print "Good tree: \n";
  print $tree->newick_format."\n";

  # a: (H, G, others)
  # b: (H#1, G, others)
  # c: (H, G#1, others)
  # d: (H#1, G#2, others)
  # e: (H#1, G#1, others)

  my $labeled_tree;
  my $labeled_tree = $self->categorize_nodes($tree,{});
  my $model_a = $self->calculate_branch_likelihood($tree,$labeled_tree);
  $self->store($model_a,'a');

  my $categories = {9606 => '1'};
  my $model_b = $self->calculate_branch_likelihood($tree,$categories);
  $self->store($model_b,'b');

  $categories = {9593 => '1'};
  my $model_c = $self->calculate_branch_likelihood($tree,$categories);
  $self->store($model_c,'c');

  $categories = {9593 => '1', 9606 => '2'};
  my $model_d = $self->calculate_branch_likelihood($tree,$categories);
  $self->store($model_d,'d');

  $categories = {9593 => '1', 9606 => '1'};
  my $model_e = $self->calculate_branch_likelihood($tree,$categories);
  $self->store($model_e,'e');

}

sub store {
  my $self = shift;
  my $model = shift;
  my $prefix = shift;

  my $key = $prefix."_";
  $self->store_tag("${key}lnL",$model->{lnL});

  my @omegas = @{$model->{omegas}};
  for (my $i=0; $i < scalar @omegas; $i++) {
    my $omega = $omegas[$i];
    $self->store_tag("${key}omega_${i}",$omega);
  }

}

sub model_string {
  my $obj = shift;
  return sprintf(" lnL: %s\n omegas: %s\n",$obj->{lnL},join(",",@{$obj->{omegas}}));
}

sub calculate_branch_likelihood {
  my $self = shift;
  my $tree = shift;
  my $taxon_categories = shift;

  # Set the root branch length to zero.
  $tree->distance_to_parent(0);

  my $tmpdir = $self->worker_temp_directory;

  my $params = $self->params();
  $params->{tree} = $tree;
  my $aln = $self->get_cdna_aln($params);
  Bio::EnsEMBL::Compara::AlignUtils->pretty_print($aln);
  delete $params->{tree};

  # Label categories in tree.
  my $labeled_tree = $self->categorize_nodes($tree,$taxon_categories);

  my $treeI = $TREE->to_treeI($labeled_tree);

  my $newick = $labeled_tree->newick_format;
  if ($newick !~ m/#/) {
    print "Only one omega category!\n";
    $params->{model} = 0;
  }

  my $obj = Bio::Greg::Codeml->branch_model_likelihood($treeI,$aln,$tmpdir,$params);

  return $obj;
}

sub categorize_nodes {
  my $self = shift;
  my $tree = shift;
  my $taxon_categories = shift;

  my $copy = $tree->copy;

  foreach my $leaf ($copy->leaves) {
    my $category = $taxon_categories->{$leaf->taxon_id} || '';
    $leaf->stable_id($leaf->stable_id."#$category") if (defined $category && $category ne '');
  }

  return $copy;
}

sub good_primate_tree {
  my $self = shift;
  my $tree = shift;

  my @good_primates = (9606,9593,9600,9544,9483);

  my @keeper_leaves = $TREE->get_leaves_for_species($tree,\@good_primates);
  my @keeper_ids = map {$_->node_id} @keeper_leaves;
  $tree = $TREE->extract_subtree_from_leaves($tree,\@keeper_ids);
  return $tree;
}


sub write_output {
  my $self = shift;

}


1;
