package Bio::Greg::Gorilla::LikelihoodTests;

use strict;
use Bio::Greg::Codeml;

use base ('Bio::Greg::Hive::Process');

my $TREE = 'Bio::EnsEMBL::Compara::TreeUtils';

sub fetch_input {
  my ($self) = @_;

  my $default_params = {
  };

  # Fetch parameters from all possible locations.
  $self->load_all_params($default_params);
}


sub run {
  my $self = shift;

  my $full_tree = $self->get_tree;
  my $tree = $self->good_primate_tree($full_tree);

  print "Base tree: \n";
  print $tree->newick_format."\n";

  # Branch models summary:
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

sub calculate_branch_likelihood {
  my $self = shift;
  my $tree = shift;
  my $taxon_categories = shift;

  # Set the branch length of the root node to zero.
  $tree->distance_to_parent(0);

  my $tmpdir = $self->worker_temp_directory;

  # A little funkiness here to make sure we get the correct alignment for our chosen sub-tree.
  my $params = $self->params();
  $params->{tree} = $tree;
  my $aln = $self->get_cdna_aln($params);
  Bio::EnsEMBL::Compara::AlignUtils->pretty_print($aln); # Debugging print the alignment.
  delete $params->{tree};

  # Label categories in tree.
  my $labeled_tree = $self->categorize_nodes($tree,$taxon_categories);

  # Turn the labeled Ensembl tree object into a BioPerl TreeI-compliant object.
  my $treeI = $TREE->to_treeI($labeled_tree);

  # Detect the case when there's only one rate category, and switch to 'model=0' for PAML.
  my $newick = $labeled_tree->newick_format;
  if ($newick !~ m/#/) {
    print "Only one omega category!\n";
    $params->{model} = 0;
  }

  # Use the Codeml.pm helper function to run Codeml.
  my $obj = Bio::Greg::Codeml->branch_model_likelihood($treeI,$aln,$tmpdir,$params);

  return $obj;
}

# Uses a mapping from ncbi_taxon_id => rate_category to create a PAML-formatted
# branch model tree for terminal branches of various species.
sub categorize_nodes {
  my $self = shift;
  my $tree = shift;
  my $taxon_categories = shift;

  # Create a copy of the tree.
  my $copy = $tree->copy;

  foreach my $leaf ($copy->leaves) {
    # If this leaf's taxon_id exists in the mapping, append that value (along with '#') to the label.
    # Otherwise, leave it alone.
    my $category = $taxon_categories->{$leaf->taxon_id} || '';
    $leaf->stable_id($leaf->stable_id."#$category") if (defined $category && $category ne '');
  }

  return $copy;
}

# Extracts the desired sub-tree with only the species of interest.
# Assumption: the previous step in the pipeline, FilterOneToOneOrthologs.pm,
# should have left us only with trees containing perfect 1-to-1 orthology
# in these species.
sub good_primate_tree {
  my $self = shift;
  my $tree = shift;

  # NCBI taxon IDs of desired species.
  my @good_primates = (9606,9593,9600,9544,9483);

  # This method from Bio::EnsEMBL::Compara::TreeUtils extracts ALL the leaves for
  # a given list of taxon IDs.
  my @keeper_leaves = $TREE->get_leaves_for_species($tree,\@good_primates);

  # Turn our Bio::EnsEMBL::Compara::Member objects into a list of (database-specific) node_ids.
  my @keeper_ids = map {$_->node_id} @keeper_leaves;

  # Helper method to generate a nice sub-tree from a list of node_ids.
  $tree = $TREE->extract_subtree_from_leaves($tree,\@keeper_ids);
  return $tree;
}

# Simple string format for debugging.
sub model_string {
  my $obj = shift;
  return sprintf(" lnL: %s\n omegas: %s\n",$obj->{lnL},join(",",@{$obj->{omegas}}));
}


# Store the lnL and omegas from the model object as tree tags.
# These will be collected into the stats_lnl table when the
# CollectLikelihoodStats.pm module is run.
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

1;
