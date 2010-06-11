package Bio::Greg::Gorilla::TopologyTests;

use strict;
use Bio::Greg::Codeml;
use Bio::Greg::Baseml;

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
  my $tree = $self->primate_tree($full_tree);

  print "Good tree: \n";
  print $tree->newick_format."\n";

  # a) ((H,C),G)
  # b) ((H,G),C)
  # c) ((C,G),H)
  my $newick;

  $newick = "(((H,C),G),O);";
  my $obj_a = $self->calculate_tree_likelihood($tree,$newick);
  $self->store($obj_a,'a');

  $newick = "(((H,G),C),O);";
  my $obj_b = $self->calculate_tree_likelihood($tree,$newick);
  $self->store($obj_b,'b');

  $newick = "(((G,C),H),O);";
  my $obj_c = $self->calculate_tree_likelihood($tree,$newick);
  $self->store($obj_c,'c');
  
}

sub store {
  my $self = shift;
  my $model = shift;
  my $prefix = shift;

  my $key = $prefix."_";
  $self->store_tag("${key}lnL",$model->{lnL});

#  my @omegas = @{$model->{omegas}};
#  for (my $i=0; $i < scalar @omegas; $i++) {
#    my $omega = $omegas[$i];
#    $self->store_tag("${key}omega_${i}",$omega);
#  }
}

sub tax_to_char {
  my $self = shift;
  my $tax = shift;

  my $taxon_to_letter = {
    9606 => 'H',
    9598 => 'C',
    9593 => 'G',
    9600 => 'O'
  };

  return $taxon_to_letter->{$tax};
}

sub calculate_tree_likelihood {
  my $self = shift;
  my $tree = shift;
  my $test_topology_newick = shift;

  my $tmpdir = "/tmp/pamltest";
  mkdir($tmpdir);
#  my $tmpdir = $self->worker_temp_directory;

  # Get the alignment.
  my $params = $self->params();
  $params->{tree} = $tree;
  my $aln = $self->get_cdna_aln($params);
  Bio::EnsEMBL::Compara::AlignUtils->pretty_print($aln,{length=>200});
  delete $params->{tree};

  # Translate the tree and alignment sequences to single-letters.
  my $map = {};
  foreach my $leaf ($tree->leaves) {
    $map->{$leaf->stable_id} = $self->tax_to_char($leaf->taxon_id);
  }
  $aln = Bio::EnsEMBL::Compara::AlignUtils->translate_ids($aln,$map);

  # Parse the tree into BioPerl object.
  my $treeI = $TREE->treeI_from_newick($test_topology_newick);

  Bio::EnsEMBL::Compara::AlignUtils->pretty_print($aln);
  my $obj = Bio::Greg::Baseml->dna_model_likelihood($treeI,$aln,$tmpdir,$params);

  return $obj;
}

sub primate_tree {
  my $self = shift;
  my $tree = shift;

  my @good_primates = (9606,9598,9593,9600);

  my @keeper_leaves = $TREE->get_leaves_for_species($tree,\@good_primates);
  my @keeper_ids = map {$_->node_id} @keeper_leaves;
  $tree = $TREE->extract_subtree_from_leaves($tree,\@keeper_ids);
  return $tree;
}


1;
