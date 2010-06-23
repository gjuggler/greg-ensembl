package Bio::Greg::Hive::FilterOneToOneOrthologs;

use strict;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::Process');

my $TREE = 'Bio::EnsEMBL::Compara::TreeUtils';

sub fetch_input {
  my $self = shift;

  ### DEFAULT PARAMETERS ###
  my $defaults = {
    one_to_one_taxon_list => '9606'
  };
  ##########################

  $self->load_all_params($defaults);
  
  my $tree = $self->get_tree();
  $self->param('tree',$tree);
}

sub run {
  my $self = shift;

  my $tree = $self->param('tree');
  my $one_to_one_taxon_list = $self->param('one_to_one_taxon_list');

  print $tree->newick_format."\n";

  my @species_list = split(',',$one_to_one_taxon_list);
  print "@species_list\n";
  my @keeper_leaves = $TREE->get_leaves_for_species($tree,\@species_list);
  my @keeper_ids = map {$_->node_id} @keeper_leaves;
  my $pruned_tree = $TREE->extract_subtree_from_leaves($tree,\@keeper_ids);
  print $pruned_tree->newick_format."\n";
  print map {$_->stable_id." "} @keeper_leaves;
  print "\n";

  # Test whether it's a good 1-1-1 orthology.
  my $is_good_tree = 1;
  my %keeper_hash;
  map {$keeper_hash{$_->taxon_id}=1} @keeper_leaves;
  map {$is_good_tree = 0 if (!defined $keeper_hash{$_})} @species_list;
  $is_good_tree = 0 if ($#keeper_leaves != $#species_list);
  if (!$is_good_tree) {
    print " -> Bad tree! Not one-to-one orthology: Doing nothing!\n";
    $self->autoflow_inputjob(0);
    return;
  } else {
    print " -> Good tree! Flowing output job...\n";
  }
}

sub write_output {
  
}

1;
