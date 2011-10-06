package Bio::Greg::AlnScore::AlnScore;

use strict;
use Bio::EnsEMBL::Hive::Process;
use Bio::Greg::Hive::Process;
use File::Path;
use File::Copy;
use File::Basename;
use FreezeThaw qw(freeze thaw cmpStr safeFreeze cmpStrHard);
use List::Util 'shuffle';

use base (
  'Bio::Greg::Hive::Process',
  'Bio::Greg::Slrsim::LoadTrees',
  'Bio::Greg::Hive::PhyloSim',
  'Bio::Greg::Hive::Align',
  'Bio::Greg::Hive::AlignmentScores',
  );

sub param_defaults {
  my $self = shift;

  my $phylosim = Bio::Greg::Hive::PhyloSim::param_defaults;
  my $align = Bio::Greg::Hive::Align::param_defaults;
  my $alignment_scores = Bio::Greg::Hive::AlignmentScores::param_defaults;

  return $self->replace($phylosim, $align, $alignment_scores, {
    force_recalc => 0
                        });
}

sub fetch_input {
  my $self = shift;

  $self->load_all_params();
  $self->create_table_from_params( $self->compara_dba, 'aln',
                                   $self->_aln_table_structure );
  $self->dbc->disconnect_when_inactive(1);
}

sub run {
  my $self = shift;

  $self->load_simulation_params;

  #  - Load & rescale tree
  my $tree_file = $self->param('tree');
  my $mpl = $self->param('mpl');
  my $indel = $self->param('indel');

  $self->param('phylosim_insertrate', $indel/2);
  $self->param('phylosim_deleterate', $indel/2);

  $tree_file = Bio::Greg::EslrUtils->baseDirectory."/projects/alnscore/".$tree_file;
  print "$tree_file\n";
  my $tree = Bio::EnsEMBL::Compara::TreeUtils->from_file($tree_file);
  $tree->distance_to_parent(0);

  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);

  $treeI->root->scale_mean_path_to($mpl);
  #$tree = Bio::EnsEMBL::Compara::TreeUtils->scale_mean_to( $tree, $mpl );

  print $treeI->ascii(1,1,0);
  print $treeI->root->newick_format."\n";

  #  - Simulate alignment
  my $defaults = $self->scheme_b;
  delete $defaults->{aligner};
  $self->set_params($defaults);

  $self->param('orig_mpl', $self->param('mpl'));
  $self->param('job_id', $self->input_job->dbID);

  my $sim_obj = $self->_simulate_alignment($treeI);
  my $true_aln = $sim_obj->{aln};
  my $true_pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($true_aln);

  $self->pretty_print($true_pep_aln, {full => 1});
  
  # Choose sizes at which to align / subset / score the subsets.
  #my @sizes = (2, 4, 6, 10, 16, 24, 32);
  my @sizes = (2, 6, 16, 24, 32);

  my $whole_aln = $self->_align($treeI, $true_aln, $true_pep_aln, 'all');

  foreach my $size (@sizes) {
    my ($sub_treeI, $sub_true_aln, $sub_true_pep_aln) = $self->_get_subsets($treeI, $true_aln, $true_pep_aln, $size);

    # Infer sub-alignment
    print "  inferring sub-alignment\n";
    my $sub_aln_inferred = $self->_align($sub_treeI, $sub_true_aln, $sub_true_pep_aln, 'n'.$size);

    # Extract sub-alignment
    print "  extracting sub-alignment\n";
    my $sub_aln_subsetted = $self->_subset_aln($whole_aln, $sub_treeI);

    # Score sub-alignments
    print $sub_treeI->ascii(1,1,0);

    print "  scoring inferred\n";
    my $scores_inf = $self->_score_aln($sub_treeI, $sub_true_aln, $sub_aln_inferred, 'aln');

    print "  scoring extracted\n";
    my $scores_sub = $self->_score_aln($sub_treeI, $sub_true_aln, $sub_aln_subsetted, 'sub');

    my $params = $self->replace($self->params, $scores_inf, $scores_sub);
    $self->store_params_in_table($self->dbc, 'aln', $params);
  }
}

sub _simulate_alignment {
  my $self = shift;
  my $treeI = shift;

  my $tree = Bio::EnsEMBL::Compara::TreeUtils->from_treeI($treeI);

  my $out_f = $self->_save_file('sim', 'perlobj');
  my $out_file = $out_f->{full_file};
  $self->param('sim_file', $out_file);

  if (!-e $out_file || $self->param('force_recalc')) {
    print "  running simulation\n";
    my $obj = $self->simulate_alignment($tree);
    $self->frz($out_file, $obj);
  }

  print "  loading simulated aln from file [$out_file]\n";
  my $obj = $self->thw($out_file);
  return $obj;
}

sub _align {
  my $self = shift;
  my $treeI = shift;
  my $true_aln = shift;
  my $true_pep_aln = shift;
  my $lbl = shift;

  my $tree = Bio::EnsEMBL::Compara::TreeUtils->from_treeI($treeI);

  my $out_f = $self->_save_file($lbl, 'fasta');
  my $out_file = $out_f->{full_file};

  if (!-e $out_file || $self->param('force_recalc')) {
    print "  running alignment for $lbl\n";
    # Calls Bio::Greg::Hive::Align->align
    my $inferred_aln = $self->align($tree, $true_aln, $true_pep_aln);
    Bio::EnsEMBL::Compara::AlignUtils->to_file($inferred_aln, $out_file);
  }

  print "  loading inferred aln from file [$out_file]\n";
  my $inferred_aln = Bio::EnsEMBL::Compara::AlignUtils->from_file($out_file);
  return $inferred_aln;
}

sub _get_subsets {
  my $self = shift;
  my $tree = shift;
  my $true_aln = shift;
  my $true_pep_aln = shift;
  my $n_keepers = shift;

  my $sub_tree_f = $self->_save_file('sub_'.$n_keepers, 'nh');
  my $sub_tree_file = $sub_tree_f->{full_file};

  if (!-e $sub_tree_file || $self->param('force_recalc')) {
    my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
    my @leaves = $treeI->leaves;
    my @shuffled_leaves = shuffle(@leaves);
    my @first_n = @shuffled_leaves[0..($n_keepers-1)];

    my $sub_treeI = $treeI->root->slice_with_internals(@first_n);
    #print $sub_treeI->ascii(1,1,0);
    $sub_treeI->contract_linear_paths(1);

    if ($sub_treeI->child_count == 1) {
      # We want to make the bifurcating node the new root, and give all left-over
      # branch length to the two next children.
      my $root = $sub_treeI;
      my $root_bl = $root->branch_length;
      my ($new_root) = $sub_treeI->children;
      $new_root->lengthen($root_bl);
      $sub_treeI = new Bio::Tree::Tree(-root => $new_root);
    }

    open(OUT, ">$sub_tree_file");
    print OUT $sub_treeI->to_newick."\n";
    close(OUT);
    #my $sub_tree = Bio::EnsEMBL::Compara::TreeUtils->extract_subtree_from_leaf_objects($tree, \@first_n);
    #Bio::EnsEMBL::Compara::TreeUtils->to_file($sub_tree, $sub_tree_file);
  }

  my $sub_treeI = Bio::Tree::Tree->from_file($sub_tree_file);
  
  my $sub_tree = Bio::EnsEMBL::Compara::TreeUtils->from_file($sub_tree_file);
  my $sub_aln = Bio::EnsEMBL::Compara::ComparaUtils->restrict_aln_to_tree($true_aln, $sub_tree);
  my $sub_pep_aln = Bio::EnsEMBL::Compara::ComparaUtils->restrict_aln_to_tree($true_pep_aln, $sub_tree);
  
  $sub_aln = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns($sub_aln);
  $sub_pep_aln = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns($sub_pep_aln);

  return ($sub_treeI, $sub_aln, $sub_pep_aln);
}

sub _subset_aln {
  my $self = shift;
  my $whole_aln = shift;
  my $sub_tree = shift;

  my $sub_aln = Bio::EnsEMBL::Compara::ComparaUtils->restrict_aln_to_tree($whole_aln, $sub_tree);
  $sub_aln = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns($sub_aln);
  return $sub_aln;
}

sub _score_aln {
  my $self = shift;
  my $treeI = shift;
  my $true_aln = shift;
  my $aln = shift;
  my $prefix = shift;

  my $true_pep_aln = $self->_tx_aln($true_aln);
  my $pep_aln = $self->_tx_aln($aln);

  $pep_aln = Bio::EnsEMBL::Compara::AlignUtils->sort_by_tree($pep_aln, $treeI);
  $self->pretty_print($pep_aln);

  print "  tcs\n";
  my $tcs = Bio::EnsEMBL::Compara::AlignUtils->total_column_score($true_pep_aln, $pep_aln);
  print "  sps\n";
  my $sps = Bio::EnsEMBL::Compara::AlignUtils->sum_of_pairs_score($true_pep_aln, $pep_aln);

  print "  other\n";
  my $aln_scores_calc = Bio::EnsEMBL::Compara::AlignUtils->_correct_subtree_calc($treeI, $true_pep_aln, $pep_aln);
  my $cbl = $aln_scores_calc->{complete_match_bl} / $aln_scores_calc->{aligned_bl};
  my $pbl = $aln_scores_calc->{partial_bl} / $aln_scores_calc->{aligned_bl};

  my $total_bl = $treeI->root->total_branch_length;
  my $mpl = $treeI->root->mean_path_length;
  my $mean_bl = $treeI->root->mean_branch_length;

  my $params = {
    $prefix.'_tcs' => $tcs,
    $prefix.'_sps' => $sps,
    $prefix.'_cbl' => $cbl,
    $prefix.'_pbl' => $pbl,
    n_seqs => scalar($treeI->leaves),
    mpl => sprintf("%.3f", $mpl),
    mean_bl => $mean_bl,
    total_bl => $total_bl
  };
  return $params;
}

sub _save_file {
  my $self = shift;
  my $filename_base = shift;
  my $ext = shift;

  my $rep = $self->param('replicate');
  my $id = $self->param('label');

  my $filename = "${id}_${filename_base}_${rep}";

  my $file_params = {
    id => $id,
    filename => $filename,
    extension => $ext,
    subfolder => 'data',
  };

  my $file_obj = $self->save_file($file_params);
  return $file_obj;
}


sub _aln_table_structure {
  my $self = shift;

  return {
    job_id => 'int',
    data_id => 'int',
    replicate => 'int',

    tree => 'char16',
    aligner => 'char16',
    n_seqs => 'int',
    mpl => 'float',
    orig_mpl => 'float',

    total_bl => 'float',
    mean_bl => 'float',
    
    sub_sps => 'float',
    aln_sps => 'float',
    sub_tcs => 'float',
    aln_tcs => 'float',
    sub_pbl => 'float',
    aln_pbl => 'float',
    sub_cbl => 'float',
    aln_cbl => 'float',

    unique_keys => 'data_id,replicate,n_seqs'

  };
}

1;
