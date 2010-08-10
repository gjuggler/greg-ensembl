package Bio::Greg::Gorilla::LikelihoodTests;

use strict;
use Bio::Greg::Codeml;

use base ( 'Bio::Greg::Hive::Process', 'Bio::Greg::StatsCollectionUtils',
  'Bio::Greg::Hive::Align' );

my $TREE = 'Bio::EnsEMBL::Compara::TreeUtils';

sub fetch_input {
  my ($self) = @_;

  my $default_params = {
    output_table            => 'stats_branch',
    chimp_output_table      => 'stats_chimp_branch',
    alignment_output        => 1,
    alignment_output_folder => 'branch_alns'
  };

  # Fetch parameters from all possible locations.
  $self->load_all_params($default_params);

  # Create tables if necessary.
  $self->create_table_from_params( $self->compara_dba, $self->param('output_table'),
    $self->get_gene_stats_def );
  $self->create_table_from_params( $self->compara_dba, $self->param('chimp_output_table'),
    $self->get_gene_stats_def );
}

sub run {
  my $self = shift;

  my $full_tree = $self->get_tree;
  my $tree_no_chimp = $self->good_primate_tree( $full_tree, [ 9606, 9593, 9600, 9544, 9483 ] );
  my $full_primate_tree =
    $self->good_primate_tree( $full_tree, [ 9606, 9593, 9598, 9600, 9544, 9483 ] );
  my $tree_no_human = $self->good_primate_tree( $full_tree, [ 9593, 9598, 9600, 9544, 9483 ] );

  print "full primates: \n";
  print $full_primate_tree->newick_format . "\n";
  print "no chimp: \n";
  print $tree_no_chimp->newick_format . "\n";
  print "no human: \n";
  print $tree_no_human->newick_format . "\n";

  my $full_primate_aln = $self->get_filtered_aln($full_primate_tree);

  # Branch models summary:
  # a: (H, G, others)
  # b: (H#1, G, others)
  # c: (H, G#1, others)
  # d: (H#1, G#2, others)
  # e: (H#1, G#1, others)
  # f: (((H,C),G)#1, others)
  # g: (((H,C),G)$1, others)
  # h: (((H,C),G), others)

  my $categories;
  my $table;
  my $tx_id;
  my $tree;
  my $aln;

  foreach my $non_gorilla ( 'chimpanzee', 'human' ) {
    print "NON-GORILLA SPECIES: $non_gorilla\n";

    if ( $non_gorilla eq 'human' ) {

      # Analyze with human and gorilla.
      $tx_id = 9606;
      $tree  = $tree_no_chimp;
      my ($chimp_member) = grep { $_->taxon_id == 9598 } $full_primate_tree->leaves;
      my $aln_no_chimp = Bio::EnsEMBL::Compara::AlignUtils->remove_seq_from_aln( $full_primate_aln,
        $chimp_member->stable_id );
      $aln   = $aln_no_chimp;
      $table = $self->param('output_table');
    } else {

      # Analyze with chimp and gorilla.
      $tx_id = 9598;
      $tree  = $tree_no_human;
      my ($human_member) = grep { $_->taxon_id == 9606 } $full_primate_tree->leaves;
      my $aln_no_human = Bio::EnsEMBL::Compara::AlignUtils->remove_seq_from_aln( $full_primate_aln,
        $human_member->stable_id );
      $aln   = $aln_no_human;
      $table = $self->param('chimp_output_table');
    }

    Bio::EnsEMBL::Compara::AlignUtils->pretty_print($aln);    # Debugging print the alignment.

    my $model_a = $self->calculate_branch_likelihood( $tree, $aln, {} );
    $self->store( $model_a, $table, 'a' );

    $categories = { $tx_id => '1' };
    my $model_b = $self->calculate_branch_likelihood( $tree, $aln, $categories );
    $self->store( $model_b, $table, 'b' );

    $categories = { 9593 => '1' };
    my $model_c = $self->calculate_branch_likelihood( $tree, $aln, $categories );
    $self->store( $model_c, $table, 'c' );

    $categories = { 9593 => '1', $tx_id => '2' };
    my $model_d = $self->calculate_branch_likelihood( $tree, $aln, $categories );
    $self->store( $model_d, $table, 'd' );

    $categories = { 9593 => '1', $tx_id => '1' };
    my $model_e = $self->calculate_branch_likelihood( $tree, $aln, $categories );
    $self->store( $model_e, $table, 'e' );

    $categories = { 9604 => '1' };
    my $model_f =
      $self->calculate_branch_likelihood( $full_primate_tree, $full_primate_aln, $categories );
    $self->store( $model_f, $table, 'f' );

    $categories = { 9604 => '$1' };
    my $model_g =
      $self->calculate_branch_likelihood( $full_primate_tree, $full_primate_aln, $categories );
    $self->store( $model_g, $table, 'g' );

    $categories = {};
    my $model_h =
      $self->calculate_branch_likelihood( $full_primate_tree, $full_primate_aln, $categories );
    $self->store( $model_h, $table, 'h' );

    $self->_get_peptides( $full_primate_tree, 9606, 'human' );
    $self->_get_peptides( $full_primate_tree, 9593, 'gorilla' );
    $self->_get_peptides( $full_primate_tree, 9598, 'chimpanzee' );
    $self->collect_values($tree);

    my $lambda_mammals  = $self->dawg_lambda($full_tree);
    my $lambda_primates = $self->dawg_lambda($full_primate_tree);
    $self->param( 'lambda_mammals',  $lambda_mammals );
    $self->param( 'lambda_primates', $lambda_primates );

    $self->store_params_in_table( $self->db_handle(1), $table, $self->params );
  }

}

sub _get_peptides {
  my $self     = shift;
  my $tree     = shift;
  my $taxon_id = shift;
  my $prefix   = shift;

  my @proteins = grep { $_->taxon_id == $taxon_id } $tree->leaves;
  if ( scalar @proteins > 0 ) {
    my $member = $proteins[0];
    $self->param( $prefix . '_protein', $member->stable_id );
    $self->param( $prefix . '_gene',    $member->gene_member->stable_id );
  } else {
    $self->param( $prefix . '_protein', undef );
    $self->param( $prefix . '_gene',    undef );
  }
}

sub collect_values {
  my $self = shift;
  my $tree = shift;

  # Collect gene tag values into the params hash.
  $self->param( 'tree_length',   $self->tree_length($tree) );
  $self->param( 'tree_max_path', $self->max_path($tree) );

}

sub get_filtered_aln {
  my $self = shift;
  my $tree = shift;

  my $file;

  # A little funkiness here to make sure we get the correct alignment for our chosen sub-tree.
  my $params = $self->params();
  $params->{tree} = $tree;
  my $aln      = $self->get_cdna_aln($params);
  my $orig_aln = $aln;

  #Bio::EnsEMBL::Compara::AlignUtils->pretty_print($aln); # Debugging print the alignment.
  delete $params->{tree};

  my @leaves = $tree->leaves;
  my ($member) = grep { $_->taxon_id == 9606 } @leaves;
  $member = $leaves[0] if ( !defined $member );
  my $aln_id = $member->stable_id;

  $file =
    $self->save_aln( $aln,
    { id => $aln_id, filename => $aln_id . "_orig", subfolder => 'branch_alns' } );
  $self->param( 'orig_aln_length', $aln->length );
  $self->param( 'orig_aln_file',   $file );

  # Align with Prank to try and de-align incorrectly called exons.
  my $pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($aln);
  my $prank_params = { alignment_prank_f => 1 };
  $pep_aln = $self->align_with_prank( $pep_aln, $tree, $prank_params );
  $aln = Bio::EnsEMBL::Compara::AlignUtils->peptide_to_cdna_alignment( $pep_aln, $tree );
  my $prank_aln = $aln;

  $file =
    $self->save_aln( $aln,
    { id => $aln_id, filename => $aln_id . "_realign", subfolder => 'branch_alns' } );
  $self->param( 'prank_aln_length', $aln->length );
  $self->param( 'prank_aln_file',   $file );

  # Remove 'funky' columns with triple substitutions between close pairs of species.
  my ($first_keeper)  = grep { $_->taxon_id == 9606 } $tree->leaves;
  my ($second_keeper) = grep { $_->taxon_id == 9593 } $tree->leaves;
  my ($third_keeper)  = grep { $_->taxon_id == 9598 } $tree->leaves;

  if ( $first_keeper && $second_keeper ) {
    $aln = Bio::EnsEMBL::Compara::AlignUtils->remove_triple_mutated_codons( $aln,
      $first_keeper->stable_id, $second_keeper->stable_id );
  }
  if ( $first_keeper && $third_keeper ) {
    $aln = Bio::EnsEMBL::Compara::AlignUtils->remove_triple_mutated_codons( $aln,
      $first_keeper->stable_id, $third_keeper->stable_id );
  }
  if ( $second_keeper && $third_keeper ) {
    $aln = Bio::EnsEMBL::Compara::AlignUtils->remove_triple_mutated_codons( $aln,
      $second_keeper->stable_id, $third_keeper->stable_id );
  }

  $file =
    $self->save_aln( $aln,
    { id => $aln_id, filename => $aln_id . "_final", subfolder => 'branch_alns' } );
  $self->param( 'filtered_aln_length', $aln->length );
  $self->param( 'filtered_aln_file',   $file );

  Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $orig_aln, { width => 200 } )
    ;    # Debugging print the alignment.
  Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $prank_aln, { width => 200 } )
    ;    # Debugging print the alignment.
  Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $aln, { width => 200 } )
    ;    # Debugging print the alignment.

  return $aln;
}

sub calculate_branch_likelihood {
  my $self             = shift;
  my $tree             = shift;
  my $aln              = shift;
  my $taxon_categories = shift;

  my $params = $self->params();

  # Set the branch length of the root node to zero.
  $tree->distance_to_parent(0);

  my $tmpdir = $self->worker_temp_directory;

  # Label categories in tree.
  my $labeled_tree = $self->categorize_nodes( $tree, $taxon_categories );

  print $labeled_tree->newick_format . "\n";

  # Turn the labeled Ensembl tree object into a BioPerl TreeI-compliant object.
  my $treeI = $TREE->to_treeI($labeled_tree);

  # Detect the case when there's only one rate category, and switch to 'model=0' for PAML.
  my $newick = $labeled_tree->newick_format;
  if ( $newick !~ m/(#|\$)/ ) {
    print "Only one omega category!\n";
    $params->{model} = 0;
  }

  # Use the Codeml.pm helper function to run Codeml.
  my $obj = Bio::Greg::Codeml->branch_model_likelihood( $treeI, $aln, $tmpdir, $params );

  return $obj;
}

# Uses a mapping from ncbi_taxon_id => rate_category to create a PAML-formatted
# branch model tree for terminal branches of various species.
sub categorize_nodes {
  my $self             = shift;
  my $tree             = shift;
  my $taxon_categories = shift;

  # Create a copy of the tree.
  my $copy = $tree->copy;

  foreach my $leaf ( $copy->nodes ) {

   # If this leaf's taxon_id exists in the mapping, append that value (along with '#') to the label.
   # Otherwise, leave it alone.
    my $category = $taxon_categories->{ $leaf->taxon_id } || '';
    $category = '#' . $category unless ( $category eq '' || $category =~ m/(#|\$)/ );
    if ( $leaf->is_leaf ) {
      $leaf->stable_id( $leaf->stable_id . "$category" )
        if ( defined $category && $category ne '' );
    } else {
      $leaf->name( $leaf->name . "$category" ) if ( defined $category && $category ne '' );
    }
  }

  return $copy;
}

# Extracts the desired sub-tree with only the species of interest.
# Assumption: the previous step in the pipeline, FilterOneToOneOrthologs.pm,
# should have left us only with trees containing perfect 1-to-1 orthology
# in these species.
sub good_primate_tree {
  my $self                  = shift;
  my $tree                  = shift;
  my $good_primate_arrayref = shift;

  $tree = $tree->copy;

  # NCBI taxon IDs of desired species.
  my @good_primates = @$good_primate_arrayref;

  # This method from Bio::EnsEMBL::Compara::TreeUtils extracts ALL the leaves for
  # a given list of taxon IDs.
  my @keeper_leaves = $TREE->get_leaves_for_species( $tree, \@good_primates );

  # Turn our Bio::EnsEMBL::Compara::Member objects into a list of (database-specific) node_ids.
  my @keeper_ids = map { $_->node_id } @keeper_leaves;

  # Helper method to generate a nice sub-tree from a list of node_ids.
  $tree = $TREE->extract_subtree_from_leaves( $tree, \@keeper_ids );
  return $tree;
}

# Simple string format for debugging.
sub model_string {
  my $obj = shift;
  return sprintf( " lnL: %s\n omegas: %s\n", $obj->{lnL}, join( ",", @{ $obj->{omegas} } ) );
}

# Store the lnL and omegas from the model object as tree tags.
# These will be collected into the stats_lnl table when the
# CollectLikelihoodStats.pm module is run.
sub store {
  my $self   = shift;
  my $model  = shift;
  my $table  = shift;
  my $prefix = shift;

  my $key = $prefix . "_";

  #$self->store_tag("${key}lnL",$model->{lnL});
  $self->param( "${key}lnL", $model->{lnL} );

  my @omegas = @{ $model->{omegas} };
  for ( my $i = 0 ; $i < scalar @omegas ; $i++ ) {
    my $omega = $omegas[$i];

    #$self->store_tag("${key}omega_${i}",$omega);
    $self->param( "${key}omega_${i}", $omega );
  }

}

sub get_gene_stats_def {
  my $gene_stats_def = {
    data_id      => 'int',
    node_id      => 'int',
    orig_node_id => 'string',

    human_protein      => 'string',
    human_gene         => 'string',
    gorilla_protein    => 'string',
    gorilla_gene       => 'string',
    chimpanzee_protein => 'string',
    chimpanzee_gene    => 'string',

    orig_aln_length     => 'int',
    prank_aln_length    => 'int',
    filtered_aln_length => 'int',

    orig_aln_file     => 'string',
    prank_aln_file    => 'string',
    filtered_aln_file => 'string',

    lambda_mammals  => 'float',
    lambda_primates => 'float',

    unique_keys => 'data_id,node_id'
  };
  foreach my $model ( 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h' ) {
    $gene_stats_def->{ $model . '_lnL' }     = 'float';
    $gene_stats_def->{ $model . '_omega_0' } = 'float';
    $gene_stats_def->{ $model . '_omega_1' } = 'float' unless ( $model eq 'a' );
    $gene_stats_def->{ $model . '_omega_2' } = 'float' if ( $model eq 'd' );
  }
  return $gene_stats_def;
}

1;
