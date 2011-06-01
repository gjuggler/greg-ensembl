package Bio::Greg::Gorilla::LikelihoodTests;

use strict;
use Bio::Greg::Codeml;
use File::Path;

use base (
  'Bio::Greg::Hive::Process', 'Bio::Greg::StatsCollectionUtils',
  'Bio::Greg::Hive::Align',   'Bio::Greg::Hive::CountSubstitutions'
);

my $TREE = 'Bio::EnsEMBL::Compara::TreeUtils';

sub param_defaults {
  return {
    output_table                        => 'stats_branch',
    alignment_output                    => 1,
    alignment_output_folder             => 'branch_alns_nofilter',
    aln_type                            => 'compara',
    quality_threshold                   => 30,
    fail_on_altered_tree                => 1,
    likelihood_realign_with_prank       => 0,
    likelihood_filter_substitution_runs => 1,
    species_taxon_ids                   => '9606, 9593, 9598, 9600, 9544, 9483'
  };
}

sub fetch_input {
  my ($self) = @_;

  # Fetch parameters from all possible locations.
  $self->load_all_params();

  $self->param( 'Hsap_output_table', $self->param('output_table') . '_Hsap' );
  $self->param( 'Ptro_output_table', $self->param('output_table') . '_Ptro' );

  # Create tables if necessary.
  $self->create_table_from_params( $self->compara_dba, $self->param('Hsap_output_table'),
    $self->get_gene_stats_def );
  $self->create_table_from_params( $self->compara_dba, $self->param('Ptro_output_table'),
    $self->get_gene_stats_def );

  my $actual_folder = $self->get_output_folder;

  my $aln_subdir      = $self->param('alignment_output_folder');
  my $aln_root_folder = "${actual_folder}/${aln_subdir}";
  mkpath( [$aln_root_folder] );

}

sub list_from_string {
  my $self   = shift;
  my $string = shift;
  return Bio::EnsEMBL::Compara::ComparaUtils->get_taxon_ids_from_keepers_list( $self->compara_dba,
    $string );
}

sub subtract {
  my $self      = shift;
  my $list_a    = shift;
  my @remove_us = @_;
  my $hash;
  map { $hash->{$_} = 1 } @$list_a;
  foreach my $list_b (@remove_us) {
    map { delete $hash->{$_} } @$list_b;
  }
  return [ keys %$hash ];
}

sub run {
  my $self = shift;

  my @taxon_ids = $self->list_from_string( $self->param('species_taxon_ids') );

  #print "@taxon_ids\n";
  my $all = \@taxon_ids;

  my $full_tree = $self->get_tree;
  my $full_primate_tree = $self->good_primate_tree( $full_tree, $all );

  $self->_get_peptides( $full_primate_tree, 9606, 'Hsap' );
  $self->_get_peptides( $full_primate_tree, 9593, 'Ggor' );
  $self->_get_peptides( $full_primate_tree, 9598, 'Ptro' );

  # Get GC content and stuff, using the human sequence.
  $self->collect_stats( $full_tree, 9606 );

  print "full primate tree: \n";
  print $full_primate_tree->newick_format . "\n";

  my $tree_aln_obj = $self->get_filtered_aln($full_primate_tree);
  $full_primate_tree = $tree_aln_obj->{tree};
  my $full_primate_aln = $tree_aln_obj->{aln};

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

  foreach my $non_gorilla ( 9606, 9598 ) {
    my ($non_gorilla_member) = grep { $_->taxon_id == $non_gorilla } $full_primate_tree->leaves;
    my ($gorilla_member)     = grep { $_->taxon_id == 9593 } $full_primate_tree->leaves;
    $self->param( 'gorilla_member',     $gorilla_member );
    $self->param( 'non_gorilla_member', $non_gorilla_member );

    my $gor_name = $gorilla_member->taxon->short_name;
    my $name     = $non_gorilla_member->taxon->short_name;

    $tx_id = $non_gorilla;
    my $remove_taxon = 9598;
    if ( $non_gorilla == 9598 ) {
      $remove_taxon = 9606;
    }

    $tree = $self->good_primate_tree( $full_tree, $self->subtract( $all, [$remove_taxon] ) );
    $aln = Bio::EnsEMBL::Compara::AlignUtils->remove_seq_from_aln( $full_primate_aln,
      'ens_' . $remove_taxon );
    $table = $self->param('output_table') . '_' . $name;

    my $map;
    map { $map->{ $_->node_id } = 'ens_' . $_->taxon_id } $tree->nodes;
    $tree = Bio::EnsEMBL::Compara::TreeUtils->translate_ids( $tree, $map );

    # Store the reduced tree and alignment to use in the store() method to collect subst. info.
    $self->param( 'tree',     $tree );
    $self->param( 'cdna_aln', $aln );
    $self->param( 'aa_aln',   Bio::EnsEMBL::Compara::AlignUtils->translate($aln) );

    my $model_a = $self->calculate_branch_likelihood( $tree, $aln, {} );
    $self->store( $tree, $model_a, $table, 'a' );

    $categories = { $tx_id => '1' };
    my $model_b = $self->calculate_branch_likelihood( $tree, $aln, $categories );
    $self->store( $tree, $model_b, $table, 'b' );

    $categories = { 9593 => '1' };
    my $model_c = $self->calculate_branch_likelihood( $tree, $aln, $categories );
    $self->store( $tree, $model_c, $table, 'c' );

    $categories = { 9593 => '1', $tx_id => '2' };
    my $model_d = $self->calculate_branch_likelihood( $tree, $aln, $categories );
    $self->store( $tree, $model_d, $table, 'd' );

    $categories = { 9593 => '1', $tx_id => '1' };
    my $model_e = $self->calculate_branch_likelihood( $tree, $aln, $categories );
    $self->store( $tree, $model_e, $table, 'e' );

    # Now store the FULL primate tree & align for subst. collection.
    $self->param( 'tree',     $full_primate_tree );
    $self->param( 'cdna_aln', $full_primate_aln );
    $self->param( 'aa_aln',   Bio::EnsEMBL::Compara::AlignUtils->translate($full_primate_aln) );

    $categories = { 9604 => '1' };
    my $model_f =
      $self->calculate_branch_likelihood( $full_primate_tree, $full_primate_aln, $categories );
    $self->store( $tree, $model_f, $table, 'f' );

    $categories = { 9604 => '$1' };
    my $model_g =
      $self->calculate_branch_likelihood( $full_primate_tree, $full_primate_aln, $categories );
    $self->store( $tree, $model_g, $table, 'g' );

    $categories = {};
    my $model_h =
      $self->calculate_branch_likelihood( $full_primate_tree, $full_primate_aln, $categories );
    $self->store( $tree, $model_h, $table, 'h' );

    $self->collect_values($tree);
    $self->store_params_in_table( $self->dbc, $table, $self->params );
  }

}

sub collect_stats {
  my $self     = shift;
  my $tree     = shift;
  my $taxon_id = shift;

  my ($ref_member) = grep { $_->taxon_id == $taxon_id } $tree->leaves;

  $self->param( 'name', $ref_member->get_Gene->external_name );

  my $gc = $self->gc_content($ref_member);
  $self->param( 'gc_cds', sprintf( "%.3f", $gc ) );
  my $gc3 = $self->gc3_content($ref_member);
  $self->param( 'gc_3', sprintf( "%.3f", $gc3 ) );
  my $genomic = $self->genomic_gc_content($ref_member);
  $self->param( 'gc_genomic', sprintf( "%.3f", $genomic ) );
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
    $self->param( $prefix . '_tx',      $member->get_Transcript->stable_id );
  } else {
    $self->param( $prefix . '_protein', undef );
    $self->param( $prefix . '_gene',    undef );
    $self->param( $prefix . '_tx',      undef );
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

  $self->params( 'tree_to_use', $tree );
  my $file;

  $self->param( 'reference_species', 9606 );
  my $ref_species = $self->param('reference_species');
  my @members     = $tree->leaves;
  my ($ref_member) = grep { $_->taxon_id == $self->param('reference_species') } @members;
  $ref_member = $members[0] if ( !defined $ref_member );

  my $gene_name = $ref_member->get_Gene->external_name || $ref_member->gene_member->stable_id;
  $self->param( 'gene_name', $gene_name );

  my $c_dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
    -url => 'mysql://ensadmin:ensembl@ensdb-archive:5304/ensembl_compara_58' );
  my $params = $self->params;
  my $tree_aln_obj =
    Bio::EnsEMBL::Compara::ComparaUtils->get_compara_or_genomic_aln( $c_dba, $tree, $ref_member,
    $self->params );

  if ( !defined $tree_aln_obj ) {
    $self->fail_and_die( sprintf( "Tree or alignment not good! [%s]", $self->param('aln_type') ) );
  }

  my $aln = $tree_aln_obj->{aln};
  $tree = $tree_aln_obj->{tree};
  my $extra  = $tree_aln_obj->{extra};
  my $aln_id = $gene_name;

  $self->set_params($extra);

  if ( $self->param('likelihood_filter_substitution_runs') ) {
    
    my $masked_count =
      Bio::EnsEMBL::Compara::AlignUtils->mask_high_mutation_windows( $aln, 30, 10, 5 );
    $self->param( 'filtered_mut_windows', $masked_count );
  }
  
  $file = $self->save_aln(
    $aln, {
      id        => $aln_id,
      filename  => $aln_id,
      subfolder => $self->param('alignment_output_folder')
    }
    );
  $self->param( 'aln_length', $aln->length );
  $self->param( 'aln_file',   $file->{rel_file} );

  return { tree => $tree, aln => $aln };
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

  # Turn the labeled Ensembl tree object into a BioPerl TreeI-compliant object.
  my $treeI = $TREE->to_treeI($labeled_tree);

  # Detect the case when there's only one rate category, and switch to 'model=0' for PAML.
  my $newick = $labeled_tree->newick_format;
  if ( $newick !~ m/(#|\$)/ ) {

    #warn("Only one omega category!");
    $params->{model} = 0;
  }

  $params->{aaDist} = 0;

  print "Tree to PAML: " . $newick . "\n";
  print "Aln to PAML:\n";
  $self->pretty_print( $aln, { full => 0, width => 150 } );

  #my $aa = Bio::EnsEMBL::Compara::AlignUtils->translate($aln);
  #$self->pretty_print($aa,{full => 1});

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
  my $copy = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree);

  foreach my $node ( $copy->nodes ) {

   # If this leaf's taxon_id exists in the mapping, append that value (along with '#') to the label.
   # Otherwise, leave it alone.
    my $category = $taxon_categories->{ $node->taxon_id } || '';
    if ( $category eq '' && !$node->is_leaf ) {
      $node->name('');
    }
    next if ( $category eq '' );
    $category = '#' . $category unless ( $category =~ m/(#|\$)/ );
    if ( $node->is_leaf ) {
      $node->name( $node->name . $category );
    } else {
      $node->name($category);
    }
  }

  #print "Categorized tree: ".$copy->newick_format."\n";
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
  my $tree   = shift;
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

  if ( $prefix =~ m/[ha]/i ) {

    # Get the codeml results tree.
    my $codeml_tree = $model->{tree};

    # Choose which table to store things in.
    my $subs_table = $table . '_subs';


    # Store subs from the full primate tree into an 'allsubs' table.
    if ( $prefix eq 'h' ) {
      $subs_table = $self->param('output_table') . '_allsubs';
    }
    $self->create_table_from_params( $self->hive_dba, $subs_table, $self->codon_subs_table_def );

    foreach my $node ( $codeml_tree->get_nodes ) {
      my $id = $node->id;

      # Grab the substitutions from the codeml tree node.
      my $subs_hash = $node->get_tag_values('substitutions');
      my $i         = 0;
      foreach my $sub_key ( sort { $a <=> $b } keys %$subs_hash ) {
        my $subst = $subs_hash->{$sub_key};

        #$self->hash_print($subst);
        my $obj = $self->extract_substitution_info(
          $node,
          $subst,
          $self->param('tree'),    # Previously stored these objects in the params for safe keeping.
          $self->param('aa_aln'),
          $self->param('cdna_aln'),
        );
        if ( defined $obj ) {
          $i++;
          my $params = $self->replace( $self->params, $obj );

          #$self->hash_print($params);
          $self->store_params_in_table( $self->dbc, $subs_table, $params );
        } else {

          #print "No results for subst:\n";
          #$self->hash_print($subst);
        }
      }
      print "Stored $i substitutions for $id\n" if ( $self->debug );
    }
  }

  # Get the grantham scores and N and S substitution counts.
  if ( $prefix eq 'a' ) {
    my $tree       = $model->{tree};
    my $gor_name   = $self->param('gorilla_member')->name;
    my $other_name = $self->param('non_gorilla_member')->name;

    my ($gor_node)   = grep { $_->id eq $gor_name } $tree->get_nodes;
    my ($other_node) = grep { $_->id eq $other_name } $tree->get_nodes;

    $self->store_subs_parameters( 'Ggor', $gor_node, $gor_node->get_tag_values('substitutions') );
    $self->store_subs_parameters( 'other', $other_node,
      $other_node->get_tag_values('substitutions') );
  }
}

sub store_subs_parameters {
  my $self      = shift;
  my $suffix    = shift;
  my $node      = shift;
  my $subs_hash = shift;

  my $scores = Bio::Greg::Codeml->get_grantham_score_hash;

  my $dn       = 0;
  my $ds       = 0;
  my $grantham = 0;
  foreach my $key ( sort { $a <=> $b } keys %$subs_hash ) {
    my $subst = $subs_hash->{$key};

    # This will return undefined if the substitution involves Ns or stop codons.
    my $obj = $self->extract_substitution_info(
      $node,
      $subst,
      $self->param('tree'),    # Previously stored these objects in the params for safe keeping.
      $self->param('aa_aln'),
      $self->param('cdna_aln'),
    );
    next unless ( defined $obj );

    $self->hash_print($subst);
    my $sub_key = $subst->{'aa_a'} . $subst->{'aa_b'};
    $grantham += $scores->{$sub_key};
    if ( $subst->{'aa_a'} ne $subst->{'aa_b'} ) {
      $dn += 1;
    } else {
      $ds += 1;
    }
  }

  $self->param( 'grantham_' . $suffix, $grantham );
  $self->param( 'subs_n_' . $suffix,   $dn );
  $self->param( 'subs_s_' . $suffix,   $ds );
}

sub get_gene_stats_def {
  my $gene_stats_def = {
    data_id      => 'int',
    node_id      => 'int',
    orig_node_id => 'string',

    name         => 'char32',
    Hsap_protein => 'char16',
    Hsap_gene    => 'char16',
    Hsap_tx      => 'char16',
    Ggor_protein => 'char16',
    Ggor_gene    => 'char16',
    Ggor_tx      => 'char16',

    job_id => 'int',

    aln_length => 'int',
    aln_file       => 'string',

    filtered_mut_windows => 'int',
    filtered_Hsap        => 'int',
    filtered_Ggor        => 'int',
    filtered_Ptro        => 'int',
    filtered_Ppyg        => 'int',
    filtered_Mmul        => 'int',
    filtered_Cjac        => 'int',
    filtered_Tsyr        => 'int',
    filtered_Ogar        => 'int',
    filtered_Mmur        => 'int',

    gc_cds     => 'float',
    'gc_3'     => 'float',
    gc_genomic => 'float',

    grantham_Ggor  => 'int',
    grantham_other => 'int',
    subs_n_Ggor    => 'int',
    subs_n_other   => 'int',
    subs_s_Ggor    => 'int',
    subs_s_other   => 'int',

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
