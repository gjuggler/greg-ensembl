package Bio::Greg::Mammals::SitewiseMammals;

use strict;
use Bio::Greg::Codeml;
use File::Path;

use base (
  'Bio::Greg::Hive::Process',
  'Bio::Greg::StatsCollectionUtils',
  'Bio::Greg::Hive::Align',
  'Bio::Greg::Hive::CountSubstitutions',
  'Bio::Greg::Hive::SitewiseMapper',
  'Bio::Greg::Hive::PhyloAnalysis',
  'Bio::Greg::Hive::AlignmentFilteringSteps'
);

sub param_defaults {
  return {
  };
}

sub fetch_input {
  my $self = shift;

  # Fetch parameters from all possible locations.
  $self->load_all_params();

  # Create tables if necessary.
  $self->create_table_from_params( $self->compara_dba, 'sites',
                                   $self->_sites_table_structure );

  $self->create_table_from_params( $self->compara_dba, 'sites_pos',
                                   $self->_sites_table_structure );

  $self->create_table_from_params( $self->compara_dba, 'genes',
                                   $self->_genes_table_structure );

  $self->create_table_from_params( $self->compara_dba, 'trees',
                                   $self->_trees_table_structure );

}

sub run {
  my $self = shift;

  $self->param('store_stuff', 0);

  my $tree = $self->get_tree;
  my $member = $self->_get_ref_member($tree);
  my $gene = $member->get_Gene;
  $self->param('gene_name',$gene->external_name);
  $self->param('data_id', $self->param('node_id'));
  print "Gene : ".$member->display_label."\n";

  $self->param('member', $member);
  my $gene_name = $member->get_Gene->external_name || $member->get_Gene->stable_id;
  $self->param('gene_name',$gene_name);

  # Handle the mitochondrial genetic code.
  $self->param('chr_name', ''.$member->get_Transcript->slice->seq_region_name);
  if ($self->param('chr_name') =~ m/MT/gi) {
    $self->param('bioperl_codontable_id', 2);
    $self->param('icode', 1);
  } else {
    $self->param('bioperl_codontable_id', 1);
    $self->param('icode', 0);
  }
  
  $self->put('orig_tree', Bio::EnsEMBL::Compara::TreeUtils->to_newick($tree));
  if (scalar($tree->leaves) < 10) {
    my $leaf_count = scalar($tree->leaves);
    $self->fail_and_die("small_tree", "Tree w/ $leaf_count leaves!");
  }
#  if (scalar($tree->leaves) > 80) {
#    print $tree->ascii;
#    my $leaf_count = scalar($tree->leaves);
#    $self->fail_and_die("big_tree", "Tree w/ $leaf_count leaves!");
#  }

  my $seq_length = $member->seq_length;
  if ($seq_length > 1200) {
    $self->param('aln_type', 'pagan');
  } else {
    $self->param('aln_type', 'prank_codon');
  }

  $self->dbc->disconnect_when_inactive(1);
  $self->compara_dba->disconnect_when_inactive(1);

  # Get the seq quality filtered alignments.
  my ($aln) = $self->_get_aln($tree, $member, 1);
  $self->save_hash;

  my $tree_orig = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree);
  my $orig_newick = Bio::EnsEMBL::Compara::TreeUtils->to_newick($tree_orig);

  ($tree, $aln) = $self->_get_clade($tree_orig, $aln, 'mammals');

  if (!defined $tree || scalar($tree->leaves) < 2) {
    my $leaf_count = scalar($tree->leaves);
    $self->fail_and_die("small tree", "Mammal sub-tree w/ $leaf_count leaves!");    
  }

  $self->param('dup_species_list',$self->species_with_dups($tree));
  $self->param('dup_species_count',$self->duplication_count($tree));

  my $tree_w_paralogs = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree);
  my $w_paralogs_newick = Bio::EnsEMBL::Compara::TreeUtils->to_newick($tree_w_paralogs);
  ($tree, $aln) = $self->_remove_paralogs($tree, $aln, 1);
  # Get the ref member again now that paralogs are removed.
  $member = $self->_get_ref_member($tree);
  $self->param('member', $member);
  $self->save_hash;

  $self->param('store_stuff', 1);
  my $newick = Bio::EnsEMBL::Compara::TreeUtils->to_newick($tree);
  $self->_save_trees($orig_newick, $w_paralogs_newick, $newick);
  $self->param('store_stuff', 0);

  return;

  $self->dbc->disconnect_when_inactive(1);

  # Realign with prank.
  $self->pretty_print($self->_tx_aln($aln));
  $aln = $self->_align($tree, $aln);
  $self->pretty_print($self->_tx_aln($aln));
  $self->save_hash;

  # Store realigned seqs in a table.
  #$self->_store_seqs($tree, $aln);

  # Store / mask out substitution clusters
  $self->dbc->disconnect_when_inactive(1);
  #$aln = $self->_mask_subs($tree, $aln);
  $self->save_hash;

  # Output alignment to a file.
  $self->output_alignment($tree, $aln);
  
  # Save m0-inferred substitutions.
  #my $m0_lines = $self->_save_subs($tree, $aln);

  # Run SLR.
  $self->param('store_stuff', 1);
  $self->_run_slr($tree, $aln);
  $self->save_hash;
  $self->param('store_stuff', 0);

  $self->_collect_gene_data($tree, $aln);
  $self->save_hash;
}

sub _get_aln {
  my $self = shift;
  my $tree = shift;
  my $member = shift;
  my $get_filtered = shift;


  if ($get_filtered) {
    my $existing_aln = $self->get('orig_aln');
    return $existing_aln if (defined $existing_aln && $self->param('force_recalc') != 1);
  }

  $self->param('aln_type', 'compara');
  if ($get_filtered) {
    $self->param('quality_threshold', 25);
    $self->param('store_filtered_codons', 1);
  } else {
    $self->param('quality_threshold', 0);
    $self->param('store_filtered_codons', 0);
  }
  $self->param('hive_dbc', $self->dbc);
  $self->param('process', $self);
  my $tree_aln_obj =
    Bio::EnsEMBL::Compara::ComparaUtils->get_compara_or_genomic_aln( $self->compara_dba, $tree, $member,
                                                                     $self->params );

  my $aln = $tree_aln_obj->{aln};
  if ($get_filtered) {
    $self->put('orig_aln', $aln);
  }

  return $aln;
}

sub _remove_paralogs {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $store_in_table = shift;

  return $self->remove_paralogs($tree, $aln, $store_in_table);
}

sub _save_trees {
  my $self = shift;
  my $tree_orig = shift;
  my $tree_w_paralogs = shift;
  my $tree = shift;

  my $params = $self->replace($self->params, {
    tree_orig => $tree_orig,
    tree_w_paralogs => $tree_w_paralogs,
    tree => $tree
  });

  $self->store_params_in_table($self->dbc, 'trees', $params);
}

sub _align {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  my $existing_aln = $self->get('prank_aln');
  return $existing_aln if (defined $existing_aln && $self->param('force_recalc') != 1);

  if ($self->param('aln_type') eq 'pagan') {
    $self->param('aligner', 'pagan');
  } else {
    $self->param('aligner', 'prank_codon');
  }

  my $pep_aln = $self->_tx_aln($aln);
  $aln = $self->align($tree, $aln, $pep_aln);

  $self->put('prank_aln', $aln);
  return $aln;
}

sub _mask_subs {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  print "  masking subs...\n";
  my $lines_s = 'mask_subs_lines';

  my $mask_subs_lines = $self->get($lines_s);
  if (defined $mask_subs_lines && !$self->param('force_recalc')) {
    $self->param($lines_s, $mask_subs_lines);
  }
  my $do_filtering = 0;
  my $masked_aln = $self->filter_substitution_runs($tree, $aln, $do_filtering);

  $self->put($lines_s, $self->param($lines_s));
#  $self->put('masked_aln', $masked_aln);
  return $masked_aln;
}

sub output_alignment {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  my $id = $self->param('gene_name');
  my $fldr = $self->get_output_folder . '/alns';
  my $file = $self->save_file({
    folder => $fldr,
    filename => "$id",
    id => "$id",
    extension => "fasta"
                              });
  my $filename = $file->{full_file};
  Bio::EnsEMBL::Compara::AlignUtils->to_file($aln, $filename);

  $file = $self->save_file({
    folder => $fldr,
    filename => "$id",
    id => "$id",
    extension => "nh"
                              });
  $filename = $file->{full_file};
  Bio::EnsEMBL::Compara::TreeUtils->to_file($tree, $filename);
}

sub _get_clade {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $clade = shift;

  my @primates = (9478, 9483, 9544, 9593, 9598, 9601, 9606, 30608,
  30611, 61853);
  my @glires = (9978, 9986, 10020, 10090, 10116, 10141, 43179);

  my @laur = (9365, 9615, 9646, 9685, 9739, 9796, 9823, 9913, 30538,
  42254, 59463, 132908);
  my @atlanto = (9358, 9361, 9371, 9785, 9813);

  #my @sauria = (28377, 59729, 9031, 9103);
  my @euth = (9358, 9361, 9365, 9371, 9478, 9483, 9544, 9593, 9598,
  9601, 9606, 9615, 9646, 9685, 9739, 9785, 9796, 9813, 9823, 9913,
  9978, 9986, 10020, 10090, 10116, 10141, 30538, 30608, 30611, 37347,
  42254, 43179, 59463, 61853, 132908);
  my @mamm = (9258, 9315, 9358, 9361, 9365, 9371, 9478, 9483, 9544,
  9593, 9598, 9601, 9606, 9615, 9646, 9685, 9739, 9785, 9796, 9813,
  9823, 9913, 9978, 9986, 10020, 10090, 10116, 10141, 13616, 30538,
  30608, 30611, 37347, 42254, 43179, 59463, 61853, 132908);

  my @sparse_glires = (10090, 10116, 10020, 43179, 10141);

  my @sparse_mamm = (9258, 9315, 9361, 9785, 9606, 10090, 9615);

  my @hmrd = (9606, 10090, 10116, 9615);

  my @hiq_only = (9606, 9598, 9544, 10090, 10116, 9615, 9913, 9823, 9796);
  
  my $map = {
    primates => \@primates,
    glires => \@glires,
    laurasiatheria => \@laur,
    atlantogenata => \@atlanto,
    eutheria => \@euth,
    mammals => \@mamm,
    sparse_glires => \@sparse_glires,
    sparse_mammals => \@sparse_mamm,
    hmrd => \@hmrd,
    hiq_only => \@hiq_only
  };
  my $arrayref = $map->{$clade};
  my @taxids = @$arrayref;

  my @keepers;
  foreach my $leaf ($tree->leaves) {
    my $found = grep {$leaf->taxon_id == $_} @taxids;
    if ($found) {
      push @keepers, $leaf;
      print "  keeping leaf ".$leaf->name."\n";
    }
  }

  if (scalar(@keepers) == 0) {
    return (undef, undef);
  }

  my $subtree = Bio::EnsEMBL::Compara::TreeUtils->extract_subtree_from_leaf_objects($tree, \@keepers);
  if (!$subtree) {
    return (undef, undef);
  }
  my $subaln = Bio::EnsEMBL::Compara::ComparaUtils->restrict_aln_to_tree($aln, $subtree);

  return ($subtree, $subaln);
}

sub _save_subs {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  my $lines = $self->get('mask_subs_lines');
  if ($self->param('store_stuff') != 0) {
    $self->store_subs($tree, $aln, $lines, $self->param('member'), 'subs');
  }
  return $lines;  
}

sub _store_seqs {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  return if ($self->param('store_stuff') == 0);
  
  $self->store_seqs($tree, $aln);
}

sub _run_slr {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pos_only = shift;

  $pos_only = 0 unless (defined $pos_only);

  #$self->pretty_print($self->_tx_aln($aln), {full => 1});

  $self->param('analysis_action', 'slr');

  # Collect Pfam, exon, and filter data.
  my $pep_aln = $self->_tx_aln($aln);
  my $pep_data = $self->get('pep_data');
#  if (!defined $pep_data || $self->param('force_recalc')) {
    $pep_data = {};
    my $pfam_sites = $self->collect_pfam($tree,$pep_aln);
    my $exon_sites = $self->collect_exons($tree,$pep_aln);
    $pep_data->{pfam} = $pfam_sites;
    $pep_data->{exons} = $exon_sites;
    $self->put('pep_data', $pep_data);
#  }
  
  my $clade_map = {
    primates => 1,
    glires => 2,
    laurasiatheria => 3,
    atlantogenata => 4,
    eutheria => 5,
    mammals => 6,
    sparse_glires => 7,
    sparse_mammals => 8,
    hiq_only => 9,
    hmrd => 10
  };

  my @clades = keys %$clade_map;
  my @skipped_clades;
  foreach my $clade (@clades) {
    print "Getting clade $clade\n";
    my $tree_cpy = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree);
    my ($subtree, $subaln) = $self->_get_clade($tree_cpy, $aln, $clade);
    if (!defined $subtree) {
      push @skipped_clades, $clade;
      next;
    }

    $self->pretty_print($self->_tx_aln($subaln), {full => 1});

    #print $subtree->ascii."\n";
    my $lines_key = $clade.'_slr';
    if ($pos_only) {
      $lines_key = $clade . '_slr_pos';
    }

    my $lines = $self->get($lines_key);

    my $pep_aln = $self->_tx_aln($subaln);

    if (!defined $lines || $self->param('force_recalc')) {
      $self->pretty_print($pep_aln);
      if ($pos_only) {
        $self->param('slr_positive', 1);
      } else {
        $self->param('slr_positive', 0);
      }
      $lines = $self->run_sitewise_analysis($subtree, $subaln, $pep_aln);
      if (defined $lines) {
        $self->put($lines_key, $lines);
      }
    }

    if (defined $lines) {
      my $results = $self->parse_sitewise_output_lines($subtree, $subaln, $pep_aln, $lines);
      #$self->hash_print($results);

      #if ($clade eq 'primates') {
      #  print join("\n", @$lines)."\n";
      #}

      my $param_set_id = $clade_map->{$clade};
      if ($self->param('store_stuff') != 0) {
        $self->_store_sitewise($subtree, $subaln, $pep_aln, $results, $clade, $param_set_id, $pep_data, $pos_only);
      }
      $self->save_hash;
    } else {
      push @skipped_clades, $clade;
    }
  }

  $self->param('skipped_clades', join(', ',@skipped_clades));
  $self->param('skipped_clade_count', scalar(@skipped_clades));
}

sub clade_string_map {
  my $self = shift;
  my $clade_string_map = {
    primates => 'p',
    glires => 'g',
    laurasiatheria => 'l',
    atlantogenata => 'a',
    eutheria => 'e',
    mammals => 'm',
    sparse_glires => 'sg',
    sparse_mammals => 'sm',
    hiq_only => 'h',
    hmrd => 'f'
  };
  return $clade_string_map;
}

sub _store_sitewise {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pep_aln = shift;
  my $slr_hash = shift;
  my $clade = shift;
  my $param_set_id = shift;
  my $pep_data = shift;
  my $pos_only = shift;
  
  my $ref_member = $self->param('member');
  my ($ref_seq) = grep {$_->id eq $ref_member->stable_id} $aln->each_seq;
  my ($pep_ref_seq) = grep {$_->id eq $ref_member->stable_id} $pep_aln->each_seq;
  my $slr_tree = $slr_hash->{slr_tree};
  my $slr_treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($slr_tree);

  my $clade_string_map = $self->clade_string_map;
  my $str = $clade_string_map->{$clade};
  $self->param($str.'_slr_mean_path', $slr_treeI->root->mean_path_length);
  $self->param($str.'_slr_total_length', $slr_treeI->root->total_branch_length);
  $self->param($str.'_slr_dnds', $slr_hash->{omega});
  $self->param($str.'_slr_kappa', $slr_hash->{kappa});
  $self->param($str.'_leaf_count', scalar($slr_treeI->leaves));

  # Add ungapped branch lengths to the slr hash.
  $self->add_ungapped_branch_lengths($slr_tree, $pep_aln, $slr_hash);

  $self->dbc->db_handle->{AutoCommit} = 0;

  my $n_sites = 0;
  foreach my $key (1 .. $pep_aln->length) {
    my $site_obj = $slr_hash->{$key};
    next unless (defined $site_obj);
    my $aln_position = $site_obj->{aln_position};
    my $note = $site_obj->{note};
    next if ($note eq 'single_char' || $note eq 'all_gaps');
    $n_sites++;
    my $cur_params = $self->params;

    $cur_params = $self->replace($cur_params,$site_obj);
    $cur_params->{parameter_set_id} = $param_set_id;

    if (defined $pep_ref_seq) {
      my $ref_position = $pep_ref_seq->location_from_column($aln_position);
      if (defined $ref_position && $ref_position->location_type() eq 'EXACT') {
        if ($ref_member->taxon_id == 9606) {
          my $ref_coords = $self->get_coords_from_pep_position($ref_member,$ref_position->start);
          $cur_params = $self->replace($cur_params,$ref_coords);
        }
        $cur_params->{seq_position} = $ref_position->start;
      }
    }

    if (defined $pep_data->{pfam}->{$aln_position}) {
      $cur_params = $self->replace($cur_params, $pep_data->{pfam}->{$aln_position});
    }    
    if (defined $pep_data->{exons}->{$aln_position}) {
      $cur_params = $self->replace($cur_params, $pep_data->{exons}->{$aln_position});
    }

    my $table_name = 'sites';
    $table_name = 'sites_pos' if ($pos_only);
    
    $self->store_params_in_table($self->dbc, $table_name, $cur_params);
  }

  $self->dbc->db_handle->{AutoCommit} = 1;
  $self->param($str.'_slr_sites', $n_sites);

}

sub _collect_gene_data {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $m0_lines = shift;

#  my $m0_tree = Bio::Greg::Codeml->parse_codeml_results($m0_lines);
#  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);

  my $ref_member = $self->param('member');

  # Basic details.
  my $tx = $ref_member->get_Transcript;
  eval {
    $tx = $tx->transform('chromosome');
  };
  if (defined $tx) {
    $self->param('chr_name','chr'.$tx->slice->seq_region_name);
    $self->param('chr_start',$tx->coding_region_start);
    $self->param('chr_end',$tx->coding_region_end);
    $self->param('chr_strand',$tx->strand);
  }
  $self->param('ref_taxon_id',$ref_member->taxon_id);
  $self->param('gene_name',$ref_member->get_Gene->external_name);
  $self->param('ref_gene_description',$ref_member->get_Gene->description);
  $self->param('ref_gene_id',$ref_member->get_Gene->stable_id);
  $self->param('ref_transcript_id',$ref_member->get_Transcript->stable_id);
  $self->param('ref_protein_id',$ref_member->get_Transcript->translation->stable_id);
  $self->param('ref_seq_length',$ref_member->seq_length);
  $self->param('aln_length',$aln->length / 3);
  
  # Calculate GC content.
  $self->param( 'gc_cds', sprintf( "%.3f", $self->gc_content($ref_member) ) );
  $self->param( 'gc_3', sprintf( "%.3f", $self->gc3_content($ref_member) ) );
  $self->param( 'gc_genomic', sprintf( "%.3f", $self->genomic_gc_content($ref_member) ) );

  # Tree properties.
  $self->param('tree_mean_path',$self->mean_path($tree));
  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
  $self->param('tree_total_length',$treeI->root->total_branch_length);

#  $self->param('paml_mean_path', $m0_tree->root->mean_path_length);
#  $self->param('paml_total_length', $m0_tree->root->total_branch_length);  
#  my @omegas = Bio::Greg::Codeml->extract_omegas($m0_lines);
#  $self->param('paml_dnds', $omegas[0]);

  $self->param('leaf_count',scalar($tree->leaves));

  $self->store_params_in_table($self->dbc, 'genes', $self->params);
}


sub _get_ref_member {
  my $self = shift;
  my $tree = shift;

  my $ref_member;
  my @members = $tree->leaves;

  # Try for human first.
  ($ref_member) = grep { $_->taxon_id == 9606 } @members;
  
  # Then mouse.
  ($ref_member) = grep { $_->taxon_id == 10090 } @members if (!defined $ref_member);

  # Then dog.
  ($ref_member) = grep { $_->taxon_id == 9615 } @members if (!defined $ref_member);

  # Then, take whatever we can get.
  ($ref_member) = $members[0] if (!defined $ref_member);
  
  $self->param('ref_member_id',$ref_member->stable_id);
  return $ref_member;
}

sub get {
  my $self = shift;
  my $key = shift;

  $self->load_hash;

  my $hash = $self->param('hash');
  warn("Tried to get object with key $key but doesn't exist!") unless (defined $hash->{$key});
  return $hash->{$key};
}

sub put {
  my $self = shift;
  my $key = shift;
  my $obj = shift;

  $self->load_hash;
  my $hash = $self->param('hash');
  print "Putting [$key] into hash...\n";
  $hash->{$key} = $obj;
  $self->param('hash', $hash);
}

sub save_hash {
  my $self = shift;

  if ($self->param('store_stuff') == 0) {
    return;
  }

  $self->load_hash;
  my $hash = $self->param('hash');
  print "Saving hash...\n";
  $self->hash_print($hash);
  my $full_file = $self->_hash_file;
  $self->frz($full_file, $hash);
}

sub load_hash {
  my $self = shift;

  if ($self->param('hash_loaded') == 1) {
    return $self->param('hash');
  }

  my $full_file = $self->_hash_file;
  my $hash;
  if (-e $full_file) {
    $hash = $self->thw($full_file);
  } else {
    warn("Hash file not found -- creating blank hash!");
    $hash = {};
  }
  $self->param('hash', $hash);
  print "Loaded hash:\n";
  $self->hash_print($hash);
  $self->param('hash_loaded', 1);
}

sub _hash_file {
  my $self = shift;
  my $id = $self->param('gene_name');
  my $file = $self->save_file({
    filename => "$id",
    id => "$id",
    extension => "perlobj"
                              });
  my $full_file = $file->{full_file};
  return $full_file;
}

sub _sites_table_structure {
  my $self = shift;

  my $structure = {
    # IDs.
    data_id => 'int',
    parameter_set_id => 'tinyint',

    # Genomic locations on the reference sequence
    chr_name => 'char8',
    chr_start => 'int',
    chr_end => 'int',

    # Pfam annotations.
    pfam_domain => 'char16',
    pfam_position => 'int',
    pfam_score => 'int',

    # Exon annotation.
    exon_position => 'char8',
    splice_distance => 'int',

    # Basic SLR-derived attributes.
    aln_position => 'int',
    seq_position => 'int',
    ncod => 'int',
    nongap_bl => 'float',           # This is calculated using StatsCollectionUtils' method add_ungapped_branch_lengths
    omega => 'float',
    omega_lower => 'float',
    omega_upper => 'float',
    lrt_stat => 'float',
    type => 'char16',
    note => 'char16',
    random => 'char16',

    unique_keys => 'data_id,parameter_set_id,aln_position'
  };

  return $structure;
}

sub _trees_table_structure {
  my $self = shift;

  my $structure = {
    # IDs.
    data_id => 'int',
    job_id => 'int',
    tree_orig => 'string',
    tree_w_paralogs => 'string',
    tree => 'string',
    unique_keys => 'data_id'
  };
  return $structure;
}

sub _genes_table_structure {
  my $self = shift;

  my $structure = {
    # IDs.
    data_id => 'int',
    job_id => 'int',

    # Human gene attributes.
    chr_name => 'string',
    chr_start =>'int',
    chr_end => 'int',
    chr_strand => 'tinyint',
    gene_name => 'char32',
    ref_gene_description => 'string',
    ref_taxon_id => 'int',
    ref_transcript_id => 'char32',
    ref_gene_id => 'char32',
    ref_protein_id => 'char32',

    # Alignment properties.
    aln_length => 'int',
    ref_seq_length => 'int',
    
    masked_nucs => 'int',
    masked_ids => 'string',

    # General gene / tree attributes.
    gc_cds => 'float',
    'gc_3' => 'float',
    gc_genomic => 'float',
    tree_mean_path => 'float',
    tree_total_length => 'float',
    paml_mean_path => 'float',
    paml_total_length => 'float',
    paml_dnds => 'float',
    leaf_count => 'int',

    dup_species_count => 'int',
    dup_species_list => 'string',

    skipped_clades => 'string',
    skipped_clade_count => 'int',
    data_prefix => 'char8',

    unique_keys => 'data_id'
  };

  my $map = $self->clade_string_map;
  foreach my $key (keys %$map) {
    my $clade = $map->{$key};
    $structure->{$clade.'_leaf_count'} = 'int';
    $structure->{$clade.'_slr_sites'} = 'int';
    $structure->{$clade.'_slr_dnds'} = 'float';
    $structure->{$clade.'_slr_kappa'} = 'float';
    $structure->{$clade.'_slr_mean_path'} = 'float';
    $structure->{$clade.'_slr_total_length'} = 'float';    
  }

  return $structure;
}

1;
