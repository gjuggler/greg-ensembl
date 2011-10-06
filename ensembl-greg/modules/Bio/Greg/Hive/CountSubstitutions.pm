package Bio::Greg::Hive::CountSubstitutions;

use strict;
use Bio::Greg::StatsCollectionUtils;

use base ( 'Bio::Greg::Hive::Process', 'Bio::Greg::StatsCollectionUtils' );

my $TREE    = "Bio::EnsEMBL::Compara::TreeUtils";
my $ALN     = "Bio::EnsEMBL::Compara::AlignUtils";
my $COMPARA = "Bio::EnsEMBL::Compara::ComparaUtils";

sub codon_subs_table_def {
  my $self = shift;
  return {
    data_id => 'int',
    gene_name      => 'char32',

    codeml_node_number => 'int',
    codeml_node_id => 'char32',   # Node ID retrieved from Codeml.

    taxon_id       => 'int',
    chr_name            => 'char8',    # Genomic coordinates.
    chr_start      => 'int',
    chr_strand     => 'int',
    seq_pos        => 'int',      # Position on the sequence in AA coordinates.
    aln_pos        => 'int',
    cdna_pos       => 'int',

    codon_from         => 'char4',
    codon_to           => 'char4',
    aa_from            => 'char1',
    aa_to              => 'char1',
    grantham_score     => 'float',
    confidence         => 'float',    # Codeml's confidence score.

    codon_cpg => 'tinyint',           # 1 if 'codon_from' overlapped *any* CpG sites.

    # These fields describe the single nucleotide mutation.
    nuc_from  => 'char1',
    nuc_to    => 'char1',
    mut_count => 'tinyint',           # Number of mutations within the codon.
    mut_pos   => 'tinyint',           # position in the codon: 1, 2, or 3
    mut_cpg   => 'tinyint',           # 0 = no cpg overlap; 1 = C-position; 2 = G-position.
    mut_rev_cpg =>
      'tinyint',    # Same as above, but gives the mutation's cpg position on the opposite strand.
    mut_syn  => 'tinyint',    # synonymous
    mut_nsyn => 'tinyint',    # nonsynonymous
    mut_sw   => 'tinyint',    # strong-to-weak
    mut_ws   => 'tinyint',    # weak-to-strong
    mut_ts   => 'tinyint',    # transition
    mut_tv   => 'tinyint',    # transversion
    mut_ils  => 'tinyint',   # likely ILS site

    pattern => 'char8',
    nongap_count => 'int',

    unique_keys => 'data_id,taxon_id,aln_pos',
  };
}

sub run_m0 {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  my $m0_params = {
    model => 0,
    fix_blength => 0,
    method => 0,
    fix_omega => 0,
    omega => 0.3,
    getSE => 1,
    Small_Diff => 1e-7,
  };

  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);

  my $lines;
  my $res = Bio::Greg::Codeml->branch_model_likelihood( $treeI, $aln, $self->worker_temp_directory, $m0_params );
  $lines = $res->{lines};
  
  return $lines;
}

sub load_subs {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $lines = shift;
  my $ref_member = shift;
  my $table = shift;

  my $subs_obj = $self->store_subs($tree, $aln, $lines, $ref_member, undef);
  return $subs_obj;
}

sub store_subs {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $lines = shift;
  my $ref_member = shift;
  my $table = shift;
  

  if (defined $table) {
    $self->create_table_from_params( $self->compara_dba, $table,
                                     $self->codon_subs_table_def );
  }

  if (!defined $ref_member) {
    my ($ref_member) = grep {$_->can('taxon_id') && $_->taxon_id == 9606} $tree->leaves;
    if (defined $ref_member) {
      my $mba = $self->compara_dba->get_MemberAdaptor;
      $ref_member = $mba->fetch_by_source_stable_id(undef, $ref_member->name);
      $ref_member->name($ref_member->stable_id);
    }
  }

  my $pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($aln);
  my $codeml_tree = Bio::Greg::Codeml->parse_codeml_results($lines);

  my @unwrapped_subs;

  foreach my $node ( $codeml_tree->get_nodes ) {
    my $id = $node->id;
    
    # Grab the substitutions from the codeml tree node.
    my $subs_hash = $node->get_tag_values('substitutions');

    my $i         = 0;
    my $n_ws = 0;
    my $n_sw = 0;
    my $n_ns = 0;
    my $n_s = 0;

    foreach my $sub_key ( sort { $a <=> $b } keys %$subs_hash ) {
      my $subst = $subs_hash->{$sub_key};
      
      #$self->hash_print($subst);
      my $obj = $self->extract_substitution_info(
        $node,
        $subst,
        $tree,
        $pep_aln,
        $aln,
        $ref_member
        );
      if ( defined $obj ) {
        $i++;

        $n_ws++ if ($obj->{mut_ws} == 1);
        $n_sw++ if ($obj->{mut_sw} == 1);
        $n_ns++ if ($obj->{mut_nsyn} == 1);
        $n_s++ if ($obj->{mut_syn} == 1);

        my $params = $self->replace( $self->params, $obj );

        if (defined $table) {
          $self->store_params_in_table( $self->dbc, $table, $params );
        }
        
        push @unwrapped_subs, $obj;
      } else {
        #print "No results for subst:\n";
        #$self->hash_print($subst);
      }
    }

    my $verb = 'Stored';
    $verb = 'Loaded' if (!defined $table);
    printf("${verb} %3s substitutions for %15.15s [%3s ns, %3s s]\n",
             $i, $id, $n_ns, $n_s) if ( $self->debug );
  }

  $self->param('subs', \@unwrapped_subs);
  return \@unwrapped_subs;
}

sub mask_substitution_runs {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $cache_file = shift;

  print $tree->ascii(0,0,1)."\n";

  $aln = Bio::EnsEMBL::Compara::AlignUtils->copy_aln($aln);
  my $tree_copy = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree);
  my $aln_copy = Bio::EnsEMBL::Compara::AlignUtils->copy_aln($aln);

  my $lines;
  if (defined $self->param('mask_subs_lines')) {
    print "  loading mask subs results from array\n";
    $lines = $self->param('mask_subs_lines');
  } elsif (-e $cache_file) {
    print "  loading mask subs results from file\n";
    $lines = $self->thw($cache_file);
  } else {
    print "  running m0...\n";
    $lines = $self->run_m0($tree_copy, $aln_copy);
    if (defined $cache_file) {
      $self->frz($cache_file, $lines);
    } else {
      $self->param('mask_subs_lines', $lines);
    }
  }
  
  my $m0_tree = Bio::Greg::Codeml->parse_codeml_results($lines);
  
  # Store substitutions.
  my $subs_obj = $self->load_subs($tree_copy, $aln_copy, $lines);
  my @subs = @{$subs_obj};

  my $muts;
  foreach my $s (@subs) {
    my $key = $s->{leaves_beneath};
    $muts->{$key} = [] if (!$muts->{$key});
    push @{$muts->{$key}}, $s;
  }

  my $pep_aln = $self->_tx_aln($aln);
  my $to_mask;
  foreach my $seq ($pep_aln->each_seq) {
    if (defined $muts->{$seq->id}) {
      my $dna_seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id($aln,$seq->id);
      my $dna_str = $dna_seq->seq;
      my $str = $seq->seq;
      my @cur_muts = @{$muts->{$seq->id}};

      my $win_width = 15;
      for (my $i=1; $i < $pep_aln->length - $win_width; $i++) {
        my $w_lo = $i;
        my $w_hi = $i + $win_width;
        
        my @ns_subs = grep {
          ($_->{aln_pos} >= $w_lo &&
           $_->{aln_pos} < $w_hi &&
           $_->{mut_nsyn} == 1)
        } @cur_muts;

        my $win_seq = substr($str,$i-1,$win_width);
        my $dna_win_seq = substr($dna_str,($i-1)*3,$win_width*3);

        # Find node in paml tree.
        my $m0_node = $m0_tree->find($seq->id);
        die "No m0 node!" unless ($m0_node);
        my $bl = $m0_node->branch_length;

        # Clip to a minimum branch length.
        my $min_bl = 0.05;
        $bl = $min_bl if ($bl < $min_bl);
        
        # Calc the mean # of n-syn subs per codon per unit branch length
        my $n_ns = scalar(@ns_subs);

        # Count non-synonymous codons with two or three mutations twice.
        $n_ns += scalar(grep { $_->{mut_count} > 1 } @ns_subs);

        my $nongap_str = $win_seq;
        $nongap_str =~ s/[-n]//g;
        my $n_nongap_codons = int(length($nongap_str));
        $n_nongap_codons = 1 if ($n_nongap_codons < 1);

        my $ns_per_codon = $n_ns / $n_nongap_codons;
        my $ns_per_codon_bl = $ns_per_codon / $bl;

        my $max_ns_per_codon_bl = 10;
        
        # Decrease the filtering threshold if the current window overlaps a gap.
        if ($win_seq =~ m/-/g) {
          $max_ns_per_codon_bl = 5;
        }
        # Decrease the threshold even more if the current window overlaps already-masked sequence.
        if ($win_seq =~ m/x/gi) {
          $max_ns_per_codon_bl = 5;
        }
        
        if ($win_seq =~ m/x/gi && $win_seq =~ m/-/g) {
          $max_ns_per_codon_bl = 2.5;
        }

        #print $seq->id."  " . $win_seq . "  " . $ns_per_codon_bl."\n";
        
        if ($ns_per_codon_bl > $max_ns_per_codon_bl
          ) {
          print $w_lo."  ". $seq->id."  ".$bl."  ". $ns_per_codon_bl."  ".$ns_per_codon."  ".$n_ns."\n";
          print "  ".$win_seq."\n";
          print "  ".$dna_win_seq."\n";
          foreach my $j ($w_lo .. $w_hi) {
            $to_mask->{$seq->id." ".$j} = 1;
          }
        }
      }
    }
  }

  my $n_masked = 0;
  my $masked_ids;
  foreach my $key (sort keys %$to_mask) {
    my ($id, $aln_pos) = split(' ', $key);    
    $masked_ids->{$id} = 1;
    my $seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id($aln,$id);
    my $cdna_lo = ($aln_pos-1)*3 + 0;
    my $cdna_hi = ($aln_pos-1)*3 + 2;
    
    my $orig_str = $seq->seq;
    my ($str, $n) = Bio::EnsEMBL::Compara::AlignUtils->mask_string($orig_str,$cdna_lo,$cdna_hi);
    $seq->seq($str);
    $n_masked += $n;
  }
  
  $self->store_param('masked_nucs', $n_masked);
  $self->store_param('masked_ids', join(', ', sort keys(%{$masked_ids})));
  
  return $aln;
}

sub _tx_aln {
  my $self = shift;
  my $aln = shift;
  return Bio::EnsEMBL::Compara::AlignUtils->translate($aln);
}

sub _grantham_matrix {
  my $self = shift;

  if (!defined $self->{_grantham}) {
    $self->{_grantham} = Bio::Greg::Codeml->get_grantham_score_hash;
  }
  return $self->{_grantham};
}

sub extract_substitution_info {
  my $self     = shift;
  my $node = shift;
  my $sub_obj  = shift;
  my $tree     = shift;
  my $aln      = shift;
  my $cdna_aln = shift;
  my $ref_member = shift;

  # sub_obj has the following: pos, codon_a, codon_b, aa_a, aa_b, confidence
  my $id = $node->id;
  my $sub_obj_id           = $sub_obj->{id};
  my $aln_pos      = $sub_obj->{pos};
  my $cdna_aln_pos = ( $aln_pos - 1 ) * 3 + 1;
  my $codon_a      = $sub_obj->{codon_a};
  my $codon_b      = $sub_obj->{codon_b};
  my $aa_a         = $sub_obj->{aa_a};
  my $aa_b         = $sub_obj->{aa_b};
  my $confidence   = $sub_obj->{confidence};

  return undef if ( $aa_a eq '*' || $aa_b eq '*' );
  return undef if ( $codon_a =~ m/[nx]/i || $codon_b =~ m/[nx]/i );

  my $nuc_left        = '';
  my $nuc_right       = '';
  my $codon_context_a = '';
  my $codon_context_b = '';

  my $final_params = {};

  # Get the list of leaves beneath this node.
  my $leaves_beneath = '';
  if ($node->is_Leaf) {
    $leaves_beneath = $node->id;
  } else {
    my @leaves;
    foreach my $sub_node ($node->get_all_Descendents) {
      push @leaves,$sub_node if ($sub_node->is_Leaf);
    }
    @leaves = sort {$a->id cmp $b->id} @leaves;
    $leaves_beneath = join(',',map {$_->id} @leaves);
  }

  # Store the genomic position in the reference member.
  if (defined $ref_member) {
    my $ref_seq = $aln->get_seq_by_id($ref_member->name);
    if ($ref_seq) {
      my $location = $ref_seq->location_from_column($aln_pos);
      if ($location && $location->location_type ne 'IN-BETWEEN') {
        eval {
          my $seq_pos = $location->start;
          my $coords_obj = $self->get_coords_from_pep_position($ref_member, $seq_pos);
          $final_params = $self->replace($final_params, $coords_obj);
        };
      }
    }
  }
  
  # Get the taxon ID of the node (if relevant)
  if ($node->is_Leaf) {
    my ($member) = grep {$_->name eq $node->id} $tree->leaves;
    if ($member && $member->can('taxon_id')) {
      $final_params->{taxon_id} = $member->taxon_id;
    } else {
      print "No member found for leaf node!\n";
    }
  } else {
    # Extract the taxon ID from the leaves beneath.
    my $h = "ENSP0";
    my $c = "ENSPTRP0";
    my $g = "ENSGGOP0";
    my $o = "ENSPPYP0";
    my $gibbon = "ENSNLEP0";
    my $macaque = "ENSMMUP0";
    my $marmoset = "ENSCJAP0";
    my $tarsier = "ENSTYSP0";
    my $mouse_lemur = "ENSMICP0";
    my $bushbaby = "ENSOGAP0";
    my $regex_to_tx = {
      "$h" => 9606,
      "$c" => 9598,
      "$g" => 9593,
      "$o" =>  9601,
      "$gibbon" => 61853,
      "$macaque" => 9544,
      "$marmoset" => 9483,
      "$tarsier" => 9478,
      "$mouse_lemur" => 30608,
      "$bushbaby" => 30611
    };
    my $tx_ids;
    foreach my $ensp (keys %$regex_to_tx) {
      my $taxon_id = $regex_to_tx->{$ensp};
      if ($leaves_beneath =~ m/${ensp}\d+/) {
        $tx_ids->{$taxon_id} = 1;
      }
    }

    my $taxon_id = 0;
    # Human-chimp: pan-homo artificial clade
    $taxon_id = 1234 if ($tx_ids->{9606} && $tx_ids->{9598});
    # Human-gorilla: Homininae
    $taxon_id = 207598 if ($tx_ids->{9606} && $tx_ids->{9593});
    # Human-orang: Hominidae (great apes)
    $taxon_id = 9604 if ($tx_ids->{9606} && $tx_ids->{9601});
    # Human-gibbon: Hominoidea (apes)
    $taxon_id = 314295 if ($tx_ids->{9606} && $tx_ids->{61853});
    # Human-macaque: Haplorrhini
    $taxon_id = 376913 if ($tx_ids->{9606} && $tx_ids->{9544});
    # Human-marmoset: Catarrhini
    $taxon_id = 9526 if ($tx_ids->{9606} && $tx_ids->{9483});
    # Human-tarsier: artificial clade
    $taxon_id = 10001 if ($tx_ids->{9606} && $tx_ids->{9478});
    # Mouse lemur and galago internal artificial clade
    $taxon_id = 10002 if ($tx_ids->{30608} && $tx_ids->{30611});
    # Human-{mouse_lemur or galago}
    $taxon_id = 10003 if ($tx_ids->{9606} && ($tx_ids->{30608} || $tx_ids->{30611}) );

    $final_params->{taxon_id} = $taxon_id if (defined $taxon_id);
  }

  # Synon/nonsynon.
  my $mut_syn  = 0;
  my $mut_nsyn = 0;
  $mut_syn = 1 if ( $aa_a eq $aa_b );
  $mut_nsyn = 1 if ( $aa_a ne $aa_b );
  
  # Isolate the nucleotide difference.
  my $mut_count = 0;
  my $mut_pos   = -1;
  my $nuc_a;
  my $nuc_b;
  foreach my $i ( 1, 2, 3 ) {
    my $tmp_nuc_a = substr( $codon_a, $i - 1, 1 );
    my $tmp_nuc_b = substr( $codon_b, $i - 1, 1 );
    if ( $tmp_nuc_a ne $tmp_nuc_b ) {
      $mut_pos = $i;
      $mut_count++;
      $nuc_a = $tmp_nuc_a;
      $nuc_b = $tmp_nuc_b;
    }
  }
  my $both_nucs = $nuc_a . $nuc_b;
  if ( $mut_count > 1 ) {

    # If we have multiple substitutions, we'll only use the latter substitution
    # to collect the rest of the stats.
  }

  # Weak-strong, transition-transversion.
  my $mut_ws = 0;
  my $mut_sw = 0;
  my $mut_ts = 0;
  my $mut_tv = 0;
  $mut_ws = 1 if ( $nuc_a =~ m/[at]/i && $nuc_b =~ m/[gc]/i );
  $mut_sw = 1 if ( $nuc_a =~ m/[gc]/i && $nuc_b =~ m/[at]/i );
  $mut_ts = 1 if ( $both_nucs =~ m/(ag|ga|ct|tc)/i );
  $mut_tv = 1 if ( $both_nucs =~ m/(ac|ca|at|ta|gc|cg|gt|tg)/i );

  # CpG overlaps.
  my $mut_cpg = 0;
  $mut_cpg = 1
    if ( ( $mut_pos == 1 && $codon_a =~ m/cg./i )
    || ( $mut_pos == 2 && $codon_a =~ m/.cg/i ) );
  $mut_cpg = 2
    if ( ( $mut_pos == 2 && $codon_a =~ m/cg./i )
    || ( $mut_pos == 3 && $codon_a =~ m/.cg/i ) );
  if ( $nuc_left ne '' ) {

    # Look for CpG overlaps on the boundary, too.
    $mut_cpg = 2 if ( ( $mut_pos == 1 && $codon_context_a =~ m/cg.../i ) );
  }
  if ( $nuc_right ne '' ) {
    $mut_cpg = 1 if ( ( $mut_pos == 3 && $codon_context_a =~ m/...cg/i ) );
  }
  my $codon_cpg = 0;
  if ( $codon_context_a ne '' ) {
    $codon_cpg = 1 if ( $codon_context_a =~ m/cg/i );
  } else {
    $codon_cpg = 1 if ( $codon_a =~ m/cg/i );
  }

  # Look on opposite strand to assign mut_rev_cpg and codon_rev_cpg
  my $mut_rev_cpg = 0;

  $mut_rev_cpg = 2
    if ( ( $mut_pos == 1 && $codon_a =~ m/cg./i )
    || ( $mut_pos == 2 && $codon_a =~ m/.cg/i ) );
  $mut_rev_cpg = 1
    if ( ( $mut_pos == 2 && $codon_a =~ m/cg./i )
    || ( $mut_pos == 3 && $codon_a =~ m/.cg/i ) );

  if ( $nuc_right ne '' ) {
    $mut_rev_cpg = 2 if ( $mut_pos == 3 && $codon_context_a =~ m/...cg/i );
  }

  if ( $nuc_left ne '' ) {
    $mut_rev_cpg = 1 if ( $mut_pos == 1 && $codon_context_a =~ m/cg.../i );
  }

  my $grantham = $self->_grantham_matrix;
  my $aa_key = $sub_obj->{aa_a} . $sub_obj->{aa_b};
  my $grantham_score = $grantham->{$aa_key};

  # Get the HCG pattern relative to orang.
  my $pattern = $self->_get_sub_pattern($tree, $cdna_aln, $aln, $aln_pos);

  my $mut_ils = 0;
  $mut_ils = 1 if ($pattern eq '101' || $pattern eq '011');

  # Count non-gap sites.
  my $nongaps = Bio::EnsEMBL::Compara::AlignUtils->get_nongaps_at_column( $aln, $aln_pos );

  my $base_params = {
    codeml_node_id => $id,
    codeml_node_number => $sub_obj_id,
    leaves_beneath => $leaves_beneath,
    aln_pos    => $aln_pos,
    
    codon_from         => $sub_obj->{codon_a},
    codon_to           => $sub_obj->{codon_b},
    aa_from            => $sub_obj->{aa_a},
    aa_to              => $sub_obj->{aa_b},
    grantham_score     => $grantham_score,
    confidence         => $sub_obj->{confidence},

    codon_cpg => $codon_cpg,

    nuc_from    => $nuc_a,
    nuc_to      => $nuc_b,
    mut_count   => $mut_count,
    mut_pos     => $mut_pos,
    mut_ws      => $mut_ws,
    mut_sw      => $mut_sw,
    mut_ts      => $mut_ts,
    mut_tv      => $mut_tv,
    mut_cpg     => $mut_cpg,
    mut_rev_cpg => $mut_rev_cpg,
    mut_syn     => $mut_syn,
    mut_nsyn    => $mut_nsyn,
    mut_ils     => $mut_ils,

    pattern => $pattern,
    nongap_count => $nongaps,
  };

  $final_params = $self->replace($final_params,$base_params);
  return $final_params;
}

sub _get_sub_pattern {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pep_aln = shift;
  my $pos = shift;

  my $cdna_aln_pos = ( $pos - 1 ) * 3 + 1;

  my $map;
  foreach my $leaf ($tree->leaves) {
    if ($leaf->can('taxon_id')) {
      $map->{$leaf->name} = $leaf->taxon_id;
    }
  }
  $pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate_ids($pep_aln, $map);
  $aln = Bio::EnsEMBL::Compara::AlignUtils->translate_ids($aln, $map);

  my @seqs = $pep_aln->each_seq;
  my @cdna_seqs = $aln->each_seq;

  my ($h_seq) = grep {$_->id == 9606} @seqs;
  my ($c_seq) = grep {$_->id == 9598} @seqs;
  my ($g_seq) = grep {$_->id == 9593} @seqs;
  my ($o_seq) = grep {$_->id == 9600} @seqs;

  foreach my $seq ($h_seq, $c_seq, $g_seq, $o_seq) {
    return undef if (!defined $seq);
  }

  my $h_str = $h_seq->seq;
  my $c_str = $c_seq->seq;
  my $g_str = $g_seq->seq;
  my $o_str = $o_seq->seq;

  my ($h_cdna_seq) = grep {$_->id == 9606} @cdna_seqs;
  my ($c_cdna_seq) = grep {$_->id == 9598} @cdna_seqs;
  my ($g_cdna_seq) = grep {$_->id == 9593} @cdna_seqs;
  my ($o_cdna_seq) = grep {$_->id == 9600} @cdna_seqs;
  my $h_cdna_str = $h_cdna_seq->seq;
  my $c_cdna_str = $c_cdna_seq->seq;
  my $g_cdna_str = $g_cdna_seq->seq;
  my $o_cdna_str = $o_cdna_seq->seq;

  my $h_codon = substr($h_cdna_str, $cdna_aln_pos-1, 3);
  my $c_codon = substr($c_cdna_str, $cdna_aln_pos-1, 3);
  my $g_codon = substr($g_cdna_str, $cdna_aln_pos-1, 3);
  my $o_codon = substr($o_cdna_str, $cdna_aln_pos-1, 3);

  foreach my $codon ($h_codon, $c_codon, $g_codon, $o_codon) {
    return undef if ($codon =~ m/[-n]/i);
  }

  my $j = 1;
  my $map;
  my $o = 1;
  my $h = 0;
  my $c = 0;
  my $g = 0;
  $map->{$o_codon} = $j++;
  $h = $map->{$h_codon} || $j++;
  $map->{$h_codon} = $h;
  $c = $map->{$c_codon} || $j++;
  $map->{$c_codon} = $c;
  $g = $map->{$g_codon} || $j++;
  $map->{$g_codon} = $g;
  
  my $pattern = join('', map {$_ - 1} ($h, $c, $g));
  
  return $pattern;
}

1;
