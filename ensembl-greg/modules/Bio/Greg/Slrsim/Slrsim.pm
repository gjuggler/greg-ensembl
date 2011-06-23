package Bio::Greg::Slrsim::Slrsim;

use strict;
use Bio::EnsEMBL::Hive::Process;
use Bio::Greg::Hive::Process;
use File::Path;
use File::Copy;
use File::Basename;
use FreezeThaw qw(freeze thaw cmpStr safeFreeze cmpStrHard);

use base (
  'Bio::Greg::StatsCollectionUtils',
  'Bio::Greg::Hive::Process',
  'Bio::Greg::Hive::PhyloSim',
  'Bio::Greg::Hive::Align',
  'Bio::Greg::Hive::AlignmentScores',
  'Bio::Greg::Hive::PhyloAnalysis',
  );

sub param_defaults {
  my $self = shift;

  my $phylosim = Bio::Greg::Hive::PhyloSim::param_defaults;
  my $align = Bio::Greg::Hive::Align::param_defaults;
  my $alignment_scores = Bio::Greg::Hive::AlignmentScores::param_defaults;
  my $phylo_analysis = Bio::Greg::Hive::PhyloAnalysis::param_defaults;

  return $self->replace($phylosim, $align, $alignment_scores, $phylo_analysis, {
    genes_table                        => 'genes',
    sites_table                        => 'sites',
    force_recalc => 0
                        });
}

sub fetch_input {
  my $self = shift;

  $self->load_all_params();

  # Load params from the params txt file.
  my $params_f = $self->_save_file('params', 'txt');
  print "  loading params from file [". $params_f->{full_file}."]\n";
  my $params = $self->thw($params_f->{full_file});

  #$self->hash_print($self->params);
  #$self->hash_print($params);

  foreach my $key (keys %$params) {
    if ($params->{$key} ne $self->param($key)) {
      $self->param($key, $params->{$key});
    }
  }

  $self->create_table_from_params( $self->compara_dba, $self->param('sites_table'),
                                   $self->_sites_table_structure );
  
  $self->create_table_from_params( $self->compara_dba, $self->param('genes_table'),
                                   $self->_genes_table_structure );

  $self->dbc->disconnect_when_inactive(1);

#  print $self->get_output_folder."\n";
#  my $data_file = $self->get_output_folder . "/data.tar.gz";
#  $self->param('data_tarball', $data_file) if (-e $data_file);  
}

sub run {
  my $self = shift;

#  $self->param('force_recalc', 1);
#  $self->param('filter', 'none');
#  $self->param('maximum_mask_fraction', 0.6);
#  $self->param('alignment_score_filtering', 1);

  my $tree = $self->_get_tree;
  print $tree->ascii."\n";
  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);

  # Output the tree to file.
  my $tree_f = $self->_save_file('tree', 'nh');
  my $out_file = $tree_f->{full_file};
  $self->param('tree_file', $tree_f->{rel_file});
  Bio::EnsEMBL::Compara::TreeUtils->to_file($tree, $out_file);

  # First, simulate the true alignment.
  my $sim_obj = $self->_simulate_alignment($tree);
  my $true_aln = $sim_obj->{aln};
  my $true_sitewise_hash = $sim_obj->{omegas};
  my $true_pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($true_aln);

  $self->param('true_aln', $true_aln);
  $self->param('true_pep_aln', $true_pep_aln);
  
  # Now, align the sequences.
  my $inferred_aln = $self->_align($tree, $true_aln, $true_pep_aln);
  my $inferred_pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($inferred_aln);

  # Mask the alignment if needed.
  my $masked_aln = $self->_mask_alignment($tree, $inferred_aln, $inferred_pep_aln);
  my $masked_pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($masked_aln);

  # Run the sitewise analysis.
  my $sitewise_hash = $self->_run_sitewise($tree, $masked_aln, $masked_pep_aln);

  # Collect and store results.
  $self->_collect_and_store_results($tree, $treeI,
                                    $masked_aln, $masked_pep_aln, $sitewise_hash, 
                                    $true_aln, $true_pep_aln, $true_sitewise_hash,
                                    $inferred_aln, $inferred_pep_aln
    );

}

sub _get_tree {
  my $self = shift;

  my $tree_file = $self->param('slrsim_tree_file');
  my $mpl = $self->param('slrsim_tree_mean_path');

  $tree_file = Bio::Greg::EslrUtils->baseDirectory."/projects/slrsim/trees/".$tree_file;
  print "$tree_file\n";

  my $tree = Bio::EnsEMBL::Compara::TreeUtils->from_file($tree_file);
  $tree = Bio::EnsEMBL::Compara::TreeUtils->scale_mean_to( $tree, $mpl );
  return $tree;
}

sub _simulate_alignment {
  my $self = shift;
  my $tree = shift;

  my $out_f = $self->_save_file('sim', 'perlobj');
  my $out_file = $out_f->{full_file};

  print "OUTFILE: $out_file\n";

  $self->param('sim_file', $out_f->{rel_file});

  # Copy the cached aln over if it makes sense.
  my $cached_sim = $self->_cached_file('sim.perlobj');
  $self->param('cached_sim', $cached_sim);
  if (!-e $out_file && -e $cached_sim) {
    print "  copying cached simulation...\n";
    copy($cached_sim, $out_file) or die("Error copying cached alignment!");
  }
  
  if (!-e $out_file || $self->param('force_recalc')) {
    print "  running simulation\n";
    my $obj = $self->simulate_alignment($tree);
    $self->frz($out_file, $obj);
  } else {
    print "  loading simulated aln from file [$out_file]\n";
  }

  if (!-e $cached_sim) {
    print "  copying simulation to cache..\n";
    mkpath(dirname($cached_sim));
    copy($out_file, $cached_sim) or die("Error copying cached simulation!");
  }
  print "FILE: $out_file\n";

  my $obj = $self->thw($out_file);
  return $obj;
}

sub _align {
  my $self = shift;
  my $tree = shift;
  my $true_aln = shift;
  my $true_pep_aln = shift;

  $self->pretty_print($true_aln);
  $self->pretty_print($true_pep_aln);

  my $out_f = $self->_save_file('inferred_aln', 'fasta');
  my $out_file = $out_f->{full_file};
  $self->param('aln_file', $out_f->{rel_file});

  # Copy the cached aln over if it makes sense.
  my $cached_aln = $self->_cached_file('aln.fasta');
  $self->param('cached_aln', $cached_aln);
  if (!-e $out_file && -e $cached_aln) {
    print "  copying cached alignment...\n";
    copy($cached_aln, $out_file) or die("Error copying cached alignment!");
  }
  
  if (!-e $out_file || $self->param('force_recalc')) {
    print "  running alignment\n";
    # Calls Bio::Greg::Hive::Align->align
    my $inferred_aln = $self->align($tree, $true_aln, $true_pep_aln);
    Bio::EnsEMBL::Compara::AlignUtils->to_file($inferred_aln, $out_file);
  } else {
    print "  loading inferred aln from file [$out_file]\n";
  }

  if (!-e $cached_aln) {
    print "  copying alignment to cache..\n";
    mkpath(dirname($cached_aln));
    copy($out_file, $cached_aln) or die("Error copying cached alignment!");
  }

  my $inferred_aln = Bio::EnsEMBL::Compara::AlignUtils->from_file($out_file);
  return $inferred_aln;
}

sub _cached_file {
  my $self = shift;
  my $filename = shift;

  # Caching: common aligned files for each indel/mpl/rep
  my $aln_id = $self->param('slrsim_label').'_'.$self->param('slrsim_rep');
  my $f = $self->get_output_folder . '/common';
  my $rep = $self->param('slrsim_rep');
  my $cache_aln_dir = $f . "/$rep";
  mkpath($cache_aln_dir);
  my $cache_aln_f = $cache_aln_dir."/".$aln_id."_".$filename;

  return $cache_aln_f;
}

sub _mask_alignment {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pep_aln = shift;

  my $scores_f = $self->_save_file('aln_scores', 'perlobj');
  my $scores_file = $scores_f->{full_file};
  $self->param('aln_scores_file', $scores_f->{rel_file});
  if (!-e $scores_file || $self->param('force_recalc')) {
    print "  calculating alignment scores\n";
    my $scores_hash = $self->_get_alignment_scores($tree, $aln, $pep_aln);
    $self->hash_print($scores_hash);
    open(OUT,">$scores_file");
    print OUT freeze($scores_hash);
    close(OUT);
  } else {
    print "  loading scores from file [$scores_file]\n";
  }
  open(IN, $scores_file);
  my @lines = <IN>;
  close(IN);
  my ($scores_hash) = thaw(join('',@lines));

  my $out_f = $self->_save_file('masked_aln', 'fasta');
  my $out_file = $out_f->{full_file};
  $self->param('masked_aln_file', $out_f->{rel_file});

  print "Pep aln:\n";
  $self->pretty_print($pep_aln);

  if (!-e $out_file || $self->param('force_recalc')) {
    print "  masking alignment\n";
    # Calls Bio::Greg::Hive::AlignmentScores->mask_alignment
    my $masked_aln = $self->mask_alignment($tree, $aln, $pep_aln, $scores_hash);
    Bio::EnsEMBL::Compara::AlignUtils->to_file($masked_aln, $out_file);
  } else {
    print "  loading masked aln from file [$out_file]\n";
  }

  my $masked_aln = Bio::EnsEMBL::Compara::AlignUtils->from_file($out_file);

  print "Masked:\n";
  my $tx_masked = Bio::EnsEMBL::Compara::AlignUtils->translate($masked_aln);
  $self->pretty_print($tx_masked, {full => 1});

  return $masked_aln;
}


sub _run_sitewise {
  my $self = shift;
  my $tree_orig = shift;
  my $aln_orig = shift;
  my $pep_aln = shift;

  my $tree = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree_orig);
  my $aln = Bio::EnsEMBL::Compara::AlignUtils->copy_aln($aln_orig);

  my $leaf_map;
  my $seq_map;
  map {$leaf_map->{$_->name} = 's_'.$_->name;} $tree->leaves;
  map {$seq_map->{$_->id} = 's_'.$_->id} $aln->each_seq;
  $tree = Bio::EnsEMBL::Compara::TreeUtils->translate_ids($tree, $leaf_map);
  $aln = Bio::EnsEMBL::Compara::AlignUtils->translate_ids($aln, $seq_map);

  my $out_suffix = 'sitewise';
  if ($self->param('analysis_action') =~ m/paml/i) {
    my $model;
    my $action = $self->param('analysis_action');
    if ($action =~ m/2/g) {
      $model = 'M2';
    } elsif ($action =~ m/8/g) {
      $model = 'M8';
    } else {
      die("Unrecognized analysis_action parameter for PAML");
    }
    print "PAML MODEL: $model\n";
    if ($model) {
      $out_suffix .= '_'.$model;
    }

    my $null_f = $self->_save_file('lrt_null', 'out');
    my $null_file = $null_f->{full_file};
    $self->param('LRT null file', $null_f->{rel_file});
    my $alt_f = $self->_save_file('lrt_alt', 'out');
    my $alt_file = $alt_f->{full_file};
    $self->param('LRT alt file', $alt_f->{rel_file});

    if (!-e $null_file || !-e $alt_file || $self->param('force_recalc')) {
      my $null_lines;
      my $alt_lines;

      my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);

      if ($model =~ m/8/) {
        # M7 vs M8 LRT
        my $params = {
          NSsites => 8
        };
        my $res = Bio::Greg::Codeml->branch_model_likelihood($treeI,$aln, $self->worker_temp_directory, $params);
        $alt_lines = $res->{lines};
        
        my $params = {
          NSsites => 7
        };
        my $res = Bio::Greg::Codeml->branch_model_likelihood($treeI,$aln, $self->worker_temp_directory, $params);
        $null_lines = $res->{lines};
      } else {
        # M1 vs M2 LRT

        my $params = {
          NSsites => 2
        };
        my $res = Bio::Greg::Codeml->branch_model_likelihood($treeI,$aln, $self->worker_temp_directory, $params);
        $alt_lines = $res->{lines};
        
        my $params = {
          NSsites => 1
        };
        my $res = Bio::Greg::Codeml->branch_model_likelihood($treeI,$aln, $self->worker_temp_directory, $params);
        $null_lines = $res->{lines};
      } 

      print "  Outputting Null and Alt models!\n";
      open(OUT, ">$null_file");
      print OUT join("", @{$null_lines});
      close(OUT);
      open(OUT, ">$alt_file");
      print OUT join("", @{$alt_lines});
      close(OUT);
    }

    print "  loading LRTs from file...\n";
    print "  $null_file\n";
    open(IN, $null_file);
    my @null_lines = <IN>;
    close(IN);
    print "  $alt_file\n";
    open(IN, $alt_file);
    my @alt_lines = <IN>;
    close(IN);

    my $null_lnl = Bio::Greg::Codeml->extract_lnL(\@null_lines);
    my $alt_lnl = Bio::Greg::Codeml->extract_lnL(\@alt_lines);
    
    $self->param('lnl_null', $null_lnl);
    $self->param('lnl_alt', $alt_lnl);
  }

  $self->param('force_recalc', 1);

  my $out_f = $self->_save_file($out_suffix, 'out');
  my $out_file = $out_f->{full_file};
  $self->param('sitewise_file', $out_f->{rel_file});

  if (!Bio::EnsEMBL::Compara::AlignUtils->has_any_data($aln)) {
    print "No data!!\n";
    return {};
  }

  print "$out_file\n";
  if (!-e $out_file || $self->param('force_recalc')) {
    print("  running sitewise analysis\n");
    my $output_lines = $self->run_sitewise_analysis($tree, $aln, $pep_aln);
    open(OUT, ">$out_file");
    print OUT join("", @{$output_lines});
    close(OUT);
  }

  $self->param('force_recalc', 0);

  print("  loading sitewise results [$out_file]\n");
  my $sitewise_results = $self->parse_sitewise_file($tree, $aln, $pep_aln, $out_file);
  return $sitewise_results;
}

sub _collect_and_store_results {
  my $self = shift;
  my $tree = shift;
  my $treeI = shift;
  my $aln = shift;
  my $pep_aln = shift;
  my $sitewise_hash = shift;
  my $true_aln = shift;
  my $true_pep_aln = shift;
  my $true_sitewise_hash = shift;
  my $unfiltered_aln = shift;
  my $unfiltered_pep_aln = shift;

  $self->pretty_print($pep_aln);
  $self->pretty_print($true_pep_aln);

  # Get the sequence to act as a reference in site-wise value comparisons.
  my $reference_id = $self->param('slrsim_ref') || 'human';
  my @seqs = $pep_aln->each_seq;
  my ($ref_seq) = grep { $_->id eq $reference_id } @seqs;
  $ref_seq = $seqs[0] unless ( defined $ref_seq );
  my $seq_str = $ref_seq->seq;
  $seq_str =~ s/-//g;
  my $id = $ref_seq->id;

  #$self->dbc->do('LOCK TABLES sites WRITE');
  #$self->dbc->do('UNLOCK TABLES');
  
  my $aln_scores_calc;
  my @aln_entropies;
  my @aln_aligned;
  my @aln_match;
#  my $aln_scores_calc = Bio::EnsEMBL::Compara::AlignUtils->_correct_subtree_calc($treeI, $true_pep_aln, $pep_aln);
#  my @aln_entropies = Bio::EnsEMBL::Compara::AlignUtils->column_entropies($aln);
#  my @aln_aligned = @{$aln_scores_calc->{aligned_branchlengths}};
#  my @aln_match = @{$aln_scores_calc->{match_branchlengths}};

  foreach my $seq_position ( 1 .. length($seq_str) ) {
    my $true_column = $true_pep_aln->column_from_residue_number( $id, $seq_position );
    my $aln_column = $pep_aln->column_from_residue_number( $id, $seq_position );
    my $true_hash_item = $true_sitewise_hash->{$true_column};
    my $true_dnds = $true_hash_item->{omega};
    my $aln_hash_item = $sitewise_hash->{$aln_column};    
    my $aln_dnds = $aln_hash_item->{omega};
    
    $true_dnds = 0;
    $aln_dnds = 0;
    if ( !( defined $aln_dnds && defined $true_dnds ) ) {
      if ( !defined $true_dnds ) {
        my $str = sprintf("No true dnds! aln:%s true:%s\n", $aln_column, $true_column);
        die($str);
      } elsif ( !defined $aln_dnds ) {
        die(sprintf("  no aln dnds! aln:%s true:%s\n", $aln_column, $true_column));
      }
    }    

    my $obj = {};
    $obj->{true_type}      = $true_sitewise_hash->{$true_column}->{'type'}      || '';
    $obj->{true_dnds}      = $true_sitewise_hash->{$true_column}->{'omega'} || '';

    if (!defined $aln_hash_item) {
      $obj->{aln_dnds} = 0;
      $obj->{lrt_stat} = 0;
    } else {
      $obj->{aln_type}       = $sitewise_hash->{$aln_column}->{'type'}        || '';
      $obj->{aln_dnds}       = $sitewise_hash->{$aln_column}->{'omega'} || '';
      $obj->{aln_note}       = $sitewise_hash->{$aln_column}->{'note'}        || '';
      $obj->{lrt_stat}        = $sitewise_hash->{$aln_column}->{'lrt_stat'}    || '';
      $obj->{aln_dnds_lower} = $sitewise_hash->{$aln_column}->{'omega_lower'} || '';
      $obj->{aln_dnds_upper} = $sitewise_hash->{$aln_column}->{'omega_upper'} || '';
    }

    if (defined @aln_entropies) {
#    $obj->{entropy} = $aln_entropies[$aln_column-1];
    }
    if (defined @aln_aligned) {
#    $obj->{bl_aligned} = $aln_aligned[$aln_column-1];    
    }
    if (defined @aln_match) {
#    $obj->{bl_match} = $aln_match[$aln_column-1];
    }

    $obj->{aln_position} = $aln_column;
    $obj->{seq_position} = $seq_position;

    # Important: delete lrt_stat value if this column is all_gaps or single_char
    if ($obj->{aln_note} eq 'all_gaps' || $obj->{aln_note} eq 'single_char') {
      delete $obj->{lrt_stat};
    }
    
    my $params = $self->params;
    $params = $self->replace($params, $obj);
    printf "%d %d %f %f\n", $true_column, $aln_column, $obj->{true_dnds}, $obj->{lrt_stat} if ($seq_position < 20 || $seq_position > length($seq_str)-20);
    $self->store_params_in_table( $self->dbc, $self->param('sites_table'), $params );
    sleep(0.2);
  }

  # Collect some alignment-wide calculations.
  # Alignment.
  my $aln_length = $pep_aln->length;
  $self->param('aln_length', $aln_length);

  # Mean branch lengths.
  if (defined $aln_scores_calc) {
    my $aligned_per_site = $aln_scores_calc->{aligned_bl} / $aln_length;
    my $match_per_site = $aln_scores_calc->{match_bl} / $aln_length;
    my $mismatch_per_site = $aln_scores_calc->{mismatch_bl} / $aln_length;
    $self->param('mean_bl_aligned', $aligned_per_site);
    $self->param('mean_bl_match', $match_per_site);
    $self->param('mean_bl_mismatch', $mismatch_per_site);
  }

  my $do_alignment_stats = 0;
  if ($do_alignment_stats) {
    my $sps = Bio::EnsEMBL::Compara::AlignUtils->sum_of_pairs_score($true_pep_aln, $pep_aln);
    $self->param('sum_of_pairs_score', $sps);
    my $tcs = Bio::EnsEMBL::Compara::AlignUtils->total_column_score($true_pep_aln, $pep_aln);
    $self->param('total_column_score', $tcs);
    my $match_bl_score = $aln_scores_calc->{match_bl} / $aln_scores_calc->{aligned_bl};
    $self->param('match_bl_score', $match_bl_score);
    my $mismatch_bl_score = $aln_scores_calc->{mismatch_bl} / $aln_scores_calc->{aligned_bl};
    $self->param('mismatch_bl_score', $mismatch_bl_score);
#  my $lambda = Bio::EnsEMBL::Compara::AlignUtils->dawg_lambda($pep_aln, $tree, $self->params, $self->worker_temp_directory);
#  $self->param('lambda', $lambda);

    # Column entropies.
    my $ce_aln = Bio::EnsEMBL::Compara::AlignUtils->average_column_entropy($aln);
    $self->param('mean_entropy', $ce_aln);

  }

  # Masked fraction.
  my $unmasked_n = Bio::EnsEMBL::Compara::AlignUtils->count_residues($unfiltered_pep_aln);
  my $masked_n = Bio::EnsEMBL::Compara::AlignUtils->count_residues($pep_aln);
  $self->param('unmasked_residue_count', $unmasked_n);
  $self->param('residue_count', $masked_n);

  # Tree length.
  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
  my $tree_length = $treeI->total_branch_length;
  $self->param('tree_total_length', $tree_length);

  my $mean_path = Bio::EnsEMBL::Compara::TreeUtils->mean_path($tree);
  $self->param('tree_mean_path', $mean_path);

  $self->param('job_id', $self->job_id);
  $self->store_params_in_table( $self->dbc, $self->param('genes_table'), $self->params);
}

sub _save_file {
  my $self = shift;
  my $filename_base = shift;
  my $ext = shift;

  my $rep = $self->param('slrsim_rep');
  my $id = $self->param('slrsim_label');

  my $filename = "${id}_${filename_base}_${rep}";

  my $file_params = {
    id => $id,
    filename => $filename,
    extension => $ext,
    subfolder => 'data',
  };

  my $file_obj = $self->save_file($file_params);

  if ($self->param('data_tarball') && !-e $file_obj->{full_file}) {
    # Load the file from the tarball into our temp directory.
    #my $tarball = $self->param('data_tarball');
    #my $tmp = $self->worker_temp_directory;
    #my $rel_f = $file_obj->{rel_file};

    #my $cmd = qq^tar -zxvf $tarball $rel_f^;
    #print $cmd."\n";
    #$file_obj->{full_file} = $tmp . $filename;
  }

  return $file_obj;
}

sub _genes_table_structure {
  my $structure = {
    node_id => 'int',
    job_id => 'int',

    data_prefix => 'char4',

    slrsim_label                        => 'char64',
    slrsim_analysis_name                => 'string',

    aligner => 'string',
    filter => 'string',
    analysis_action                     => 'string',
    alignment_score_threshold           => 'float',

    slrsim_rep         => 'int',
    slrsim_tree_file   => 'string',
    slrsim_ref         => 'string',
    slrsim_tree_mean_path => 'float',

    phylosim_seq_length         => 'int',
    phylosim_omega_distribution => 'string',
    phylosim_meanlog            => 'float',
    phylosim_sdlog              => 'float',
    phylosim_insertrate  => 'float',
    phylosim_deleterate  => 'float',
    phylosim_insertmodel => 'string',
    phylosim_deletemodel => 'string',

    aln_length => 'int',
    mean_bl_aligned => 'float',
    mean_bl_match => 'float',
    mean_bl_mismatch => 'float',

    lnl_null => 'float',
    lnl_alt => 'float',

    sum_of_pairs_score => 'float',
    total_column_score => 'float',
    match_bl_score => 'float',
    mismatch_bl_score => 'float',
    lambda => 'float',
    mean_entropy  => 'float',

    maximum_mask_fraction => 'float',
    unmasked_residue_count    => 'int',
    residue_count    => 'int',

    tree_total_length => 'float',
    tree_mean_path => 'float',

    unique_keys              => 'node_id,slrsim_rep'
  };
  return $structure;
}

sub _sites_table_structure {
  my $structure = {
    node_id          => 'int',
    slrsim_rep => 'int',

    # Site-wise stuff.
    aln_position => 'int',
    seq_position => 'int',

    bl_match => 'float',
    bl_aligned => 'float',
    entropy => 'float',

    true_dnds    => 'float',
    true_type    => 'char16',

    aln_dnds       => 'float',
    aln_type       => 'char16',
    aln_note       => 'char16',
    lrt_stat => 'float',

    unique_keys => 'node_id,slrsim_rep,aln_position'
  };
  return $structure;
}


1;
