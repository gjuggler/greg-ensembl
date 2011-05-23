package Bio::Greg::Gorilla::BranchTests;

use strict;
use File::Path;
use Bio::Greg::Codeml;
use FreezeThaw qw(freeze thaw cmpStr safeFreeze cmpStrHard);
use Data::Types qw(:all);

use base (
  'Bio::Greg::StatsCollectionUtils',
  'Bio::Greg::Hive::Process',
  'Bio::Greg::Hive::Align',
  'Bio::Greg::Hive::AlignmentScores',
  'Bio::Greg::Hive::CountSubstitutions',
  'Bio::Greg::Hive::PhyloAnalysis',
);

sub param_defaults {
  return {
  };
}

sub fetch_input {
  my ($self) = @_;

  # Fetch parameters from all possible locations.
  $self->load_all_params();

  $self->{_genes_structure} = $self->genes_table;

  # Create tables if necessary.
  $self->create_table_from_params( $self->compara_dba, 'subs',
    $self->codon_subs_table_def );
}

sub run {
  my $self = shift;

  my $gene_id = $self->param('gene_id');
  my $mba = $self->compara_dba->get_MemberAdaptor;
  my $member = $mba->fetch_by_source_stable_id(undef, $gene_id);
  my $gene = $member->get_Gene;
  $self->param('gene_name', $gene->external_name);
  $self->param('data_id', $member->dbID);

  my $ortholog_tree = Bio::EnsEMBL::Compara::ComparaUtils->get_one_to_one_ortholog_tree($self->compara_dba, $member);

  my ($tree, $aln) = $self->_get_tree_and_alignment($ortholog_tree);
  $self->pretty_print($aln);
  print $tree->ascii."\n";

  # Fail if we don't have any outgroup primates.
#  my @outgroups = grep {
#    $_->taxon_id == 9483 ||
#      $_->taxon_id == 9478 ||
#      $_->taxon_id == 9600 ||
#      $_->taxon_id == 9544 ||
#      $_->taxon_id == 30611
#  } $tree->leaves;
#  if (scalar(@outgroups) == 0) {
#    $self->fail("any_outgroup", "No non-HCG primate outgroup -- FAIL!");
#  }
  
  # Fail if we don't have at least h-g-c-o sequences
  my $req = {
    9606 => 0,
    9598 => 0,
    9593 => 0,
    9600 => 0,
    9544 => 0,
    9483 => 0
  };
  map {
    $req->{$_->taxon_id} = 1;
  } $tree->leaves;
  my @missing = grep {$req->{$_} == 0} keys %$req;
  if (scalar(@missing) > 0) {
    my $missing_string = join(', ', sort {$a <=> $b} @missing);
    $self->fail("one2one", "Missing $missing_string");
  }


  

  my $params = {
    keep_species => join(',', map {$_->taxon_id} $tree->leaves)
  };
  my $genome_tree = Bio::EnsEMBL::Compara::ComparaUtils->get_genome_tree_subset($self->compara_dba, $params);

  # Match up genome tree names to gene tree IDs.
  my $id_hash;
  map {$id_hash->{$_->taxon_id} = $_->name} $ortholog_tree->leaves;
  foreach my $node ($genome_tree->nodes) {
    if ($node->is_leaf) {
      $node->name($id_hash->{$node->taxon_id});
    } else {
      $node->name('');
    }
  }

  $self->dbc->disconnect_when_inactive(1);
  $self->_run_slr($ortholog_tree, $aln);
  $self->param('aln_use_type', 'genomic_primates');
  $aln = $self->_mask_sub_runs($genome_tree, $aln);
  $self->_get_sub_patterns_from_aln($tree, $aln);
  $aln = $self->_mask_ils($genome_tree, $aln);
  $self->_test_for_aln_coverage($genome_tree, $aln);

  $self->_run_tests($genome_tree, $aln);
  $self->_plot_subs($genome_tree, $aln);
  $self->_collect_stats($genome_tree, $aln);
}

sub write_output {
  my $self = shift;

  $self->create_table_from_params( $self->dbc, 'genes', $self->{_genes_structure});
  $self->store_params_in_table($self->dbc, 'genes', $self->params);
}

sub _tx_aln {
  my $self = shift;
  my $aln = shift;
  return Bio::EnsEMBL::Compara::AlignUtils->translate($aln);
}

sub _get_tree_and_alignment {
  my $self = shift;
  my $tree = shift;

  #my @aln_types = ('compara', 'genomic_primates', 'genomic_mammals', 'genomic_all');
  my @aln_types = ('genomic_primates');
  my $use_type = 'genomic_primates';

  my $use_file;

  foreach my $type (@aln_types) {
    my $tree_copy = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree);

    my $f = $self->_save_file('aln_'.$type, 'fasta');
    my $file = $f->{full_file};
    $self->store_param('aln_'.$type.'_file', $f->{rel_file});

    if ($type eq $use_type) {
      $use_file = $f;
    }

    if (!-e $file || $self->param('force_recalc')) {
      my $ref_species = 9606;
      my @members     = $tree_copy->leaves;
      my ($ref_member) = grep { $_->taxon_id == $ref_species } @members;
      die("No ref member!") unless (defined $ref_member);

      # Fetch the alignment.
      my $params = $self->params;
      $params->{aln_type} = $type;

      $params->{quality_threshold} = 30;
      my $tree_aln_obj =
        Bio::EnsEMBL::Compara::ComparaUtils->get_compara_or_genomic_aln( 
          $self->compara_dba, $tree_copy, $ref_member, $params
        );
      if ($tree_aln_obj == -1) {
        $self->fail("alignment", "No tree_aln object from ComparaUtils -- probably just a single human sequence in the EPO alignments");
      }
      my $aln = $tree_aln_obj->{aln};
      my $aln_tree = $tree_aln_obj->{tree};
      
      # Realign compara alignments w/ PRANK...
#      if ($type eq 'compara') {
#        $params->{alignment_pagan_codon_model} = 1;
#        $aln = $self->align_with_pagan($aln, $aln_tree, $params);
#      }

      Bio::EnsEMBL::Compara::AlignUtils->to_file($aln, $file);
      #Bio::EnsEMBL::Compara::AlignUtils->to_file($raw_aln, $raw_file);
    }
  }

  my $tree_f = $self->_save_file('tree', 'nh');
  my $tree_file = $tree_f->{full_file};
  $self->store_param('tree_file', $tree_f->{rel_file});
  if (!-e $tree_file || $self->param('force_recalc')) {
    Bio::EnsEMBL::Compara::TreeUtils->to_file($tree, $tree_file);
  }
  my $tree_from_file = Bio::EnsEMBL::Compara::TreeUtils->from_file($tree_file);
  my $aln = Bio::EnsEMBL::Compara::AlignUtils->from_file($use_file->{full_file});

  # Get sub-tree of primate sequences.
  $tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_clade(
    $self->compara_dba, $tree, "Primates"
    );
  
  # Restrict the tree & alignment to each other.
  $tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_aln($tree, $aln);
  $aln = Bio::EnsEMBL::Compara::ComparaUtils->restrict_aln_to_tree($aln, $tree);

  #$self->param('aln_use_type', $use_type);

  return ($tree, $aln);
}

sub _run_slr {
  my $self = shift;
  my $tree_orig = shift;
  my $aln_orig = shift;

  $self->param('analysis_action', 'slr');
  my $tree = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree_orig);
  my $ref_member = $self->_ref_member($tree);

  my $f = $self->_save_file('aln_compara', 'fasta');
  my $file = $f->{full_file};
  if (!-e $file || $self->param('force_recalc')) {
    # Get and save the Compara alignment.
    my $params = $self->params;
    $params->{aln_type} = 'compara';
    $params->{quality_threshold} = 0;
    my $tree_aln_obj =
    Bio::EnsEMBL::Compara::ComparaUtils->get_compara_or_genomic_aln( 
      $self->compara_dba, $tree_orig, $ref_member, $params
    );
    my $aln = $tree_aln_obj->{aln};  
    Bio::EnsEMBL::Compara::AlignUtils->to_file($aln, $file);
  }  
  my $aln = Bio::EnsEMBL::Compara::AlignUtils->from_file($file);
  
  $tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_clade(
    $self->compara_dba, $tree, "Mammalia"
    );
  
  # Restrict the tree & alignment to each other.
  $tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_aln($tree, $aln);
  $aln = Bio::EnsEMBL::Compara::ComparaUtils->restrict_aln_to_tree($aln, $tree);

  print $tree->ascii."\n";
  $self->pretty_print($aln);

  my $tree_copy = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree);
  foreach my $node ($tree_copy->nodes) {
    my $bl = $node->distance_to_parent;
    $bl = 0.01 if ($bl < 0.01);
    $bl = 1 if ($bl > 1);
    $node->distance_to_parent($bl);
  }

  my $pep_aln = $self->_tx_aln($aln);
  my $slr_f = $self->_save_file('slr', 'out');
  my $slr_file = $slr_f->{full_file};
  if (!-e $slr_file || $self->param('force_recalc')) {
    print "  running SLR\n";
    my $output_lines = $self->run_sitewise_analysis($tree_copy,$aln, $pep_aln);
    $self->frz($slr_file, $output_lines);
  }

  print "  loading SLR results from file\n";
  my $output_lines = $self->thw($slr_file);
  #print join("\n", @$output_lines)."\n";
  my $results = $self->parse_sitewise_output_lines($tree_copy, $aln, $pep_aln, $output_lines);
  $self->store_param('slr_dnds', 0.0 + $results->{omega});

  # Store sitewise results into 'sites' table.
  my $params = $self->params;
  $params->{omega_table} = 'sites';
  $self->store_sitewise($tree_copy, $pep_aln, $results, $params);

  $self->param('slr_results', $results);
}

sub _mask_ils {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  my $to_mask;
  my @subs = @{$self->param('subs')};
  
  my $muts;
  foreach my $s (@subs) {
    if ($s->{mut_ils} == 1) {
      my $aln_pos = $s->{aln_pos};
      $to_mask->{$aln_pos + 0} = 1;
    }
  }

  my $n_masked = 0;
  foreach my $key (sort keys %$to_mask) {
    my $aln_pos = $key;
    foreach my $seq ($aln->each_seq) {
      next unless ($seq->id =~ m/(ENSP0|ENSGGOP0|ENSPTRP0)/gi);
      my $cdna_lo = ($aln_pos-1)*3 + 0;
      my $cdna_hi = ($aln_pos-1)*3 + 2;
    
      my $orig_str = $seq->seq;
      my ($str, $n) = Bio::EnsEMBL::Compara::AlignUtils->mask_string($orig_str,$cdna_lo,$cdna_hi);
      $seq->seq($str);
    }
    $n_masked += 1;
  }
  
  $self->store_param('masked_ils', $n_masked);

  # Write the modified alignment to file.
  my $f = $self->_save_file('aln_masked_ils', 'fasta');
  Bio::EnsEMBL::Compara::AlignUtils->to_file($aln, $f->{full_file});
  
  print "  wrote mutation-masked alignment\n";
  print "  (N masked columns: $n_masked)\n";

  return $aln;  
}

sub _mask_sub_runs {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  my $store_unmasked_file = 1;
  if ($store_unmasked_file) {
    my $f = $self->_save_file('aln_unmasked', 'fasta');
    Bio::EnsEMBL::Compara::AlignUtils->to_file($aln, $f->{full_file});
  }

  my $tree_copy = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree);
  my $aln_copy = Bio::EnsEMBL::Compara::AlignUtils->copy_aln($aln);

  my $m0_f = $self->_save_file('mask_subs_results', 'perlobj');
  if (!-e $m0_f->{full_file} || $self->param('force_recalc')) {
    print "  calculating mask subs results...\n";
    my $lines = $self->run_m0($tree_copy, $aln_copy);
    my @lines = @{$lines};
    $self->frz($m0_f->{full_file}, $lines);
  }
  print "  loading mask subs results from file\n";
  my $lines = $self->thw($m0_f->{full_file});
  my $m0_tree = Bio::Greg::Codeml->parse_codeml_results($lines);

  # Store substitutions.
  $self->store_subs($tree_copy, $aln_copy, $lines);

  my @subs = @{$self->param('subs')};
  
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

  # Write the modified alignment to file.
  my $f = $self->_save_file('aln_masked_subs', 'fasta');
  Bio::EnsEMBL::Compara::AlignUtils->to_file($aln, $f->{full_file});
  
  my $data_id = $self->param('data_id');
  $self->dbc->do("delete from subs where data_id=${data_id};");
  
  print "  wrote mutation-masked alignment\n";
  print "  (N masked nucs: $n_masked)\n";

  return $aln;
}

sub _test_for_aln_coverage {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  print "  testing for alignment coverage\n";
  my $no_coverage;
  foreach my $seq ($aln->each_seq) {
    my $str = $seq->seq;

    my $aln_len = $aln->length;
    
    $str =~ s/[-n]//gi;
    if (length($str) / $aln_len <= 0.8) {
      print "Coverage: " . $seq->id . "\n";
      $no_coverage->{$seq->id} = length($str) / $aln_len;
    }
  }

  if (scalar(keys %$no_coverage) > 0) {
    my @coverages = map { $_ .": ". sprintf("%.3f",$no_coverage->{$_})} keys %$no_coverage;
    my $cov_string = join('; ', @coverages);
    $self->store_param('poor_coverage', $cov_string);
  }
}

sub _run_tests {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  my $tree_copy = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree);
  my $aln_copy = Bio::EnsEMBL::Compara::AlignUtils->copy_aln($aln);
  
  my $m0_f = $self->_save_file('paml_m0_results', 'perlobj');
  if (!-e $m0_f->{full_file} || $self->param('force_recalc')) {
    print "  calculating m0 results...\n";
    my $lines = $self->run_m0($tree_copy, $aln_copy);
    my @lines = @{$lines};
    $self->frz($m0_f->{full_file}, $lines);
  }
  print "  loading m0 results from file\n";
  my $lines = $self->thw($m0_f->{full_file});
  # Store substitutions.
  $self->store_subs($tree_copy, $aln_copy, $lines);

  # Store m0 omega.
  my @omegas = Bio::Greg::Codeml->extract_omegas($lines);
  $self->store_param('m0_dnds', $omegas[0]);
  # Apply m0 branch lengths to our tree.
  my $m0_tree = Bio::Greg::Codeml->parse_codeml_results($lines);
  Bio::EnsEMBL::Compara::TreeUtils->transfer_branchlengths($m0_tree, $tree);
  $tree->branch_length(0.05); # Set small b.l. on the root node.

  my $branch_f = $self->_save_file('paml_branch_results', 'perlobj');
  my $branch_file = $branch_f->{full_file};
  if (!-e $branch_file || $self->param('force_recalc')) {  
    print "  calculating branch results...\n";

    # We'll have two sets of branch results: all, noH, and noC.
    my $branch_results;
    
    my $config = {
      0 => [],
      1 => ['H'],
      2 => ['C'],

      3 => ['H', 'CH'],
      4 => ['C', 'CH'],
      5 => ['G'],

      6 => ['CGH', 'CG', 'HG'],
      7 => ['CH', 'G', 'H', 'C'],

      8 => ['H', 'C'],
      9 => ['G', 'H', 'CH'],
      10 => ['C', 'G', 'CH']
    };

    my $cur_tree;
    my $cur_aln;

    my $noH;
    my $noC;

    my ($human_leaf) = grep {$_->taxon_id==9606} $tree_copy->leaves;
    my ($chimp_leaf) = grep {$_->taxon_id==9598} $tree_copy->leaves;

    # NoH = no human sequence.
    if ($human_leaf) {
      ($cur_tree, $cur_aln) = 
        Bio::EnsEMBL::Compara::ComparaUtils->remove_from_tree_and_aln($tree_copy, $aln_copy, $human_leaf->name);
      print "  running noH\n";
      $noH = $self->run_branches($cur_tree, $cur_aln, $config);
    }
    
    # NoC = no chimp sequence.
    if ($chimp_leaf) {
      ($cur_tree, $cur_aln) = 
        Bio::EnsEMBL::Compara::ComparaUtils->remove_from_tree_and_aln($tree_copy, $aln_copy, $chimp_leaf->name);
      print "  running noC\n";
      $noC = $self->run_branches($cur_tree, $cur_aln, $config);
    }
    
    # all = all sequences.
    print "  running all\n";
    my $all = $self->run_branches($tree_copy, $aln_copy, $config);
    
    my $obj = {
      all => $all
    };
    if ($noH) {
      $obj->{noH} = $noH;
    }
    if ($noC) {
      $obj->{noC} = $noC;
    };
    
    $self->frz($branch_file, $obj);
  }
  print "  loading branch results file\n";
  my $branch_results = $self->thw($branch_file);
  $self->parse_results($branch_results, $tree, $aln);  

  $self->param('run_branchsite', 1);

  if ($self->param('run_branchsite')) {

    my $sites_f = $self->_save_file('paml_branchsites_results', 'perlobj');
    my $sites_file = $sites_f->{full_file};
    $self->param('branchsites_results_file', $sites_f->{rel_file});  
    if (!-e $sites_file || $self->param('force_recalc')) {  
      
      my $config = {
        11 => ['H'],
        12 => ['C'],
        15 => ['G'],
        16 => ['CGH', 'CG', 'HG'],
        17 => ['CH', 'G', 'H', 'C'],
        18 => ['CGHO', 'GHO', 'CGO'], # The 2nd and 3rd items are for the no-C and no-H datasets.
      };
      
      my $cur_tree;
      my $cur_aln;
      
      # all = all sequences.
      print "  running all\n";
      my $all = $self->run_branchsites($tree_copy, $aln_copy, $config);
      
      my $noH;
      my $noC;

      # NoH = no human sequence.
      my ($human_leaf) = grep {$_->taxon_id==9606} $tree_copy->leaves;
      my ($chimp_leaf) = grep {$_->taxon_id==9598} $tree_copy->leaves;

#      if ($human_leaf) {
#        ($cur_tree, $cur_aln) = 
#          Bio::EnsEMBL::Compara::ComparaUtils->remove_from_tree_and_aln($tree_copy, $aln_copy, $human_leaf->name);
#        print "  running noH\n";
#        $noH = $self->run_branchsites($cur_tree, $cur_aln, $config);
#      }
      
#      # NoC = no chimp sequence.
#      if ($chimp_leaf) {
#        ($cur_tree, $cur_aln) = 
#          Bio::EnsEMBL::Compara::ComparaUtils->remove_from_tree_and_aln($tree_copy, $aln_copy, $chimp_leaf->name);
#        print "  running noC\n";
#        $noC = $self->run_branchsites($cur_tree, $cur_aln, $config);
#      }
      
      my $obj = {
        all => $all
      };
      if ($noH) {
        $obj->{noH} = $noH;
      }
      if ($noC) {
        $obj->{noC} = $noC;
      };
      
      $self->frz($sites_file, $obj);

    }
    print "  loading branch-sites results file\n";
    my $branchsites_results = $self->thw($sites_file);
    $self->parse_results($branchsites_results, $tree, $aln);
  }
  
}

sub run_branches {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $config = shift;

  my $lines_hash;

  foreach my $key (sort {$a <=> $b} keys %$config) {
    my $fg_branches = $config->{$key};

    print "  running $key...\n";
    my $copy = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree);
    my $lines = $self->run_branch($copy, $aln, $fg_branches, $key);
    $lines_hash->{$key} = $lines;
  }
  return $lines_hash;
}

sub run_branchsites {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $config = shift;

  my $lines_hash;
  foreach my $key (sort {$a <=> $b} keys %$config) {
    my $fg_branches = $config->{$key};

    print "  running $key...\n";
    my $copy = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree);
    my $results = $self->run_branchsite($copy, $aln, $fg_branches, $key);
    $lines_hash->{$key} = $results;
  }
  return $lines_hash;
}

sub parse_results {
  my $self = shift;
  my $results = shift;
  my $tree = shift;
  my $aln = shift;

  foreach my $major_key (keys %$results) {
    # Major key will be noC, noH, and all
    my $result_set = $results->{$major_key};
    foreach my $minor_key (keys %$result_set) {
      # Minor key will be 0, 1, 2, etc. for each test.
      my $lines = $result_set->{$minor_key};

      if (ref($lines) =~ m/hash/gi) {
        # We've got a branchsites model on our hands.
        my $alt = $lines->{alt};
        my $null = $lines->{null};
        my $lnL = Bio::Greg::Codeml->extract_lnL($alt);
        $self->store_param("${major_key}_${minor_key}_alt_lnL", Bio::Greg::Codeml->extract_lnL($alt));
        $self->store_param("${major_key}_${minor_key}_null_lnL", Bio::Greg::Codeml->extract_lnL($null));
      } else {
        my $lnL = Bio::Greg::Codeml->extract_lnL($lines);
        my @omegas = Bio::Greg::Codeml->extract_omegas($lines);
        
        $self->store_param("${major_key}_${minor_key}_lnL", $lnL);
        for ( my $i = 0 ; $i < scalar @omegas ; $i++ ) {
          my $omega = $omegas[$i];
          $self->store_param("${major_key}_${minor_key}_w_${i}", $omega);
        }
      }
    }
  }
}

sub run_branch {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $fg_branches = shift;
  my $test_num = shift;

  my $key = "$test_num";

  my ($labeled_tree, $labeled_aln) = $self->label_branches($tree, $aln, $fg_branches);
  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($labeled_tree);
  print $treeI->ascii;
  print $treeI->newick."\n";
  my $res = Bio::Greg::Codeml->branch_model_likelihood( $treeI, $labeled_aln, $self->worker_temp_directory, $self->params );
  return $res->{lines};
}

sub run_branchsite {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $fg_branches = shift;
  my $test_num = shift;

  my $key = "$test_num";

  my ($labeled_tree, $labeled_aln) = $self->label_branches($tree, $aln, $fg_branches);
  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($labeled_tree);
  print $treeI->ascii;

  # All branch-sites models should have some labeled branches. If not, then we're missing
  # a needed species from the tree -- return undef for the whole test.
  my $newick = $treeI->to_newick;
  my $has_fg_branch = 0;
  foreach my $node ($treeI->nodes) {
    my $id = $node->id;
    $has_fg_branch = 1 if ($id =~ m/[\#\$]/);
  }
  if (!$has_fg_branch) {
    print "NO FG BRANCH\n";
    return {
      alt => [''],
      null => ['']
    };
  }

  # Null model.
  my $null = $self->params;
  $null->{model} = 2;
  $null->{NSsites} = 2;
  $null->{fix_omega} = 1;
  $null->{omega} = 1;

  # Alternative model.
  my $alt = $self->params;
  $alt->{model} = 2;
  $alt->{NSsites} = 2;
  $alt->{fix_omega} = 0;
  $alt->{omega} = 0.3;

  print "  running alt\n";

  my @starting_omegas = (0.3, 0.8, 1.5);
  my $res_alt;
  my $starting_omega;
  do {
    eval {
      $starting_omega = shift @starting_omegas;
      $alt->{omega} = $starting_omega;
      print "Trying with omega=$starting_omega\n";
      $res_alt = Bio::Greg::Codeml->branch_model_likelihood( $treeI, $labeled_aln, $self->worker_temp_directory, $alt );
    } or do {
      print "Failed codeml with omega=$starting_omega\n";
      next;
    };
  } while (!defined $res_alt && scalar(@starting_omegas) > 0);

  if (!defined $res_alt) {
    return {
      alt => [''],
      null => ['']
    };
  }

  print "  running null\n";
  my $res_null = Bio::Greg::Codeml->branch_model_likelihood( $treeI, $labeled_aln, $self->worker_temp_directory, $null );

  return {
    alt => $res_alt->{lines},
    null => $res_null->{lines}
  };
}


sub label_branches {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $fg_branches = shift;

  my $copy = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree);

  my $char_to_taxon = {
    'H' => 9606,
    9606 => 'H',
    'C' => 9598,
    9598 => 'C',
    'G' => 9593,
    9593 => 'G',
    'O' => 9600,
    9600 => 'O',
    'M' => 9544,
    9544 => 'M',
    'X' => 9483,
    9483 => 'X'
  };

  foreach my $fg_branch (@$fg_branches) {
    foreach my $node ($copy->nodes) {
      my @leaves = $node->leaves;
      my @leaf_strings = map {
        if ($char_to_taxon->{$_->taxon_id}) {
          $char_to_taxon->{$_->taxon_id};
        }
      } @leaves;
      my $leaf_string = join('',sort {$a cmp $b} @leaf_strings);

      if ($leaf_string eq $fg_branch) {
        $node->name($node->name . '#1');
      }
    }
  }

  return ($copy, $aln);
}

sub _plot_aln_f {
  my $self = shift;
  my $filename = shift;

  my $tmp = $self->worker_temp_directory;

  my $f = $self->_save_file($filename, 'fasta');
  my $pdf_f = $self->_save_file($filename.'_nuc', 'pdf');
  my $tx_pdf_f = $self->_save_file($filename.'_aa', 'pdf');
  my $pdf_file = $pdf_f->{full_file};
  my $tx_pdf_file = $tx_pdf_f->{full_file};

  if (-e $pdf_file && -e $tx_pdf_file && !$self->param('force_recalc')) {
    return;
  }
  if (!-e $f->{full_file}) {
    print "No alignment for file ".$f->{full_file}."!\n";
  }

  my $aln = Bio::EnsEMBL::Compara::AlignUtils->from_file($f->{full_file});
  my $tx_aln = $self->_tx_aln($aln);

  Bio::EnsEMBL::Compara::AlignUtils->to_file($aln, "$tmp/nuc.fasta");  
  Bio::EnsEMBL::Compara::AlignUtils->to_file($tx_aln, "$tmp/aa.fasta");  

  my $rcmd = qq^
library(ggplot2)
library(phylosim)

sim <- PhyloSim()
readAlignment(sim, "$tmp/nuc.fasta")
pdf(file="${pdf_file}", width=30, height=6)
  plotAlignment(sim, axis.text.size=5, aln.plot.chars=F)
dev.off()
readAlignment(sim, "$tmp/aa.fasta")
pdf(file="${tx_pdf_file}", width=30, height=6)
  plotAlignment(sim, axis.text.size=5, aln.plot.chars=T)
dev.off()
^;
  Bio::Greg::EslrUtils->run_r($rcmd);  
}

sub _plot_subs {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  $aln = Bio::EnsEMBL::Compara::AlignUtils->copy_aln($aln);

  my $plot_tree = 1;
  my $plot_alignment = 1;
  my $plot_all_aln_types = 0;
  my $plot_histogram = 0;
  my $plot_scatter = 0;

  if ($plot_all_aln_types) {
    my $use_type = $self->param('aln_use_type');
    $self->param('aln_use_type', 'genomic_all'); # Temporarily cause the default base directory to be used.
    $self->_plot_aln_f('aln_compara');
    $self->_plot_aln_f('aln_genomic_all');
    $self->_plot_aln_f('aln_genomic_primates');
    $self->_plot_aln_f('aln_genomic_mammals');
    $self->param('aln_use_type', $use_type);
  }

  my $plot_aln_f = $self->_save_file('plot_aln', 'pdf');
  my $plot_aln_file = $plot_aln_f->{full_file};
  $self->param('plot_aln_file', $plot_aln_f->{rel_file});  

  my $plot_tree_f = $self->_save_file('plot_tree', 'pdf');
  my $plot_tree_file = $plot_tree_f->{full_file};
  $self->param('plot_tree_file', $plot_tree_f->{rel_file});  

  my $plot_hist_f = $self->_save_file('plot_hist', 'pdf');
  my $plot_xy_f = $self->_save_file('plot_xy', 'pdf');
  my $plot_hist_file = $plot_hist_f->{full_file};
  my $plot_xy_file = $plot_xy_f->{full_file};

  if (!-e $plot_aln_file || $self->param('force_recalc')) {  
    
    my @features;
    my @mutations;
    my @dnds;
    
    my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);

    my $id_to_species;
    map {
      $id_to_species->{$_->name} = $_->ensembl_alias_name;
    } $tree->leaves;

    my $tree_out = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
    $tree_out->translate_ids($id_to_species);
    my $nh_str = $tree_out->root->as_text('newick');
    
    map {
      $_->branch_length(1);
    } $treeI->nodes;

    my $subs_obj = $self->param('subs');
    my @subs = @{$subs_obj};

    my $muts;
    foreach my $s (@subs) {
      my $key = $s->{leaves_beneath};
      $muts->{$key} = [] if (!$muts->{$key});
      push @{$muts->{$key}}, $s;
    }

    my $pos = {
      'Human' => 1,
      'Chimpanzee' => 2,
      'Gorilla' => 3,
      'Orangutan' => 4,
      'Macaque' => 5,
      'Marmoset' => 6,
      'Tarsier' => 7,
      'Mouse Lemur' => 8,
      'Bushbaby' => 9
    };

    my @keys = keys %$muts;
    @keys = sort {
      my @bs = split(',', $b);
      my @as = split(',', $a);
      my $max_b = 0;
      my $max_a = 0;
      foreach my $cur (@as) {
        my $num = $pos->{$id_to_species->{$cur}};
        $max_a = $num if ($num > $max_a);
      }
      foreach my $cur (@bs) {
        my $num = $pos->{$id_to_species->{$cur}};
        $max_b = $num if ($num > $max_b);
      }
      
      return (
        $max_b <=> $max_a ||
        scalar(@bs) <=> scalar(@as)
        );
    } @keys;
    
    foreach my $key (@keys) {
      # Take all mutations for a given branch, split the syn ones in half and
      # 'engulf' the nsyn's w/ the syn's.
      my @cur_branch_muts = @{$muts->{$key}};
      my @syn_rows = grep {$_->{mut_nsyn} == 0} @cur_branch_muts;
      my @nsyn_rows = grep {$_->{mut_nsyn} == 1} @cur_branch_muts;

      my @rows = sort {
        $a->{mut_nsyn} <=> $b->{mut_nsyn} ||
          $a->{aln_pos} <=> $b->{aln_pos}
#          $a->{grantham_score} <=> $b->{grantham_score}
      } @cur_branch_muts;
      
      foreach my $row (@rows) {
        # Using the 'leaves_beneath' string from the table,
        # parse out the taxon IDs of the leaves beneath,
        # and match up to the node from our TreeI object
        my $leaves = $row->{leaves_beneath};
        
        my @ids = split(",", $leaves);
        
        my $org_string = join("/", map {$id_to_species->{$_}} @ids);
        
        my @id_nodes = map {$treeI->find($_)} @ids;
        my $tree_node;
        if (scalar(@id_nodes) > 1) {
          $tree_node = $treeI->lca(@id_nodes);
          if (!defined $tree_node) {
            print "@id_nodes\n";
          }
        } else {
          $tree_node = $id_nodes[0];
          if (!defined $tree_node) {
            print $treeI->ascii;
            print "@id_nodes\n";
            print "@ids\n";
          }
        }
        
        # Create a new internal node representing this mutation event in the tree.
        my $new_node = new $tree_node;
        $new_node->id('m_'.$row->{aln_pos});
        $new_node->sort_order($tree_node->sort_order);
        $new_node->set_tag_value('nsyn', $row->{mut_nsyn});

        my $gs = $row->{grantham_score};
        if ($row->{mut_nsyn} == 0) {
          $gs = 0;
        } else {
          $gs = $gs;
          $gs = 10 if ($gs < 10);
        }
        $row->{grantham_score} = $gs;

        $new_node->set_tag_value('pos', $row->{aln_pos});
        $new_node->set_tag_value('grantham', $row->{grantham_score});
        $new_node->set_tag_value('from', $row->{nuc_from});
        $new_node->set_tag_value('to', $row->{nuc_to});
        
        $tree_node->split_branch_with_node($new_node, 1);
        $new_node->branch_length(1);
        $tree_node->branch_length(0);
        
        my $tax_name = $tree_node->id;
        if ($id_to_species->{$tree_node->id}) {
          $tax_name = $id_to_species->{$tree_node->id};
        }
        
        $self->hash_print($row);

        my $f = Bio::SeqFeature::Generic->new(
          -start => $row->{aln_pos},
          -end => $row->{aln_pos}+1,
          -score => $row->{mut_nsyn},
          -source => $org_string." subs"
          );
        push @features, $f;

        $row->{org_string} = $org_string;
        push @mutations, $row;
      }
    }
    
    map {
      $_->id($id_to_species->{$_->id});
    } $treeI->leaves;
    
    map {
      $_->id($id_to_species->{$_->id});
    } $aln->each_seq;
    
    my $nhx_str = $treeI->root->as_text('nhx');
    my $gene_name = $self->param('gene_name');
    
    my $tmp = $self->worker_temp_directory;
    
    my $tree_subs_file = "$tmp/${gene_name}_subs.nhx";
    my $tree_file = "$tmp/${gene_name}.nh";
    my $aln_file = "$tmp/${gene_name}.fasta";
    my $gff_file = "$tmp/${gene_name}.gff";
    
    open(OUT, ">$tree_subs_file");
    print OUT $nhx_str."\n";
    close(OUT);
    
    Bio::EnsEMBL::Compara::TreeUtils->to_file($treeI, $tree_file);      
    Bio::EnsEMBL::Compara::AlignUtils->to_file($self->_tx_aln($aln), $aln_file);  
        
    # Collect SLR dnds features.
    my $pep_aln = $self->_tx_aln($aln);
    my $results = $self->param('slr_results');
    foreach my $key (1 .. $pep_aln->length) {
      my $obj = $results->{$key};
      next unless (defined $obj->{aln_position});
      next if ($obj->{signed_lrt} eq '-inf');
      my $f = Bio::SeqFeature::Generic->new(
        -start => $obj->{aln_position},
        -end => $obj->{aln_position}+1,
        -score => $obj->{signed_lrt},
        -source => "SLR Mammals"
        );
      push @features, $f;
      $obj->{aln_pos} = $obj->{aln_position};
      push @dnds, $obj;
    }

    $self->output_hash_tsv(\@mutations, "$tmp/muts.txt");
    $self->output_hash_tsv(\@dnds, "$tmp/slr.txt");

    my $gff = Bio::Tools::GFF->new(
      -file => ">$gff_file",
      -gff_version => 3
      );
    foreach my $f (@features) {
      $f->seq_id('chr1');
      $f->strand('+');
      $gff->write_feature($f);
    }
    $gff->close;

    my $gene_name = $self->param('gene_name');
    
    my $rcmd = qq^
library(ggplot2)

# OK, let's do the histograms.
slr <- read.table("$tmp/slr.txt", header=T, sep="\t")
muts <- read.table("$tmp/muts.txt", header=T, sep="\t")
merged <- merge(slr, muts, by=c('aln_pos'))
print(str(merged))

# X-Y scatterplot of dnds vs subst. grantham
if (${plot_scatter} == 1) {
  ns.muts <- subset(merged, mut_nsyn == 1)
  ns.muts[, 'is.gorilla'] <- as.factor(ns.muts[, 'org_string'] == 'Gorilla')
  p <- ggplot(ns.muts, aes(x=signed_lrt, y=grantham_score, colour=is.gorilla, size=is.gorilla, shape=is.gorilla))
  p <- p + theme_bw()
  p <- p + geom_point()
  p <- p + scale_size_manual(values=c(1, 3))
  p <- p + scale_x_continuous("Conservation in Mammals (SLR)")
  p <- p + scale_y_continuous("Grantham score")
  p  <- p + scale_shape(solid=TRUE)
  p <- p + opts(title="Grantham score vs Mammalian Conservation")
  pdf(file="${plot_xy_file}", width=10, height=5)
  print(p)
  dev.off()
}

if (${plot_histogram} == 1) {
  source("~/src/greg-ensembl/scripts/conservation_histogram.R")
  merged[, 'org_string'] <- as.character(merged[, 'org_string'])

  all.sites <- merged[!duplicated(merged[, 'aln_pos']),]
  hgc.subs <- subset(merged, org_string %in% c('Gorilla', 'Chimpanzee', 'Human', 'Human/Chimpanzee', 'Orangutan'))

  #hc.subs <- subset(merged, org_string %in% c('Human/Chimpanzee'))
  #if (nrow(hc.subs) > 0) {
  #  hc.subs[, 'org_string'] <- 'Human'
  #  hgc.subs <- rbind(hgc.subs, hc.subs)
  #  hc.subs[, 'org_string'] <- 'Chimpanzee'
  #  hgc.subs <- rbind(hgc.subs, hc.subs)
  #}

  p <- do.histogram(
    all.sites, 
    hgc.subs, 
    hist.column='signed_lrt', 
    types.column='org_string', 
    types.label="Lineage"
  )
  p <- p + scale_size(name="Grantham score", to=c(0.2, 2))
  p <- p + theme_bw()
  p <- p + scale_x_continuous(name="Conservation in Mammals (SLR statistic)")
  p <- p + opts(
    title=paste("Selection in Mammals for ", "$gene_name")
  )
  pdf(file="${plot_hist_file}", width=10, height=6)
  print(p)
  dev.off()
}

library(phylosim)
source("~/src/greg-ensembl/projects/phylosim/PhyloSimPlots.R")

sim <- PhyloSim()
readAlignment(sim, "${aln_file}")
phylo <- read.nhx("$nh_str")
sim\$.phylo <- phylo
if (${plot_alignment} == 1) {
  pdf(file="${plot_aln_file}", width=24, height=8)
  plotAlignment(sim, aln.gff.file="${gff_file}", axis.text.size=5, aln.plot.chars=T)
  dev.off()
}

if (${plot_tree} == 1) {
  sim <- PhyloSim()
  phylo <- read.nhx("$nhx_str")
  sim\$.phylo <- phylo
  pdf(file="${plot_tree_file}")
  xlab <- "Number of inferred substitutions since LCA"
  obj <- plotTree(sim, tree.do.plot=F, color.by='pos', size.by='grantham', tree.xlab=xlab, tree.xlim.expand=c(0.05, 0.5))
  p <- obj[['grob']]
  p <- p + theme_bw()
  p <- p + scale_colour_gradientn("Alignment Position", colour=rainbow(5))
  #p <- p + scale_colour_manual("Substitution Type",
  #  values=c("gray", "red"),
  #  breaks=c(0, 1),
  #  labels=c('Synonymous', 'Nonsynonymous')
  #)
  p <- p + scale_size_continuous(name="Grantham score", to=c(1, 7))
  #p <- p + scale_size_continuous("Grantham Score")
  p <- p + opts(title="${gene_name}")
  print(p)
  dev.off()
}

^;
    Bio::Greg::EslrUtils->run_r($rcmd);  
  }
}

sub _ref_member {
  my $self = shift;
  my $tree = shift;

  my ($ref_member) = grep {$_->taxon_id == 9606} $tree->leaves;
  my $mba = $self->compara_dba->get_MemberAdaptor;
  $ref_member = $mba->fetch_by_source_stable_id(undef, $ref_member->name);
  $ref_member->name($ref_member->stable_id);

  return $ref_member;
}

sub _collect_stats {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  my $member = $self->_ref_member($tree);

  # IDs
  $self->store_param('tx_id', $member->get_Transcript->stable_id);
  $self->store_param('gene_id', $member->get_Gene->stable_id);
  $self->store_param('protein_id', $member->stable_id);
  $self->store_param('gene_name', $member->get_Gene->external_name);

  # Alignment stats.
  my $pep_aln = $self->_tx_aln($aln);
  $self->store_param('aln_length', $pep_aln->length);

  # GC content and whatnot.
  $self->store_param('gc_cds', $self->gc_content($member));
  $self->store_param('gc_3', $self->gc3_content($member));
  $self->store_param('gc_genomic', $self->genomic_gc_content($member));
  
}

sub _get_sub_patterns_from_aln {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  my $map;
  foreach my $leaf ($tree->leaves) {
    $map->{$leaf->name} = $leaf->taxon_id;
  }
  $aln = Bio::EnsEMBL::Compara::AlignUtils->translate_ids($aln, $map);

  # Substitution patterns.
  my $pattern_hash;
  foreach my $i (0,1) {
    foreach my $j (0,1) {
      foreach my $k (0,1) {
        $pattern_hash->{$i.$j.$k} = 0;
      }
    }
  }

  my @seqs = $aln->each_seq;

  my ($h_seq) = grep {$_->id == 9606} @seqs;
  my ($c_seq) = grep {$_->id == 9598} @seqs;
  my ($g_seq) = grep {$_->id == 9593} @seqs;
  my ($o_seq) = grep {$_->id == 9600} @seqs;
  my $h_str = $h_seq->seq;
  my $c_str = $c_seq->seq;
  my $g_str = $g_seq->seq;
  my $o_str = $o_seq->seq;

  my @col_array;
  my $h = '';
  my $c = '';
  my $g = '';
  my $o = '';
  foreach my $i (1 .. $aln->length) {
    $h = substr($h_str, $i-1, 1);
    $c = substr($c_str, $i-1, 1);
    $g = substr($g_str, $i-1, 1);
    $o = substr($o_str, $i-1, 1);
    my $char = '-';
    next if ($h eq $char || $g eq $char || $c eq $char || $o eq $char);
    $char = 'N';
    next if ($h eq $char || $g eq $char || $c eq $char || $o eq $char);
    $char = '.';
    next if ($h eq $char || $g eq $char || $c eq $char || $o eq $char);
    $char = '';
    next if ($h eq $char || $g eq $char || $c eq $char || $o eq $char);
    
    my $pattern = join('', map {if ($_ eq $o) {0;} else {1;}} ($h, $c, $g));
    $pattern_hash->{$pattern}++;
  }
  foreach my $key (keys %$pattern_hash) {
    $self->store_param('ptrn_'.$key, $pattern_hash->{$key});
  }
}


sub _save_file {
  my $self = shift;
  my $filename_base = shift;
  my $ext = shift;

  my $id = $self->param('gene_name');

  my $filename = "${id}_${filename_base}";

  my $subfolder = 'data';
  if ($self->param('aln_use_type') ne 'genomic_all') {
    $subfolder = 'data/' . $self->param('aln_use_type');
  }

  my $file_params = {
    id => $id,
    filename => $filename,
    extension => $ext,
    subfolder => $subfolder,
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

  if (!defined $self->param('data_prefix')) {
    $self->store_param('data_prefix', $file_obj->{hash_folder});
  }

  return $file_obj;
}

sub fail {
  my $self = shift;
  my $reason_type = shift;
  my $reason_msg = shift;

  $self->fail_and_die($reason_type, $reason_msg);
}

sub genes_table {
  my $self = shift;

  return {
    job_id => 'int',
    unique_keys => 'data_id'
  };
}

sub sites_table {
  my $self = shift;
}

sub thw {
  my $self = shift;
  my $file = shift;

  open(IN, $file);
  my @lines = <IN>;
  close(IN);
  my ($obj) = thaw(join('',@lines));
  return $obj;
}

sub frz {
  my $self = shift;
  my $file = shift;
  my $obj = shift;

  open(OUT,">$file");
  print OUT freeze($obj);
  close(OUT);
}

1;
