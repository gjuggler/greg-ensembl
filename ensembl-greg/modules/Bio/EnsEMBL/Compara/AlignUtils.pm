package Bio::EnsEMBL::Compara::AlignUtils;

use strict;
use Bio::AlignIO;
use Bio::EnsEMBL::Compara::LocalMember;
use Bio::Tools::CodonTable;
use File::Path;
use Cwd;

#
# A grab bag of useful methods for tree manipulations.
#

my $TREE = "Bio::EnsEMBL::Compara::TreeUtils";
my $ALN = "Bio::EnsEMBL::Compara::AlignUtils";
my $COMPARA = "Bio::EnsEMBL::Compara::ComparaUtils";


sub translate_chars {
  my $class = shift;
  my $aln = shift;
  my $map = shift;

  my $aln_copy = $class->copy_aln($aln);

  foreach  my $seq ($aln_copy->each_seq) {
    my $str = $seq->seq;

    foreach my $key (keys %$map) {
      my $val = $map->{$key};
      $str =~ s/$key/$val/g;
    }
    $seq->seq($str);
  }
  
  return $aln_copy;
}

sub translate_ensembl {
  my $class = shift;
  my $aln = shift;

  my $map = {
    "ENSP0.*" => 'Human',
    "ENSPTRP0.*" => 'Chimpanzee',
    "ENSGGOP.*" => 'Gorilla',
    "ENSPPYP0.*" => 'Orangutan',
    "ENSMMUP0.*" => 'Macaque',
    "ENSCJAP0.*" => 'Marmoset',
    "ENSTSYP0.*" => 'Tarsier',
    "ENSMICP0.*" => 'MouseLemur',
    "ENSOGAP0.*" => 'Bushbaby',
  };

  return $class->translate_ids($aln, $map);
}

sub translate_ids {
  my $class = shift;
  my $aln = shift;
  my $map = shift;
  my $params = shift;
  
  my $ensure_unique = 1;
  $ensure_unique = $params->{ensure_unique} if (defined $params->{ensure_unique});

  $aln->set_displayname_flat;
  my $new_aln = $aln->new;

  my $used_ids;

  foreach my $seq ($aln->each_seq) {
    my $id = $seq->id;
    my $nse = $seq->get_nse;

    my $new_id;
    if (defined $id && defined $map->{$id}) {
      $new_id = $map->{$id};
    } else {
      # Treat all map entries as regular expressions.
      foreach my $map_key (keys %$map) {
        if ($id =~ m/$map_key/) {
          $new_id = $map->{$map_key};
        }
      }
    }

    if ($new_id) {
      if ($ensure_unique) {
        while ($used_ids->{$new_id}) {
          $new_id =~ m/_(\d+)$/;
          my $num = $1;
          #print "ID in use: [$new_id]\n";
          $new_id =~ s/_\d+$//;
          my $new_int = $num + 1;
          $new_id .= "_$new_int";
          #print "going to use [$new_id]\n";
        }
      }

      $used_ids->{$new_id} = 1;
      $id = $new_id;
    }

    $new_aln->add_seq(Bio::LocatableSeq->new(-seq => $seq->seq,
					     -id => $id));
  }
  return $new_aln;
}


sub hmmbuild {
  my $class = shift;
  my $aln = shift;
  my $tree = shift;
  my $temp_dir = shift;
  # Uses Hmmer: http://selab.janelia.org/software/hmmer3/3.0/hmmer-3.0.tar.gz

}

sub has_any_data {
  my $class = shift;
  my $aln = shift;

  my $has_data = 0;
  foreach my $seq ($aln->each_seq) {
    $has_data = 1 if ($seq->seq =~ m/[^NX]/gi);
  }
  return $has_data;
}

sub dawg_lambda {
  # Setup steps:
  # 1) Download lambda.py from the DAWG package and put it somewhere in your PATH.
  # http://scit.us/projects/files/dawg/Releases/dawg-release.tar.gz
  my $class = shift;
  my $aln = shift;
  my $tree = shift;
  my $params = shift;
  my $temp_dir = shift;
  if (!defined $temp_dir) {
    $temp_dir = '/tmp/dawg';
    rmtree([$temp_dir]);
    mkpath([$temp_dir]);
  }
  $temp_dir =~ s/\/$//;

  my $aln_out = "$temp_dir/aln.fa";
  my $tree_out = "$temp_dir/tree.nh";
  my $lambda_out = "$temp_dir/lambda.out";

  $class->to_file($aln,$aln_out);
  Bio::EnsEMBL::Compara::TreeUtils->to_file($tree,$tree_out);
  
  my $cmd = "lambda.pl $tree_out $aln_out > $lambda_out";
  my $ret = system("$cmd");

  open(IN,"$lambda_out");
  my @lines = <IN>;
  close(IN);
  
  my $lambda = -1;
  foreach my $line (@lines) {
    if ($line =~ m/Lambda = (.*)/) {
      $lambda = $1;
    }
  }

  #rmtree([$temp_dir]);
  return $lambda;
}


sub indelign{
  # Setup steps:
  # 1) Install GSL from ftp://ftp.gnu.org/gnu/gsl/gsl-1.13.tar.gz
  #    Note: this can be a big pain in the ass, linking the GSL libraries to indelign. Be warned!
  # 2) Install Indelign from http://europa.cs.uiuc.edu/indelign/Indelign.tar.gz
  # 3) Make sure the executable is in your PATH.
  my $class = shift;
  my $sa = shift; # SimpleAlign of DNA/codon sequences.
  my $tree = shift; # Note: Indelign requires a ROOTED tree as input!!!
  my $params = shift;
  my $temp_dir = shift;
  if (!defined $temp_dir) {
    $temp_dir = '/tmp/indelign';
    rmtree([$temp_dir]);
    mkpath([$temp_dir]);
  }

  $temp_dir =~ s/\/$//;

  my $aln_orig = "$temp_dir/aln.fa";
  my $control_file = "$temp_dir/input.txt";
  my $out_aln = "$temp_dir/out_aln.fa";
  my $out_anno = "$temp_dir/out_anno.fa";
  my $ancestor_anno = "$temp_dir/AncAnnotation.txt";
  my $anchor_stat = "/homes/greg/src/Indelign-2.0.3/samples/AnchorStat.txt";
  if ($ENV{'USER'} =~ /gj1/) {
    $anchor_stat = "/nfs/users/nfs_g/gj1/bin/AnchorStat.txt";
  }
  Bio::EnsEMBL::Registry->set_disconnect_when_inactive(1);
  
  my $stdout = "$temp_dir/out.txt";
  my $stderr = "$temp_dir/err";

  $class->to_file($sa,$aln_orig);
  print $tree->newick_format()."\n";
  my $orig_dir = cwd();
  chdir($temp_dir);

  # Call the FindAnchor program.
  my $num_seqs = scalar($sa->each_seq);
  my $ld_path = "/homes/greg/lib";
  my $rc = system("export LD_LIBRARY_PATH=$ld_path; FindAnchor $aln_orig $num_seqs $anchor_stat");
  die "FindAnchor error!" if ($rc);

  # Write the input file.
  # Annoyingly, the parameters have to be in a certain order for Indelign to work (!@$!@$).
  # Use a numerical prefix to sort, and then remove it before writing the control file.
  my $indelign_params = {
    '0_NumSeq' => $num_seqs,
    '1_Tree' => Bio::EnsEMBL::Compara::TreeUtils->to_newick($tree),
    '2_Out1' => $out_aln,
    '3_Out2' => $out_anno,
    '4_LenDist' => 'PL',
    '5_DiffLenDist' => 'false',
    '6_NotRE' => 'false',
#    '7_prA' => 0.25,
#    '8_prC' => 0.25,
#    '9_prG' => 0.25,
#    '10_prT' => 0.25,
  };
  open(OUT,">$control_file");
  foreach my $param (sort {$a <=> $b} keys %{$indelign_params}) {
    my $val = $indelign_params->{$param};
    $param =~ s/(\d+)_//g;
    printf OUT "%s = %s\n",$param,$val;
  }
  close(OUT);

  # Call the Indelign2 program.
  my $cmd = "export LD_LIBRARY_PATH=$ld_path; Indelign2 $aln_orig ${aln_orig}-anchor1 ${aln_orig}-anchor2 $control_file $anchor_stat -A 1>$stdout 2>$stderr";
  print "$cmd\n";
  my $rc = system($cmd);
  die "Indelign error!" if ($rc != 0);

  # TODO: collect the output into site-wise, cumulative insertion and deletion counts.
  my @insertions = (0) x $sa->length;
  my @deletions = (0) x $sa->length;

  my @leaf_lines = read_lines($out_anno);
  my @anc_lines = read_lines($ancestor_anno);
  my @std_lines = read_lines($stderr);

  my $ins_rate = 0;
  my $del_rate = 0;
  foreach my $line (@std_lines) {
    if ($line =~ m/constIns = (\d+(\.\d+)?)/i) {
      $ins_rate = $1;
    }
    if ($line =~ m/constDel = (\d+(\.\d+)?)/i) {
      $del_rate = $1;
    }
  }

  foreach my $line (@leaf_lines,@anc_lines) {
    if ($line =~ m/(\d+)\s+(\d+)\s+(Insertion|Deletion)/) {
      my $start = $1;
      my $stop = $2;
      my $in_del = $3;
      $start -= 1;
      $stop -= 1;
      if ($in_del eq 'Insertion') {
        map {$insertions[$_] += 1} ($start..$stop);
      } else {
        map {$deletions[$_] += 1} ($start..$stop);
      }
    }
  }

  chdir $orig_dir;
  return (\@insertions,\@deletions,$ins_rate,$del_rate);
}

sub read_lines {
  my $file = shift;
  open(IN,$file);
  my @lines = <IN>;
  close(IN);
  return @lines;
}

# Calculate the average column entropy (ACE) of an alignment.
sub column_entropies {
  my $class = shift;
  my $sa = shift; # SimpleAlign object, codon alignment.

# Column entropy = -sum(Pk*log(Pk)), where Pk is the proportion of k in a column
  my $column_entropy_sum;
  my @column_entropy_array;
  foreach my $pos (1 .. $sa->length/3) {
    my @col_array = $class->get_column_array($sa,$pos,3);
#    print join(" ",@col_array)."\n";

    # Mask out gaps and masked codons.
    @col_array = grep {$_ !~ m/(NNN|---)/i} @col_array;

    # Collect the codon proportions into a hash.
    my $n = scalar(@col_array);

    my %codon_hash;
    foreach my $codon (@col_array) {
      if (!defined $codon_hash{$codon}) {
        $codon_hash{$codon} = 1/$n;
      } else {
        $codon_hash{$codon} = $codon_hash{$codon} + 1/$n;
      }
    }
    
    # Add up the entropy for each codon.
    my $column_entropy;
    foreach my $cdn (keys %codon_hash) {
      my $Pk = $codon_hash{$cdn};
      $column_entropy = $column_entropy + ($Pk * log($Pk));
    }
    $column_entropy = -$column_entropy;
    push @column_entropy_array, $column_entropy;

    # Keep a running sum of column-wise CEs.
    #$column_entropy_sum = $column_entropy_sum + $column_entropy;
  }
  # Return the mean CE for this alignment.
  #my $ace = $column_entropy_sum / $sa->length;
  return @column_entropy_array;
}

sub average_column_entropy {
  my $class = shift;
  my $sa = shift;

  # Get the array of column entropies.
  my @ce_array = $class->column_entropies($sa);

  my $sum;
  map {$sum = $sum + $_} @ce_array;
  return $sum / scalar(@ce_array);
}


sub total_column_score {
  my $class = shift;
  my $true_aln = shift;
  my $test_aln = shift;

#  $class->pretty_print($true_aln,{length=>500});
#  $class->pretty_print($test_aln,{length=>500});

  # Index the alignments by sequence residue #.
  my $true_obj = $class->to_arrayrefs($true_aln);
  my $test_obj = $class->to_arrayrefs($test_aln);

  # Alternative method, allowing us to count columns with filtered sites in the test alignment as being correct.
  my $correct_col_count = 0;
  my $true_col_count = $true_aln->length;
  TRUE_COL: foreach my $i (1..$true_aln->length) {
    my $has_match_in_test_aln = 0;

    TEST_COL: foreach my $j (1..$test_aln->length) {      
      my @test_resnums = ();
      my @true_resnums = ();
      my $good_match_count = 0;
      ID: foreach my $seq ($true_aln->each_seq) {
        my $id = $seq->id;
        my $true_resnum = $true_obj->{$id}->[$i-1];
        my $test_resnum = $test_obj->{$id}->[$j-1];

        $good_match_count++ if ($true_resnum eq $test_resnum && $true_resnum ne '-');

        push @test_resnums,$test_resnum;
        push @true_resnums,$true_resnum;

        if ($true_resnum ne $test_resnum && $test_resnum ne 'X') {
          next TEST_COL;
        }
      }
      if ($good_match_count == 0) {
        #next TEST_COL;
      }
      $has_match_in_test_aln = 1;
      last TEST_COL;
    }

    if ($has_match_in_test_aln) {
      $correct_col_count++;
    }
  }
    
  return ($correct_col_count / $true_col_count);
}

sub sum_of_pairs_score {
  my $class = shift;
  my $true_aln = shift;
  my $test_aln = shift;

  # Index the alignments by sequence residue #.
  my $true_obj = $class->to_arrayrefs($true_aln);
  my $test_obj = $class->to_arrayrefs($test_aln);

  # Put all aligned pairs into the hashtable.
  my $true_pairs = $class->store_pairs($true_aln,$true_obj);
  my $test_pairs = $class->store_pairs($test_aln,$test_obj);

  my $true_pair_count = scalar(keys %$true_pairs);
  my $test_pair_count = scalar(keys %$test_pairs);
  
  my $correct_pair_count = 0;
  foreach my $key (keys %$test_pairs) {
    $correct_pair_count++ if ($true_pairs->{$key});
  }
  
  return ($correct_pair_count/$test_pair_count);
}


sub total_aligned_bl {
  my $class = shift;
  my $tree = shift;
  my $ref = shift;
  my $aln = shift;

  my $obj = $class->_correct_subtree_calc($tree, $ref, $aln);
  return $obj->{total_aligned_bl};
}

# Returns the total corectly aligned sites*branchlength
sub total_correct_bl {
  my $class = shift;
  my $tree = shift;
  my $ref = shift;
  my $aln = shift;

  my $obj = $class->_correct_subtree_calc($tree, $ref, $aln);
  return $obj->{total_correct_bl};
}

sub total_incorrect_bl {
  my $class = shift;
  my $tree = shift;
  my $ref = shift;
  my $aln = shift;

  my $obj = $class->_correct_subtree_calc($tree, $ref, $aln);
  return $obj->{total_incorrect_bl};
}

sub correct_subtree_score {
  my $class = shift;
  my $tree = shift;
  my $ref = shift;
  my $aln = shift;

  my $obj = $class->_correct_subtree_calc($tree, $ref, $aln);
  return $obj->{correct_subtree_score};
}

sub incorrect_subtree_score {
  my $class = shift;
  my $tree = shift;
  my $ref = shift;
  my $aln = shift;

  my $obj = $class->_correct_subtree_calc($tree, $ref, $aln);
  return $obj->{incorrect_subtree_score};
}

sub _correct_subtree_calc {
  my $class = shift;
  my $treeI = shift;
  my $ref = shift;
  my $aln = shift;

  my $total_bl = $treeI->total_branch_length;

  my $ref_obj = $class->to_arrayrefs($ref);
  my $aln_obj = $class->to_arrayrefs($aln);
  my $ref_pairs = $class->store_pairs($ref, $ref_obj);
  my $aln_pairs = $class->store_pairs($aln, $aln_obj);

  my @id_list = map { $_->id } $ref->each_seq;
  my $tree_bl_hash;
  my $tree_strings_hash;
  my $scores_hash;
  map { $scores_hash->{$_} = [] } @id_list;

  my $correct_sum;
  my $incorrect_sum;
  my $aligned_sum;
  my $match_sum;
  my $mismatch_sum;
  my $unambiguous_sum;
  my $partial_sum;

  my $correct_sum;
  my $nongap_sum;

  my @correct_bls;
  my @incorrect_bls;
  my @aligned_bls;
  my @match_bls;
  my @mismatch_bls;

  sub get_subtree_residue_string {
    my $aln = shift;
    my $node = shift;
    my $i = shift;

    my @ids = map {$_->id} $node->leaves;
    my @residues = map {$class->get_residue($aln,$_,$i)} @ids;
    my $string = join('',@residues);
    return $string;
  }

  foreach my $i ( 1 .. $aln->length ) {
    my @nongap_ids_at_pos = grep { $aln_obj->{$_}->[$i] !~ m/[-X]/i } @id_list;

    my $aligned_bl = 0;
    my $nongap_node_strings;
    my @node_strings = $class->get_subtree_node_strings($treeI, \@nongap_ids_at_pos, $tree_strings_hash);
    if (scalar(@node_strings) == 1) {
      $aligned_bl = 0;
    } else {
      map {$nongap_node_strings->{$_} = 1} @node_strings;
      foreach my $node ($treeI->nodes) {
        next if ($node->is_leaf);

        my @children = $node->children; # Require bifurcation.
        #next if (scalar(@children) == 1);

        if (scalar(@children) != 2) {
          my $str = "Not two children for node!";
          $str .= "\n";
          $str .= $node->ascii;
          $str .= "\n";
          $str .= join(',', @children)."\n";
          die($str);
        }
        my $a = $children[0];
        my $b = $children[1];
        my $a_str = get_subtree_residue_string($aln, $a, $i);
        my $b_str = get_subtree_residue_string($aln, $b, $i);
        if ($a_str =~ m/^[-X]+$/ || $b_str =~ m/^[-X]+$/) {
          # All gaps or masked residues on one side or the other -- don't add to aligned bl
        } else {
          $aligned_bl += $node->children_branch_length;
        }
      }
    }

    my $seen_aligned_ids;
    #my $seen_misaligned_ids;
    my $correctly_aligned_hash; 
    my $aligned_cluster_count = 0;
    
    #print STDERR get_column_string($aln, $i) . "\n";
#    my $column_score_string = '';
    foreach my $this_seq_id (@id_list) {
      my $this_res_num = $aln_obj->{$this_seq_id}->[$i];

      # Current residue is a gap!
      if ( $this_res_num eq '-' ) {
        $scores_hash->{$this_seq_id}->[$i] = -1;
        #$column_score_string .= '-';
        next;
      }

      foreach my $other_seq_id (@id_list) {
        next if ( $this_seq_id eq $other_seq_id );
        my $other_res_num = $aln_obj->{$other_seq_id}->[$i];
        if ($other_res_num eq '-') {
          next;
        }
        my $pair_string = join( '_', $this_seq_id, $this_res_num, $other_seq_id, $other_res_num );
        my ($seq_a, $seq_b) = sort {$a cmp $b} ($this_seq_id,$other_seq_id);
        if ( defined $ref_pairs->{$pair_string} && $ref_pairs->{$pair_string} == 1 ) {
          # Collect correctly-aligned pairs
          if (!defined $seen_aligned_ids->{$this_seq_id}) {
            $aligned_cluster_count++;
            $seen_aligned_ids->{$this_seq_id} = $aligned_cluster_count;
            $seen_aligned_ids->{$other_seq_id} = $aligned_cluster_count; 
          } elsif (!defined $seen_aligned_ids->{$other_seq_id}) {
            $seen_aligned_ids->{$other_seq_id} = $aligned_cluster_count; 
          }
        } else {
          # Collect incorrectly-aligned pairs
          #$seen_misaligned_ids->{$this_seq_id} = 1;
          #$seen_misaligned_ids->{$other_seq_id} = 1;
        }
      }
    }

    # Use a match / mismatch classification approach.
    # The two branches below a given internal node are assigned to the
    # match state if any aligned pair can be found passing through this
    # node.
    my $match_bl = 0;
    my $mismatch_bl = 0;
    my $unambiguously_match_bl = 0;
    my $partial_match_bl = 0;

    foreach my $node ($treeI->nodes) {
      next if ($node->is_leaf);
      my @children = $node->children;

      next if (scalar(@children) == 1);

      die ("Not two children for node ".$node->id) unless (scalar(@children) == 2);

      my $a = $children[0];
      my $b = $children[1];
      
      my $hash;
      my $found_shared_cluster = 0;
      my $any_a_nongap = 0;
      my $any_b_nongap = 0;
      my $found_any_mismatch = 0;

      my $n_pairs = 0;

      my $pair_hash;
      my $pair_count = 0;
      my $correct_pair_count = 0;
      foreach my $x ($a->leaves) {
        my $x_chr = $aln_obj->{$x->id}->[$i];
        foreach my $y ($b->leaves) {
          my $y_chr = $aln_obj->{$y->id}->[$i];
          next if ($x_chr =~ m/[-X]/i);
          next if ($y_chr =~ m/[-X]/i);
          next if (defined $pair_hash->{$x->id.$y->id});
          $pair_count++;
          my $x_cluster = $seen_aligned_ids->{$x->id};
          my $y_cluster = $seen_aligned_ids->{$y->id};
          if (defined $x_cluster && defined $y_cluster && $x_cluster == $y_cluster) {
            $correct_pair_count++;
          }
          $pair_hash->{$y->id.$x->id} = 1;
          $pair_hash->{$x->id.$y->id} = 1;
        }
      }
      #print $i."  ". $correct_pair_count . " / " . $pair_count . "\n";

      if ($pair_count > 0) {
        my $child_bl = $node->children_branch_length;
        $partial_match_bl += ($correct_pair_count / $pair_count) * $child_bl;

        $unambiguously_match_bl += $child_bl if ($correct_pair_count == $pair_count);
        $mismatch_bl += $child_bl if ($correct_pair_count == 0);
        $match_bl += $child_bl if ($correct_pair_count > 0);
      }
    }

    my $do_extra_tests = 0;
    if ($do_extra_tests) {
      my $score = 0;
      if ($aligned_bl > 0) {
        $score = $match_bl / $aligned_bl;
      }
      
      my $column_string = $class->get_column_string($aln, $i);
      #print $column_string."\n";
      foreach my $k (1 .. $aligned_cluster_count) {
        my $cluster_string = '';
        foreach my $seq ($aln->each_seq) {
          my $res = $class->get_residue($aln, $seq->id, $i);
          if ($res =~ m/[-X]/i) {
            $cluster_string .= ' ';
            next;
          }
          my $char = ' ';
          my $cluster = $seen_aligned_ids->{$seq->id};
          $char = $cluster if (defined $cluster && $cluster == $k);
          $cluster_string .= $char;
        }
        #print STDERR "$cluster_string\n";
      }
    }

    #my $total = $treeI->total_branch_length;
    #print "$aligned_bl $match_bl $mismatch_bl $total\n";

    # Sanity check: we should never have more correct BL than we have nongap BL.
    #die if ($score > 1);

    #push @match_bls, $match_bl;
    #push @mismatch_bls, $mismatch_bl;
    #push @aligned_bls, $aligned_bl;

    $match_sum += $match_bl;
    $mismatch_sum += $mismatch_bl;
    $aligned_sum += $aligned_bl;
    $unambiguous_sum += $unambiguously_match_bl;
    $partial_sum += $partial_match_bl;
  }

  my $obj = {
#    match_branchlengths => \@match_bls,
#    mismatch_branchlengths => \@mismatch_bls,
#    aligned_branchlengths => \@aligned_bls,

    match_bl => $match_sum,
    mismatch_bl => $mismatch_sum,
    aligned_bl => $aligned_sum,
    complete_match_bl => $unambiguous_sum,
    partial_bl => $partial_sum
  };

  return $obj;
}

sub nongap_branch_lengths {
  my $class = shift;
  my $tree = shift;
  my $aln = shift;

  my $obj = $class->to_arrayrefs($aln);
  my @id_list = map { $_->id } $aln->each_seq;

  my $tree_bl_hash;
  my @branchlengths;
  foreach my $i ( 1 .. $aln->length ) {
    my @nongap_ids_at_pos = grep { $obj->{$_}->[$i] ne '-' } @id_list;
    my $nongap_bl = $class->_subtree_bl( $tree, \@nongap_ids_at_pos, $tree_bl_hash );

    push @branchlengths, $nongap_bl;
  }

  return @branchlengths;
}

sub _subtree_bl {
  my $class    = shift;
  my $tree    = shift;
  my $seq_ids = shift;
  my $bl_hash = shift;

  my @id_array = @$seq_ids;
  @id_array = sort {$a cmp $b} @id_array;
  my $key = join('_',@id_array);
  my $existing_value = $bl_hash->{$key};
  return $existing_value if (defined $existing_value);
    
  my $subtree = Bio::EnsEMBL::Compara::TreeUtils->extract_subtree_from_leaves($tree,\@id_array, 1);

  my $total = Bio::EnsEMBL::Compara::TreeUtils->total_distance($subtree);
  if (scalar($subtree->leaves) == 1) {
    $total = 0;
  }
  $bl_hash->{$key} = $total;
  return $total;
}


sub get_residue {
  my $class = shift;
  my $aln = shift;
  my $seq_id = shift;
  my $pos = shift;

  my @seqs = $aln->each_seq;
  my ($seq) = grep {$_->id eq $seq_id} @seqs;
  return $seq->subseq($pos,$pos);
}


sub get_subtree_node_strings {
  my $class = shift;
  my $tree    = shift;
  my $seq_ids = shift;
  my $bl_hash = shift;

  # Optimization: check the branch-length hash, to see if we've already
  # calculated the total branch length for the given list of IDs.
  my @id_array = @$seq_ids;
  @id_array = sort {$a cmp $b} @id_array;
  my $key = join('_',@id_array);
  my $existing_value = $bl_hash->{$key};
  return @{$existing_value} if (defined $existing_value);

  # Get the minimum spanning subtree from a list of IDs
  my $subtree = $tree->slice_by_ids(@id_array);
  #print $subtree->ascii;

  if (scalar($subtree->leaves) == 1) {
    my @leaves = $subtree->leaves;
    my @node_strings = ($leaves[0]->id);
    $bl_hash->{$key} = \@node_strings;
    return @node_strings;
  }

  my @node_strings;
  foreach my $node ($subtree->nodes) {
    push @node_strings, $node->enclosed_leaves_string;
  }
  
  $bl_hash->{$key} = \@node_strings;
  return @node_strings;
}


sub store_pairs {
  my $class = shift;
  my $aln = shift;
  my $obj = shift;

  my %pairs;
  my $pair_string;
  foreach my $i (1..$aln->length) {
    foreach my $key_a (keys %$obj) {
      
      my $resnum_a = $obj->{$key_a}->[$i];
      next if ($resnum_a =~ /[-X]/i); # Ignore pairs with gaps or filtered sites.

      foreach my $key_b (keys %$obj) {
        next if ($key_a eq $key_b);
        my $resnum_b = $obj->{$key_b}->[$i];
        next if ($resnum_b =~ /[-X]/i); # Ignore pairs with gaps for filtered sites.

        $pair_string = join('_',$key_a,$resnum_a,$key_b,$resnum_b);
        #print "$pair_string\n";
        $pairs{$pair_string} = 1;
        $pair_string = join('_',$key_b,$resnum_b,$key_a,$resnum_a);
        $pairs{$pair_string} = 1;
      }
    }
  }

  return \%pairs;
}

sub check_against_possibly_filtered_resnum_list {
  my $class = shift;
  my @resnum_a = shift;
  my @resnum_b = shift;
}

# Encode an alignment as an array of arrayrefs of residue numbers.
sub to_arrayrefs {
  my $class = shift;
  my $aln = shift;

  my $seq_objs;
  foreach my $seq ($aln->each_seq) {
    my @sites = ('-');
    foreach my $i (1..$aln->length) {
      my $range = $seq->location_from_column($i);
      my $aa = $seq->subseq($i,$i);
      if ($aa =~ m/X/i) {
        push @sites, 'X';
      } elsif (!defined $range) {
        push @sites, '-';
      } elsif ($range->location_type eq 'EXACT') {
        push @sites, $range->start
      } else {
        push @sites, '-';
      }
    }
    $seq_objs->{$seq->id} = \@sites;
  }
  return $seq_objs;
}


sub to_aln {
  my $class = shift;
  my $alnF = shift;

  if (ref $alnF) {
    return $alnF;
  } elsif (-e $alnF) {
    return $class->from_file($alnF);
  } else {
    return $class->from_string($alnF);
  }
}

sub from_file {
  my $class = shift;
  my $alnF = shift;
  my $in = Bio::AlignIO->new(-file => $alnF,
			     -format => "fasta");
  my $sa = $in->next_aln();
  $in->close();
  return $sa;
}

# Returns a SimpleAlign from a String.
sub from_string {
  my $class = shift;
  my $alnF = shift;

  open(my $fake_fh, "+<", \$alnF);
  my $in = Bio::AlignIO->new(-fh => $fake_fh,
			  -format => "fasta");
  my $sa = $in->next_aln();
  $in->close();
  return $sa;
}

sub dump_ungapped_seqs {
  my $class = shift;
  my $aln = shift;
  my $alnF = shift;

  open(OUT,">$alnF");
  foreach my $seq ($aln->each_seq) {
    my $seq_str = $seq->seq;
    $seq_str =~ s/-//g;
    print OUT ">".$seq->id."\n";
    print OUT $seq_str."\n";
  }
  close(OUT);
  return $alnF;
}  

sub to_file {
  my $class = shift;
  my $aln = shift;
  my $alnF = shift;

  open(OUT,">$alnF");
  foreach my $seq ($aln->each_seq) {
    print OUT ">".$seq->id."\n";
    print OUT $seq->seq."\n";
  }
  close(OUT);
  return $alnF;
}

# Counts the number of sequences that are neither gapped nor masked at the given column.
# @created GJ 2009-01-09
# @bugfix GJ 2009-03-20
sub get_nongaps_at_column {
  my $class = shift;
  my $aln = shift;
  my $pos = shift;
  
  my $nongap_count = 0;
  my $str = '';
  foreach my $seq ($aln->each_seq) {
    next if ($pos > $seq->length);
    my $residue = $seq->subseq($pos,$pos);
    #print "$pos $residue\n";
    #$str .= $residue;
    $nongap_count++ if ($residue !~ /[-]/i);
  }
  #print Spos." ". $nongap_count." ".$str."\n";

  return $nongap_count;
}

# Counts the number of sequences that are neither gapped nor masked at the given column.
# @created GJ 2009-01-09
# @bugfix GJ 2009-03-20
sub count_residues_at_column {
  my $class = shift;
  my $aln = shift;
  my $pos = shift;
  my $char = shift;

  $char = "-x" unless (defined $char);

  my $char_count = 0;
  foreach my $seq ($aln->each_seq) {
    next if ($pos >= $seq->length);
    my $residue = $seq->subseq($pos,$pos);
    $char_count++ if ($residue =~ m/[$char]/i);
  }
  return $char_count;
}

sub count_residues {
  my $class = shift;
  my $aln = shift;
  my $cdna_option = shift;

  my $n = 0;
  foreach my $seq ($aln->each_seq) {
    my $str = $seq->seq;
    if ($cdna_option) {
      $str =~ s/[-n]//gi;
    } else {
      $str =~ s/[-x]//gi
    }
    
    $n += length($str);
  }
  return $n;
}

sub get_ungapped_branchlength {
  my $class = shift;
  my $aln = shift;
  my $in_tree = shift;
  my $pos = shift;

  my $tree = $in_tree->copy;

  #print " -> $pos\n";
  my @keep_node_ids = ();
  foreach my $seq ($aln->each_seq) {
    #print " ID: ". $seq->id."\n";
    my $residue = $seq->subseq($pos,$pos);
    if ($residue !~ /[-x]/i) {
      my $leaf = $tree->find_leaf_by_name($seq->id);
      push @keep_node_ids, $leaf->node_id;
    }
  }

  my $subtree = Bio::EnsEMBL::Compara::TreeUtils->extract_subtree_from_leaves($tree,\@keep_node_ids);
  my $total = Bio::EnsEMBL::Compara::TreeUtils->total_distance($subtree);
  #print " -> $total\n";
  return $total;
}

sub get_column_array {
  my $class = shift;
  my $aln = shift;
  my $pos = shift;
  my $codon_width = shift;

  $codon_width = 1 unless (defined $codon_width);

  my $lo_i = ($pos-1)*$codon_width + 1;
  my $hi_i = ($pos)*$codon_width;

  my @chars = ();
  foreach my $seq ($aln->each_seq) {
    my $residue = $seq->subseq($lo_i,$hi_i);
    push @chars,$residue;
  }
  return @chars;
}

sub get_column_string {
  my $class = shift;
  my $aln = shift;
  my $pos = shift;

  my $str = "";
  foreach my $seq ($aln->each_seq) {
    my $residue = $seq->subseq($pos,$pos);
    $str .= $residue;
  }
  return $str;
}


# Retrieves a sequence from a Bio::AlignI object based on the given ID.
sub get_seq_with_id {
  my $class = shift;
  my $aln = shift;
  my $id = shift;
  
  foreach my $seq ($aln->each_seq) {
    return $seq if ($seq->id eq $id);
  }
  
  #warn("aln_get_seq_with_id: No sequence with id $id found!\n");
  return undef;
}

sub translate {
  my $class = shift;
  my $aln = shift;
  my $params = shift;

  my $codon_table_id = 1;
  $codon_table_id = $params->{bioperl_codontable_id} if (defined $params->{bioperl_codontable_id});

  my $sa = $aln->new;
  foreach my $seq ($aln->each_seq) {
    my $tx = $seq->translate(-codontable_id => $codon_table_id);
    die "Translation for ".$seq->id." contains stop codon!\n" if ($tx =~ m/\*/);
    my $tmp = $tx->seq;
    $tmp =~ s/-//g;
    $tx->end(length($tmp));
    $sa->add_seq($tx);
  }
  return $sa;
}


# Generates an Ensembl cigar line (MD-line) for each sequence in the alignment.
# Returns a hashref of sequences indexed by IDs.
sub cigar_lines {
  my $class = shift;
  my $aln = shift;
  
  my $lines;
  foreach my $pos (1..$aln->no_sequences) {   # $aln->no_sequences... what an awful method name!
    my $seq = $aln->get_seq_by_pos($pos);
    $lines->{$seq->display_id} = $COMPARA->cigar_line($seq->seq);
  }
  return $lines;
}

sub filter_stop_codons {
  my $class = shift;
  my $aln = shift;
  my $params = shift;

  my $new_aln = $aln->new;
  foreach my $seq ($aln->each_seq) {
    $new_aln->add_seq($class->_filter_seq_stops($seq, $params));
  }

  return $new_aln;
}

sub filter_frameshifting_indels {
  my $class = shift;
  my $aln = shift;
  my $ref_seq_id = shift;

  my $ref_seq = $class->get_seq_with_id($aln,$ref_seq_id);
  my $str = $ref_seq->seq;
  my $str_nogaps = $str;
  $str_nogaps =~ s/-//g;
  
  # 1) Remove columns where the reference sequence has a gap w/ length != multiple-of-three
  my @remove_me;

  $_ = $str;
  while (m/(-+)/g) {
    my $match = $1;
    my $len = length($match);

#    next if ($len % 3 == 0); # Multiples of 3 are OK!

    my $end = pos($_);
    my $start = $end - $len;
    my $match_substr = substr($_, $start-1, $len+2);
    #print "$match $match_substr\n";

    # TODO: do some better test for "good" indels here...

    # OK, so we want to remove some columns...
    foreach my $i ($start .. ($end-1)) {
      push @remove_me, $i+1;
    }
  }

  $aln = $class->remove_columns($aln, \@remove_me);

  # 2) Go through the columns by threes, and turn to gaps any non-complete codons (i.e. gap
  #    in one or two positions) or stop codons into all gaps.
  foreach my $seq ($aln->each_seq) {
    my $id = $seq->id;
    my $str = $seq->seq;
    for (my $i=1; $i <= $aln->length - 2; $i+= 3) {
      my $codon = substr($str, $i-1, 3);
      if (($codon =~ m/(tag|tga|taa)/gi) or
          ($codon ne '---' && $codon =~ m/-/g)) {
        #print "$id $i $codon\n";
        substr($str, $i-1, 3) = '---';
      }
    }
    $seq->seq($str);
  }

  return $aln;
}

sub ensure_multiple_of_three {
  my $class = shift;
  my $aln = shift;

  if ($aln->length % 3 != 0) {
    my $diff = $aln->length % 3;
    print "Alignment not the right length! diff[$diff]\n";
    $aln = $aln->slice(1,$aln->length - $diff);
  }
  return $aln;
}

sub _filter_seq_stops {
  my $class = shift;
  my $seq = shift;
  my $params = shift;

  my $be_less_stringent = $params->{be_less_stringent};
  $be_less_stringent = 0 unless (defined $be_less_stringent);

  my $codon_table_id = 1;
  $codon_table_id = $params->{bioperl_codontable_id} if (defined $params->{bioperl_codontable_id});

  for (my $i=1; $i <= $seq->length-2; $i+= 3) {
    my $seq_str = $seq->seq;

    my $codon_str = $seq->subseq($i, $i+2);
    if ($codon_str eq '---') {
      next;
    }

    my $has_gaps = 0;
    $has_gaps = 1 if ($codon_str =~ m/-/i);

    my $aa = new Bio::PrimarySeq(-seq => $seq->subseq($i,$i+2),
      -codontable_id => $codon_table_id)->translate->seq;
    
    #substr($seq_str,$i-1,3) = '---' if (substr($seq_str,$i-1,3) =~ m/(tag|tga|taa)/gi);
    #print "$aa\n";
    my $is_stop;
    if ($be_less_stringent) {
      $is_stop = 1 if ($be_less_stringent && $codon_str =~ m/(tag|tga|taa)/gi); # Only filter out full-on stop codons.
    } else {
      $is_stop = 1 if ($aa =~ m/([x\*])/gi); # More stringent -- anything that doesn't translate correctly.
    }
    #print " Masking stop $1 at $i\n" if ($is_stop);

    if ($is_stop && !$has_gaps) {
      substr($seq_str,$i-1,3) = 'NNN';
    } elsif ($is_stop && $has_gaps) {
      substr($seq_str,$i-1,3) = '---';
    }
    #print "  Masking stop $1 at $i\n"  if (substr($seq_str,$i-1,3) =~ m/(tag|tga|taa)/gi);
    $seq->seq($seq_str);
  }

  return $seq;
}

sub seq_index {
  my $class = shift;
  my $aln = shift;
  my $id = shift;

  my $i=0;
  foreach my $seq ($aln->each_seq) {
    if ($seq->id eq $id) {
      return $i;
    }
    $i++;
  }
  return -1;
}

sub get_nongap_indices {
  my $class = shift;
  my $aln = shift;
  my $id = shift;

  my $ref_seq = $class->get_seq_with_id($aln, $id);
  warn("Seq [$id] not found while flattening alignment!") unless (defined $ref_seq);

  my $seq = $ref_seq->seq;
  my @nongap_columns;
  while ($seq =~ m/[^-]/g) {
    push @nongap_columns, (pos($seq));
  }
  return @nongap_columns;
}
  
=head2 flatten_to_sequence
 Title     : flatten_to_sequence
 Usage     : $aln->flatten_to_sequence($my_favorite_sequence_index)
 Function  : Flattens the alignment (removes all insertions) relative to the sequence at the given index.
 Returns   : a new Bio::SimpleAlign object, flattened appropriately.
 Argument  : The index of the sequence against which to flatten the alignment. (Use get_pos_by_id to find a sequence index by name)
=cut
  
sub flatten_to_sequence {
  my $class = shift;
  my $aln = shift;
  my $id = shift;
  
  my $ref_seq = $aln->get_seq_by_id($id);
  warn("Seq [$id] not found while flattening alignment!") unless (defined $ref_seq);
  my $display_id = $ref_seq->display_id;
  my $seq_str = $ref_seq->seq;
    
  # Go through the reference sequence and remove gapped columns in the new alignment.
  my $gap_char = $aln->gap_char();
  
  #FIXME: This is less efficient than it would be if we grouped the deletions into gaps. Sue me.
  my @remove_cols;
  while ($seq_str =~ m/[$gap_char]/g) {
    #print pos($seq_str)."\n";
    my @start_end = (pos($seq_str)-1,pos($seq_str)-1);
    push @remove_cols, \@start_end;
  }
  my $new_aln = $aln->_remove_columns_by_num(\@remove_cols);
  return $new_aln;
}

# Remove columns included in the arrayref. The first alignment column is 1.
sub remove_columns {
  my $class = shift;
  my $aln = shift;
  my $column_number_arrayref = shift;

  my @column_numbers = sort @{$column_number_arrayref};

  my @remove_cols;
  foreach my $column (@column_numbers) {
    my @start_end = ($column-1, $column-1);
    push @remove_cols, \@start_end;
  }
  my $new_aln = $aln->_remove_columns_by_num(\@remove_cols);
  return $new_aln;
}

sub combine_alns {
  my $class = shift;
  my @alns = @_;

  my $seq_hash;

  my $genomic_coords = {};

  foreach my $aln (@alns) {
    # Collect all species from all alignments into a hash.
    map {$seq_hash->{$_->id} = [] if (!defined $seq_hash->{$_->id});} $aln->each_seq;

    # Copy over the new genomic coords.
    $genomic_coords = Bio::EnsEMBL::Compara::ComparaUtils->replace_params( $genomic_coords, $aln->annotation->{_genomic_coords} );
  }

  foreach my $aln (@alns) {
    # Fill in each sequence in the hash with the current alignment, or gaps if that species is currently missing.
    foreach my $id (keys %$seq_hash) {
      my @chars = @{$seq_hash->{$id}};
      my ($seq) = grep {$_->id eq $id} $aln->each_seq;
      #my $seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id($aln,$id);
      if (defined $seq) {
	push @chars, split(//,$seq->seq);
      } else {
	#print "Missing seq for $id from alignment!\n";
	push @chars, split(//,'-' x $aln->length);
      }
      $seq_hash->{$id} = \@chars;
    }
  }

  # Make a final alignment.
  my $aln = new Bio::SimpleAlign;
  foreach my $id (keys %$seq_hash) {
    my $seq = join('',@{$seq_hash->{$id}});
    $aln->add_seq(Bio::LocatableSeq->new(-seq => $seq,
					 -id => $id));
  }

  $aln->annotation->{_genomic_coords} = $genomic_coords;
  return $aln;
}

sub contains_sequence {
  my $class = shift;
  my $aln = shift;
  my $seq_str = shift;
  my $params = shift;

  print $params->{bioperl_codontable_id}."\n";

  my $codontable_id = 1;
  $codontable_id = $params->{bioperl_codontable_id} if (defined $params->{bioperl_codontable_id});

  my $ref_aa = new Bio::PrimarySeq(-seq=>$seq_str)->translate(-codontable_id => $codontable_id)->seq;
  $ref_aa =~ s/[\*X]//g;
  
  print "REF AA: $ref_aa\n";

  foreach my $seq ($aln->each_seq) {
    my $compare_str = $seq->seq;
    my $nogaps = $compare_str;
    $nogaps =~ s/-//g;
    my $nogaps_aa = new Bio::PrimarySeq(-seq=>$nogaps)->translate( -codontable_id => $codontable_id)->seq;
    $nogaps_aa =~ s/[\*X]//g;

    print "NOGAPS: $nogaps_aa\n";

    #print $nogaps_aa." ".$seq->id."\n";
    #print $ref_aa."\n";
    
    return 1 if ($nogaps_aa eq $ref_aa);

    # Try taking one aa off of the nogaps sequence.
    my $nogaps_shorter = substr($nogaps_aa, 0, length($nogaps_aa)-1);
    return 1 if ($nogaps_shorter eq $ref_aa);
  }
  return 0;
}


# GJ 2009-03-17
# Returns a map with key = (new column #) and val = (old column #)
sub remove_blank_columns {
  my $class = shift;
  my $aln = shift;

  my $new_to_old = {};
  my $old_to_new = {};
  my @cols_to_remove = ();
  my $j = 1;
  for (my $i=1; $i <= $aln->length; $i++) {
    my $col = $class->get_column_string($aln,$i);
    if ($col =~ m/^[-]+$/) { # Regex for a blank line.
      # I guess the remove_columns_by_num method is zero-based? GJ 2009-03-19
      push @cols_to_remove, [$i-1,$i-1];
      $old_to_new->{$i} = -1;
    } else {
      $new_to_old->{$j} = $i;
      $old_to_new->{$i} = $j;
      $j++;
    }
  }
  $aln = $aln->_remove_columns_by_num(\@cols_to_remove);
  $aln->{_old_to_new} = $old_to_new;
  $aln->{_new_to_old} = $new_to_old;
  return $aln;
}

# Removes blank columns in threes, and returns a codon-wise column mapping.
sub remove_blank_columns_in_threes {
  my $class = shift;
  my $aln = shift;

  my @cols_to_remove = ();
  for (my $i=1; $i <= $aln->length-2; $i+= 3) {
    my $col = $class->get_column_string($aln,$i);
    my $col2 = $class->get_column_string($aln,$i+1);
    my $col3 = $class->get_column_string($aln,$i+2);
    if ($col =~ m/^[-]+$/ && $col2 =~ m/^[-]+$/ && $col3 =~ m/^[-]+$/ ) { # Regex for a blank line.
      push @cols_to_remove, [$i-1,$i-1+2];
    }
  }
  $aln = $aln->_remove_columns_by_num(\@cols_to_remove);
  return $aln;
}

# Removes blank columns in threes, and returns a codon-wise column mapping.
sub remove_gappy_columns_in_threes {
  my $class = shift;
  my $aln = shift;

  my @cols_to_remove = ();
  my $n = 0;
  for (my $i=1; $i <= $aln->length-2; $i+= 3) {
    my $col = $class->get_column_string($aln,$i);
    my $col2 = $class->get_column_string($aln,$i+1);
    my $col3 = $class->get_column_string($aln,$i+2);
    if ($col =~ m/[-]/ || $col2 =~ m/[-]/ || $col3 =~ m/[-]/ ) { # Regex for any gaps in the codon..
      push @cols_to_remove, [$i-1,$i-1+2];
      $n += 3;
    }
  }
  print "Removing $n gappy columns...\n";
  
  $aln = $aln->_remove_columns_by_num(\@cols_to_remove);
  return $aln;
}

sub remove_regex_columns_in_threes {
  my $class = shift;
  my $aln = shift;
  my $regex = shift;

  my $cols_hash;
  my $n = 0;
  foreach my $seq ($aln->each_seq) {
    my $seq_str = $seq->seq;

    print $seq->id."\n";
    my @codons = unpack('(A3)*', $seq_str);
    print scalar(@codons)."\n";
    my $i=0;
    foreach my $codon (@codons) {
      if ($codon =~ m/$regex/) {
        $cols_hash->{$i*3} = 1;
        $cols_hash->{$i*3+1} = 1;
        $cols_hash->{$i*3+2} = 1;
      }
      $i += 1;
    }
  }

  my @cols_to_remove = keys %$cols_hash;
  my $n = scalar(@cols_to_remove);
  print "Removing $n regex columns...\n";

  $aln = $class->quick_remove_columns($aln, \@cols_to_remove);
  return $aln;
}

sub quick_remove_columns {
  my $class = shift;
  my $aln = shift;
  my $cols_arrayref = shift;

  my @cols = @$cols_arrayref;
  @cols = sort {$b <=> $a} @cols;

  my $keep_hash;
  map {$keep_hash->{$_-1} = 1} 1 .. $aln->length;

  foreach my $col (@cols) {
    delete $keep_hash->{$col};
  }
  
  my @keepers = sort {$a <=> $b} keys %$keep_hash;

  my $new_aln = new $aln;
  foreach my $seq ($aln->each_seq) {
    my $str = $seq->seq;
    my @arr = unpack("(A1)*", $str);
    @arr = @arr[@keepers];

    my $new_seq = new $seq;
    $new_seq->id($seq->id);
    $new_seq->seq(join('', @arr));
    $new_aln->add_seq($new_seq);
  }
  return $new_aln;
}

sub remove_funky_stretches {
  my $class = shift;
  my $aln = shift;
  my $id1 = shift;
  my $id2 = shift;
  my $max_allowable_run_of_substitutions = shift;

  my $seq1 = $class->get_seq_with_id($aln,$id1);
  my $seq2 = $class->get_seq_with_id($aln,$id2);

  return $aln if (!defined $seq1 || !defined $seq2);

  # Add a space to the seq in the beginning to make it one-based like Bio::Seq.
  my $str1 = $seq1->seq;
  my $str2 = $seq2->seq;

  my $current_substitution_run = 0;
  my $current_substitution_run_codon_start = -1;
  my @cols_to_remove;
  my $n = 0;
  for (my $i=1; $i <= $aln->length-2; $i+= 3) {
    for (my $j=$i; $j <= $i+2; $j++) {
      my $c1 = substr($str1,$j-1,1);
      my $c2 = substr($str2,$j-1,1);

      if ($c1 ne $c2 && $c1 !~ m/[-xn]/gi && $c2 !~ m/[-xn]/gi) {
	$current_substitution_run_codon_start = $i if ($current_substitution_run_codon_start == -1);
	$current_substitution_run++;
      } else {
	if ($current_substitution_run >= $max_allowable_run_of_substitutions) {
	  # Stretch the lower and upper bounds to the nearest codon boundary.
	  my $lo = $current_substitution_run_codon_start;
	  my $hi = $j - 1;
	  #print "$lo -> $hi\n";
	  my $overhang = ($hi-$lo+1) % 3;
	  $hi += 3 - $overhang if ($overhang);
	  my $len = $hi - $lo + 1;
	  
	  push @cols_to_remove, [$lo-1,$hi-1];
	  #print "$lo - $hi\n";
	  $n += $len;
	  $current_substitution_run_codon_start = -1;
	  $current_substitution_run = 0;
        } else {
	  # Reset the counters to zero.
  	  $current_substitution_run_codon_start = -1;
  	  $current_substitution_run = 0;
        }
      }
    }
  }  

  print "[$id1 - $id2] No clustered muts!\n" if ($n == 0);
  print "[$id1 - $id2] Removing $n clustered muts...\n" if ($n > 0);
  $aln = $aln->_remove_columns_by_num(\@cols_to_remove);
  return $aln;
}


sub mask_high_mutation_windows {
  my $class = shift;
  my $aln = shift;
  my $window_size = shift;
  my $window_step = shift;
  my $max_mutations_in_window = shift;

  my $num_sites_masked = 0;

  for (my $i=1; $i <= $aln->length-($window_size-1); $i+= $window_step) {
    my $lo = $i;
    my $hi = $i+$window_size;
    $hi = $aln->length if ($hi > $aln->length);

    # Count up differences between seq_a and all other seqs.
    foreach my $seq_a ($aln->each_seq) {
      my $diff_count = 0;

      foreach my $seq_b ($aln->each_seq) {
        next if ($seq_b == $seq_a);

        my $n_diffs = $class->aln_string_diff($seq_a->subseq($lo,$hi),$seq_b->subseq($lo,$hi));
        $diff_count += $n_diffs;
      }

      # If the differences are too many, mask out all of $seq_a in this window.
      if ($diff_count / scalar($aln->each_seq) > $max_mutations_in_window) {
        printf "(%d %d) %s N_DIFFS: %.3f\n",$lo,$hi,$seq_a->id,$diff_count / scalar($aln->each_seq);
        my $orig_seq = $seq_a->seq;
        my ($seq,$n_masked) = $class->mask_string($orig_seq,$lo-1,$hi-1);
        $num_sites_masked += $n_masked;
        $seq_a->seq($seq);
      }
    }
  }
  return $num_sites_masked;
}

sub mask_string {
  my $class = shift;
  my $string = shift;
  my $lo = shift;
  my $hi = shift;

  my $n = 0;
  for (my $i=$lo; $i <= $hi; $i++) {
    if (substr($string,$i,1) ne '-') {
      substr($string,$i,1,'N');
      $n++;
    }
  }
  return ($string,$n);
}

sub aln_string_diff {
  my $class = shift;
  my $string1 = shift;
  my $string2 = shift;

  my $n_diffs = 0;
  for (my $i=0; $i < length($string1); $i++) {
    my $char1 = substr($string1,$i,1);
    my $char2 = substr($string2,$i,1);
    
    if ($char1 ne '-' && $char2 ne '-' && $char1 ne 'N' && $char2 ne 'N' && $char1 ne $char2) {
      $n_diffs++;
    }
  }
  return $n_diffs;
}


# Removes codons where all three nucleotide positions differ between two species.
sub remove_triple_mutated_codons {
  my $class = shift;
  my $aln = shift;
  my $id1 = shift;
  my $id2 = shift;

  my $seq1 = $class->get_seq_with_id($aln,$id1);
  my $seq2 = $class->get_seq_with_id($aln,$id2);

  my @removed_column_strings = ();
  my @cols_to_remove = ();
  my $n=0;
  for (my $i=1; $i <= $aln->length-2; $i+= 3) {
    my $bad_pos_count = 0;
    for (my $j=$i; $j <= $i+2; $j++) {
      my $char1 = $seq1->subseq($j,$j);
      my $char2 = $seq2->subseq($j,$j);
      if ($char1 ne $char2 && $char1 ne '-' && $char2 ne '-') {
	$bad_pos_count++;
      }
    }
    if ($bad_pos_count >= 2) {
      push @cols_to_remove, [$i-1,$i-1+2];
      push @removed_column_strings, $class->get_column_string($aln,$i);
      push @removed_column_strings, $class->get_column_string($aln,$i+1);
      push @removed_column_strings, $class->get_column_string($aln,$i+2);
      push @removed_column_strings, "";
      $n++;
    }
  }
  print "Removing $n funky codons...\n";
  #print join("\n",@removed_column_strings);
  $aln = $aln->_remove_columns_by_num(\@cols_to_remove);
  return $aln;
}


sub map_alignment_position {
  my $self = shift;
  my $aln_a = shift;
  my $column_in_a = shift;
  my $aln_b = shift;

  # Assumption: alignment A has at least one non-gap sequence at the given position.

  foreach my $seq ($aln_a->each_seq) {
    if ($seq->subseq($column_in_a,$column_in_a) =~ m/[-]/) {
      next;
    } else {
      my $location = $seq->location_from_column($column_in_a);
      if (defined $location && $location->location_type() eq 'EXACT') {
	my $seq_position = $location->start;

	# So we found the sequence position. Now, look for that residue in alignment B.
	my $seq_in_b = $aln_b->get_seq_by_id($seq->id);
	if (!defined $seq_in_b) {
	  next;
        } else {
	  my $column_in_b = $seq_in_b->column_from_residue_number($seq_position);
	  return $column_in_b;
        }
      }
    }
  }
  return undef;
}

=head2 sort_by_tree
  
  Title     : sort_by_tree
  Usage     : Bio::EnsEMBL::Compara::AlnUtils->sort_by_tree($aln,$tree);
  Function  : Given a Bio::AlignI object and a Bio::Tree::TreeI object, creates a new AlignI sorted by the tree.
  Returns   : a new Bio::SimpleAlign object, sorted appropriately.
  Argument  : ($aln,$tree) - a Bio::AlignI object, and a Bio::TreeI or Bio::EnsEMBL::Compara::ProteinTree object

=cut
sub sort_by_tree {
  my $class = shift;
  my $aln = shift;
  my $tree = shift;
  
  if (! ref $tree) {
    # Check that the file exists, and try to load it.
    my $infile = $tree;
    $aln->throw("Tree file not found: $infile") unless (-e $infile);
    my $temp = Bio::TreeIO->new('-file' => $infile,
				'-format' => 'newick');
    $tree = $temp->next_tree;
    $temp->close();
  } else {
    # If we're given an object hashref, we assume it's a TreeI object.
    if ($tree->isa('Bio::EnsEMBL::Compara::NestedSet')) {
      $tree = $TREE->to_treeI($tree);
    } elsif (! $tree->isa('Bio::Tree::TreeI')) {
      # Try converting ProteinTree to TreeI.
      warn "sort_by_tree: given a hashref that is not a TreeI object!\n";
    }
  }
  
  my $new_aln = $aln->new;
  foreach my $leaf ($tree->get_leaf_nodes) {
    my $name = $leaf->id;
    #print "$name\n";
    my $seq = $class->get_seq_with_id($aln,$name);
    warn("Seq with id [$name] not found in aln!\n") unless (defined $seq);
    $new_aln->add_seq($seq);
  }
  return $new_aln;
}

# Takes zero-based column coordinates to keep.
sub mask_columns {
  my $class = shift;
  my $aln = shift;
  my $columns_ref = shift;
  my $mask_char = shift;

  my %keeper_sites;
  map {$keeper_sites{int($_)+1}=1} @{$columns_ref};

  foreach my $i (1 .. $aln->length) {
    if (!exists $keeper_sites{$i}) {
      foreach my $seq ($aln->each_seq) {
	my $seq_str = $seq->seq;
	substr($seq_str,$i-1,1) = $mask_char unless (substr($seq_str,$i-1,1) eq '-');
        $seq->seq($seq_str);
      }
    }
  }
  return $aln;
}

sub prank_filter {
  my $class = shift;
  my $aln = shift;
  my $tree = shift;
  my $params = shift;

  my $threshold = $params->{'prank_filtering_threshold'} || 5;
  my $char = $params->{'prank_mask_character'} || 'X';

  my ($scores,$blocks) = $class->get_prank_filter_matrices($aln,$tree,$params);
  my $filtered_aln = $class->mask_below_score($aln,$threshold,$scores,$char);

  return $filtered_aln;
}


sub filter_sites {
  my $class = shift;
  my $aln = shift;
  my $seq_name = shift;
  my $sites_ref = shift;
  my $mask_char = shift;

  my @sites = keys %{$sites_ref};

  foreach my $seq ($aln->each_seq) {
    my $label = $seq->display_id;
    if ($label eq $seq_name) {
      # We're at the right sequences now. Filter out the sites!
      my $seq_str = $seq->seq;
      foreach my $site (@sites) {
        substr($seq_str,$site-1,1) = $mask_char unless (substr($seq_str,$site-1,1) eq '-');
        $seq->seq($seq_str);
      }
      last;
    }
  }
  return $aln;
}

sub gblocks_filter {
  my $class = shift;
  my $aln = shift;
  my $tree = shift;
  my $params = shift;

  my $gblocks_char = $params->{'gblocks_mask_character'} || "U";

  my $node_id = $tree->node_id;
  my $dbc = $tree->adaptor->db->dbc;
  my $sth = $dbc->prepare("select aln_start,aln_end FROM gblocks_aln WHERE node_id=$node_id;");
  $sth->execute();

  my %keeper_sites;
  while (my @row = $sth->fetchrow_array) {
    my $start = $row[0];
    my $end = $row[1];

    foreach my $i ($start .. $end) {
      $keeper_sites{$i} = 1;
    }
  }

  foreach my $i (1 .. $aln->length) {
    if (!defined $keeper_sites{$i}) {
      # Mask out this column with "U"s
      foreach my $seq ($aln->each_seq) {
        my $seq_str = $seq->seq;
        substr($seq_str,$i-1,1) = $gblocks_char unless (substr($seq_str,$i-1,1) eq '-');
        $seq->seq($seq_str);
      }
    }
  }
  return $aln;
}

=head2 mask_below_score
  Title     : mask_below_score
  Usage     : $aln->mask_below_score($threshold, $score_hashref)
  Function  : Masks out nucleotides (with N) or amino acids (with X) in the alignment that are LESS THAN
             the given threshold score. So, with a threshold of 3, residues with score 3 are kept, but score 2 is turned to X or N.
             This was mainly created to allow thresholding based on T_Coffee's score_ascii output.
 Returns   : a new Bio::SimpleAlign object with the masking applied. Xs are applied to protein alphabets, and Ns to nucleotides.
 Argument  : $score_hashref should be a hash_ref with keys as the sequence labels, and values
             as alignment-shaped strings with quality scores instead of residues.
=cut

sub mask_below_score {
  my $class = shift;
  my $aln = shift;
  my $threshold = shift;
  my $score_hashref = shift;
  my $mask_char = shift;

  if ($threshold > 9) {
    warn "WARNING: Bio::SimpleAlign::mask_below_score: threshold $threshold capped to the maximum of 9.\n";
    $threshold = 9;
  }
  
  #  Alignment string
  #    ABC---DEF
  #  Score string
  #    995---339

  my $alphabet = 'protein'; # Default to a protein alphabet.
  # Create a new shell SimpleAlign object.
  my $new_aln = $aln->new;

  # Go through each sequence in the alignment and mask accordingly.
  foreach my $seq ( $aln->each_seq ) {
    $alphabet = $seq->alphabet;
    my $label = $seq->display_id;
    my $score_string = $score_hashref->{$label};

    die("No score string found for sequence [$label]!") unless (defined $score_string);

    $score_string =~ s/[^\d-]/9/g;   # Convert non-digits and non-dashes into 9s.
    # (The above is necessary because t_coffee leaves leftover letters in the sequence score strings)
    $score_string =~ s/[-]/9/g;   # Convert dashes to 9s, because we will never mask those out.
    #print $score_string."\n";
    my @score_array = split(//,$score_string);
    
    my @mask_array = map {$_ >= int($threshold) ? 1 : 0 } @score_array;  # Map the score array to an array of 0s or 1s
    #print join "", @mask_array,"\n";
    my @seq_array = split(//,$seq->seq);
    foreach my $i (0..scalar(@seq_array)-1) {
      my $char = $seq_array[$i];
      if ($mask_array[$i] == 0 && $alphabet eq 'protein') {
        next if ($char eq 'X' || $char eq 'x');
        if ($char =~ /[a-z]/) {
          $seq_array[$i] = lc($mask_char);
        } else {
          $seq_array[$i] = uc($mask_char);
        }
      } elsif ($mask_array[$i] == 0) {
        if ($char =~ /[a-z]/) {
          $seq_array[$i] = lc($mask_char);
        } else {
          $seq_array[$i] = uc($mask_char);
        }
      }
    }
    my $new_str = join "", @seq_array;

    # We've got the new alignment sequence, just put it into the new SimpleAlign.
    my $tmp = $new_str;
    $tmp =~ s/-//g;
    
    my $new_seq = new Bio::LocatableSeq(-seq => $new_str, -id => $label, -end => length($tmp));
    $new_aln->add_seq($new_seq);
  }

  return $new_aln;
}


# GJ 2009-03-01 : Combines a sequence string and cigar line to return an alignment string.
sub combine_seq_cigar {
  my $class = shift;
  my $sequence = shift;		# Sequence string. (not bioseq!)
  my $cigar_line = shift;	# CIGAR line, in Ensembl format.

#  print $sequence."\n";
#  print $cigar_line . "\n";

  $cigar_line =~ s/([MD])/$1 /g;
  
  my @cigar_segments = split " ",$cigar_line;
  my $alignment_string = "";
  my $seq_start = 0;
  foreach my $segment (@cigar_segments) {
    if ($segment =~ /^(\d*)D$/) {
      my $length = $1;
      $length = 1 if ($length eq "");
      $alignment_string .= "-" x $length;
    } elsif ($segment =~ /^(\d*)M$/) {
      my $length = $1;
      $length = 1 if ($length eq "");
      $alignment_string .= substr($sequence,$seq_start,$length);
      $seq_start += $length;
    }
  }
#  print $alignment_string."\n";
  return $alignment_string;
}


# Masks out the entire sequence with ID $seq_id from the AlignI object $aln.
sub mask_seq_from_aln {
  my $class = shift;
  my $aln = shift;
  my $seq_id = shift;
  my $params = shift; # Optional params hashref.
  my $mask_char = shift;

  my $mask_char_from_params = $params->{alignment_quality_mask_character} || "X";

  $mask_char = $mask_char_from_params unless ($mask_char);

  my $seq = $class->get_seq_with_id($aln,$seq_id);
  my $seq_str = $seq->seq;
  
  printf "MASKING ENTIRE SEQUENCE %s: %.120s\n",$seq_id,$seq_str;
  #print "MASKING ENTIRE SEQUENCE $seq_id: ".$seq_str."\n";

  if ($seq->alphabet eq 'protein') {
    $seq_str =~ s/[^-]/$mask_char/g;
  } else {
    $seq_str =~ s/[^-]/N/g;
  }
  #print "new str: $seq_str\n";
  $seq->seq($seq_str);
  return $aln;
}

sub remove_seq_from_aln {
  my $class = shift;
  my $aln = shift;
  my $seq_id = shift;

  my $new_aln = $aln->new;

  foreach my $seq ($aln->each_seq) {
    next if ($seq->id eq $seq_id);
    my $new_seq = new Bio::LocatableSeq(-seq => $seq->seq, -id => $seq->id);
    $new_aln->add_seq($new_seq) if ($seq->id ne $seq_id);
  }

  return $new_aln;
}

sub remove_empty_seqs {
  my $class = shift;
  my $aln = shift;

  my @empty_seqs;
  foreach my $seq ($aln->each_seq) {
    my $str = $seq->seq;
    #print "Testing for empty:\n";
    #print $str."\n";
    if ($str !~ m/[^-X]/g) {
      push @empty_seqs, $seq;
    }
  }

  my $n = scalar(@empty_seqs);
  if ($n > 0) {
    print "  found $n empty seqs!\n";
  }

  foreach my $seq (@empty_seqs) {
    $aln = $class->remove_seq_from_aln($aln, $seq->id);
  }

  return $aln;
}

sub pretty_print {
  my $class = shift;
  my $aln = shift;
  my $params = shift;

  my $full = $params->{'full'} || 0;
  my $length = $params->{'length'} || $params->{width} || 150;

  $full = 1 if (!$full && $length > $aln->length + 10);
  
  if ($full) {
    my $num_slices = int($aln->length / $length);
    
    for (my $i=0; $i <= $num_slices; $i++) {
      my $start = $i*$length+1;
      my $end = ($i+1)*$length;

      my $top_string = "";
      for (my $j=0; $j < $length && ($start+$j) < $aln->length; $j++) {
	if ($j % 10 == 0) {
	  my $pos_string = "".($start+$j);
	  $top_string .= $pos_string;
	  $j += length($pos_string)-1;
	} else {
	  $top_string .= " ";
	}
      }
      printf(" %-20.20s  %s\n","",$top_string);

      foreach my $seq ($aln->each_seq) {
        $end = $seq->length if ($end > $seq->length);
        last if ($start > $seq->length || $end > $seq->length);
	my $seq_str = $seq->subseq($start,$end);
	printf(">%-20.20s  %s\n",$seq->id,$seq_str);
      }
    }
    return;
  }
  
  $length = ($length - 20)/2;

  my $aln_length = ($aln->length)." cols";
  my $empty_str = " "x$length;
  printf(" %-20.20s  %s   [ %10s ]   %s\n",
	 "",
	 $empty_str,
	 $aln_length,
	 $empty_str
	 );
	 
  foreach my $seq ($aln->each_seq) {
    my $seq_str = $seq->seq;
    my $end = substr($seq_str,-$length);
    my $start = substr($seq_str,0,$length);
    my $mid_str = "";
    if (length($seq_str) > 2*$length) {
      my $mid = substr($seq_str,$length,-$length);
      my $bases = $mid;
      $bases =~ s/-//g;
      my $gaps = $mid;
      $gaps =~ s/[^-]//g;
      my $percent = length($gaps)/length($mid) * 100;
      $mid_str = sprintf("(%1.1d%% gaps)",
			 $percent
			 );
    }

    printf(">%-20.20s  %s ... %10s ... %s\n",$seq->id,$start,$mid_str,$end);
  }
}


sub output {
  my $class = shift;
  my $aln = shift;
  my $file = shift;

  open(OUT,">$file");
  foreach my $seq ($aln->each_seq) {
    my $seq_str = $seq->seq;
    printf OUT ">%s\n%s\n",$seq->id,$seq_str;
  }
  close(OUT);
}


sub get_genomic_align {
  my $class = shift;
  my $tree = shift;
  my $params = shift;

  my $ref_taxon_id = 9606;
  if ($params->{ref_taxon_id}) {
    $ref_taxon_id = $params->{ref_taxon_id};
  }

  my ($ref_member) = grep {$_->taxon_id == 9606} $tree->leaves;
  my $ref_gene = $ref_member->get_Gene;
  my $ref_tx = $ref_member->transcript;
  my @exons = @{$ref_tx->get_all_translateable_Exons};

  foreach my $exon (@exons) {
    my $pep_exon = $exon->peptide($ref_tx);

    my $slice = $exon->slice;
    $exon = $exon->transfer($ref_tx->slice);
    my $dna_seq = $exon->seq->seq;
  }

}

# Given a peptide alignment and CDNA alignment, threads the CDNA sequence through the peptide
# alignment so both alignments are in 'sync'. Returns a new copy of the altered CDNA alignment.
sub apply_peptide_alignment_to_cdna_alignment {
  my $class = shift;
  my $pep_aln = shift;
  my $cdna_aln = shift;

  my $new_cdna = $cdna_aln->new;

  foreach my $seq ($pep_aln->each_seq) {
    my $cdna_seq = $class->get_seq_with_id($cdna_aln, $seq->id);
    die("No CDNA seq found!") unless (defined $cdna_seq);

    my $pep_string = $seq->seq;

    my $cdna_str = $cdna_seq->seq;
    $cdna_str =~ s/-//g; # Remove gaps from cdna string.
    
    # Thread the unaligned CDNA through the aligned pep.
    my $cdna_aln_string = $class->cdna_alignment_string($cdna_str,$pep_string);
    $cdna_aln_string =~ s/ //g;

    my $new_seq = new Bio::LocatableSeq(-seq => $cdna_aln_string, -id => $seq->id, -start => 1, -end => length($cdna_str));
    $new_cdna->add_seq($new_seq);
  }

  return $new_cdna;
}

# Thread cdna sequences through a given peptide alignment and ProteinTree.
sub peptide_to_cdna_alignment {
  my $class = shift;
  my $aln = shift;
  my $tree = shift;

  my $cdna_aln = $aln->new;

  my @leaves = $tree->leaves;
  foreach my $seq ($aln->each_seq) {
    my $pep_string = $seq->seq;
    my ($member) = grep {$_->name eq $seq->id} @leaves;
    die ("No member found for ".$seq->id."!") unless (defined $member);

    my $cdna_seq = $member->sequence_cds;

    my $cdna_aln_string = $class->cdna_alignment_string($cdna_seq,$pep_string);
    $cdna_aln_string =~ s/ //g;
    my $new_seq = new Bio::LocatableSeq(-seq => $cdna_aln_string, -id => $seq->id);
    $cdna_aln->add_seq($new_seq);
  }
  return $cdna_aln;
}

sub copy_aln {
  my $class = shift;
  my $aln = shift;
  
  my $new_aln = $aln->new;

  foreach my $seq ($aln->each_seq) {
    my $new_seq = new Bio::LocatableSeq(-seq => $seq->seq, -id => $seq->id);
    $new_aln->add_seq($new_seq);
  }
  return $new_aln;
}

# Stolen from AlignedMember.pm, wrapped into a standalone method.
sub cdna_alignment_string {
  my $class = shift;
  my $cdna = shift;  # CDNA string (no gaps)
  my $alignment_string = shift; # Peptide alignment string (with gaps)

  my $cdna_len = length($cdna);
  my $start = 0;
  my $cdna_align_string = '';
  
  #printf "%s %s\n", $self->stable_id,$self->member_id;
  # foreach my $pep (split(//, $self->alignment_string)) { # Speed up below
  foreach my $pep (unpack("A1" x length($alignment_string), $alignment_string)) {
    if($pep eq '-') {
      $cdna_align_string .= '--- ';
    } else {
      my $codon = substr($cdna, $start, 3);
      unless (length($codon) == 3) {
        # sometimes the last codon contains only 1 or 2 nucleotides.
        # making sure that it has 3 by adding as many Ns as necessary
        $codon .= 'N' x (3 - length($codon));
      }
      #print $pep." $codon ";
      if ($codon =~ m/(tga|tag|taa)/ig) {
        #print "$codon ";
        print " -> Stop codon found in alignment string! Masking...\n";
        $codon = 'N' x length($codon);
      }
      $cdna_align_string .= $codon . ' ';
      $start += 3;
    }
  }
  return $cdna_align_string;
}

sub has_stop_codon {
  my $class = shift;
  my $aln = shift;
  my $params = shift;

  my $codonTableId = 1;
  $codonTableId = $params->{bioperl_codontable_id} if (defined $params->{bioperl_codontable_id});

  my @seqs = $aln->each_seq;
  my @seq_strings = map {$_->seq} @seqs;

  my $code = Bio::Tools::CodonTable->new( -codontable_id => $codonTableId);

  for (my $i=0; $i < $aln->length-2; $i+= 3) {
    foreach my $str (@seq_strings) {
      my $codon = substr($str,$i,3);
      #print $codon;
      next if ($codon =~ m/^-+$/);

      my $ps = new Bio::PrimarySeq(-seq => $codon);
      my $aa = $ps->translate()->seq;
      #print "$i $codon $aa\n";
      if ($code->is_ter_codon($codon)) {
        print "Terminal codon: $codon  $i\n";
        return {
                aln_pos => $i
               };
      }
    }
  }  
  return 0;
}


sub get_substitution_runs {
  my $class = shift;
  my $aln = shift;

  my ($seq1,$seq2) = $aln->each_seq;
  my $str1 = $seq1->seq;
  my $str2 = $seq2->seq;

  my $sub_run_length = 0;
  my @sub_runs;
  for (my $i=0; $i <= $aln->length; $i++) {
    my $c1 = substr($str1,$i-1,1);
    my $c2 = substr($str2,$i-1,1);

    if ($c1 ne $c2 && $c1 !~ m/[-nx\.]/i && $c2 !~ m/[-nx\.]/i) {
      $sub_run_length++;
    } elsif ($sub_run_length > 0) {
      my $loc1 = $seq1->location_from_column($i);
      my $seq1_pos = $loc1->start if (defined $loc1);
      my $loc2 = $seq2->location_from_column($i);
      my $seq2_pos = $loc2->start if (defined $loc2);

      my $run_obj = {
                     ref_pos => $seq1_pos,
                     other_pos => $seq2_pos,
                     length => $sub_run_length,
                     aln_pos => $i
                    };
      push @sub_runs, $run_obj;
      $sub_run_length = 0;
    }
  }
  return @sub_runs;
}

sub keep_by_id {
  my $class = shift;
  my $aln = shift;
  my $seq_id_arrayref = shift;

  my $new_aln = $aln->new;
  my @ids = @{$seq_id_arrayref};
  foreach my $id (@ids) {
    my $seq = $class->get_seq_with_id($aln,$id);
    if (defined $seq) {
      $new_aln->add_seq($seq);
    }
  }
  return $new_aln;
}

1;# Keep perl happy
