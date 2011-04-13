#!/usr/bin/perl -w

use strict;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::TreeIO;
use Getopt::Long;


my $string = qq^(((((((((((((((((Human:0.003731,Chimp:0.005501):0.013010,Gorilla:0.010000)
:0.010000,Orangutan:0.020000):0.020000,Rhesus:0.031571):0.010000,
Marmoset:0.010000):0.020000,Tarsier:0.100000):0.050000,
(Mouse_lemur:0.084110,Bushbaby:0.145437):0.033956):0.020000,
TreeShrew:0.203975):0.020000,(((((Mouse:0.104920,Rat:0.109421):0.020000,
Kangaroo_rat:0.200000):0.050000,Guinea_Pig:0.150000):0.050000,
Squirrel:0.150000):0.020000,(Rabbit:0.100000,Pika:0.200000)
:0.050000):0.020000):0.020000,(((Alpaca:0.100000,
(Dolphin:0.120000,Cow:0.162368):0.050000):0.050000,
((Horse:0.150000,(Cat:0.098674,Dog:0.114682):0.002783):0.050000,
(Microbat:0.142600,Megabat:0.142460):0.003381):0.007170):0.030000,
(Hedgehog:0.279121,Shrew:0.309867):0.023929):0.040000):0.010000,
(((Elephant:0.110021,Rock_hyrax:0.150000):0.003000,Tenrec:0.265218):0.066339,
(Armadillo:0.178799,Sloth:0.200000):0.090000):0.053170)Mammalia:0.213469,
Opossum:0.320721):0.088647,Platypus:0.488110):0.118797,
((Chicken:0.230000,Zebra_finch:0.160000):0.160000,Lizard:0.513962)
:0.093688):0.151358,X_tropicalis:0.778272):0.174596,
(((Tetraodon:0.203933,Fugu:0.239587):0.203949,
(Stickleback:0.314162,Medaka:0.501915):0.055354):0.346008,
Zebrafish:0.730028):0.174596):0.000000,Lamprey:0.100000); ^;
my $tree = Bio::Tree::Tree->from_string($string);
print $tree->ascii;

print $tree->total_branch_length."\n";
my $slice = $tree->slice_by_ids('Human', 'Chimp', 'Gorilla', 'Orangutan', 'Rhesus', 'Marmoset', 'Tarsier', 'Mouse_lemur', 'Bushbaby');
my @children = $slice->children;
$slice = $children[0];
print $slice->ascii;
print $slice->total_branch_length."\n";
print $slice->max_distance_to_leaf."\n";

my $mammals = $tree->find('Mammalia');
print $mammals->total_branch_length."\n";

my ($tree_f, $ref_f, $aln_f);

GetOptions(
  'tree=s' => \$tree_f,
  'ref=s' => \$ref_f,
  'aln=s' => \$aln_f
  );

my $in;

# Load alignments.
$in = Bio::AlignIO->new(-file => $ref_f);
my $ref = $in->next_aln();
$in->close();
$in = Bio::AlignIO->new(-file => $aln_f);
my $aln = $in->next_aln();
$in->close();

# Load tree.
$in = Bio::TreeIO->new(-file => $tree_f);
my $tree = $in->next_tree();
$in->close();

# Exit if we have non-unique IDs or the IDs don't match up.
my (%ref_ids, %tree_ids, %aln_ids);
map {$tree_ids{$_->id} = 1} $tree->leaves;
die ("Non-unique IDs in the tree!") if (scalar(keys %tree_ids) != scalar($tree->leaves));
map {$ref_ids{$_->id} = 1} $ref->each_seq;
die ("Non-unique IDs in the ref!") if (scalar(keys %ref_ids) != scalar($ref->each_seq));
map {$aln_ids{$_->id} = 1} $aln->each_seq;
die ("Non-unique IDs in the aln!") if (scalar(keys %aln_ids) != scalar($aln->each_seq));

die ("Non equal number of sequences!") if (scalar($aln->each_seq) != scalar($ref->each_seq) or
                                           scalar($aln->each_seq) != scalar($tree->leaves));

my $css = get_correct_subtree_score($tree, $ref, $aln);
printf "CSS %.4f\n", $css;
my $tcs = total_column_score($ref, $aln);
printf "TCS %.4f\n", $tcs;
my $sps = sum_of_pairs_score($ref, $aln);
printf "SPS %.4f\n", $sps;

sub get_correct_subtree_score {
  my $tree = shift;
  my $ref = shift;
  my $aln = shift;

  my $total_bl = $tree->total_branch_length;

  my $ref_obj = to_arrayrefs($ref);
  my $aln_obj = to_arrayrefs($aln);
  my $ref_pairs = store_pairs($ref, $ref_obj);
  my $aln_pairs = store_pairs($aln, $aln_obj);

  my @id_list = map { $_->id } $ref->each_seq;
  my $tree_bl_hash;
  my $tree_strings_hash;
  my $scores_hash;
  map { $scores_hash->{$_} = [] } @id_list;

  my $correct_sum;
  my $nongap_sum;

  foreach my $i ( 1 .. $aln->length ) {
    my @nongap_ids_at_pos = grep { $aln_obj->{$_}->[$i] ne '-' } @id_list;

    #my $total_nongap_bl = get_subtree_bl( $tree, \@nongap_ids_at_pos, $tree_bl_hash );
    my $total_nongap_bl = 0;
    my $nongap_node_strings;
    my @node_strings = get_subtree_node_strings($tree, \@nongap_ids_at_pos, $tree_strings_hash);
    if (scalar(@node_strings) == 1) {
      $total_nongap_bl = 0;
    } else {
      map {$nongap_node_strings->{$_} = 1} @node_strings;
      foreach my $node ($tree->nodes) {
        if (defined $nongap_node_strings->{$node->enclosed_leaves_string}) {
          $total_nongap_bl += $node->branch_length;
        }
      }
    }

    my $seen_aligned_ids;
    my $correctly_aligned_hash; 
    my $aligned_cluster_count = 0;
    
    #print STDERR get_column_string($aln, $i) . "\n";
    my $column_score_string = '';
    foreach my $this_seq_id (@id_list) {
      my $this_res_num = $aln_obj->{$this_seq_id}->[$i];

      # Current residue is a gap!
      if ( $this_res_num eq '-' ) {
        $scores_hash->{$this_seq_id}->[$i] = -1;
        $column_score_string .= '-';
        next;
      }

      my @correctly_aligned_ids = ($this_seq_id);
      foreach my $other_seq_id (@id_list) {
        next if ( $this_seq_id eq $other_seq_id );
        my $other_res_num = $aln_obj->{$other_seq_id}->[$i];
        if ($other_res_num eq '-') {
          next;
        }
        my $pair_string = join( '_', $this_seq_id, $this_res_num, $other_seq_id, $other_res_num );
        my ($seq_a, $seq_b) = sort {$a cmp $b} ($this_seq_id,$other_seq_id);
        if ( defined $ref_pairs->{$pair_string} && $ref_pairs->{$pair_string} == 1 ) {
          if (!defined $seen_aligned_ids->{$this_seq_id}) {
            $aligned_cluster_count++;
            $seen_aligned_ids->{$this_seq_id} = $aligned_cluster_count;
            $seen_aligned_ids->{$other_seq_id} = $aligned_cluster_count; 
          } elsif (!defined $seen_aligned_ids->{$other_seq_id}) {
            $seen_aligned_ids->{$other_seq_id} = $aligned_cluster_count; 
          }
        }
      }
    }

    # For each correctly-aligned cluster, get the subtree branch length.
    my $correct_node_strings;
    my $total_correct_bl = 0;
    foreach my $k (1 .. $aligned_cluster_count) {
      my (@ids) = grep {$seen_aligned_ids->{$_} == $k} keys %$seen_aligned_ids;
      my @node_strings = get_subtree_node_strings($tree, \@ids, $tree_strings_hash);
      map {$correct_node_strings->{$_} = 1} @node_strings;
    }
    foreach my $node ($tree->nodes) {
      if (defined $correct_node_strings->{$node->enclosed_leaves_string}) {
        $total_correct_bl += $node->branch_length;
      }
    }

    my $score = 0;
    if ($total_nongap_bl > 0) {
      $score = $total_correct_bl / $total_nongap_bl;
    }

    foreach my $k (1 .. $aligned_cluster_count) {
      my $cluster_string = '';
      foreach my $seq ($aln->each_seq) {
        my $res = get_residue($aln, $seq->id, $i);
        if ($res eq '-') {
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

    #printf STDERR "  %.2f $total_correct_bl $total_nongap_bl \n", $score;

    # Sanity check: we should never have more correct BL than we have nongap BL.
    die if ($score > 1);

    $correct_sum += $total_correct_bl;
    $nongap_sum += $total_nongap_bl;
    #print STDERR "\n";
  }

  my $total_score = $correct_sum / $nongap_sum;
  return $total_score;
}

sub get_subtree_bl {
  my $tree    = shift;
  my $seq_ids = shift;
  my $bl_hash = shift;

  # Optimization: check the branch-length hash, to see if we've already
  # calculated the total branch length for the given list of IDs.
  my @id_array = @$seq_ids;
  @id_array = sort {$a cmp $b} @id_array;
  my $key = join('_',@id_array);
  my $existing_value = $bl_hash->{$key};
  return $existing_value if (defined $existing_value);

  # Get the minimum spanning subtree from a list of IDs
  my $subtree = $tree->slice_by_ids(@id_array);

  my $cur_total = $subtree->total_branch_length;
  $cur_total = $cur_total - $subtree->branch_length;
  if (scalar($subtree->leaves) == 1) {
    $cur_total = 0;
  }

  $bl_hash->{$key} = $cur_total;
  return $cur_total;
}

sub get_subtree_node_strings {
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

  my @node_strings;
  foreach my $node ($subtree->nodes) {
    push @node_strings, $node->enclosed_leaves_string;
  }
  
  $bl_hash->{$key} = \@node_strings;
  return @node_strings;
}


# Encode an alignment as an array of arrayrefs of residue numbers.
sub to_arrayrefs {
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

sub store_pairs {
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

sub get_residue {
  my $aln = shift;
  my $seq_id = shift;
  my $pos = shift;

  my @seqs = $aln->each_seq;
  my ($seq) = grep {$_->id eq $seq_id} @seqs;
  return $seq->subseq($pos,$pos);
}

sub get_column_string {
  my $aln = shift;
  my $pos = shift;

  my $str = "";
  foreach my $seq ($aln->each_seq) {
    my $residue = $seq->subseq($pos,$pos);
    $str .= $residue;
  }
  return $str;
}

sub total_column_score {
  my $true_aln = shift;
  my $test_aln = shift;

  # Index the alignments by sequence residue #.
  my $true_obj = to_arrayrefs($true_aln);
  my $test_obj = to_arrayrefs($test_aln);

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
  my $true_aln = shift;
  my $test_aln = shift;

  # Index the alignments by sequence residue #.
  my $true_obj = to_arrayrefs($true_aln);
  my $test_obj = to_arrayrefs($test_aln);

  # Put all aligned pairs into the hashtable.
  my $true_pairs = store_pairs($true_aln,$true_obj);
  my $test_pairs = store_pairs($test_aln,$test_obj);

  my $true_pair_count = scalar(keys %$true_pairs);
  my $test_pair_count = scalar(keys %$test_pairs);
  
  my $correct_pair_count = 0;
  foreach my $key (keys %$test_pairs) {
    $correct_pair_count++ if ($true_pairs->{$key});
  }
  
  return ($correct_pair_count/$test_pair_count);
}
