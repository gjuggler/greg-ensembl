package Bio::Greg::Hive::AlignmentFilteringSteps;

use strict;
use Bio::Greg::StatsCollectionUtils;
use Bio::Align::DNAStatistics;

use base ( 'Bio::Greg::Hive::Process',
           'Bio::Greg::StatsCollectionUtils' ,
           'Bio::Greg::Hive::CountSubstitutions' 
);

sub remove_paralogs {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $store_results = shift;
  
  my $tax_id_hash;
  foreach my $leaf ($tree->leaves) {
    my $tx_id = $leaf->taxon_id;
    if (!defined $tax_id_hash->{$tx_id}) {
      $tax_id_hash->{$tx_id} = 0;
    } else {
      $tax_id_hash->{$tx_id} = 1;
    }
  }

  my $paralog_count = 0;
  map {$paralog_count += $tax_id_hash->{$_}} keys %$tax_id_hash;
  if ($store_results) {
    $self->store_param('species_with_paralogs', $paralog_count);
  }

  # Calculate distances between each sequence. Turn Ns to gaps so they don't mess up
  # the calculations.
  my $map = {
    'N' => '-'
  };
  my $aln_copy = Bio::EnsEMBL::Compara::AlignUtils->translate_chars($aln, $map);
  my $stats = Bio::Align::DNAStatistics->new();
  my $jcmatrix = $stats->distance(-align => $aln_copy, 
                                  -method => 'D_JukesCantor');

  my $total_avg_dist = 0;
  my $seen_hash;
  my $n_dists = 0;
  my $tmp_total = 0;
  foreach my $l_i ($tree->leaves) {
    foreach my $l_j ($tree->leaves) {
      my $id_i = $l_i->stable_id;
      my $id_j = $l_j->stable_id;
      next if ($id_i eq $id_j);
      next if ($seen_hash->{$id_i.$id_j});
      $seen_hash->{$id_i.$id_j} = 1;
      $seen_hash->{$id_j.$id_i} = 1;
      my $dist = $jcmatrix->get_entry($id_i, $id_j);
      next if ($dist == -1 || $dist >= 5);
      $tmp_total += $dist;
      $n_dists += 1;
    }
  }
  $total_avg_dist = $tmp_total / $n_dists;

  my $dbh;
  if ($store_results) {
    $self->create_table_from_params( $self->dbc, 'paralogs',
                                     $self->_paralog_table_structure );
    $dbh = $self->db_handle;
    $dbh->begin_work;
    
    print "  storing paralogs...\n";
  }

  my $n = 0;
  my $tot_length = 0;
  foreach my $leaf ($tree->leaves) {
    my $str = $leaf->sequence;
    $str =~ s/N-//gi;  
    $tot_length += length($str);
    $n++;
  }
  my $mean_length = $tot_length / $n;  
  
  foreach my $taxon_id (keys %$tax_id_hash) {
    if ($tax_id_hash->{$taxon_id} > 0) {
      my @entries = grep {$_->taxon_id == $taxon_id} $tree->leaves;
      @entries = sort {$a->seq_length <=> $b->seq_length || 
                         $a->stable_id cmp $b->stable_id} @entries;

      my $min = $entries[0];
      my $mean_distances;

      foreach my $entry (@entries) {
        my $this_id = $entry->stable_id;
        my @dists;
        my $n_dists = 0;
        foreach my $other_entry (@entries) {
          next if ($other_entry == $entry);
          my $that_id = $other_entry->stable_id;
          my $dist = $jcmatrix->get_entry($that_id, $this_id);
          if (defined $dist && $dist < 5 && $dist > 0) {
            push @dists, $dist;
            $n_dists++;
          } 
        }
        my $mean_dist = -1;
        $mean_dist = ( List::Util::sum(@dists)/$n_dists) if ($n_dists > 0);
        $mean_distances->{$entry->stable_id} = $mean_dist;

        my $min_seq = $min->sequence;
        $min_seq =~ s/N-//gi;
        my $min_seq_length = length($min_seq);
        $min_seq_length = 1 if ($min_seq_length == 0);
        my $cur_seq = $entry->sequence;
        $cur_seq =~ s/N-//gi;
        my $cur_seq_length = length($cur_seq);
        $cur_seq_length = 1 if ($cur_seq_length == 0);
        my $min_ratio = $min_seq_length / $mean_length;
        my $cur_ratio = $cur_seq_length / $mean_length;

        if ($min_ratio < 0.5 && $cur_ratio > 0.5) {
          # Skip the distance comparison if the current length is >50% and
          # the "min" entry's length is <50%
          $min = $entry;
          next;
        } elsif ($min_ratio > 0.5 && $cur_ratio < 0.5) {
          next;
        }

        my $min_dist = $mean_distances->{$min->stable_id};
        if ($mean_dist == -1 && $min_dist == -1) {
          # If both seqs failed to give a distance, compare lengths.
          $min = $entry if ($cur_ratio > $min_ratio);
        } elsif ($min_dist == -1 && $mean_dist > -1) {
          $min = $entry;
        } elsif ($min_dist > -1 && $mean_dist > -1) {
          my $dist_ratio = $mean_dist / $min_dist;
          if ($dist_ratio < 1) {
            $min = $entry;
          } elsif ($dist_ratio == 1) {
            # If both seqs have equal distances, keep based on length.
            $min = $entry if ($cur_ratio > $min_ratio);
          }
        }
      }

      my $kept_seq = $min->sequence;
      $kept_seq =~ s/-//g;
      my $kept_len = length($kept_seq);
      my $kept_len_ratio = $kept_len / $mean_length;
      my $kept_dist = $mean_distances->{$min->stable_id};

      my $n_paralogs = scalar(@entries);

      foreach my $entry (@entries) {
          # Remove all but the closest paralog from the alignment.
        if ($entry->stable_id ne $min->stable_id) {
          my $cur_id = $entry->stable_id;
          print "Removing $cur_id\n";

          my $ungap_seq = $entry->sequence;
          $ungap_seq =~ s/-//g;
          my $cur_len = length($ungap_seq);
          my $cur_len_ratio = $cur_len / $mean_length;

          my $p = $self->replace($self->params,{
            taxon_id => $entry->taxon_id,
            removed_id => $entry->stable_id,
            removed_len => $cur_len,
            removed_len_ratio => $cur_len_ratio,
            removed_dist => $mean_distances->{$entry->stable_id},
            kept_id => $min->stable_id,
            kept_len => $kept_len,
            kept_len_ratio => $kept_len_ratio,
            kept_dist => $kept_dist,
            n_paralogs => $n_paralogs,
            avg_dist => $total_avg_dist
                                 });
          if ($store_results) {
            $self->store_params_in_table($dbh, 'paralogs', $p);
          }
          $aln = Bio::EnsEMBL::Compara::AlignUtils->remove_seq_from_aln($aln, $entry->stable_id);
        }
      }
    }
  }

  if ($store_results) {
    $dbh->commit;
    print "  Done!\n";
  }

  $tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_aln( $tree, $aln );
  return ($tree, $aln);
}

sub _paralog_table_structure {
  my $self = shift;
  return {
    data_id => 'int',
    taxon_id => 'int',
    removed_id => 'char32',
    removed_len => 'int',
    removed_len_ratio => 'float',
    removed_dist => 'float',
    kept_id => 'char32',
    kept_len => 'int',
    kept_len_ratio => 'float',
    kept_dist => 'float',
    n_paralogs => 'int',
    avg_dist => 'float',
    unique_keys => 'data_id,removed_id'
  };
}

sub _filter_bad_seqs {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

#  $aln = Bio::EnsEMBL::Compara::AlignUtils->remove_seq_from_aln($aln, $entry->stable_id);
#  $tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_aln( $tree, $aln );

  return ($tree, $aln);
}

sub filter_substitution_runs {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $do_filtering = shift;

  # Call the method from CountSubstitutions.pm
  return $self->mask_substitution_runs($tree, $aln, $do_filtering);
}

1;
