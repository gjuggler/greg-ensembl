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
  $self->store_param('species_with_paralogs', $paralog_count);

  # Calculate distances between each sequence. Turn Ns to gaps so they don't mess up
  # the calculations.
  my $map = {
    'N' => '-'
  };
  my $aln_copy = Bio::EnsEMBL::Compara::AlignUtils->translate_chars($aln, $map);
  my $stats = Bio::Align::DNAStatistics->new();
  my $jcmatrix = $stats->distance(-align => $aln_copy, 
                                  -method => 'D_JukesCantor');

  foreach my $taxon_id (keys %$tax_id_hash) {
    if ($tax_id_hash->{$taxon_id} > 0) {
      my @entries = grep {$_->taxon_id == $taxon_id} $tree->leaves;
      @entries = sort {$a->seq_length <=> $b->seq_length || 
                         $a->stable_id cmp $b->stable_id} @entries;

      my $min_dist = 9999;
      my $min_id = $entries[0];
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
        my $mean_dist = 0;
        $mean_dist = ( List::Util::sum(@dists)/$n_dists) if ($n_dists > 0);

        $min_id = $this_id if ($mean_dist <= $min_dist && $mean_dist > 0);
        $min_dist = $mean_dist if ($mean_dist <= $min_dist && $mean_dist > 0);
      }
      print "Min dist for $taxon_id: $min_id $min_dist\n";

      foreach my $entry (@entries) {
          # Remove all but the closest paralog from the alignment.
        if ($entry->stable_id ne $min_id) {
          my $cur_id = $entry->stable_id;
          print "Removing $cur_id\n";
          $aln = Bio::EnsEMBL::Compara::AlignUtils->remove_seq_from_aln($aln, $entry->stable_id);
        }
      }
    }
  }

  $tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_aln( $tree, $aln );
  return ($tree, $aln);
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

  return $self->mask_substitution_runs($tree, $aln);
}

1;
