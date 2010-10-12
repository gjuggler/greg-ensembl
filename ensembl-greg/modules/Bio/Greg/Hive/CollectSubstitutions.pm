package Bio::Greg::Hive::CollectSubstitutions;

use strict;
use Bio::Greg::Codeml;
use File::Path;
use Bio::Greg::Gorilla::Utils;
use Bio::Align::Utilities qw(:all);

use Time::HiRes qw(sleep);
use DateTime;
use DateTime::Format::MySQL;

use base ( 'Bio::Greg::Hive::Process', 'Bio::Greg::Hive::Align',
  'Bio::Greg::StatsCollectionUtils' );

my $TREE = 'Bio::EnsEMBL::Compara::TreeUtils';

sub table_def {
  my $self = shift;

  return {
    data_id => 'int',

    ref_taxon_id     => 'int',
    ref_stable_id    => 'char32',
    ref_seq_name     => 'char32',
    ref_seq_pos_hg19 => 'int',
    ref_seq_pos_hg18 => 'int',
    ref_seq_strand   => 'char8',
    ref_nucleotide   => 'char8',
    ref_aa           => 'char8',

    other_taxon_id   => 'int',
    other_stable_id  => 'char32',
    other_seq_name   => 'char32',
    other_seq_pos    => 'int',
    other_seq_strand => 'char8',
    other_nucleotide => 'char8',
    other_aa         => 'char8',

    substitution_type => 'char32',

    unique_keys => 'ref_seq_name,ref_seq_pos'
  };
}

sub fetch_input {
  my $self           = shift;
  my $default_params = {
    table           => 'stats_subs',
    genome_taxon_id => '9606',
    other_taxon_id  => '9598',
  };

  # Fetch parameters from all possible locations.
  $self->load_all_params($default_params);
  $self->create_table_from_params( $self->hive_dba, $self->param('table'), $self->table_def );
}

sub run {
  my $self = shift;

  my $compara = $self->compara_dba;
  my $ma      = $compara->get_MemberAdaptor;

  my $stable_id      = $self->param('stable_id');
  my $other_taxon_id = $self->param('other_taxon_id');

  my $member = $ma->fetch_by_source_stable_id( undef, $stable_id );
  my ( $homology, $other_member ) =
    Bio::EnsEMBL::Compara::ComparaUtils->get_one_to_one_ortholog( $compara, $member,
    $other_taxon_id );

  $self->param( 'data_id', $member->dbID );

  #$self->new_data_id unless (defined $self->param('data_id'));

  if ( defined $other_member ) {
    print $member->stable_id . "\n";
    print $other_member->stable_id . "\n";
    $self->collect_and_store_substitutions($homology);
  } else {
    print "No ortholog for [$stable_id]!\n";
  }
}

sub collect_and_store_substitutions {
  my $self     = shift;
  my $homology = shift;

  my @genes = @{ $homology->gene_list };
  my ($ref_gene)   = grep { $_->taxon_id == $self->param('genome_taxon_id') } @genes;
  my ($other_gene) = grep { $_->taxon_id == $self->param('other_taxon_id') } @genes;

  my $ref_member   = $ref_gene->get_canonical_peptide_Member;
  my $other_member = $other_gene->get_canonical_peptide_Member;

  my $cur_params = $self->params;
  $self->hash_print($cur_params);

  my $aln_aa   = $homology->get_SimpleAlign;
  my $aln_cdna = $homology->get_SimpleAlign('cdna');

  return unless ( $aln_cdna->length > 1 );
  $self->pretty_print($aln_cdna);

  my $newick = sprintf( "(%s,%s);", $ref_member->stable_id, $other_member->stable_id );

  my $compara  = $self->compara_dba;
  my $pta      = $compara->get_ProteinTreeAdaptor;
  my $tree     = $pta->fetch_node_by_node_id( $homology->node_id );
  my @node_ids = ( $ref_member->stable_id, $other_member->stable_id );
  $tree = Bio::EnsEMBL::Compara::TreeUtils->extract_subtree_from_leaves( $tree, \@node_ids );

  # Re-scale tree if it's huge.
  if ( Bio::EnsEMBL::Compara::TreeUtils->total_distance($tree) > 5 ) {
    $tree = Bio::EnsEMBL::Compara::TreeUtils->scale_to( $tree, 5 );
  }
  print $tree->newick_format . "\n";

  ( $aln_aa, $aln_cdna ) = $self->get_prank_aln( $aln_aa, $aln_cdna, $tree );

  my @seqs = $aln_cdna->each_seq;
  my ($ref_seq)   = grep { $_->id eq $ref_member->stable_id } @seqs;
  my ($other_seq) = grep { $_->id eq $other_member->stable_id } @seqs;

  my $ref_tx   = $ref_member->get_Transcript;
  my $other_tx = $other_member->get_Transcript;

  my $ref_aa;
  my $other_aa;
  foreach my $i ( 1 .. $aln_cdna->length ) {
    if ( ( $i - 1 ) % 3 == 0 ) {

      # Store the amino acids.
      $ref_aa = $ref_seq->trunc( $i, $i + 2 )->translate->seq;
      $other_aa = $other_seq->trunc( $i, $i + 2 )->translate->seq;
      print "$ref_aa $other_aa\n";
    }
    my $ref_nuc = $ref_seq->subseq( $i, $i );
    my $other_nuc = $other_seq->subseq( $i, $i );

    #    next if ($ref_nuc =~ m/[-xn]/ig || $other_nuc =~ m/[-xn]/ig);
    #    print "$ref_nuc $other_nuc\n";

    my $loc1 = $ref_seq->location_from_column($i);
    my $loc2 = $other_seq->location_from_column($i);

    my $offset1 = $ref_tx->cdna_coding_start;
    my $offset2 = $other_tx->cdna_coding_start;

    if ( defined $loc1 && defined $loc2 ) {
      my $l1       = $loc1->start + $offset1 - 1;
      my $l2       = $loc2->start + $offset2 - 1;
      my @genomic1 = $ref_tx->cdna2genomic( $l1, $l1 );
      my @genomic2 = $other_tx->cdna2genomic( $l2, $l2 );
      if ( @genomic1 && @genomic2 ) {
        my $coord1 = $genomic1[0];
        printf( "%s %s %s\t", $ref_tx->seq_region_name, $coord1->start, $ref_nuc );
        my $coord2 = $genomic2[0];
        printf( "%s %s %s\n", $other_tx->seq_region_name, $coord2->start, $other_nuc );

        my $substitution_type = undef;
        $substitution_type = 'synonymous' if ( $ref_nuc ne $other_nuc );
        $substitution_type = 'nonsynonymous' if ( $ref_aa ne $other_aa && $ref_nuc ne $other_nuc );

        next unless defined($substitution_type);
        next if ( $ref_nuc =~ m/[-xn]/gi || $other_nuc =~ m/[-xn]/gi );
        next unless ( $coord1->isa('Bio::EnsEMBL::Mapper::Coordinate') );
        next unless ( $coord2->isa('Bio::EnsEMBL::Mapper::Coordinate') );

        my $hg18 =
          $self->get_hg18_from_hg19( $ref_tx->seq_region_name, $coord1->start, $coord1->strand );

        my $obj = {
          ref_taxon_id     => $ref_member->taxon_id,
          ref_stable_id    => $ref_member->stable_id,
          ref_seq_name     => $ref_tx->seq_region_name,
          ref_seq_strand   => $coord1->strand,
          ref_seq_pos_hg19 => $coord1->start,
          ref_seq_pos_hg18 => $hg18,
          ref_nucleotide   => $ref_nuc,
          ref_aa           => $ref_aa,

          other_taxon_id   => $other_member->taxon_id,
          other_stable_id  => $other_member->stable_id,
          other_seq_name   => $other_tx->seq_region_name,
          other_seq_pos    => $coord2->start,
          other_seq_strand => $coord2->strand,
          other_nucleotide => $other_nuc,
          other_aa         => $other_aa,

          substitution_type => $substitution_type
        };
        my $store_params = $self->replace( $cur_params, $obj );

        $self->store_params_in_table( $self->hive_dba->dbc, $self->param('table'), $store_params );
      }
    }
  }

}

sub get_prank_aln {
  my $self     = shift;
  my $pep_aln  = shift;
  my $cdna_aln = shift;
  my $tree     = shift;

  print " -> Aligning with Prank!\n";
  print $tree->newick_format . "\n";
  my $n = $cdna_aln->length;
  print "Before: $n\n";

  # Align with Prank to try and de-align incorrectly called exons.
  my $prank_params = { alignment_prank_f => 0 };
  my $new_aa = $self->align_with_prank( $pep_aln, $tree, $prank_params );
  my $new_cdna = Bio::EnsEMBL::Compara::AlignUtils->peptide_to_cdna_alignment( $new_aa, $tree );

  $n = $new_cdna->length;
  print "After $n:\n";
  $self->pretty_print( $new_cdna, { length => 200 } );

  return ( $new_aa, $new_cdna );
}
