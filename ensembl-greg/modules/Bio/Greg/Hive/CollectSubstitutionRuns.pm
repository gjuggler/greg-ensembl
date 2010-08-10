package Bio::Greg::Hive::CollectSubstitutionRuns;

use strict;
use Bio::Greg::Codeml;
use File::Path;
use Bio::Greg::Gorilla::Utils;
use Bio::Align::Utilities qw(:all);

use Time::HiRes qw(sleep);
use DateTime;
use DateTime::Format::MySQL;

use base ( 'Bio::Greg::Hive::Process', 'Bio::Greg::Hive::Align' );

my $TREE = 'Bio::EnsEMBL::Compara::TreeUtils';

sub table_def {
  my $self = shift;

  return {
    data_id => 'int',

    member_id      => 'int',
    ref_taxon_id   => 'int',
    other_taxon_id => 'int',

    length    => 'int',
    aln_pos   => 'int',
    ref_pos   => 'int',
    other_pos => 'int',

    ref_stable_id   => 'char32',
    other_stable_id => 'char32',
    ref_exon_id     => 'char32',
    other_exon_id   => 'char32',

    unique_keys => 'member_id,other_taxon_id,aln_pos'
  };
}

sub fetch_input {
  my $self           = shift;
  my $default_params = {
    ref_taxon_id    => 9606,          # set to 'fan_jobs' to fan out alignment gathering jobs.
    other_taxon_ids => '9593',
    aln_type        => 'aa',
    table           => 'stats_subs'
  };

  # Fetch parameters from all possible locations.
  $self->load_all_params($default_params);
  $self->create_table_from_params( $self->compara_dba, $self->param('table'), $self->table_def );

  $self->param( 'data_id', $self->param('member_id') );
}

sub run {
  my $self = shift;

  my $member_id = $self->param('member_id');
  my $compara   = $self->compara_dba;
  my $mba       = $compara->get_MemberAdaptor;
  my $ha        = $compara->get_HomologyAdaptor;
  my $gdba      = $compara->get_GenomeDBAdaptor;

  my $nhx =
    Bio::EnsEMBL::Compara::ComparaUtils->get_genome_tree_nhx( $compara, { labels => 'mnemonics' } );
  print "$nhx\n";

  my $member   = $mba->fetch_by_dbID($member_id);
  my $this_gdb = $gdba->fetch_by_taxon_id( $self->param('ref_taxon_id') );

  my @taxon_ids = $self->get_taxon_ids_from_keepers_list( $self->param('other_species') );

  foreach my $taxon_id (@taxon_ids) {
    my $other_gdb = $gdba->fetch_by_taxon_id($taxon_id);

    my @homologies = @{ $ha->fetch_all_by_Member_paired_species( $member, $other_gdb->name ) };
    foreach my $homology (@homologies) {
      my $desc = $homology->description;
      if ( $desc eq 'ortholog_one2one' ) {
        my @genes = @{ $homology->gene_list };

        my ($ref_member)   = grep { $_->taxon_id == $member->taxon_id } @genes;
        my ($other_member) = grep { $_->taxon_id == $taxon_id } @genes;
        $ref_member   = $ref_member->get_canonical_peptide_Member;
        $other_member = $other_member->get_canonical_peptide_Member;
        next if ( !defined $ref_member->get_Transcript );
        next if ( !defined $other_member->get_Transcript );
        my @ref_exons   = @{ $ref_member->get_Transcript->get_all_translateable_Exons() };
        my @other_exons = @{ $other_member->get_Transcript->get_all_translateable_Exons() };

        my $aln;
        if ( $self->param('aln_type') eq 'cdna' ) {

          # Get the CDNA alignment.
          print "CDNA!\n";
          $aln = $homology->get_SimpleAlign('cdna');
        } elsif ( $self->param('aln_type') eq 'epo' ) {

          # Get the EPO alignment.
          print $member->stable_id . "\n";
          my ( $genomic_cdna, $genomic_aa ) =
            Bio::EnsEMBL::Compara::ComparaUtils->genomic_aln_for_member( $compara, $member );
          $aln = $genomic_cdna;
          my @species_names = map { $_->name } ( $this_gdb, $other_gdb );
          $aln = Bio::EnsEMBL::Compara::AlignUtils->keep_by_id( $aln, \@species_names );
          next if ( scalar( $aln->each_seq ) < 2 );
          $self->pretty_print( $aln, { length => 200 } );
        } else {

          # Get the reg'ler AA alignment.
          $aln = $homology->get_SimpleAlign();
        }

        #        $aln = $aln->remove_gaps(undef,1);

        $self->pretty_print( $aln, { length => 100, full => 1 } );
        my @substitution_run_lengths =
          Bio::EnsEMBL::Compara::AlignUtils->get_substitution_runs($aln);
        foreach my $run_obj (@substitution_run_lengths) {

          #         $self->hash_print($run_obj);
          printf "%d %d\n", $run_obj->{aln_pos}, $run_obj->{length};

          # Get the exon for the given sequence position.
          my $ref_exon =
            $self->get_exon_from_position( $ref_member, $run_obj->{ref_pos}, \@ref_exons );
          my $other_exon =
            $self->get_exon_from_position( $other_member, $run_obj->{other_pos}, \@other_exons );

          my $cur_params = $self->params;
          $cur_params = $self->replace(
            $cur_params,
            $run_obj, {
              other_taxon_id => $taxon_id,
              ref_stable_id  => $member->stable_id,
            }
          );
          $cur_params->{other_exon_id} = $other_exon->stable_id if ( defined $other_exon );
          $cur_params->{ref_exon_id}   = $ref_exon->stable_id   if ( defined $ref_exon );
          $cur_params->{other_stable_id} = $other_member->stable_id;
          $cur_params->{ref_stable_id}   = $ref_member->stable_id;

          $self->store_params_in_table( $self->dbc, $self->param('table'), $cur_params );
        }
      }
    }
  }
}

sub get_exon_from_position {
  my $self          = shift;
  my $member        = shift;
  my $pos           = shift;
  my $exon_arrayref = shift;

  my $tr    = $member->get_Transcript;
  my @exons = @$exon_arrayref;

  #print "Num exons: ".scalar(@exons)."\n";
  my ($best_exon) = grep {
    my $tmp_exon = $_->transfer( $tr->slice );
    if ( !$tmp_exon ) {
      1;
    }
    my @coords = $tr->genomic2pep( $tmp_exon->start, $tmp_exon->end, $tmp_exon->strand );
    @coords = grep { $_->isa('Bio::EnsEMBL::Mapper::Coordinate') } @coords;

    if ( scalar(@coords) == 1 ) {
      my $c = $coords[0];

      #print $tmp_exon->stable_id." ". $c->start." ".$c->end." ".$pos."\n";
      1 if ( $c->start <= $pos && $c->end > $pos );
    }
    1;
  } @exons;
  return $best_exon;
}

1;
