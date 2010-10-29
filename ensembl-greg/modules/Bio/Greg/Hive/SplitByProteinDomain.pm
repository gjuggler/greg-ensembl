package Bio::Greg::Hive::SplitByProteinDomain;

use strict;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::Process');

sub fetch_input {
  my $self = shift;

  ### DEFAULT PARAMETERS ###
  my $defaults = {
    domain_definition_taxon => '',
    min_domain_length       => 50,
    min_domain_score        => 50,
  };
  ##########################

  $self->load_all_params($defaults);
}

sub run {
  my $self = shift;

  $self->autoflow_inputjob(0);

  my $tree = $self->get_tree;
  my $aln  = $self->get_aln;

  my $min_domain_score  = $self->param('min_domain_score');
  my $min_domain_length = $self->param('min_domain_length');
  my $domain_taxon      = $self->param('domain_definition_taxon');

  my $domain_claims = {};

  foreach my $member ( $tree->leaves ) {
    next if ( $domain_taxon && $member->taxon_id != $domain_taxon );

    my $tx = $member->get_Translation;
    next if ( !defined $tx );
    my @features = @{ $tx->get_all_ProteinFeatures('PFam') };
    foreach my $f (@features) {
      my $pf_id = $f->display_id;

      # Start and end in the hit (Pfam domain) coordinates.
      my $pf_lo = $f->hstart;
      my $pf_hi = $f->hend;

      # Start and end in the query (Ensembl protein) coordinates.
      my $lo = $f->start;
      my $hi = $f->end;
      $hi = $tx->length if ( $hi > $tx->length );

      my $pf_score = $f->score;

      # Find the alignment coordinates.
      my $name   = $member->stable_id;
      my $aln_lo = $aln->column_from_residue_number( $name, $lo );
      my $aln_hi = $aln->column_from_residue_number( $name, $hi );

      if ( $hi - $lo > $min_domain_length && $pf_score > $min_domain_score ) {
        my $rounded_lo = int( $aln_lo / 10 ) * 10;
        my $rounded_hi = int( $aln_hi / 10 ) * 10;
        my $pf_key     = join( ".", $pf_id, $rounded_lo, $rounded_hi );
        if ( !$self->is_domain_claimed( $domain_claims, $pf_id, $aln_lo, $aln_hi ) ) {
          $self->flow_domain( $pf_id, $aln_lo, $aln_hi );
          $self->claim_domain( $domain_claims, $pf_id, $aln_lo, $aln_hi );
        }
      }
    }
  }
}

sub claim_domain {
  my $self   = shift;
  my $hash   = shift;
  my $pf_id  = shift;
  my $aln_lo = shift;
  my $aln_hi = shift;

  foreach my $i ( $aln_lo .. $aln_hi ) {
    my $key = $pf_id . $i;
    $hash->{$key} = 1;
  }
}

sub is_domain_claimed {
  my $self   = shift;
  my $hash   = shift;
  my $pf_id  = shift;
  my $aln_lo = shift;
  my $aln_hi = shift;

  foreach my $i ( $aln_lo .. $aln_hi ) {
    my $key = $pf_id . $i;
    return 1 if ( $hash->{$key} );
  }
  return 0;
}

sub flow_domain {
  my $self      = shift;
  my $domain_id = shift;
  my $aln_lo    = shift;
  my $aln_hi    = shift;

  print "Domain: $domain_id $aln_lo $aln_hi\n";

  my $added_params = {
    domain_id       => $domain_id,
    alignment_slice => "$aln_lo-$aln_hi"
  };

  my $input_params = $self->_parse_string('input_id');
  my $output_params = { %$input_params, %$added_params };

  $self->new_data_id($output_params);

  my ($job_id) = @{ $self->dataflow_output_id( $output_params, 1 ) };
  print "  -> Fanned domain job: $job_id \n";
}

1;
