package Bio::Greg::Hive::SplitByParameterSet;

use strict;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::Process');

sub fetch_input {
  my $self = shift;

  ### DEFAULT PARAMETERS ###
  my $defaults = {
    domain_definition_taxon => '',
    min_domain_length => 50
  };
  ##########################

  $self->load_all_params($defaults);
}

sub run {
  my $self = shift;

  my $tree = $self->get_tree;
  my $aln = $self->get_aln;

  my $min_domain_length = $self->param('min_domain_length');
  my $domain_taxon = $self->param('domain_definition_taxon');

  foreach my $member ($tree->leaves) {
    next if ($domain_taxon && $member->taxon_id != $domain_taxon);

    my $tx   = $member->get_Translation;
    next if ( !defined $tx );
    my @features = @{ $tx->get_all_ProteinFeatures('PFam') };
    foreach my $f (@features) {
      my $pf_id = $f->display_id;
      # Start and end in the hit (Pfam domain) coordinates.
      my $pf_lo = $f->hstart;       
      my $pf_hi = $f->hend;
      # Start and end in the query (Ensembl protein) coordinates.
      my $lo    = $f->start;
      my $hi    = $f->end;
      $hi = $tx->length if ( $hi > $tx->length );
      
      # Find the alignment coordinates.
      my $name = $member->stable_id;
      my $aln_lo = $aln->column_from_residue_number( $name, $lo );
      my $aln_hi = $aln->column_from_residue_number( $name, $hi );

      if ($aln_hi - $aln_lo > $min_domain_length) {
	$self->flow_domain($pf_id,$aln_lo,$aln_hi);
      }
    }
  }
}

sub flow_domain {
  my $self = shift;
  my $domain_id = shift;
  my $aln_lo = shift;
  my $aln_hi = shift;

  my $added_params = {
    domain_id => $domain_id,
    alignment_slice => "$aln_lo-$aln_hi"
  };

  my $input_params = $self->_parse_string('input_id');
  my $output_params = { %$input_params, %$added_params };

  my ($job_id) = @{$self->dataflow_output_id($output_params, 1)};
  print "  -> Fanned domain job: $job_id \n";
}

1;
