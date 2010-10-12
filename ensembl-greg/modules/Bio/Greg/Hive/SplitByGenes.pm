package Bio::Greg::Hive::SplitByGenes;

use strict;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::Process');

sub fetch_input {
  my $self = shift;

  ### DEFAULT PARAMETERS ###
  my $defaults = { 
    genome_taxon_id => 9606, 
    split_by_genes_limit => 0
  };
  ##########################

  $self->load_all_params($defaults);
}

sub run {
  my $self = shift;

  $self->autoflow_inputjob(0);

  my $taxon_id = $self->param('genome_taxon_id');

  my $compara  = $self->compara_dba;
  my $member_a = $compara->get_MemberAdaptor;

  my $lim = $self->param('split_by_genes_limit');
  my $i=0;
#  my @members = @{$member_a->fetch_all_by_source_taxon('ENSEMBLPEP',$taxon_id)};
#  print "Members: ".scalar(@members)."\n";
  foreach my $member ( @{ $member_a->fetch_all_by_source_taxon('ENSEMBLGENE', $taxon_id)} ) {
#    next unless ( $gene->biotype eq 'protein_coding' );

    last if ($lim > 0 && $i > $lim);

    $i++;
    $self->fan_gene($member);
  }
}

sub fan_gene {
  my $self = shift;
  my $member = shift;

  my $pep_member = $member->get_canonical_peptide_Member;
  return unless (defined $pep_member);

  my $added_params  = { stable_id => $pep_member->stable_id };
  my $input_params  = $self->_parse_string('input_id');
  my $output_params = { %$input_params, %$added_params };

  my ($job_id) = @{ $self->dataflow_output_id( $output_params, 1 ) };
  print "  -> Fanned member job: $job_id \n";
  $self->hash_print($output_params);
}

1;
