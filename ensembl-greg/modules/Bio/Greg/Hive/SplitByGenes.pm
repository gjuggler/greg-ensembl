package Bio::Greg::Hive::SplitByGenes;

use strict;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::Process');

sub fetch_input {
  my $self = shift;

  ### DEFAULT PARAMETERS ###
  my $defaults = { genome_taxon_id => 9606, };
  ##########################

  $self->load_all_params($defaults);
}

sub run {
  my $self = shift;

  $self->autoflow_inputjob(0);

  my $compara  = $self->compara_dba;
  my $gdba     = $compara->get_GenomeDBAdaptor;
  my $gdb      = $gdba->fetch_by_taxon_id( $self->param('genome_taxon_id') );
  my $core_dba = $gdb->db_adaptor;
  my $gene_a   = $core_dba->get_GeneAdaptor;

  foreach my $gene ( @{ $gene_a->fetch_all } ) {
    next unless ( $gene->biotype eq 'protein_coding' );

    $self->fan_gene($gene);
  }
}

sub fan_gene {
  my $self = shift;
  my $gene = shift;

  my $canonical_ts = $gene->canonical_transcript;
  print $canonical_ts->stable_id . "\n";

  my $added_params  = { transcript_stable_id => $canonical_ts->stable_id };
  my $input_params  = $self->_parse_string('input_id');
  my $output_params = { %$input_params, %$added_params };

  my ($job_id) = @{ $self->dataflow_output_id( $output_params, 1 ) };
  print "  -> Fanned gene job: $job_id \n";
  $self->hash_print($output_params);
}

1;
