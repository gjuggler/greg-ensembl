package Bio::Greg::Hive::SplitByGenomeWindows;

use strict;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::Process');

sub fetch_input {
  my $self = shift;

  ### DEFAULT PARAMETERS ###
  my $defaults = {
    genome_taxon_id       => 9606,
    genome_window_width   => 1000000,
    genome_window_overlap => 0
  };
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
  my $slice_a  = $core_dba->get_SliceAdaptor;
  my @slices   = @{ $slice_a->fetch_all('toplevel') };

  foreach my $slice (@slices) {
    $self->fan_slice_jobs($slice);
  }
}

sub fan_slice_jobs {
  my $self  = shift;
  my $slice = shift;

  my $slice_name = $slice->seq_region_name;

  my $window_size = $self->param('genome_window_width');
  my $window_step = $window_size - $self->param('genome_window_overlap');
  for ( my $i = $slice->start ; $i < $slice->end ; $i += $window_step ) {
    my $start = $i;
    my $end   = $i + $window_size;

    if ( $start > $slice->end ) {
      last;
    }
    if ( $end > $slice->end ) {
      $end = $slice->end;
    }

    my $added_params = {
      genome_slice_name  => $slice_name,
      genome_slice_start => $start,
      genome_slice_end   => $end
    };
    my $input_params = $self->_parse_string('input_id');
    my $output_params = { %$input_params, %$added_params };

    my ($job_id) = @{ $self->dataflow_output_id( $output_params, 1 ) };
    print "  -> Fanned genome window job: $job_id \n";
    $self->hash_print($output_params);
  }
}

1;
