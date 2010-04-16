package Bio::Greg::Hive::SplitBySlidingWindow;

use strict;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::Process');

sub fetch_input {
  my $self = shift;

  ### DEFAULT PARAMETERS ###
  my $defaults = {
    window_width => 100,
    window_step => 30,
  };
  ##########################

  $self->load_all_params($defaults);
}

sub run {
  my $self = shift;

  $self->autoflow_inputjob(0);

  my $tree = $self->get_tree;
  my $aln = $self->get_aln;

  my $aln_lo = 1;
  my $aln_hi = $self->param('window_width');
  while ($aln_hi <= $aln->length) {
    $self->flow_window($aln_lo,$aln_hi);

    $aln_lo += $self->param('window_step');
    $aln_hi = $aln_lo + $self->param('window_width');
  }
}

sub flow_window {
  my $self = shift;
  my $aln_lo = shift;
  my $aln_hi = shift;

  print "Window: $aln_lo $aln_hi\n";

  my $added_params = {
    alignment_slice => "$aln_lo-$aln_hi"
  };

  my $input_params = $self->_parse_string('input_id');
  my $output_params = { %$input_params, %$added_params };

  $self->new_data_id($output_params);

  my ($job_id) = @{$self->dataflow_output_id($output_params, 1)};
  print "  -> Fanned window job: $job_id \n";
}

1;
