package Bio::Greg::Hive::SplitBySlidingWindow;

use strict;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::Process');

sub fetch_input {
  my $self = shift;

  ### DEFAULT PARAMETERS ###
  my $defaults = {
    window_width => 100,
    window_step => 20,
  };
  ##########################

  $self->load_all_params($defaults);
}

sub run {
  my $self = shift;

  $self->autoflow_inputjob(0);

  my $tree = $self->get_tree;
  my $aln = $self->get_aln;
  Bio::EnsEMBL::Compara::AlignUtils->pretty_print($aln,{length => 150});

  my $width = $self->param('window_width');
  my $step = $self->param('window_step');
  print "Splitting by sliding window!\n";
  print " -> width: $width\n";
  print " -> step: $step\n";
  print " -> aln length: ".$aln->length."\n";
  my $aln_lo = 1;
  my $aln_hi = 1+$width;
  $aln_hi = $aln->length if ($aln_hi > $aln->length); # Dummy check.
  while ($aln_hi <= $aln->length) {
    $self->flow_window($aln_lo,$aln_hi);

    $aln_lo += $step;

    last if ($aln_hi == $aln->length);
    $aln_hi = $aln_lo + $width;
    if ($aln_hi > $aln->length) {
      $aln_hi = $aln->length;
      $aln_lo = $aln_hi - $width;
      last if ($aln_lo < 0); # Dummy check.
    }
  }
}

sub flow_window {
  my $self = shift;
  my $aln_lo = shift;
  my $aln_hi = shift;

  print "Window: $aln_lo $aln_hi\n";

  # See Bio::Greg::Hive::Process->_get_aln method, where the input alignments are actually
  # split up into slices.
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
