package Bio::Greg::Slrsim::LoadTree;

use strict;
use Bio::EnsEMBL::Utils::Exception;

use base ('Bio::Greg::Hive::Process');

sub param_defaults {
  return {
  };
}

sub fetch_input {
  my $self = shift;

  $self->load_all_params();
}

sub run {
  my $self = shift;

  my $num_reps = $self->param('slrsim_replicates');
  my $tree_string = $self->param('tree_string');
  my $tree = Bio::EnsEMBL::Compara::TreeUtils->from_newick($tree_string);

  my $node_id = $self->param('node_id');

  foreach my $i (1 .. $num_reps) {
    $self->param('slrsim_rep', $i);
    $self->fan_reps($node_id, $i, $self->params);
  }
}

sub fan_reps {
  my $self    = shift;
  my $node_id = shift;
  my $sim_rep = shift;
  my $params  = shift;

  print " -> Node ID: $node_id\n";
  $self->param( 'node_id', $node_id );

  $params->{'slrsim_rep'}         = $sim_rep;
  delete $params->{'slrsim_tree_lengths'};
  delete $params->{'tree_string'};
  delete $params->{'output_folder'};

  $self->hash_print($params);

  my $params_file = $params->{slrsim_label} . "_params_" . $params->{slrsim_rep};
  my $file_params = {
    id => $params->{slrsim_label},
    filename => $params_file,
    extension => 'txt',
    subfolder => 'data'
  };
  my $f = $self->save_file($file_params);
  print "Full file: ".$f->{full_file}."\n";
  $self->frz($f->{full_file}, $params);

  my $output_params = {
    node_id => $node_id,
    slrsim_rep => $params->{slrsim_rep},
    slrsim_label => $params->{slrsim_label}
  };
  #my $data_id = $self->new_data_id($output_params);

  # Create a unique ID based on the node ID, # of reps, and current rep.
  my $n_reps = $self->param('slrsim_replicates');
  my $data_id = int($node_id * $n_reps) + $params->{slrsim_rep};
  $output_params->{data_id} = $data_id;

  my ($job_id) = @{ $self->dataflow_output_id( $output_params, 1 ) };
  print "  -> Created job: $job_id \n";

}


1;
