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

  my $node_id = $self->_load_tree($tree);

  foreach my $i (1 .. $num_reps) {
    $self->param('slrsim_rep', $i);
    $self->fan_reps($tree, $node_id, $i, $self->params);
  }
}

sub _load_tree {
  my $self = shift;
  my $tree = shift;

  my $mba = $self->compara_dba->get_MemberAdaptor;
  my $pta = $self->compara_dba->get_ProteinTreeAdaptor;

  my $node = $tree->copy;

  my $ug = new Data::UUID;
  my $unique_string = $ug->create_str();

  my $params = $self->params;

  my $tree_mult = $params->{'slrsim_tree_mult'};
  if ($tree_mult) {
    Bio::EnsEMBL::Compara::TreeUtils->scale( $node, $tree_mult );
  }
  my $tree_length = $params->{'slrsim_tree_length'};
  if ($tree_length) {
    $node = Bio::EnsEMBL::Compara::TreeUtils->scale_to( $node, $tree_length );
  }
  my $tree_max_path = $params->{'slrsim_tree_max_path'};
  if ($tree_max_path) {
    print "Tree max path: $tree_max_path\n";
    $node = Bio::EnsEMBL::Compara::TreeUtils->scale_max_to( $node, $tree_max_path );
  }
  my $tree_mean_path = $params->{'slrsim_tree_mean_path'};
  if ($tree_mean_path) {
    print "Tree mean path: $tree_mean_path\n";
    $node = Bio::EnsEMBL::Compara::TreeUtils->scale_mean_to( $node, $tree_mean_path );
  }

  my $final_length = Bio::EnsEMBL::Compara::TreeUtils->total_distance($node);
  $self->param('slrsim_tree_length', $final_length);

  print "Final lenth: $final_length\n";

  # Go through each leaf and store the member objects.
  foreach my $leaf ( $node->leaves ) {
    $leaf->stable_id( $leaf->name );
    $leaf->source_name($unique_string);
    $mba->store($leaf);
    $leaf->adaptor(undef);
  }

  # Store the tree in the XYZ_node table.
  $pta->store($node);

  my $node_id = $node->node_id;

  $node->release_tree;

  return $node_id;
}

sub fan_reps {
  my $self    = shift;
  my $tree    = shift;
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
  my $data_id = $self->new_data_id($output_params);
  $output_params->{data_id} = $data_id;

  my ($job_id) = @{ $self->dataflow_output_id( $output_params, 1 ) };
  print "  -> Created job: $job_id \n";

}


1;
