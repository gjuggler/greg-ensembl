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

  my $tree_string = $self->param('tree_string');
  my $tree = Bio::EnsEMBL::Compara::TreeUtils->from_newick($tree_string);
  my $sim_rep = $self->param('slrsim_rep');

  $self->load_tree_into_database($tree,$sim_rep,$self->params);

}

sub load_tree_into_database {
  my $self    = shift;
  my $tree    = shift;
  my $sim_rep = shift;
  my $params  = shift;

  my $mba = $self->compara_dba->get_MemberAdaptor;
  my $pta = $self->compara_dba->get_ProteinTreeAdaptor;

  my $node = $tree->copy;

  my $md5 = Digest::MD5->new;
  $md5->add($params);
  $md5->add($sim_rep);
  my $unique_string = substr( $md5->hexdigest, 20 );

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

  print "Final lenth: $final_length\n";

  # Go through each leaf and store the member objects.
  foreach my $leaf ( $node->leaves ) {

    #print "leaf: ".$leaf->name."\n";
    $leaf->stable_id( $leaf->name );
    $leaf->source_name($unique_string);
    $mba->store($leaf);
    $leaf->adaptor(undef);
  }

  # Store the tree in the XYZ_node table.
  $pta->store($node);

  my $node_id = $node->node_id;
  print " -> Node ID: $node_id\n";
  $self->param( 'node_id', $node_id );

  #  printf "  > %.50s\n", $node->newick_format;

  # Store all parameters as tags.
  $params->{'slrsim_rep'}         = $sim_rep;
  $params->{'slrsim_tree_length'} = $final_length;
  delete $params->{'slrsim_tree_lengths'};
  my $sim_param_str = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
  foreach my $tag ( sort keys %{$params} ) {
    $self->store_tag( $tag, $params->{$tag} );
    sleep(0.1);
  }

  # Store the whole parameter set as a string.
  $self->store_tag( "params_slrsim", $sim_param_str );

  my $output_params = {
    node_id          => $node_id,
    parameter_set_id => $self->parameter_set_id,
    experiment_name  => $params->{'experiment_name'}
  };
  my $data_id = $self->new_data_id($output_params);
  $output_params->{data_id} = $data_id;

  my ($job_id) = @{ $self->dataflow_output_id( $output_params, 1 ) };
  print "  -> Created job: $job_id \n";

  $node->release_tree;
}


1;
