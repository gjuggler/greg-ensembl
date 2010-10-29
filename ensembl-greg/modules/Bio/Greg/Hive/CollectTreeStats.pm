package Bio::Greg::Hive::CollectTreeStats;

use strict;
use Bio::EnsEMBL::Compara::Member;
use Bio::EnsEMBL::Compara::ComparaUtils;

use base ( 'Bio::Greg::Hive::Process', 'Bio::Greg::StatsCollectionUtils' );

sub table_def {
  my $self = shift;

  return {
    data_id      => 'int',
    node_id      => 'int',
    orig_node_id => 'int',

    tree_length    => 'float',
    tree_max_path  => 'float',
    tree_mean_path => 'float',
    tree_newick    => 'string',

    seq_length_mean => 'float',
    gc_content_mean => 'float',

    human_protein => 'char32',
    gor_protein   => 'char32',

    unique_keys => 'data_id,'
  };
}

sub fetch_input {
  my $self = shift;

  ### DEFAULT PARAMETERS ###
  my $params = {};
  $params->{'table'} = 'stats_trees';
  #########################

  $self->load_all_params($params);
  $self->create_table_from_params( $self->compara_dba, $self->param('table'), $self->table_def );
}

sub run {
  my $self = shift;

  my $cur_params = $self->params;

  my ( $tree, $sa_aln, $cdna_aln ) =
    Bio::EnsEMBL::Compara::ComparaUtils->tree_aln_cdna( $self->compara_dba, $cur_params );

  $cur_params->{'tree_length'}    = $self->tree_length($tree);
  $cur_params->{'tree_max_path'}  = $self->max_path($tree);
  $cur_params->{'tree_mean_path'} = $self->mean_path($tree);
  $cur_params->{'leaf_count'}     = scalar( $tree->leaves );

  $cur_params->{'tree_newick'} = $tree->newick_format;

  $cur_params->{'seq_length_mean'} = $self->seq_length_mean($tree);
  $cur_params->{'gc_content_mean'} = $self->gc_content_mean($tree);

  # Collect human protein.
  my @human_proteins = grep { $_->taxon_id == 9606 } $tree->leaves;
  if ( scalar @human_proteins > 0 ) {
    my $member = $human_proteins[0];
    $cur_params->{'human_protein'} = $member->stable_id;
  }

  # Collect gorilla protein.
  my @gor_proteins = grep { $_->taxon_id == 9593 } $tree->leaves;
  if ( scalar @gor_proteins > 0 ) {
    my $member = $gor_proteins[0];
    $cur_params->{'gor_protein'} = $member->stable_id;
  }

  $self->store_params_in_table( $self->db_handle, $self->param('table'), $cur_params );

  $tree->release_tree;
}

1;
