package Bio::Greg::Gorilla::CollectTopologyStats;

use strict;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Hive::Process;

use base ( 'Bio::Greg::Hive::Process', 'Bio::Greg::StatsCollectionUtils' );

my $gene_stats_def = {
  data_id       => 'int',
  node_id       => 'int',
  orig_node_id  => 'string',
  human_protein => 'string',
  human_gene    => 'string',
  gor_gene      => 'string',
  gor_protein   => 'string',
  unique_keys   => 'data_id,node_id'
};

foreach my $model ( 'a', 'b', 'c' ) {
  $gene_stats_def->{ $model . '_lnL' } = 'float';
}

sub fetch_input {
  my ($self) = @_;

  my $params = { genes_table => 'stats_topology', };

  $self->load_all_params($params);

  # Create tables if necessary.
  $self->create_table_from_params( $self->compara_dba, $self->param('genes_table'),
    $gene_stats_def );
}

sub run {
  my $self = shift;

  my $node_id          = $self->param('node_id');
  my $parameter_set_id = $self->param('parameter_set_id');

  $self->get_gene_data( $node_id, $parameter_set_id );
}

sub get_gene_data {
  my $self             = shift;
  my $node_id          = shift;
  my $parameter_set_id = shift;

  my $cur_params = $self->params;

  my $pta = $self->compara_dba->get_ProteinTreeAdaptor;
  $pta->protein_tree_member("protein_tree_member");
  my ( $tree, $sa_aln, $cdna_aln );
  ( $tree, $sa_aln, $cdna_aln ) =
    Bio::EnsEMBL::Compara::ComparaUtils->tree_aln_cdna( $self->compara_dba, $cur_params );

  print $tree->newick_format . "\n";

  # Collect human protein.
  my @human_proteins = grep { $_->taxon_id == 9606 } $tree->leaves;
  my @human_genes    = map  { $_->gene_member } @human_proteins;
  if ( scalar @human_proteins > 0 ) {
    my $member = $human_proteins[0];
    $cur_params->{'human_protein'} = $member->stable_id;
    $cur_params->{'human_gene'}    = $member->gene_member->stable_id;
  }

  # Collect gorilla protein.
  my @gor_proteins = grep { $_->taxon_id == 9593 } $tree->leaves;
  my @gor_genes    = map  { $_->gene_member } @gor_proteins;
  if ( scalar @gor_proteins > 0 ) {
    my $member = $gor_proteins[0];
    $cur_params->{'gor_protein'} = $member->stable_id;
    $cur_params->{'gor_gene'}    = $member->gene_member->stable_id;
  }

  # Collect gene tag values into the params hash.
  $cur_params->{'tree_length'}   = $self->tree_length($tree);
  $cur_params->{'tree_max_path'} = $self->max_path($tree);

  # Store values in our output table.
  my $table = $cur_params->{'genes_table'};
  $self->store_params_in_table( $self->db_handle, $table, $cur_params );
  $tree->release_tree;

}

1;
