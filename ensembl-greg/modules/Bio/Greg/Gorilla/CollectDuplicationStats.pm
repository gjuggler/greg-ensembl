package Bio::Greg::Gorilla::CollectDuplicationStats;

use strict;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Hive::Process;

use Bio::Greg::Gorilla::Utils;

use base ('Bio::Greg::Hive::Process','Bio::Greg::StatsCollectionUtils');

my $gene_stats_def = {
  data_id => 'int',
  node_id => 'int',
  orig_node_id => 'string',
  human_protein => 'string',
  human_gene    => 'string',
  human_names => 'string',
  gor_gene => 'string',
  gor_protein => 'string',
  gor_names => 'string',
  name => 'string',
  desc => 'string',
  unique_keys => 'data_id,node_id'
};
foreach my $species ('H','G','C','O') {
  $gene_stats_def->{$species.'_count'} = 'int';
}

sub fetch_input {
  my ($self) = @_;

  my $params = {
    collect_duplication_species => '9606,9598,9593,9600',
    genes_table    => 'stats_dups',
  };

  $self->load_all_params($params);

  # Create tables if necessary.
  $self->create_table_from_params( $self->compara_dba, $self->param('genes_table'),
    $gene_stats_def );
}

sub run {
  my $self = shift;

  my $node_id          = $self->param('node_id');
  my $parameter_set_id = $self->param('parameter_set_id');

  $self->get_gene_data( $node_id, $parameter_set_id);
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

  my $taxon_list = $self->param('collect_duplication_species');
  my @taxon_ids = split(',',$taxon_list);

  foreach my $taxon_id (@taxon_ids) {
    my @proteins = grep { $_->taxon_id == $taxon_id } $tree->leaves;
    my $count = scalar @proteins;
    my $letter = Bio::Greg::Gorilla::Utils->taxon_letter($taxon_id);
    $letter = $taxon_id if (!defined $letter);
    $cur_params->{$letter.'_count'} = $count;
  }


  # Collect human protein.
  my @human_proteins = grep { $_->taxon_id == 9606 } $tree->leaves;
  my @human_genes    = map  { $_->gene_member } @human_proteins;
  if ( scalar @human_proteins > 0 ) {
    my $member = $human_proteins[0];
    $cur_params->{'human_protein'} = $member->stable_id;
    $cur_params->{'human_gene'}    = $member->gene_member->stable_id;
    $cur_params->{'human_desc'} = $member->get_Gene->description;
  }
  $cur_params->{'human_names'} = join(", ",map {$_->get_Gene->external_name} @human_genes);

  # Collect gorilla protein.
  my @gor_proteins = grep { $_->taxon_id == 9593 } $tree->leaves;
  my @gor_genes    = map  { $_->gene_member } @gor_proteins;
  if ( scalar @gor_proteins > 0) {
    my $member = $gor_proteins[0];
    $cur_params->{'gor_protein'} = $member->stable_id;
    $cur_params->{'gor_gene'}    = $member->gene_member->stable_id;
    $cur_params->{'gor_desc'} = $member->get_Gene->description;
  }
  $cur_params->{'gor_names'} = join(", ",map {$_->get_Gene->external_name} @gor_genes);

  $cur_params->{'description'} = $cur_params->{'human_desc'} || $cur_params->{'gor_desc'};

  my @all_proteins = (@human_proteins,@gor_proteins);
  if (scalar @all_proteins > 0) {
    $cur_params->{'name'} = $all_proteins[0]->get_Gene->external_name;
}

  # Collect gene tag values into the params hash.
  $cur_params->{'tree_length'}      = $self->tree_length($tree);
  $cur_params->{'tree_max_path'}    = $self->max_path($tree);


  # Store values in our output table.
  my $table = $cur_params->{'genes_table'};
  $self->store_params_in_table( $self->db_handle, $table, $cur_params );
  $tree->release_tree;

}


1;
