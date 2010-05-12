package Bio::Greg::Hive::CollectSitewiseStats;

use strict;
use Time::HiRes qw(sleep);

use Cwd;
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::Process;

use base ( 'Bio::Greg::Hive::Process', 'Bio::Greg::StatsCollectionUtils' );

sub fetch_input {
  my ($self) = @_;

  my $params = {
    omega_table => 'sitewise_omega',
    sites_table => 'stats_sites',
    genes_table => 'stats_genes'
  };

  $self->load_all_params($params);

  # Create tables if necessary.
  $self->create_table_from_params(
    $self->compara_dba,
    $self->param('genes_table'),
    $self->get_gene_table_structure
  );
  $self->create_table_from_params(
    $self->compara_dba,
    $self->param('sites_table'),
    $self->get_sites_table_structure
  );
}

sub run {
  my $self = shift;

  print " -> Getting gene data...\n";
  $self->get_gene_data();
  print " -> Getting sites data...\n";
  $self->get_sites_data();
}

sub get_gene_data {
  my $self = shift;

  my $gene_data = $self->data_for_gene;

#  $self->hash_print($gene_data);
  $self->store_params_in_table( $self->compara_dba, $self->param('genes_table'), $gene_data );
}

sub data_for_gene {
  my $self = shift;

  # Take a snapshot of the current params.
  my $cur_params = $self->params;

  my $pta = $self->compara_dba->get_ProteinTreeAdaptor;
  $pta->protein_tree_member("protein_tree_member");
  my ( $tree, $sa_aln, $cdna_aln );
  ( $tree, $sa_aln, $cdna_aln ) =
    Bio::EnsEMBL::Compara::ComparaUtils->tree_aln_cdna( $self->compara_dba, $cur_params );

  $cur_params->{'data_id'} = $self->data_id;
  $cur_params->{'parameter_set_id'} = $self->parameter_set_id;

  # Collect gene tag values into the params hash.
  $cur_params->{'tree_length'}    = $self->tree_length($tree);
  $cur_params->{'tree_max_path'}  = $self->max_path($tree);
  $cur_params->{'tree_mean_path'} = $self->mean_path($tree);
  $cur_params->{'leaf_count'}     = scalar( $tree->leaves );

  $cur_params->{'tree_newick'} = $tree->newick_format;

  $cur_params->{'seq_length_mean'} = $self->seq_length_mean($tree);
  $cur_params->{'gc_content_mean'} = $self->gc_content_mean($tree);

  $cur_params->{'mean_copy_count'} = $self->mean_copy_count($tree);

  $cur_params->{'duplication_count'}    = $self->duplication_count($tree);
  $cur_params->{'duplication_fraction'} = $self->duplication_fraction($tree);

  my $psc_hash = $self->get_psc_hash( $self->compara_dba->dbc, $self->params );
  $cur_params->{'sitewise_value_count'} = $self->sitewise_count($psc_hash);
  if ( scalar( keys(%$psc_hash) ) > 0 ) {
    $cur_params->{'psc_count'}      = $self->psc_count( $psc_hash, 0 );
    $cur_params->{'weak_psc_count'} = $self->psc_count( $psc_hash, 1 );
    $cur_params->{'max_lrt'}        = $self->max_lrt($psc_hash);
  }

  my $ps = $cur_params->{'parameter_set_id'};
  $cur_params->{'slr_dnds'}      = $self->param('slr_omega');
  $cur_params->{'slr_kappa'}     = $self->param('slr_kappa');
  $cur_params->{'hyphy_dnds'}    = $self->param('hyphy_dnds');
  $cur_params->{'hyphy_dnds_lo'} = $self->param('hyphy_dnds_lo');
  $cur_params->{'hyphy_dnds_hi'} = $self->param('hyphy_dnds_hi');

  return $cur_params;
}

sub get_sites_data {
  my $self = shift;

  my $tree = $self->get_tree;
  my $aln  = $tree->get_SimpleAlign;

  foreach my $aln_position ( 1 .. $aln->length ) {
    my $site_data = $self->data_for_site($aln_position);
    if ( defined $site_data ) {
      $self->store_params_in_table( $self->compara_dba, $self->param('sites_table'), $site_data );
    }
  }

}

sub data_for_site {
  my $self         = shift;
  my $aln_position = shift;

  # Get the hash object containing site-wise data.
  my $site_hash = $self->param('site_hash');
  if ( !defined $site_hash ) {
    my $site_params = $self->params;
    $site_hash = $self->get_psc_hash( $self->compara_dba->dbc, $site_params );
    $self->param( 'site_hash', $site_hash );

    my $tree = $self->get_tree;
    my $aln  = $tree->get_SimpleAlign;
    $self->param( 'aln_length', $aln->length );
  }

  my $site_data = $self->params;
  my $site      = $site_hash->{$aln_position};

  return undef unless ( defined $site );

  $site_data->{omega}       = $site->{omega};
  $site_data->{omega_lower} = $site->{omega_lower};
  $site_data->{omega_upper} = $site->{omega_upper};
  $site_data->{lrt_stat}    = $site->{lrt_stat};
  $site_data->{ncod}        = $site->{ncod};
  $site_data->{type}        = $site->{type};
  $site_data->{note}        = $site->{note};

  $site_data->{'aln_position'}          = $site->{aln_position};
  $site_data->{'aln_position_fraction'} = $site->{aln_position} / $self->param('aln_length');

  $site_data->{'chr_name'}  = $site->{chr_name};
  $site_data->{'chr_start'} = $site->{chr_start};
  $site_data->{'chr_end'}   = $site->{chr_end};

  $site_data->{'data_id'} = $self->data_id;
  $site_data->{'parameter_set_id'} = $self->parameter_set_id;

  return $site_data;
}

sub get_gene_table_structure {
  my $self = shift;

  my $structure = {
    'data_id'                 => 'int',

    'node_id' => 'int',
    'orig_node_id' => 'int',

    'parameter_set_id'   => 'tinyint',
    'parameter_set_name' => 'string',

    'tree_length'    => 'float',
    'tree_max_path'  => 'float',
    'tree_mean_path' => 'float',
    'leaf_count'     => 'int',
    'tree_newick'    => 'string',

    'seq_length_mean' => 'float',
    'gc_content_mean' => 'float',

    'duplication_count'    => 'int',
    'duplication_fraction' => 'int',
    'mean_copy_count'      => 'float',

    'sitewise_value_count' => 'int',
    'psc_count'            => 'int',
    'weak_psc_count'       => 'int',
    'max_lrt'              => 'float',

    'slr_dnds'      => 'float',
    'slr_kappa'     => 'float',
    'hyphy_dnds'    => 'float',
    'hyphy_dnds_lo' => 'float',
    'hyphy_dnds_hi' => 'float',

    unique_keys => 'data_id,parameter_set_id',
    extra_keys => 'data_id,parameter_set_id,orig_node_id'
  };

  return $structure;
}

sub get_sites_table_structure {
  my $self = shift;

  my $structure = {
    data_id               => 'int',
    node_id => 'int',
    orig_node_id => 'int',
    parameter_set_id => 'tinyint',

    omega       => 'float',
    omega_lower => 'float',
    omega_upper => 'float',
    lrt_stat    => 'float',
    ncod        => 'float',
    type        => 'string',
    note        => 'string',

    aln_position          => 'int',
    aln_position_fraction => 'float',

    unique_keys => 'data_id,parameter_set_id,aln_position',
    extra_keys => 'data_id,parameter_set_id,orig_node_id'
  };

  return $structure;
}

1;
