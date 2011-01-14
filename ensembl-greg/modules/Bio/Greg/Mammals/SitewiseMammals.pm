package Bio::Greg::Gorilla::LikelihoodTests;

use strict;
use Bio::Greg::Codeml;
use File::Path;

use base (
  'Bio::Greg::Hive::Process', 'Bio::Greg::StatsCollectionUtils',
  'Bio::Greg::Hive::Align',   'Bio::Greg::Hive::CountSubstitutions'
);

my $TREE = 'Bio::EnsEMBL::Compara::TreeUtils';

sub param_defaults {
  return {
    aln_type                            => 'compara',
    genes_table                        => 'genes',
    sites_table                        => 'sites',
    alignment_output                    => 1,
    alignment_output_folder             => 'alns',
    quality_threshold                   => 30,
  };
}

sub fetch_input {
  my $self = shift;

  # Fetch parameters from all possible locations.
  $self->load_all_params();

  # Create tables if necessary.
  $self->create_table_from_params( $self->compara_dba, $self->param('sites_table'),
                                   $self->sites_table_structure );

  $self->create_table_from_params( $self->compara_dba, $self->param('genes_table'),
                                   $self->genes_table_structure );
}

sub sites_table_structure {
  my $self = shift;

  my $structure = {
    data_id => 'int',
    node_id => 'int',
    parameter_set_id => 'tinyint',

    gene_name => 'char32',

    chr_name => 'string',
    chr_start => 'int',
    chr_end => 'int',
    
    aln_position => 'int',
    ncod => 'int',
    omega => 'float',
    omega_lower => 'float',
    omega_upper => 'float',
    lrt_stat => 'float',
    type => 'char16',
    note => 'char16'
  };

  return $structure;
}

sub gene_table_structure {
  my $self = shift;

  my $structure = {
    data_id => 'int',
    node_id => 'int',
    parameter_set_id => 'tinyint',

    gene_name => 'char32',
    
    chr_name => 'string',
    chr_start =>'int',
    chr_end => 'int',

    Hsap_protein => 'char16',
    Hsap_gene => 'char16',
    Hsap_tx => 'char16',
    
    aln_length => 'int',
    aln_file => 'string',
    filtered_sites => 'int',
    filtered_mut_windows => 'int',

    gc_cds => 'float',
    `gc_3` => 'float',
    gc_genomic => 'float',

    n_leaves => 'int',
    n_sites => 'int',
    n_pos_sites => 'int',
    pval => 'float',

    slr_file => 'string',
    tree_file => 'string',
    aln_file => 'string',

    slr_dnds => 'float',
    slr_kappa => 'float',

    unique_keys => 'data_id,node_id,parameter_set_id'
  };
}
