package Bio::Greg::Hive::CollectSubstitutions;

use strict;
use Bio::Greg::Codeml;
use File::Path;
use Bio::Greg::Gorilla::Utils;
use Bio::Align::Utilities qw(:all);

use Time::HiRes qw(sleep);
use DateTime;
use DateTime::Format::MySQL;

use base ( 'Bio::Greg::Hive::Process', 'Bio::Greg::Hive::Align' );

my $TREE = 'Bio::EnsEMBL::Compara::TreeUtils';

sub table_def {
  my $self = shift;

  return {
    data_id => 'int',

    ref_taxon_id   => 'int',
    ref_stable_id  => 'char32',
    ref_seq_name   => 'char32',
    ref_seq_pos    => 'int',
    ref_nucleotide => 'char8',

    other_taxon_id   => 'int',
    other_stable_id  => 'char32',
    other_seq_name   => 'char32',
    other_seq_pos    => 'int',
    other_nucleotide => 'char8',

    substitution_type => 'char32',

    unique_keys => 'ref_seq_name,ref_seq_pos'
  };
}

sub fetch_input {
  my $self           = shift;
  my $default_params = {
    genome_taxon_id => 9606,     # set to 'fan_jobs' to fan out alignment gathering jobs.
    other_taxon_id  => '9598',
  };

  # Fetch parameters from all possible locations.
  $self->load_all_params($default_params);
  $self->create_table_from_params( $self->compara_dba, $self->param('table'), $self->table_def );

  $self->param( 'data_id', $self->param('transcript_stable_id') );
}

sub run {
  my $self = shift;

  my $compara  = $self->compara_dba;
  my $gdba     = $compara->get_GenomeDBAdaptor;
  my $gdb      = $gdba->fetch_by_taxon_id( $self->param('genome_taxon_id') );
  my $core_dba = $gdb->db_adaptor;
  my $gene_a   = $core_dba->get_GeneAdaptor;

}
