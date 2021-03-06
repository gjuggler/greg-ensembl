package Bio::Greg::Iso::CollectIsoStats;

use strict;
use Time::HiRes qw(sleep);

use Cwd;
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::CollectSitewiseStats');

sub run {
  my $self = shift;

  my $check = $self->check_tree_aln;
  if ( $check == -1 ) {
    $self->store_tag( "no_stats", "Tree or align doesn't look good!" );
    return;
  }

  $self->SUPER::run;
}

sub get_sites_table_structure {
  my $self = shift;

  my $added_structure = { unique_keys => 'aln_position,node_id,parameter_set_id' };

  my $structure = $self->SUPER::get_sites_table_structure;
  $structure = $self->replace_params( $structure, $added_structure );
  return $structure;
}

sub get_gene_table_structure {
  my $self = shift;

  my $added_structure = {
    'orig_node_id'  => 'int',
    'human_protein' => 'string',
    'human_gene'    => 'string',

    'chr_name'   => 'string',
    'chr_start'  => 'int',
    'chr_end'    => 'int',
    'chr_strand' => 'int',
  };

  my $structure = $self->SUPER::get_gene_table_structure();
  $structure = $self->replace_params( $structure, $added_structure );
  return $structure;
}

sub data_for_gene {
  my $self = shift;

  my $data = $self->SUPER::data_for_gene();
  return undef unless ( defined $data );
  my $tree = $self->get_tree;

  # Collect human protein.
  my @human_proteins = grep { $_->taxon_id == 9606 } $tree->leaves;
  my @human_genes    = map  { $_->get_Gene } @human_proteins;
  if ( scalar @human_proteins > 0 ) {
    my $member = $human_proteins[0];
    $data->{'human_protein'}      = $member->stable_id;
    $data->{'human_gene'}         = $member->get_Gene->stable_id;
    $data->{'human_protein_list'} = join( ",", map { $_->stable_id } @human_proteins );
    $data->{'human_gene_list'}    = join( ",", map { $_->stable_id } @human_genes );
  }

  # Collect protein coords.
  if ( scalar @human_proteins > 0 ) {
    my $member    = $human_proteins[0];
    my $tscr_orig = $member->get_Transcript;
    my $tscr      = $tscr_orig->transform("chromosome");
    if ( defined $tscr ) {
      my $chr    = "chr" . $tscr->slice->seq_region_name;
      my $strand = $tscr->strand;
      my $start  = $tscr->coding_region_start;
      my $end    = $tscr->coding_region_end;
      $data->{chr_name} = $chr;
      $data->{start}    = $start;
      $data->{end}      = $end;
      $data->{strand}   = $strand;
    }
  }

  return $data;
}

sub data_for_site {
  my $self         = shift;
  my $aln_position = shift;

  $self->param( 'genome', 1 );

  my $data = $self->SUPER::data_for_site($aln_position);
  return undef unless ( defined $data );

  return $data;
}

1;
