package Bio::Greg::PrimateHIV::CollectPrimateHIVStats;

use strict;
use Time::HiRes qw(sleep);

use Cwd;
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::CollectSitewiseStats');

sub get_gene_table_structure {
  my $self = shift;
  
  my $added_structure = {
    domain_id => 'string',
    alignment_slice => 'string',
    
    human_protein => 'string',
    human_gene    => 'string',
    human_gene_list => 'string',
    human_protein_list => 'string',
    
    chr_name => 'string',
    chr_start => 'int',
    chr_end => 'int',
    chr_strand => 'int',
  };
  
  my $structure = $self->SUPER::get_gene_table_structure();
  $structure = $self->replace_params($structure,$added_structure);
  return $structure;
}

sub get_sites_table_structure {
  my $self = shift;

  my $added_structure = {
    domain_id => 'string',
    filter_value => 'int',
    
    chr_name  => 'string',
    chr_start => 'int',
    chr_end   => 'int',
  };
  
  my $structure = $self->SUPER::get_sites_table_structure;
  $structure = $self->replace_params($structure,$added_structure);
  return $structure;
}

sub data_for_gene {
  my $self = shift;
  
  my $data = $self->SUPER::data_for_gene();  
  my $tree = $self->get_tree;

  # Collect human protein.
  my @human_proteins = grep { $_->taxon_id == 9606 } $tree->leaves;
  my @human_genes    = map  { $_->get_Gene } @human_proteins;
  if ( scalar @human_proteins > 0 ) {
    my $member = $human_proteins[0];
    $data->{'human_protein'} = $member->stable_id;
    $data->{'human_gene'}    = $member->get_Gene->stable_id;
    $data->{'human_protein_list'} = join( ",", map { $_->stable_id } @human_proteins );
    $data->{'human_gene_list'} = join( ",", map { $_->stable_id } @human_genes );
  }

  # Collect protein coords.
  if ( scalar @human_proteins > 0) {
    my $member = $human_proteins[0];
    my $tscr_orig = $member->get_Transcript;
    my $tscr      = $tscr_orig->transform("chromosome");
    if ( defined $tscr ) {
      my $chr    = "chr" . $tscr->slice->seq_region_name;
      my $strand = $tscr->strand;
      my $start  = $tscr->coding_region_start;
      my $end    = $tscr->coding_region_end;
      $data->{chr_name}    = $chr;
      $data->{chr_start}  = $start;
      $data->{chr_end}    = $end;
      $data->{chr_strand} = $strand;
    }
  }

  return $data;
}

sub data_for_site {
  my $self = shift;
  my $aln_position = shift;

  $self->param('genome',1);

  my $data = $self->SUPER::data_for_site($aln_position);
  return undef unless (defined $data);

  my $tag_hash = $self->param('tag_hash');
  if (!defined $tag_hash) {
    $tag_hash = $self->get_tag_hash( $self->compara_dba->dbc, $self->params);
    $self->param('tag_hash',$tag_hash);
  }

  my $site_tags = $tag_hash->{$aln_position};
  #$self->hash_print($site_tags);
  if (defined $site_tags) {
    $data->{'filter_value'} = $site_tags->{'FILTER'};
  }
  return $data;
}

1;