package Bio::Greg::Mammals::CollectMammalsStats;

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

  print "Hey!\n";
  my $check = $self->check_tree_aln;
  if ($check == -1) {
    $self->store_tag("no_stats","Tree or align doesn't look good!");
    return;
  }

  $self->SUPER::run;
}

sub get_sites_table_structure {
  my $self = shift;

  my $added_structure = {
    filter_value => 'tinyint',
    domain => 'char32',
    
    chr_name  => 'char16',
    chr_start => 'int',
    chr_end   => 'int',

  };
  
  my $structure = $self->SUPER::get_sites_table_structure;
  $structure = $self->replace_params($structure,$added_structure);
  return $structure;
}

sub get_gene_table_structure {
  my $self = shift;
  
  my $added_structure = {
    'human_protein' => 'char32',
    'human_gene'    => 'char32',
    'human_gene_list' => 'string',
    'human_protein_list' => 'string',
    'human_gene_count' => 'smallint',

    'chr_name' => 'char16',
    'chr_start' => 'int',
    'chr_end' => 'int',
    'chr_strand' => 'tinyint',
    };
  
  my $structure = $self->SUPER::get_gene_table_structure();
  $structure = $self->replace_params($structure,$added_structure);
  return $structure;
}

sub data_for_gene {
  my $self = shift;

  $self->param('genome',1);
  $self->param('filtered',1);
  $self->param('alignment_filtering_value',3);

  my $data = $self->SUPER::data_for_gene();
  return undef unless (defined $data);
  my $tree = $self->get_tree;

  # Collect human protein.
  my @human_proteins = grep { $_->taxon_id == 9606 } $tree->leaves;
  my @human_genes    = map  { $_->gene_member } @human_proteins;
  if ( scalar @human_proteins > 0 ) {
    my $member = $human_proteins[0];
    $data->{'human_protein'} = $member->stable_id;
    $data->{'human_gene'}    = $human_genes[0]->stable_id;
    $data->{'human_protein_list'} = join( ",", map { $_->stable_id } @human_proteins );
    $data->{'human_gene_list'} = join( ",", map { $_->stable_id } @human_genes );
  }
  $data->{human_gene_count} = scalar(@human_proteins);

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
  $self->param('filtered',0);
  $self->param('alignment_filtering_value',0);

  my $data;
  my $data = $self->SUPER::data_for_site($aln_position);
  return undef unless (defined $data);

  my $tag_hash = $self->param('tag_hash');
  if (!defined $tag_hash) {
    $tag_hash = $self->get_tag_hash( $self->compara_dba->dbc, $self->params);
    $self->param('tag_hash',$tag_hash);
  }

  my $site_tags = $tag_hash->{$aln_position};
#  $self->hash_print($site_tags);
  if (defined $site_tags) {
    $data->{'filter_value'} = $site_tags->{'FILTER'};
    $data->{'domain'} = $site_tags->{'DOMAIN'};
#    print "DOMAIN: ".$site_tags->{'DOMAIN'}."\n";
  }
  return $data;
}

1;
