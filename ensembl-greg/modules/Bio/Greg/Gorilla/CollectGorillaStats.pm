package Bio::Greg::Gorilla::CollectGorillaStats;

use strict;
use Time::HiRes qw(sleep);

use Cwd;
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::Process','Bio::Greg::StatsCollectionUtils');

my $gene_stats_def = {
  parameter_set_name => 'string',
  
  tree_length      => 'float',
  tree_max_path    => 'float',
  tree_mean_path    => 'float',
  leaf_count       => 'int',
  tree_newick => 'string',

  'gor.tree.pattern' => 'string',
  'gor.synonymous' => 'int',
  'gor.nonsynonymous' => 'int',

  seq_length_mean     => 'float',
  gc_content_mean     => 'float',
  cpg_obs_exp_mean    => 'float',

  human_protein => 'string',
  human_gene    => 'string',
  gor_gene => 'string',
  gor_protein => 'string',

  chr_name => 'string',
  start => 'int',
  end => 'int',
  strand => 'int',
  
  duplication_count => 'int',
  duplication_fraction => 'int',

  psc_count            => 'int',
  weak_psc_count       => 'int',
  
  'hyphy.dnds' => 'float',
  'hyphy.dnds.lo' => 'float',
  'hyphy.dnds.hi' => 'float',
  'slr.dnds' => 'float',
  
  mean_copy_count => 'float',
};

my $site_stats_def = {
  omega       => 'float',
  omega_lower => 'float',
  omega_upper => 'float',
  lrt_stat    => 'float',
  ncod        => 'float',
  type        => 'string',
  note        => 'string',

  aln_position          => 'int',
  aln_position_fraction => 'float',

  chr_name  => 'string',
  chr_start => 'int',
  chr_end   => 'int',

  unique_keys => 'aln_position,node_id,parameter_set_id'
};

sub fetch_input {
  my ($self) = @_;

  my $params = {
    omega_table                       => 'sitewise_omega',
    collect_gor_stats_parameter_sets => 'all',
    collect_gor_stats_sites_table    => 'stats_sites',
    collect_gor_stats_genes_table    => 'stats_genes',
    collect_gor_stats_pscs_table     => 'stats_pscs'
  };

  $self->load_all_params($params);


  # Create tables if necessary.
  $self->create_table_from_params( $self->compara_dba, $self->param('collect_gor_stats_genes_table'),
    $gene_stats_def );
  $self->create_table_from_params( $self->compara_dba, $self->param('collect_gor_stats_sites_table'),
    $site_stats_def );
}

sub run {
  my $self = shift;

  my $node_id          = $self->param('node_id');
  my $parameter_set_id = $self->param('parameter_set_id');

  $self->get_gene_data( $node_id, $parameter_set_id);
  $self->get_sites_data( $node_id, $parameter_set_id);
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

  my ( $cdna_nogap, $new_to_old, $old_to_new ) =
    Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns($cdna_aln);
  my ( $sa_nogap, $new_to_old, $old_to_new ) =
    Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns($sa_aln);

  #Bio::EnsEMBL::Compara::AlignUtils->pretty_print($sa_nogap);

  # Collect gene tag values into the params hash.
  $cur_params->{'tree_length'}      = $self->tree_length($tree);
  $cur_params->{'tree_max_path'}    = $self->max_path($tree);
  $cur_params->{'tree_mean_path'}   = $self->mean_path($tree);
  $cur_params->{'leaf_count'}       = scalar( $tree->leaves );

  $cur_params->{'tree_newick'} = $tree->newick_format;

  # Stuff coming from the counts_genes table.
  my $format = "SELECT %s FROM counts_genes WHERE node_id=$node_id";
  $cur_params->{'gor.tree.pattern'} = $self->mysql_getval($tree,sprintf $format,'tree_pattern');
  $cur_params->{'gor.synonymous'} = $self->mysql_getval($tree,sprintf $format,'synon');
  $cur_params->{'gor.nonsynonymous'} = $self->mysql_getval($tree,sprintf $format,'nonsynon');
  

  $cur_params->{'seq_length_mean'}  = $self->seq_length_mean($tree);
  $cur_params->{'gc_content_mean'}  = $self->gc_content_mean($tree);
  $cur_params->{'cpg_obs_exp_mean'} = $self->cpg_obs_exp_mean($cdna_aln);

  # Collect human protein.
  my @human_proteins = grep { $_->taxon_id == 9606 } $tree->leaves;
  my @human_genes    = map  { $_->get_Gene } @human_proteins;
  if ( scalar @human_proteins > 0 ) {
    my $member = $human_proteins[0];
    $cur_params->{'human_protein'} = $member->stable_id;
    $cur_params->{'human_gene'}    = $member->get_Gene->stable_id;
  }

  # Collect gorilla protein.
  my @gor_proteins = grep { $_->taxon_id == 9593 } $tree->leaves;
  my @gor_genes    = map  { $_->get_Gene } @gor_proteins;
  if ( scalar @gor_proteins > 0) {
    my $member = $gor_proteins[0];
    $cur_params->{'gor_protein'} = $member->stable_id;
    $cur_params->{'gor_gene'}    = $member->get_Gene->stable_id;
  }


  # Collect gorilla protein coords.
  if ( scalar @gor_proteins > 0) {
    my $member = $gor_proteins[0];
    my $tscr_orig = $member->get_Transcript;
    my $tscr      = $tscr_orig->transform("chromosome");
    if ( defined $tscr ) {
      my $chr    = "chr" . $tscr->slice->seq_region_name;
      my $strand = $tscr->strand;
      my $start  = $tscr->coding_region_start;
      my $end    = $tscr->coding_region_end;
      $cur_params->{chr_name}    = $chr;
      $cur_params->{start}  = $start;
      $cur_params->{end}    = $end;
      $cur_params->{strand} = $strand;
    }
  }

  $cur_params->{'mean_copy_count'} = $self->mean_copy_count($tree);

  $cur_params->{'duplication_count'} = $self->duplication_count($tree);
  $cur_params->{'duplication_fraction'} = $self->duplication_fraction($tree);

  my $psc_hash = $self->get_psc_hash( $self->compara_dba->dbc, $cur_params );
  if ( scalar( keys(%$psc_hash) ) > 0 ) {
    $cur_params->{'psc_count'}      = $self->psc_count( $psc_hash, 0 );
    $cur_params->{'weak_psc_count'} = $self->psc_count( $psc_hash, 1 );
  }

  my $ps = $cur_params->{'parameter_set_id'};

  $cur_params->{'slr.dnds'} = $cur_params->{ 'slr_omega_' . $ps };
  $cur_params->{'hyphy.dnds'} = $cur_params->{ 'hyphy_omega_' . $ps };
  $cur_params->{'hyphy.dnds.lo'} = $cur_params->{ 'hyphy_omega_lo_' . $ps };
  $cur_params->{'hyphy.dnds.hi'} = $cur_params->{ 'hyphy_omega_hi_' . $ps };

  # Store values in our output table.
  my $table = $cur_params->{'collect_gor_stats_genes_table'};
  $self->store_params_in_table( $self->compara_dba, $table, $cur_params );

  $tree->release_tree;
}

sub get_sites_data {
  my $self             = shift;
  my $node_id          = shift;
  my $parameter_set_id = shift;

  my $cur_params = $self->params;
  Bio::EnsEMBL::Compara::ComparaUtils->hash_print($cur_params);

  my $pta = $self->compara_dba->get_ProteinTreeAdaptor;
  $pta->protein_tree_member("protein_tree_member");
  my ( $tree, $sa_aln, $cdna_aln );
  eval {
    ( $tree, $sa_aln, $cdna_aln ) =
      Bio::EnsEMBL::Compara::ComparaUtils->tree_aln_cdna( $self->compara_dba, $cur_params );
  };
  if ($@) {
    print "Sites data collection error!\n";
    return;
  }

  my $psc_params = $self->replace_params( $cur_params, { 'genome' => 1 } );
  my $psc_hash = $self->get_psc_hash( $self->compara_dba->dbc, $psc_params );
  my $tag_hash = $self->get_tag_hash( $self->compara_dba->dbc, $psc_params );

  my $aln_length = $sa_aln->length;
  foreach my $aln_position ( sort { $a <=> $b } keys %$psc_hash ) {
    my $site = $psc_hash->{$aln_position};

    $cur_params->{omega}       = $site->{omega};
    $cur_params->{omega_lower} = $site->{omega_lower};
    $cur_params->{omega_upper} = $site->{omega_upper};
    $cur_params->{lrt_stat}    = $site->{lrt_stat};
    $cur_params->{ncod}        = $site->{ncod};
    $cur_params->{type}        = $site->{type};
    $cur_params->{note}        = $site->{note};

    $cur_params->{'aln_position'}          = $site->{aln_position};
    $cur_params->{'aln_position_fraction'} = $site->{aln_position} / $aln_length;

    $cur_params->{'chr_name'}  = $site->{chr_name};
    $cur_params->{'chr_start'} = $site->{chr_start};
    $cur_params->{'chr_end'}   = $site->{chr_end};

    # Store values in our output table.
    my $table = $cur_params->{'collect_gor_stats_sites_table'};
    $self->store_params_in_table( $self->compara_dba, $table, $cur_params );
  }
  $tree->release_tree;
}

1;
