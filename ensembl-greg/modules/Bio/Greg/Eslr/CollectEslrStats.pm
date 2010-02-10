package Bio::Greg::Eslr::CollectEslrStats;

use strict;
use Time::HiRes qw(sleep);

use Cwd;
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::Process;

use Bio::Greg::ProcessUtils;
use Bio::Greg::StatsCollectionUtils;

our @ISA = qw(Bio::EnsEMBL::Hive::Process Bio::Greg::ProcessUtils);

my $utils = "Bio::Greg::StatsCollectionUtils";

my $dba;
my $pta;

my $tree;
my $params;

my $results;

my $gene_stats_def = {
  tree_length           => 'float',
  tree_max_path         => 'float',
  tree_mean_path        => 'float',
  tree_max_branch       => 'float',
  tree_mean_branch      => 'float',
  leaf_count            => 'int',

  column_entropy_mean   => 'float',
  seq_length_mean       => 'float',
  gc_content_mean       => 'float',
    
  human_gene_count      => 'int',
  human_gene_list       => 'string',

  human_gene            => 'string',
  human_chr             => 'string',
  human_start           => 'int',
  human_end             => 'int',
  human_strand          => 'int',

  duplication_count     => 'int',
  duplication_fraction  => 'float',

  psc_count             => 'int',
  weak_psc_count        => 'int',
  omega_mean            => 'float',
  omega_mean_excl_pscs  => 'float',
  omega_m0              => 'float',
  kappa                 => 'float',
  slr_lnl               => 'float',

# Coming soon...
#  ins_rate              => 'float',
#  del_rate              => 'float'
  };

my $site_stats_def = {
  omega                 => 'float',
  omega_lower           => 'float',
  omega_upper           => 'float',
  lrt_stat              => 'float',
  ncod                  => 'float',
  type                  => 'string',
  note                  => 'string',

  aln_position          => 'int',
  aln_position_fraction => 'float',

  chr_name              => 'string',
  chr_start             => 'int',
  chr_end               => 'int',


# Should we integrate these into the table now, or let this be done in R?
#  column_entropy        => 'float',
#  indel_count           => 'int',
#  pfam_domain           => 'string',
#  sec_structure         => 'string',
#  solvent_acc           => 'int',
#  exon_type             => 'string',
#  splice_dist           => 'int',

  unique_keys            => 'aln_position,node_id,parameter_set_id'  
};


sub fetch_input {
  my ($self) = @_;

  $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-DBCONN=>$self->db->dbc);
  $pta = $dba->get_ProteinTreeAdaptor;

  $params = {
    omega_table => 'sitewise_omega',
    collect_eslr_stats_parameter_sets => 'all',
    collect_eslr_stats_sites_table => 'stats_sites',
    collect_eslr_stats_genes_table => 'stats_genes',
    collect_eslr_stats_pscs_table => 'stats_pscs'
  };

  my $p_params = $self->get_params($self->parameters);
  my $i_params = $self->get_params($self->input_id);
  my $node_id = $i_params->{'protein_tree_id'} || $i_params->{'node_id'};
  my $t_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_tree_tags($dba,$node_id);

  $params = $self->replace_params($params,$p_params,$i_params,$t_params);
  Bio::EnsEMBL::Compara::ComparaUtils->hash_print($params);

  # Create tables if necessary.
  $self->create_table_from_params($dba,$params->{'collect_eslr_stats_genes_table'},$gene_stats_def);
  $self->create_table_from_params($dba,$params->{'collect_eslr_stats_sites_table'},$site_stats_def);

}

sub run {
  my $self = shift;

  my $node_id = $params->{'node_id'};
  my $param_set_string = $params->{'collect_eslr_stats_parameter_sets'};

  my @param_sets;
  if ($param_set_string eq 'all') {
    my $query = qq^select distinct(parameter_set_id) FROM parameter_set order by parameter_set_id;^;
    @param_sets = @{$dba->dbc->db_handle->selectcol_arrayref($query)};
  } else {
    @param_sets = split(",",$param_set_string);
  }

  foreach my $ps_id (@param_sets) {
#    $self->get_gene_data($node_id,$ps_id);
  }

  foreach my $ps_id (@param_sets) {
    $self->get_sites_data($node_id,$ps_id);
  }

}

sub get_gene_data {
  my $self = shift;
  my $node_id = shift;
  my $parameter_set_id = shift;

  my $param_set_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_param_set($dba->dbc,$parameter_set_id);
  my $cur_params = $self->replace_params($params,$param_set_params);
  Bio::EnsEMBL::Compara::ComparaUtils->hash_print($cur_params);

  $pta->protein_tree_member("protein_tree_member");
  
  my ($tree,$sa_aln,$cdna_aln);
  eval {
    ($tree,$sa_aln,$cdna_aln) = Bio::EnsEMBL::Compara::ComparaUtils->tree_aln_cdna($dba,$cur_params);
  };
  return if ($@);

  print $tree->newick_format."\n";

  my ($cdna_nogap,$new_to_old,$old_to_new) = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns($cdna_aln);
  my ($sa_nogap,$new_to_old,$old_to_new) = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns($sa_aln);

  #Bio::EnsEMBL::Compara::AlignUtils->pretty_print($sa_nogap);
  
  # Collect gene tag values into the params hash.
  $cur_params->{'tree_length'} = $utils->tree_length($tree);
  $cur_params->{'tree_max_path'} = $utils->max_path($tree);
  $cur_params->{'tree_mean_path'} = $utils->mean_path($tree);
  $cur_params->{'tree_max_branch'} = $utils->max_branch($tree);
  $cur_params->{'tree_mean_branch'} = $utils->mean_branch($tree);
  $cur_params->{'leaf_count'} = scalar($tree->leaves);

  $cur_params->{'column_entropy_mean'} = sprintf("%.3f",Bio::EnsEMBL::Compara::AlignUtils->average_column_entropy($cdna_nogap));
  $cur_params->{'seq_length_mean'} = $utils->seq_length_mean($tree);
  $cur_params->{'gc_content_mean'} = $utils->gc_content_mean($tree);

  my @hum_gen = grep {$_->taxon_id==9606} $tree->leaves;
  $cur_params->{'human_gene_count'} = scalar(@hum_gen);
  if (scalar @hum_gen > 0) {
    my $first_human_gene = $hum_gen[0];
    $cur_params->{'human_gene_list'} = join(",",map {$_->stable_id} @hum_gen);
    $cur_params->{'human_gene'} = $first_human_gene->stable_id;

    my $tscr = $first_human_gene->get_Transcript;
    $tscr = $tscr->transform("chromosome");
    if (defined $tscr) {
      my $chr = "chr".$tscr->slice->seq_region_name;
      my $strand = $tscr->slice->strand;
      my $start = $tscr->coding_region_start;
      my $end = $tscr->coding_region_end;
      $cur_params->{human_chr} = $chr;
      $cur_params->{human_strand} = $strand;
      $cur_params->{human_start} = $start;
      $cur_params->{human_end} = $end;
    }
  }

  $cur_params->{'duplication_count'} = $utils->mysql_getval($tree,"SELECT num_dups_under_node($node_id)");
  $cur_params->{'duplication_fraction'} = sprintf "%.3f", $utils->mysql_getval($tree,"SELECT num_dups_under_node($node_id)/node_count($node_id)");

  my $psc_hash = $utils->get_psc_hash($dba->dbc,$cur_params);
  $cur_params->{'psc_count'} = $utils->psc_count($psc_hash,0);
  $cur_params->{'weak_psc_count'} = $utils->psc_count($psc_hash,1);
  $cur_params->{'omega_mean'} = $utils->omega_average($psc_hash);
  $cur_params->{'omega_mean_excl_pscs'} = $utils->omega_average_exclude_pscs($psc_hash);
 
  my $ps = $cur_params->{'parameter_set_id'};
  $cur_params->{'omega_m0'} = $cur_params->{'slr_omega_'.$ps};
  $cur_params->{'kappa'} = $cur_params->{'slr_kappa_'.$ps};
  $cur_params->{'slr_lnl'} = $cur_params->{'slr_lnL_'.$ps};

# Indel rates aren't yet implemented (need to speed up Indelign first.)
#  my ($ins,$del,$ins_rate,$del_rate) = Bio::EnsEMBL::Compara::AlignUtils->indelign($sa_nogap,$tree,$cur_params);
#  $cur_params->{'indelign_ins_rate'} = $cur_params->{'indelign_ins_rate_'.$ps};
#  $cur_params->{'indelign_del_rate'} = $cur_params->{'indelign_del_rate_'.$ps};

  # Store values in our output table.
  my $table = $cur_params->{'collect_eslr_stats_genes_table'};
  $self->store_params_in_table($dba,$table,$cur_params);
}


sub get_sites_data {
  my $self = shift;
  my $node_id = shift;
  my $parameter_set_id = shift;

  my $param_set_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_param_set($dba->dbc,$parameter_set_id);
  my $cur_params = $self->replace_params($params,$param_set_params);
  Bio::EnsEMBL::Compara::ComparaUtils->hash_print($cur_params);

  $pta->protein_tree_member("protein_tree_member");

  my ($tree,$sa_aln,$cdna_aln);
  eval {
    ($tree,$sa_aln,$cdna_aln) = Bio::EnsEMBL::Compara::ComparaUtils->tree_aln_cdna($dba,$cur_params);
  };
  return if ($@);
  
  my $psc_params = $self->replace_params($cur_params,{genome => 1});
  my $psc_hash = $utils->get_psc_hash($dba->dbc,$psc_params);

  my $aln_length = $sa_aln->length;

  foreach my $key (sort {$a <=> $b} keys %$psc_hash) {
    my $site = $psc_hash->{$key};
    
    $cur_params->{omega} = $site->{omega};
    $cur_params->{omega_lower} = $site->{omega_lower};
    $cur_params->{omega_upper} = $site->{omega_upper};
    $cur_params->{lrt_stat} = $site->{lrt_stat};
    $cur_params->{ncod} = $site->{ncod};
    $cur_params->{type} = $site->{type};
    $cur_params->{note} = $site->{note};
    
    $cur_params->{'aln_position'} = $site->{aln_position};
    $cur_params->{'aln_position_fraction'} = $site->{aln_position} / $aln_length;

    $cur_params->{'chr_name'} = $site->{chr_name};
    $cur_params->{'chr_start'} = $site->{chr_start};
    $cur_params->{'chr_end'} = $site->{chr_end};

    # Store values in our output table.
    my $table = $cur_params->{'collect_eslr_stats_sites_table'};
    $self->store_params_in_table($dba,$table,$cur_params);
  }
}

1;
