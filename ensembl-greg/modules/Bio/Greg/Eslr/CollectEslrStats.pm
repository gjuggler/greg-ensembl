package Bio::Greg::Eslr::CollectEslrStats;

use strict;
use Time::HiRes qw(sleep);

use Cwd;
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::Process;

use Bio::Greg::ProcessUtils;

our @ISA = qw(Bio::EnsEMBL::Hive::Process Bio::Greg::ProcessUtils);

my $dba;
my $pta;

my $tree;
my $params;

my $results;

my $gene_stats_def = {
  tree_length_total     => 'float',
  tree_length_max       => 'float',
  tree_length_avg       => 'float',
  num_leaves            => 'int',
  avg_column_entropy    => 'float',
  avg_seq_length        => 'float',
    
  num_human_genes       => 'int',
  human_gene            => 'string',
  human_genes           => 'string',
  human_chr             => 'string',
  human_str             => 'string'

  duplications          => 'int',
  duplication_fraction  => 'float',

  num_pscs              => 'int',
  num_pscs_weak         => 'int',
  mean_omega            => 'float',

  omega                 => 'float',
  kappa                 => 'float',
  slr_lnl               => 'float',

  ins_rate              => 'float',
  del_rate              => 'float'
  };

my $site_stats_def = {
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
    $self->get_gene_data($node_id,$ps_id);
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
  $cur_params->{'tree_length_total'} = total_distance($tree);
  $cur_params->{'tree_length_max'} = max_distance($tree);
  $cur_params->{'tree_length_avg'} = avg_distance($tree);
  $cur_params->{'num_leaves'} = scalar($tree->leaves);
  $cur_params->{'avg_column_entropy'} = sprintf("%.3f",Bio::EnsEMBL::Compara::AlignUtils->average_column_entropy($cdna_nogap));
  $cur_params->{'avg_seq_length'} = avg_seq_length($tree);

  my @hum_gen = grep {$_->taxon_id==9606} $tree->leaves;
  $cur_params->{'num_human_genes'} = scalar(@hum_gen);
  if (scalar @hum_gen > 0) {
    $cur_params->{'human_genes'} = join(",",map {$_->stable_id} @hum_gen);
    my $first_human_gene = $hum_gen[0];
    $cur_params->{'human_gene'} = $first_human_gene->stable_id;

    my $tscr = $first_human_gene->get_Transcript;
    $tscr = $tscr->transform("chromosome");
    if (defined $tscr) {
      my $chr = "chr".$tscr->slice->seq_region_name;
    }

  }

  $cur_params->{'duplications'} = mysql_getval($tree,"SELECT num_dups_under_node($node_id)");
  $cur_params->{'duplication_fraction'} = sprintf "%.3f", mysql_getval($tree,"SELECT num_dups_under_node($node_id)/node_count($node_id)");

  my $psc_hash = get_psc_hash($dba->dbc,$cur_params);
  $cur_params->{'num_pscs'} = psc_count($psc_hash,0);
  $cur_params->{'num_pscs_weak'} = psc_count($psc_hash,1);
  $cur_params->{'mean_omega'} = omega_average($psc_hash);

  my $ps = $cur_params->{'parameter_set_id'};
  $cur_params->{'omega'} = $cur_params->{'slr_omega_'.$ps};
  $cur_params->{'kappa'} = $cur_params->{'slr_kappa_'.$ps};
  $cur_params->{'lnl'} = $cur_params->{'slr_lnL_'.$ps};

#  my ($ins,$del,$ins_rate,$del_rate) = Bio::EnsEMBL::Compara::AlignUtils->indelign($sa_nogap,$tree,$cur_params);
  $cur_params->{'indelign_ins_rate'} = $cur_params->{'indelign_ins_rate_'.$ps};
  $cur_params->{'indelign_del_rate'} = $cur_params->{'indelign_del_rate_'.$ps};

  # Store values in our output table.
  my $table = $cur_params->{'collect_eslr_stats_genes_table'};
  $self->store_params_in_table($dba,$table,$cur_params);
}


sub avg_seq_length {
  my $tree = shift;
  
  my $seq_len = 0;
  map {$seq_len += $_->seq_length} $tree->leaves;
  return sprintf "%.3f", $seq_len / scalar($tree->leaves);
}

sub avg_distance {
  my $tree = shift;
  my $dist = 0;
  map {$dist += dist_to_root($_);} $tree->leaves;
  $dist = $dist / scalar($tree->leaves);
  return sprintf "%.3f", $dist;
}

sub dist_to_root {
  my $leaf = shift;

  my $d = $leaf->distance_to_parent;
  my $p = $leaf->parent;
  while ($p) {
    $d += $p->distance_to_parent;
    $p = $p->parent;
  }
  return $d;
}

sub max_distance {
  my $tree = shift;
  my $max = $tree->max_distance;
  return sprintf "%.3f", $max;
}

sub total_distance {
  my $tree = shift;
  my $dist = 0;
  map {$dist += $_->distance_to_parent} $tree->nodes;
  return sprintf "%.3f", $dist;
}

sub mysql_getval {
  my $tree = shift;
  my $cmd = shift;

  my $dbc = $tree->adaptor->dbc;
  my $sth = $dbc->prepare($cmd);
  $sth->execute();
  my $val = @{$sth->fetchrow_arrayref()}[0];
  $val = 'NA' unless (defined $val);
  $sth->finish;
  return $val;
}

sub get_psc_hash {
  my $dbc = shift;
  my $params = shift;

  my $table = $params->{'omega_table'};
  my $pset = $params->{'parameter_set_id'};
  my $node_id = $params->{'node_id'};
  my $cmd = qq^SELECT aln_position,omega,omega_lower,omega_upper,lrt_stat,ncod,type,note 
    FROM $table WHERE parameter_set_id=$pset and node_id=$node_id
    AND ncod > 4 AND note != 'random' AND omega_upper > omega^;
  print "$cmd\n";
  my $sth = $dbc->prepare($cmd);
  $sth->execute;
  my $obj = $sth->fetchall_hashref('aln_position');
  $sth->finish;
  return $obj;
}

sub psc_count {
  my $psc_hash = shift;
  my $include_weak_pscs = shift;

  my @obj_array = map {$psc_hash->{$_}} keys %$psc_hash;
  my @psc_objs;
  if ($include_weak_pscs) {
    @psc_objs = grep {$_->{type} =~ /positive[1234]/} @obj_array;
  } else {  
    @psc_objs = grep {$_->{type} =~ /positive[34]/} @obj_array;
  }

  return scalar(@psc_objs);
}


sub omega_average {
  my $hash = shift;

  my @obj_array = map {$hash->{$_}} keys %$hash;
  
  my $omega_total = 0;
  map {$omega_total += $_->{omega}} @obj_array;

  return 'NA' if (scalar @obj_array == 0);
  return sprintf "%.3f", $omega_total/scalar(@obj_array);
}

1;
