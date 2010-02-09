package Bio::Greg::Slrsim::CollectStats;

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

my $slrsim_stats_def = {
  aln_position       => 'int',

  slrsim_scheme_name => 'string',
  alignment_name     => 'string',
  filtering_name     => 'string',
  sitewise_name      => 'string',
  species_name       => 'string',

  slrsim_rep         => 'int',
  slrsim_file        => 'string',
  slrsim_ref         => 'string',
  slrsim_tree_length => 'float',
  
  phylosim_simulation_program => 'string',
  phylosim_seq_length => 'int',
  phylosim_omega_distribution => 'string',
  phylosim_ins_rate  => 'float',
  phylosim_del_rate  => 'float',

  parameter_set_name => 'string',

  tree_length        => 'float',
  tree_max_branch    => 'float',
  tree_mean_branch  => 'float',
  tree_max_path      => 'float',
  tree_mean_path     => 'float',
  num_leaves         => 'int',

  true_dnds          => 'float',
  true_type          => 'string',
  true_e             => 'float',
  aln_dnds           => 'float',
  aln_type           => 'string',
  aln_e              => 'float',
  ncod               => 'int',
  lrt                => 'float',
  paml_lrt           => 'float',
  slr_lnL            => 'float',
  slr_omega          => 'float',

  sum_of_pairs_score => 'float',
  total_column_score => 'float',
  mean_column_entropy => 'float',

  unique_keys        => 'aln_position,node_id,parameter_set_id'
};

sub fetch_input {
  my ($self) = @_;

  $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-DBCONN=>$self->db->dbc);
  $pta = $dba->get_ProteinTreeAdaptor;

  $params = {
    collect_slrsim_stats_parameter_sets => 'all',
    collect_slrsim_stats_table => 'stats_slrsim'
  };

  my $p_params = $self->get_params($self->parameters);
  my $i_params = $self->get_params($self->input_id);
  my $node_id = $i_params->{'protein_tree_id'} || $i_params->{'node_id'};
  my $t_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_tree_tags($dba,$node_id);

  $params = $self->replace_params($params,$p_params,$i_params,$t_params);
  Bio::EnsEMBL::Compara::ComparaUtils->hash_print($params);

  # Create tables if necessary.
  $self->create_table_from_params($dba,$params->{'collect_slrsim_stats_table'},$slrsim_stats_def);
}

sub run {
  my $self = shift;

  my $node_id = $params->{'node_id'};
  my $param_set_string = $params->{collect_slrsim_stats_parameter_sets};

  print "PSS:".$param_set_string."\n";

  my @param_sets;
  if ($param_set_string eq 'all') {
    my $query = qq^select distinct(parameter_set_id) FROM parameter_set order by parameter_set_id;^;
    @param_sets = @{$dba->dbc->db_handle->selectcol_arrayref($query)};
  } else {
    @param_sets = split(",",$param_set_string);
  }

  foreach my $ps_id (@param_sets) {
    $self->get_data_for_node($node_id,$ps_id);
  }
}

sub get_data_for_node {
  my $self = shift;
  my $node_id = shift;
  my $parameter_set_id = shift;

  my $param_set_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_param_set($dba->dbc,$parameter_set_id);
  my $cur_params = $self->replace_params($params,$param_set_params);
  Bio::EnsEMBL::Compara::ComparaUtils->hash_print($cur_params);

  my $sa_true;
  my $sa_aln;
  my @true_entropies;
  my @aln_entropies;
  my $sum_of_pairs_score;
  my $total_column_score;
  my $mean_column_entropy;
  my $tree;
  eval {
    $pta->protein_tree_member("protein_tree_member");
    $tree = $pta->fetch_node_by_node_id($node_id);
    $sa_true = $tree->get_SimpleAlign();
    my $cdna_aln = $tree->get_SimpleAlign(-cdna => 1);
    @true_entropies = Bio::EnsEMBL::Compara::AlignUtils->column_entropies($cdna_aln);

    $pta->protein_tree_member($cur_params->{'alignment_table'});
    $tree = $pta->fetch_node_by_node_id($node_id);
    $sa_aln = $tree->get_SimpleAlign();
    $cdna_aln = $tree->get_SimpleAlign(-cdna => 1);
    @aln_entropies = Bio::EnsEMBL::Compara::AlignUtils->column_entropies($cdna_aln);

    $sum_of_pairs_score = Bio::EnsEMBL::Compara::AlignUtils->sum_of_pairs_score($sa_true,$sa_aln);
    $total_column_score = Bio::EnsEMBL::Compara::AlignUtils->total_column_score($sa_true,$sa_aln);
    $mean_column_entropy = Bio::EnsEMBL::Compara::AlignUtils->total_column_score($sa_true,$sa_aln);
  };  
  die("Hold up: ".$@) if ($@);
  return if (!$sa_true || !$sa_aln);
  
  # Get the sequence to act as a reference in site-wise value comparisons.
  my $reference_id = '';
  if (defined $cur_params->{'slrsim_ref'}) {
    $reference_id = $cur_params->{'slrsim_ref'};
  }
  my @seqs = $sa_true->each_seq;
  my ($ref_seq) = grep {$_->id eq $reference_id} @seqs;
  #die ("Reference was defined in params but not found in aln!") if ($reference_id ne '' && !defined $ref_seq);
  $ref_seq = $seqs[0] if (!defined $ref_seq);
  my $ref_name = $ref_seq->id;
  my $str = $ref_seq->seq;
  my $nogaps = $str;
  $nogaps =~ s/-//g;

  # Calculate branch length stats.
  $cur_params->{tree_length} = $utils->tree_length($tree);
  $cur_params->{tree_max_branch} = $utils->max_branch($tree);
  $cur_params->{tree_mean_branch} = $utils->mean_branch($tree);
  $cur_params->{tree_max_path} = $utils->max_path($tree);
  $cur_params->{tree_mean_path} = $utils->mean_path($tree);
  $cur_params->{num_leaves} = scalar($tree->leaves);

  $cur_params->{sum_of_pairs_score} = $sum_of_pairs_score;
  $cur_params->{total_column_score} = $total_column_score;
  $cur_params->{mean_column_entropy} = $mean_column_entropy;

  # Get all the site-wise data from the omega table.
  my $aln_table_name = $cur_params->{'omega_table'};
  my $sth1 = $pta->prepare("SELECT aln_position,omega,type FROM sitewise_omega WHERE node_id=?;");
  my $sth2 = $pta->prepare("SELECT aln_position,omega,type,note,ncod,lrt_stat FROM $aln_table_name WHERE node_id=? AND parameter_set_id=?;");
  $sth1->execute($node_id);
  $sth2->execute($node_id,$parameter_set_id);
  my $true_omegas = $sth1->fetchall_hashref('aln_position');
  my $aln_omegas = $sth2->fetchall_hashref('aln_position');

  for (my $i=1; $i <= length($nogaps); $i++) {
    my $obj;
    my $true_col = $sa_true->column_from_residue_number($ref_name,$i);
    my $aln_col = $sa_aln->column_from_residue_number($ref_name,$i);

    $obj->{aln_position} = $aln_col;
    $obj->{true_dnds} = $true_omegas->{$true_col}->{'omega'};
    $obj->{aln_dnds} = $aln_omegas->{$aln_col}->{'omega'};
    if (!($obj->{aln_dnds} && $obj->{true_dnds})) {
      if ($cur_params->{'sitewise_action'} eq '') {
        # Do nothing.
        $obj->{aln_dnds} = 0;
        $obj->{true_dnds} = 0;
      } else {
        printf " =>Skipping! aln:%s  %s  true:%s  %s\n",$aln_col,$obj->{aln},$true_col,$obj->{true};
        next;
      }
    }
    $obj->{aln_type} = $aln_omegas->{$aln_col}->{'type'} || '';
    $obj->{true_type} = $true_omegas->{$true_col}->{'type'} || '';
    $obj->{aln_note} = $aln_omegas->{$aln_col}->{'note'} || '';
    $obj->{ncod} = $aln_omegas->{$aln_col}->{'ncod'} || 0;
    $obj->{true_e} = sprintf "%.3f", $true_entropies[$true_col];
    $obj->{aln_e} = sprintf "%.3f", $aln_entropies[$aln_col];
    $obj->{lrt} = $aln_omegas->{$aln_col}->{'lrt_stat'} || 99;

    # Store PAML LRTs if relevant.
    foreach my $tag (keys %$cur_params) {
      if ($tag =~ m/paml lrt/i) {
        $obj->{paml_lrt} = $cur_params->{$tag};
      }
    }

    # Store values in our output table.
    $obj = $self->replace_params($obj,$cur_params);
    my $table = $cur_params->{'collect_slrsim_stats_table'};
    $self->store_params_in_table($dba,$table,$obj);
  }
}



1;
