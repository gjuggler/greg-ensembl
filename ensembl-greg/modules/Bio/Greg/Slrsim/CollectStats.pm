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
  phylosim_insertrate  => 'float',
  phylosim_deleterate  => 'float',
  phylosim_insertmodel  => 'string',
  phylosim_deletemodel  => 'string',
  phylosim_domains     => 'string',

  parameter_set_name => 'string',

  tree_length        => 'float',
  tree_length_slr    => 'float',
  tree_max_branch    => 'float',
  tree_mean_branch   => 'float',
  tree_max_path      => 'float',
  tree_mean_path     => 'float',
  leaf_count         => 'int',

  true_dnds          => 'float',
  true_type          => 'string',
  true_entropy       => 'float',
  aln_dnds           => 'float',
  aln_type           => 'string',
  aln_entropy        => 'float',
  ungapped_branch_length => 'float',
  ncod               => 'int',
  lrt                => 'float',

  gene_lrt_paml      => 'float',
  gene_lnl_slr       => 'float',
  gene_omega_slr  => 'float',
  gene_kappa_slr     => 'float',

  sum_of_pairs_score => 'float',
  total_column_score => 'float',

  column_entropy_mean_true => 'float',
  column_entropy_mean_aln  => 'float',

  seq_length_mean          => 'float',
  gc_content_mean          => 'float',
  site_count               => 'float',
  unfiltered_site_count    => 'float',
  unfiltered_site_fraction => 'float',
  alignment_score_threshold => 'float',

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
  my $cdna_true;
  my $cdna_aln;
  my @true_entropies;
  my @aln_entropies;
  my $sum_of_pairs_score;
  my $total_column_score;
  my $tree;
  eval {
    my $true_aln_params = $self->replace_params($cur_params,{alignment_table => 'protein_tree_member', alignment_score_filtering => 0});
    ($tree,$sa_true,$cdna_true) = Bio::EnsEMBL::Compara::ComparaUtils->tree_aln_cdna($dba,$true_aln_params);
 
    $true_aln_params = $self->replace_params($cur_params,{alignment_score_filtering => 0});
    ($tree,$sa_aln,$cdna_aln) = Bio::EnsEMBL::Compara::ComparaUtils->tree_aln_cdna($dba,$cur_params);

    Bio::EnsEMBL::Compara::AlignUtils->pretty_print($sa_true,{length=>200});
    Bio::EnsEMBL::Compara::AlignUtils->pretty_print($sa_aln,{length=>200});
    
    @true_entropies = Bio::EnsEMBL::Compara::AlignUtils->column_entropies($cdna_true);
    @aln_entropies = Bio::EnsEMBL::Compara::AlignUtils->column_entropies($cdna_aln);
    $sum_of_pairs_score = Bio::EnsEMBL::Compara::AlignUtils->sum_of_pairs_score($sa_true,$sa_aln);
    $total_column_score = Bio::EnsEMBL::Compara::AlignUtils->total_column_score($sa_true,$sa_aln);
    print "sps: $sum_of_pairs_score  tcs: $total_column_score\n";
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
#  die ("Reference was defined in params but not found in aln!") if ($reference_id ne '' && !defined $ref_seq);
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
  $cur_params->{leaf_count} = scalar($tree->leaves);

  $cur_params->{sum_of_pairs_score} = $sum_of_pairs_score;
  $cur_params->{total_column_score} = $total_column_score;

  $cur_params->{column_entropy_mean_true} = Bio::EnsEMBL::Compara::AlignUtils->average_column_entropy($sa_true);
  $cur_params->{column_entropy_mean_aln} = Bio::EnsEMBL::Compara::AlignUtils->average_column_entropy($sa_aln);

  $cur_params->{seq_length_mean} = $utils->seq_length_mean($tree);
  $cur_params->{gc_content_mean} = $utils->gc_content_mean($cdna_true);

  $cur_params->{site_count} = $utils->site_count($sa_aln);
  $cur_params->{unfiltered_site_count} = $utils->unfiltered_site_count($sa_aln);
  $cur_params->{unfiltered_site_fraction} = $utils->unfiltered_site_count($sa_aln) / $utils->site_count($sa_aln);

  # Get all the site-wise data from the omega table.
  my $aln_table_name = $cur_params->{'omega_table'};
  my $sth1 = $pta->prepare("SELECT aln_position,omega,type,note,ncod,lrt_stat FROM sitewise_omega WHERE node_id=?;");
  my $sth2 = $pta->prepare("SELECT aln_position,omega,type,note,ncod,lrt_stat FROM $aln_table_name WHERE node_id=? AND parameter_set_id=?;");
  $sth1->execute($node_id);
  $sth2->execute($node_id,$parameter_set_id);
  my $true_omegas = $sth1->fetchall_hashref('aln_position');
  my $aln_omegas = $sth2->fetchall_hashref('aln_position');

  # Store PAML LRTs if relevant.
  foreach my $tag (keys %$cur_params) {
#    print "$tag\n";
    if ($tag =~ m/paml lrt/i) {
      $cur_params->{gene_lrt_paml} = $cur_params->{$tag};
    }
  }
  
  my $ps = $cur_params->{'parameter_set_id'};
  $cur_params->{'gene_omega_slr'} = $cur_params->{'slr_omega_'.$ps};
  $cur_params->{'gene_kappa_slr'} = $cur_params->{'slr_kappa_'.$ps};
  $cur_params->{'gene_lnl_slr'} = $cur_params->{'slr_lnL_'.$ps};

  # Get the SLR-inferred tree.
  my $newick = $cur_params->{'slr_tree_'.$ps};
  my $slr_tree = Bio::EnsEMBL::Compara::TreeUtils->from_newick($newick);
  $cur_params->{'tree_length_slr'} = $utils->tree_length($slr_tree);


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
        if (!$obj->{true_dnds}) {
          printf " =>Skipping! aln:%s  %s  true:%s  %s\n",$aln_col,$obj->{aln},$true_col,$obj->{true};
          next;
        } elsif (!$obj->{aln_dnds}) {
          # Continue; this should be counted as a false-negative in our results.
        }
      }
    }
    my @array = Bio::EnsEMBL::Compara::AlignUtils->get_column_array($sa_aln,$aln_col);
    #print join("",@array)." ".$obj->{aln_dnds}." ".$aln_omegas->{$aln_col}->{'note'}."\n";

    $obj->{aln_type} = $aln_omegas->{$aln_col}->{'type'} || '';
    $obj->{true_type} = $true_omegas->{$true_col}->{'type'} || '';
    $obj->{aln_note} = $aln_omegas->{$aln_col}->{'note'} || '';
    $obj->{ncod} = $aln_omegas->{$aln_col}->{'ncod'} || 0;
    $obj->{true_entropy} = sprintf("%.3f",$true_entropies[$true_col]||0);
    $obj->{aln_entropy} = sprintf("%.3f", $aln_entropies[$aln_col]||0);
    $obj->{lrt} = $aln_omegas->{$aln_col}->{'lrt_stat'} || 99;

    $obj->{ungapped_branch_length} = Bio::EnsEMBL::Compara::AlignUtils->get_ungapped_branchlength($sa_aln,$tree,$aln_col);

    # Store values in our output table.
    $obj = $self->replace_params($obj,$cur_params);
    my $table = $cur_params->{'collect_slrsim_stats_table'};
    $self->store_params_in_table($dba,$table,$obj);
  }
}



1;
