package Bio::Greg::Slrsim::CollectStats;

use strict;
use Time::HiRes qw(sleep);

use Cwd;
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::NestedSet;
use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Hive::Process;

our @ISA = qw(Bio::EnsEMBL::Hive::Process);

my $dba;
my $pta;

my $tree;
my $params;

my $results;

sub fetch_input {
  my ($self) = @_;

  $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-DBCONN=>$self->db->dbc);
  $pta = $dba->get_ProteinTreeAdaptor;

  $params = {
    alignment_table => 'aln_mcoffee',
    parameter_sets => '2,3',
  };

  $params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_string($params,$self->parameters);
  $params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_string($params,$self->input_id);
}

sub run {
  my $self = shift;

  my $node_id = $params->{'node_id'};

  my @parameter_sets = split(',',$params->{'parameter_sets'});
  foreach my $ps_id (@parameter_sets) {
    $self->get_data_for_node($node_id,$ps_id);
  }

}

sub get_data_for_node {
  my $self = shift;
  my $node_id = shift;
  my $parameter_set_id = shift;

  my $sa_true;
  my $sa_aln;
  my @true_entropies;
  my @aln_entropies;
  my $tree;
  eval {
    $pta->protein_tree_member("protein_tree_member");
    $tree = $pta->fetch_node_by_node_id($node_id);
    $sa_true = $tree->get_SimpleAlign();
    my $cdna_aln = $tree->get_SimpleAlign(-cdna => 1);
    @true_entropies = Bio::EnsEMBL::Compara::AlignUtils->column_entropies($cdna_aln);

    $pta->protein_tree_member("aln_mcoffee");
    $tree = $pta->fetch_node_by_node_id($node_id);
    $sa_aln = $tree->get_SimpleAlign();
    $cdna_aln = $tree->get_SimpleAlign(-cdna => 1);
    @aln_entropies = Bio::EnsEMBL::Compara::AlignUtils->column_entropies($cdna_aln);
  };
  die if ($!);
  return if (!$sa_true || !$sa_aln);
  
  my $sim_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_tag($tree,"params_slrsim");
  Bio::EnsEMBL::Compara::ComparaUtils->hash_print($sim_params);
  my $param_set_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_param_set($tree->adaptor,$parameter_set_id);

  # These two come from alternative sources.
  $sim_params->{'sim_rep'} = $tree->get_tagvalue('sim_rep');
  $sim_params->{'parameter_set_name'} = $param_set_params->{'parameter_set_name'};

  # Put sim_length, tree_mult into tree_length
  $sim_params->{'tree_length'} = $sim_params->{'sim_length'} || $sim_params->{'tree_mult'} || '';

  # Load simulation parameters into our output array.
  my @sim_tags = qw(sim_file simulation_program ins_rate del_rate omega_distribution
     seq_length sim_name sim_rep parameter_set_name tree_length);
  my @sim_vals = map {
    if (defined $sim_params->{$_}) {
      $sim_params->{$_};
    } else {
      'NA';
    }
  } @sim_tags;
  
  # Get the sequence to act as a reference in site-wise value comparisons.
  my $reference_id = '';
  if (defined $sim_params->{'sim_ref'}) {
    $reference_id = $sim_params->{'sim_ref'};
  }

  my @seqs = $sa_true->each_seq;
  my ($ref_seq) = grep {$_->id eq $reference_id} @seqs;
  die ("Reference was defined in params but not found in aln!") if ($reference_id ne '' && !defined $ref_seq);
  $ref_seq = $seqs[0] if (!defined $ref_seq);
  my $ref_name = $ref_seq->id;

  my $str = $ref_seq->seq;
  my $nogaps = $str;
  $nogaps =~ s/-//g;

  my $aln_table_name = $param_set_params->{'output_table'};
  my $sth1 = $pta->prepare("SELECT aln_position,omega,type FROM sitewise_omega WHERE node_id=?;");
  my $sth2 = $pta->prepare("SELECT aln_position,omega,type,note,ncod,lrt_stat FROM $aln_table_name WHERE node_id=? AND parameter_set_id=?;");
  $sth1->execute($node_id);
  $sth2->execute($node_id,$parameter_set_id);
  
  my $true_omegas = $sth1->fetchall_hashref('aln_position');
  my $aln_omegas = $sth2->fetchall_hashref('aln_position');

  my $table = 'sitewise_stats';
  my $sth = $pta->prepare("REPLACE INTO $table (node_id, parameter_set_id, aln_position, stats) VALuES (?,?,?,?)");

  # Put the tab-delimited header string into the meta table.
  my $header = join("\t",@sim_tags,
             qw(node_id parameter_set true aln true_type aln_type aln_note ncod true_e aln_e lrt))."\n";
  my $header_sth = $pta->prepare("REPLACE INTO meta VALUES (123,1,?,?)");
  $header_sth->execute('slrsim_stats_header',$header);
  $header_sth->finish;

  for (my $i=1; $i <= length($nogaps); $i++) {
    my $true_col = $sa_true->column_from_residue_number($ref_name,$i);
    my $aln_col = $sa_aln->column_from_residue_number($ref_name,$i);

    my $true = $true_omegas->{$true_col}->{'omega'};
    my $aln = $aln_omegas->{$aln_col}->{'omega'};
    if (!($aln && $true)) {
      print " =>Skipping! aln:$aln_col $aln  true:$true_col $true\n";
      next;
    }
    my $aln_type = $aln_omegas->{$aln_col}->{'type'} || '';
    my $true_type = $true_omegas->{$aln_col}->{'type'} || '';
    my $aln_note = $aln_omegas->{$aln_col}->{'note'} || '';
    my $ncod = $aln_omegas->{$aln_col}->{'ncod'} || 0;
    my $true_e = sprintf "%.3f", $true_entropies[$true];
    my $aln_e = sprintf "%.3f", $aln_entropies[$aln];
    my $lrt = $aln_omegas->{$aln_col}->{'lrt_stat'} || 99;
    
    my @vals = (@sim_vals,$node_id,$parameter_set_id,$true,$aln,$true_type,$aln_type,$aln_note,$ncod,$true_e,$aln_e,$lrt);
    my $str = join("\t",@vals);
    print $str."\n";
    $sth->execute($node_id,$parameter_set_id,$aln_col,$str);
  }

  $sth->finish;
}

1;
