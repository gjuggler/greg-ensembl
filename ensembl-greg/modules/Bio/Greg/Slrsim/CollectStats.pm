package Bio::Greg::Slrsim::CollectStats;

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

my @sim_tags_to_collect = 
  qw(
     node_id

     slrsim_scheme_name
     alignment_name
     filtering_name
     species_name
     sitewise_name

     slrsim_rep
     slrsim_file
     slrsim_ref
     slrsim_tree_length
     phylosim_simulation_program
     phylosim_seq_length
     phylosim_omega_distribution
     phylosim_ins_rate
     phylosim_del_rate
     parameter_set_name
     );


sub fetch_input {
  my ($self) = @_;

  $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-DBCONN=>$self->db->dbc);
  $pta = $dba->get_ProteinTreeAdaptor;

  $params = {
    collect_stats_parameter_sets => 'all',
    sitewise_stats_table => 'sitewise_stats'
  };

  my $p_params = $self->get_params($self->parameters);
  my $i_params = $self->get_params($self->input_id);
  my $node_id = $i_params->{'protein_tree_id'} || $i_params->{'node_id'};
  my $t_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_tree_tags($dba,$node_id);

  $params = $self->replace_params($params,$p_params,$i_params,$t_params);
  Bio::EnsEMBL::Compara::ComparaUtils->hash_print($params);
}

sub run {
  my $self = shift;

  my $node_id = $params->{'node_id'};
  my $param_set_string = $self->get('parameter_sets');

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

sub get {
  my $self = shift;
  my $key = shift;
  return $params->{'collect_stats_'.$key};
}

sub get_data_for_node {
  my $self = shift;
  my $node_id = shift;
  my $parameter_set_id = shift;

  my $param_set_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_param_set($dba->dbc,$parameter_set_id);
  my $new_params = $self->replace_params($params,$param_set_params);
  Bio::EnsEMBL::Compara::ComparaUtils->hash_print($new_params);

  my $sa_true;
  my $sa_aln;
  my @true_entropies;
  my @aln_entropies;
  my $sum_of_pairs_score;
  my $total_column_score;
  my $tree;
  eval {
    $pta->protein_tree_member("protein_tree_member");
    $tree = $pta->fetch_node_by_node_id($node_id);
    $sa_true = $tree->get_SimpleAlign();
    my $cdna_aln = $tree->get_SimpleAlign(-cdna => 1);
    @true_entropies = Bio::EnsEMBL::Compara::AlignUtils->column_entropies($cdna_aln);

    $pta->protein_tree_member($new_params->{'alignment_table'});
    $tree = $pta->fetch_node_by_node_id($node_id);
    $sa_aln = $tree->get_SimpleAlign();
    $cdna_aln = $tree->get_SimpleAlign(-cdna => 1);
    @aln_entropies = Bio::EnsEMBL::Compara::AlignUtils->column_entropies($cdna_aln);

    $sum_of_pairs_score = Bio::EnsEMBL::Compara::AlignUtils->sum_of_pairs_score($sa_true,$sa_aln);
    $total_column_score = Bio::EnsEMBL::Compara::AlignUtils->total_column_score($sa_true,$sa_aln);
  };  
  die("Hold up: ".$@) if ($@);
  return if (!$sa_true || !$sa_aln);
  
  # Load desired parameters into our output array.
  my @sim_tags = @sim_tags_to_collect;
  my @sim_vals = map {
    if (defined $new_params->{$_}) {
      $new_params->{$_};
    } else {
      'NA';
    }
  } @sim_tags;
  my @sitewise_tags = qw(true aln true_type aln_type aln_note ncod true_e aln_e lrt);
  
  # Get the sequence to act as a reference in site-wise value comparisons.
  my $reference_id = '';
  if (defined $new_params->{'slrsim_ref'}) {
    $reference_id = $new_params->{'slrsim_ref'};
  }

  my @seqs = $sa_true->each_seq;
  my ($ref_seq) = grep {$_->id eq $reference_id} @seqs;
#  die ("Reference was defined in params but not found in aln!") if ($reference_id ne '' && !defined $ref_seq);
  $ref_seq = $seqs[0] if (!defined $ref_seq);
  my $ref_name = $ref_seq->id;

  my $str = $ref_seq->seq;
  my $nogaps = $str;
  $nogaps =~ s/-//g;

  my $aln_table_name = $new_params->{'omega_table'};
  my $sth1 = $pta->prepare("SELECT aln_position,omega,type FROM sitewise_omega WHERE node_id=?;");
  my $sth2 = $pta->prepare("SELECT aln_position,omega,type,note,ncod,lrt_stat FROM $aln_table_name WHERE node_id=? AND parameter_set_id=?;");
  $sth1->execute($node_id);
  $sth2->execute($node_id,$parameter_set_id);
  
  my $true_omegas = $sth1->fetchall_hashref('aln_position');
  my $aln_omegas = $sth2->fetchall_hashref('aln_position');

  my $table = 'sitewise_stats';
  my $sth = $pta->prepare("REPLACE INTO $table (node_id, parameter_set_id, aln_position, stats) VALuES (?,?,?,?)");

  # Put the tab-delimited header string into the META table.
  my $header = join("\t",@sim_tags,
                    @sitewise_tags)."\n";
  my $header_sth = $pta->prepare("REPLACE INTO meta VALUES (123,1,?,?)");
  $header_sth->execute('slrsim_stats_header',$header);
  $header_sth->finish;

  for (my $i=1; $i <= length($nogaps); $i++) {
    my $obj;
    my $true_col = $sa_true->column_from_residue_number($ref_name,$i);
    my $aln_col = $sa_aln->column_from_residue_number($ref_name,$i);

    $obj->{true} = $true_omegas->{$true_col}->{'omega'};
    $obj->{aln} = $aln_omegas->{$aln_col}->{'omega'};
    if (!($obj->{aln} && $obj->{true})) {
      printf " =>Skipping! aln:%s  %s  true:%s  %s\n",$aln_col,$obj->{aln},$true_col,$obj->{true};
      next;
    }
    $obj->{aln_type} = $aln_omegas->{$aln_col}->{'type'} || '';
    $obj->{true_type} = $true_omegas->{$true_col}->{'type'} || '';
    $obj->{aln_note} = $aln_omegas->{$aln_col}->{'note'} || '';
    $obj->{ncod} = $aln_omegas->{$aln_col}->{'ncod'} || 0;
    $obj->{true_e} = sprintf "%.3f", $true_entropies[$true_col];
    $obj->{aln_e} = sprintf "%.3f", $aln_entropies[$aln_col];
    $obj->{lrt} = $aln_omegas->{$aln_col}->{'lrt_stat'} || 99;
    $obj->{sum_of_pairs_score} = $sum_of_pairs_score;
    $obj->{total_column_score} = $total_column_score;
    
    my @sitewise_vals = map {
      if (defined $obj->{$_}) {
        $obj->{$_};
      } else {
        'NA';
      }
    } @sitewise_tags;
    my $str = join("\t",@sim_vals,@sitewise_vals);
    print $str."\n";
    $sth->execute($node_id,$parameter_set_id,$aln_col,$str);
  }

  $sth->finish;
}

1;
