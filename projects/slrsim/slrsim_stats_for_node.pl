#!/usr/bin/env perl

use warnings;
use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::AlignUtils;
use Bio::Greg::ComparaLite::HiveUtils;
use Bio::Greg::EslrUtils;
use Bio::EnsEMBL::Hive::DBSQL::DBAdaptor;
use File::Path;
use File::Basename;
use Cwd;

my ($url,$clean) = undef;
GetOptions('url=s' => \$url,
	   'clean' => \$clean
	   );
$url = 'mysql://slrsim:slrsim@mysql-greg.ebi.ac.uk:4134/gj1_slrsim_1' if (!$url);
my $project_base = getcwd();

# Load the adaptors.
Bio::EnsEMBL::Registry->no_version_check(1);
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-url => $url);
my $dbc = $dba->dbc;
my $dbh = $dbc->db_handle;
my $pta = $dba->get_ProteinTreeAdaptor();

get_data_for_node($ARGV[0],$ARGV[1]);

sub get_data_for_node {
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

    Bio::EnsEMBL::Compara::AlignUtils->indelign($cdna_aln,$tree,{});

    exit(0);
  };
  return if (!$sa_true || !$sa_aln);

  

  my $sim_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_tag($tree,"sim_params");
  my $param_set_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_param_set($tree->adaptor,$parameter_set_id);

  # Get the sequence to act as a reference in site-wise value comparisons.
  my $reference_id = '';
  if (defined $sim_params->{'reference_id'}) {
    $reference_id = $sim_params->{'reference_id'};
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
  my $sth1 = $dbh->prepare("SELECT aln_position,omega FROM sitewise_aln WHERE node_id=?;");
  my $sth2 = $dbh->prepare("SELECT aln_position,omega,type,note,ncod FROM $aln_table_name WHERE node_id=? AND parameter_set_id=?;");
  $sth1->execute($node_id);
  $sth2->execute($node_id,$parameter_set_id);
  
  my $true_omegas = $sth1->fetchall_hashref('aln_position');
  my $aln_omegas = $sth2->fetchall_hashref('aln_position');

  print join("\t",qw(ref_seq true aln aln_type aln_note ncod true_e aln_e))."\n";
  for (my $i=1; $i <= length($nogaps); $i++) {
    my $true_col = $sa_true->column_from_residue_number($ref_name,$i);
    my $aln_col = $sa_aln->column_from_residue_number($ref_name,$i);

#    print "${true_col} ${aln_col}\n";

    my $true = $true_omegas->{$true_col}->{'omega'};
    my $aln = $aln_omegas->{$aln_col}->{'omega'};
    next unless ($aln && $true);
    my $aln_type = $aln_omegas->{$aln_col}->{'type'};
    my $aln_note = $aln_omegas->{$aln_col}->{'note'};
    my $ncod = $aln_omegas->{$aln_col}->{'ncod'};
    my $true_e = sprintf "%.3f", $true_entropies[$true];
    my $aln_e = sprintf "%.3f", $aln_entropies[$aln];

    print join("\t",($ref_name,$true,$aln,$aln_type,$aln_note,$ncod,$true_e,$aln_e))."\n";
  }
}
