#!/usr/bin/env perl

use warnings;
use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;
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
$url = 'mysql://greg:TMOqp3now@mysql-greg.ebi.ac.uk:4134/gj1_slrsim_1' if (!$url);

my ($mysql_args,$database,$project_base) = undef;
$mysql_args = Bio::Greg::ComparaLite::HiveUtils->hive_url_to_mysql_args($url);
$database = Bio::Greg::ComparaLite::HiveUtils->hive_url_to_hashref($url)->{'database'};
$project_base = getcwd();

# Load the adaptors.
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-url => $url);
my $dbc = $dba->dbc;
my $dbh = $dbc->db_handle;
my $nsa = $dba->get_NestedSetAdaptor();
my $mba = $dba->get_MemberAdaptor();
my $pta = $dba->get_ProteinTreeAdaptor();

my $GENE_HAS = 'gene_has';
my $GENE_INFER = 'gene_infer';
my $GENE_INFER_WEAK = 'gene_infer_weak';
my $TREE_LENGTH = 'tree_length';
my $G_PP = 'gene_++';
my $G_PN = 'gene_+-';
my $G_NP = 'gene_-+';
my $G_NN = 'gene_--';
my $PP = 'count_++';
my $PN = 'count_+-';
my $NP = 'count_-+';
my $NN = 'count_--';

# Get an array of the simulation sets from the protein_tree_tag table.
my $cmd = qq^SELECT distinct value FROM protein_tree_tag where tag="sim_set"^;
my @sim_sets = Bio::Greg::EslrUtils->mysql_array($dbc,$cmd);
@sim_sets = sort @sim_sets;

my @parameter_sets = (2,5,6);
#print "@sim_sets\n";
my $i=0;

my $sth1 = $dbh->prepare("SELECT omega FROM sitewise_aln WHERE node_id=? AND aln_position=?;");
my $sth2 = $dbh->prepare("SELECT omega,type,note FROM omega_mc WHERE node_id=? AND aln_position=? AND parameter_set_id=?;");

foreach my $sim_set (@sim_sets) {
  my @node_ids = Bio::Greg::EslrUtils->mysql_array($dbc,"SELECT distinct node_id FROM protein_tree_tag WHERE tag='sim_set' AND value='$sim_set';");
  $sim_set =~ s/sim_set_//;
  $sim_set =~ s/_[0-9]+//;
  
  foreach my $parameter_set (@parameter_sets) {
    $i++;
    my $cmd = qq^SELECT parameter_value FROM parameter_set where parameter_set_id=$parameter_set AND parameter_name='name'^;
    my @arr = Bio::Greg::EslrUtils->mysql_array($dbc,$cmd);
    my $name = $arr[0];

    my $s_pp = 0;
    my $s_pn = 0;
    my $s_np = 0;
    my $s_nn = 0;
    
    my $g_pp = 0;
    my $g_pn = 0;
    my $g_np = 0;
    my $g_nn = 0;
    my $n_nodes = scalar(@node_ids);
    my $tree_length;
    foreach my $node_id (@node_ids) {
      my $stats = get_stats_for_node($node_id,$parameter_set);
      $g_pp += $stats->{$G_PP};
      $g_pn += $stats->{$G_PN};
      $g_np += $stats->{$G_NP};
      $g_nn += $stats->{$G_NN};
      
      $s_pp += $stats->{$PP};
      $s_pn += $stats->{$PN};
      $s_np += $stats->{$NP};
      $s_nn += $stats->{$NN};
      
      $tree_length = $stats->{$TREE_LENGTH};
    }
    
    my $n = scalar(@node_ids);
    
#    printf "%s %s (%d reps)\n",$sim_set,$name,scalar(@node_ids);
    my $s_acc = 0;
    my $s_pow = 0;
    $s_acc = $s_pp/($s_pp+$s_np) if ($s_pp + $s_np > 0);
    $s_pow = $s_pp/($s_pp+$s_pn) if ($s_pp + $s_pn > 0);
#    printf "  SITES Acc:%.3f  Pow:%.3f \n", $s_acc, $s_pow;
    my $g_acc = 0;
    my $g_pow = 0;
    $g_acc = $g_pp/($g_pp+$g_np) if ($g_pp + $g_np > 0);
    $g_pow = $g_pp/($g_pp+$g_pn) if ($g_pp + $g_pn > 0);
#    printf "  GENES Acc:%.3f  Pow:%.3f \n", $g_acc, $g_pow;
    my @values = ($s_acc,$s_pow,$g_acc,$g_pow);
    my @value_s = map {sprintf("%.3f",$_)} @values;
    print join("\t",(
                 "sim_set",
                 "name",
                 "reps",
                 "tree_length",
                 "s_acc",
                 "s_pow",
                 "g_acc",
                 "g_pow"
               ))."\n" if ($i==1);
    print join("\t",(
                  $sim_set,
                 $name,
                 $n_nodes,
                 $tree_length,
                 @value_s
                ))."\n";
  }
}

$sth1->finish;
$sth2->finish;

sub get_stats_for_node {
  my $node_id = shift;
  my $parameter_set_id = shift;
  
  my $sa_true;
  my $sa_aln;
  my $tree;
  eval {
    $pta->protein_tree_member("protein_tree_member");
    $tree = $pta->fetch_node_by_node_id($node_id);
    $sa_true = $tree->get_SimpleAlign();
    #Bio::EnsEMBL::Compara::AlignUtils->pretty_print($sa_true,{length => 200});
    
    $pta->protein_tree_member("aln_mcoffee");
    $tree = $pta->fetch_node_by_node_id($node_id);
    $sa_aln = $tree->get_SimpleAlign();
    #Bio::EnsEMBL::Compara::AlignUtils->pretty_print($sa_aln,{length => 200});
  };
  my @seqs = $sa_true->each_seq;
#  my @seq_f = grep {$_->id eq 'human'} @seqs;
#  my $seq = $seq_f[0];
#  $seq = $seqs[0] if (!$seq);
  my $seq= $seqs[0];
  my $name = $seq->id;
  my $str = $seq->seq;
#  print $name."\n";
  my $nogaps = $str;
  $nogaps =~ s/-//g;

  my $stats;
  $stats->{$GENE_HAS} = 0;
  $stats->{$GENE_INFER} = 0;
  $stats->{$GENE_INFER_WEAK} = 0;
  $stats->{$PP} = 0;
  $stats->{$PN} = 0;
  $stats->{$NP} = 0;
  $stats->{$NN} = 0;
  $stats->{$G_PP} = 0;
  $stats->{$G_PN} = 0;
  $stats->{$G_NP} = 0;
  $stats->{$G_NN} = 0;

  return $stats if (!$sa_true || !$sa_aln);


  for (my $i=1; $i <= length($nogaps); $i++) {
    my $true_col = $sa_true->column_from_residue_number($name,$i);
    my $aln_col = $sa_aln->column_from_residue_number($name,$i);

    #my $cmd = qq^SELECT omega FROM sitewise_aln WHERE node_id=$node_id AND aln_position=$true_col;^;
    $sth1->execute($node_id,$true_col);
    my @arr = $sth1->fetchrow_array();
    my $omg = $arr[0];
    $stats->{$GENE_HAS} = 1 if ($omg > 1);

    # Get the inferred omega for the same site in the same protein.
    #@arr = $dbh->selectrow_array(qq^SELECT omega,type,note FROM omega_mc WHERE node_id=$node_id AND aln_position=$aln_col AND parameter_set_id=$parameter_set_id;^);
    $sth2->execute($node_id,$aln_col,$parameter_set_id);
    my @arr2 = $sth2->fetchrow_array();
    my ($omg2,$type2,$note2) = @arr2;
    if (scalar(@arr2) != 0) {
#      $stats->{$GENE_INFER_WEAK} = 1 if ($type2 =~ m/positive[1234]/);
      $stats->{$GENE_INFER} = 1 if ($type2 =~ m/positive[34]/);
      
      # Compare to this site's true omega.
      my $method = 'confident';
      if ($method eq 'estimate') {
        $stats->{$PP}++ if ($omg > 1 && $omg2 > 1);
        $stats->{$PN}++ if ($omg > 1 && $omg2 <= 1);
        $stats->{$NP}++ if ($omg <= 1 && $omg2 > 1);
        $stats->{$NN}++ if ($omg <= 1 && $omg2 <= 1);
      } elsif ($method eq 'confident') {
        $stats->{$PP}++ if ($omg > 1 && $type2 =~ m/positive[1234]/);
        $stats->{$PN}++ if ($omg > 1 && $type2 !~ m/positive[1234]/); 
        $stats->{$NP}++ if ($omg <= 1 && $type2 =~ m/positive[1234]/);
        $stats->{$NN}++ if ($omg <= 1 && $type2 !~ m/positive[1234]/);
      }
    }
  }

#  $sth1->finish;
#  $sth2->finish;

  $stats->{$G_PP} = 1 if ($stats->{$GENE_HAS} && $stats->{$GENE_INFER});
  $stats->{$G_PN} = 1 if ($stats->{$GENE_HAS} && !$stats->{$GENE_INFER});
  $stats->{$G_NP} = 1 if (!$stats->{$GENE_HAS} && $stats->{$GENE_INFER});
  $stats->{$G_NN} = 1 if (!$stats->{$GENE_HAS} && !$stats->{$GENE_INFER});

  $stats->{$TREE_LENGTH} = Bio::EnsEMBL::Compara::TreeUtils->total_distance($tree);

  return $stats;
}
