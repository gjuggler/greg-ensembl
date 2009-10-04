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
$url = 'mysql://greg:TMOqp3now@mysql-greg.ebi.ac.uk:4134/gj1_slrsim_1';

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


# Get an array of the simulation sets from the protein_tree_tag table.
my $cmd = qq^SELECT distinct value FROM protein_tree_tag where tag="sim_set"^;
my @sim_sets = Bio::Greg::EslrUtils->mysql_array($dbc,$cmd);

foreach my $sim_set (@sim_sets) {
  my @node_ids = Bio::Greg::EslrUtils->mysql_array($dbc,"SELECT distinct node_id FROM protein_tree_tag WHERE tag='sim_set' AND value='$sim_set';");
  
  my $pos_sel = 0;
  foreach my $node_id (@node_ids) {
    my $stats = genewise_stats_for_node($node_id);
    $pos_sel++ if ($stats->{'has_pos_sel'});
  }

  my $n = scalar(@node_ids);
  print "$sim_set  pos-sel:$pos_sel total:$n\n";
  print $pos_sel/$n."\n";
}


sub genewise_stats_for_node {
  my $node_id = shift;
  
  $pta->protein_tree_member("protein_tree_member");
  my $tree = $pta->fetch_node_by_node_id($node_id);
  my $sa_true = $tree->get_SimpleAlign();
  #Bio::EnsEMBL::Compara::AlignUtils->pretty_print($sa_true,{length => 200});

  #$pta->protein_tree_member("aln_mcoffee");
  #$tree = $pta->fetch_node_by_node_id($node_id);
  #my $sa_aln = $tree->get_SimpleAlign();
  #Bio::EnsEMBL::Compara::AlignUtils->pretty_print($sa_aln,{length => 200});

  my @seqs = $sa_true->each_seq;
  my $seq = $seqs[0];
  my $name = $seq->id;
  my $str = $seq->seq;
  my $nogaps = $str;
  $nogaps =~ s/-//g;

  my $stats;
  $stats->{'has_pos_sel'} = 0;
  for (my $i=1; $i <= length($nogaps); $i++) {
    my $true_col = $sa_true->column_from_residue_number($name,$i);
    #my $aln_col = $sa_aln->column_from_residue_number($name,$i);

    my $cmd = qq^SELECT omega FROM sitewise_aln WHERE node_id=$node_id AND aln_position=$true_col;^;
    my @arr = $dbh->selectrow_array($cmd);
    my $omg = $arr[0];
    $stats->{'has_pos_sel'} = 1 if ($omg > 1);

    #@arr = $dbh->selectrow_array(qq^SELECT omega,type,note FROM omega_tr WHERE node_id=$node_id AND aln_position=$true_col;^);
    #my ($omg2,$type2,$note2) = @arr;
    #if (scalar(@arr) != 0) {
#      $stats->{'true_'.$type2}++;
    #}

    #@arr = $dbh->selectrow_array(qq^SELECT omega,type,note FROM omega_mc WHERE node_id=$node_id AND aln_position=$aln_col;^);
    #my ($omg3,$type3,$note3) = @arr;
    #if (scalar(@arr) != 0) {
#      $stats->{'aln_'.$type3}++;
    #  next;
    #}
  }

  #print "$node_id\n";
  #Bio::EnsEMBL::Compara::ComparaUtils->hash_print($stats);
  return $stats;
}
