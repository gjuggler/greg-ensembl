#!/usr/bin/env perl

use warnings;
use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
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
my $pta = $dba->get_ProteinTreeAdaptor();

get_data_for_node($ARGV[0],$ARGV[1]);

sub get_data_for_node {
  my $node_id = shift;
  my $parameter_set_id = shift;
  
  my $sa_true;
  my $sa_aln;
  my $tree;
  eval {
    $pta->protein_tree_member("protein_tree_member");
    $tree = $pta->fetch_node_by_node_id($node_id);
    $sa_true = $tree->get_SimpleAlign();
    
    $pta->protein_tree_member("aln_mcoffee");
    $tree = $pta->fetch_node_by_node_id($node_id);
    $sa_aln = $tree->get_SimpleAlign();
  };
  return if (!$sa_true || !$sa_aln);

  my @seqs = $sa_true->each_seq;
  my $seq= $seqs[0];
  my $name = $seq->id;
  my $str = $seq->seq;
  my $nogaps = $str;
  $nogaps =~ s/-//g;

  my $sth1 = $dbh->prepare("SELECT aln_position,omega FROM sitewise_aln WHERE node_id=?;");
  my $sth2 = $dbh->prepare("SELECT aln_position,omega,type,note FROM omega_mc WHERE node_id=? AND parameter_set_id=?;");
  $sth1->execute($node_id);
  $sth2->execute($node_id,$parameter_set_id);
  
  my $true_omegas = $sth1->fetchall_hashref('aln_position');
  my $aln_omegas = $sth2->fetchall_hashref('aln_position');

  print "true\taln\n";
  for (my $i=1; $i <= length($nogaps); $i++) {
    my $true_col = $sa_true->column_from_residue_number($name,$i);
    my $aln_col = $sa_aln->column_from_residue_number($name,$i);

    my $true = $true_omegas->{$true_col}->{'omega'};
    my $aln = $aln_omegas->{$aln_col}->{'omega'};

    print "$true\t$aln\n";
  }
}
