#!/usr/bin/env perl

use warnings;
use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;
use File::Path;
use Bio::Greg::EslrPlots;
use Bio::Greg::EslrUtils;

Bio::EnsEMBL::Registry->no_version_check(1);

my $url;
my $id;
my $file = "/homes/greg/public_html/plot.pdf";
GetOptions('url=s' => \$url,
           'id=s' => \$id,
           'file=s' => \$file,
	   );

my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-url => $url);
my $dbc = $dba->dbc;
my $dbh = $dbc->db_handle;
my $pta = $dba->get_ProteinTreeAdaptor();
my $mba = $dba->get_MemberAdaptor();

die("No ID given!") unless ($id);

sub node_from_stable_id {
  my $stable_id = shift;

  return $stable_id if ($stable_id =~ m/$[0-9]+^/);

  print "STABLE ID: $stable_id\n";
  
  my $ext_member = Bio::Greg::EslrUtils->find_member_by_external_id($dba,$stable_id);
  $stable_id = $ext_member->stable_id if (defined $ext_member);

  my $member = $dba->get_MemberAdaptor->fetch_by_source_stable_id(undef,$stable_id);
  $member = $member->get_canonical_peptide_Member;
  print $member->stable_id."\n";
  if (defined $member) {
    my $member_id = $member->dbID;

    my $cmd = qq^
SELECT ptt.node_id from protein_tree_tag ptt, protein_tree_node ptn1, protein_tree_node ptn2, protein_tree_member ptm
WHERE ptm.member_id=$member_id AND ptn2.node_id=ptm.node_id
AND ptt.node_id=ptn1.node_id AND ptt.tag LIKE "%slr%"
AND ptn2.left_index BETWEEN ptn1.left_index AND ptn1.right_index limit 1;
^;
    print $cmd."\n";
    my $sth = $dba->dbc->prepare($cmd);
    $sth->execute;
    my $node_id = $sth->fetchrow_array();
    $sth->finish;

#    my $tree = $pta->fetch_by_Member_root_id($member);
#    if (defined $tree) {
#      my $node_id = $tree->node_id;
    if (defined $node_id) {
      print "node: $node_id\n";

      my $sql = qq^SELECT data_id from stats_genes where node_id=$node_id;^;
      print $sql."\n";
      my $sth = $dba->dbc->prepare($sql);
      $sth->execute;
      if (my @row = $sth->fetchrow_array) {
        my $data_id = $row[0];
        print "DATA ID: $data_id\n";
        return {data_id => $data_id,
                node_id => $node_id};
      }
    }
  }
  return undef;
}

my $params = node_from_stable_id($id);

die("No node ID found!") unless ($params);

my $plot_params = {
  dba => $dba,
  parameter_set_id => 1,
  remove_blank_columns => 1
};
my $final_params = Bio::EnsEMBL::Compara::ComparaUtils->replace_params($plot_params,$params);
Bio::Greg::EslrPlots->plotTreeWithOmegas($file,$final_params);
