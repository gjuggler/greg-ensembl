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
GetOptions('url=s' => \$url,
           'id=s' => \$id
	   );

my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-url => $url);
my $dbc = $dba->dbc;
my $dbh = $dbc->db_handle;
my $pta = $dba->get_ProteinTreeAdaptor();
my $mba = $dba->get_MemberAdaptor();

die("No ID given!") unless ($id);

sub node_from_stable_id {
  my $stable_id = shift;

  my $member = Bio::Greg::EslrUtils->find_member_by_external_id($dba,$stable_id);
  if (defined $member) {
    $member = $member->get_canonical_peptide_Member;

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
    
    print "node: $node_id\n";
    return $node_id;
  }
  return undef;
}

my $node_id = node_from_stable_id($id);

die("No node ID found!") unless ($node_id);

my $plot_params = {
  dba => $dba,
  node_id => $node_id,
  parameter_set_id => 2,
  remove_blank_columns => 1
};
my $tree = Bio::EnsEMBL::Compara::ComparaUtils->get_tree_for_comparative_analysis($dba,$plot_params);
die("No tree found!") unless ($tree);
Bio::Greg::EslrPlots->plotTreeWithOmegas("~/public_html/test.pdf",$plot_params,$tree);
