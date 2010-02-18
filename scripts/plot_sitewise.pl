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
  my $node_id = '';
  if ($stable_id =~ /ENSP/i) {
    my $sth = $dba->dbc->prepare("select sitewise_node_for_peptide(\"${stable_id}\",1);");
    $sth->execute();
    $node_id = @{$sth->fetchrow_arrayref()}[0];
  } elsif ($stable_id =~ /ENSG/i) {
    my $sth = $dba->dbc->prepare("select sitewise_node_for_gene(\"${stable_id}\",1);");
    $sth->execute();
    $node_id = @{$sth->fetchrow_arrayref()}[0];
  }
  return $node_id;
}

my $node_id;
if ($id =~ m/EN/gi) {
  $node_id = node_from_stable_id($id);
} else {
  $node_id = $id;
}

die("No node ID found!") unless ($node_id);

my $plot_params = {
  dba => $dba,
  node_id => $node_id,
};
my $tree = $pta->fetch_node_by_node_id($node_id);
die("No tree found!") unless ($tree);
Bio::Greg::EslrPlots->plotTreeWithOmegas("~/public_html/test.pdf",$plot_params,$tree);

$plot_params = {
  dba => $dba,
  node_id => $node_id,
};
$tree = $pta->fetch_tree_at_node_id($node_id);
die("No tree found!") unless ($tree);
Bio::Greg::EslrPlots->plotTreeWithOmegas("~/public_html/test2.pdf",$plot_params,$tree);

