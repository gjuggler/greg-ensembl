#!/usr/bin/env perl

use warnings;
use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::Greg::EslrUtils;
use File::Path;
use File::Basename;

my ($url) = 'mysql://ensadmin:ensembl@ens-research/gj1_nodesets_57';
GetOptions('url=s' => \$url);

my $clean = 1;
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-url => $url);
my $dbc = $dba->dbc;
my $pta = $dba->get_ProteinTreeAdaptor;
my $mba = $dba->get_MemberAdaptor;

sub test_dist {
  my $stable_id = shift;

  my $member = $mba->fetch_by_source_stable_id(undef, $stable_id);
  my $tree = $pta->fetch_by_Member_root_id($member);
  my $dist = Bio::EnsEMBL::Compara::ComparaUtils->protein_tree_taxonomy_distance($tree);
}

#my $tree = Bio::EnsEMBL::Compara::TreeUtils->from_newick("(a:1,(b:1,c:1):1);");
#my $tree = Bio::EnsEMBL::Compara::NestedSet->new_from_TreeI($treeI);
#print $tree->newick_format."\n";
#foreach my $node ($tree->nodes) {
#  $node->re_root if (defined $node->parent && !$node->is_leaf);
#  last;
#}
#print $tree->newick_format."\n";
#$tree = $tree->unroot;
#print $tree->newick_format."\n";

test_dist("ENSG00000139618");
test_dist("ENSG00000007372");
