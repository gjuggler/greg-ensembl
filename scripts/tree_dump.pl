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
my $tree_output;
my $aln_output;
GetOptions('url=s' => \$url,
           'id=s' => \$id,
           'tree=s' => \$tree_output,
           'aln=s' => \$aln_output
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

my $params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_tree_tags($dba,$node_id);
$params->{node_id} = $node_id;
Bio::EnsEMBL::Compara::ComparaUtils->hash_print($params);
my ($tree,$sa) = Bio::EnsEMBL::Compara::ComparaUtils->get_tree_and_alignment($dba,$params);
    
# Output the tree.
open(OUT,">$tree_output");
my $tree2 = $tree->copy;
$tree2 = $tree2->minimize_tree;
print OUT $tree2->newick_format()."\n";
close(OUT);

# Output the alignment.
open(OUT,">$aln_output");
foreach my $seq ($sa->each_seq) {
  print OUT ">".$seq->id."\n";
  print OUT $seq->seq."\n";
}
close(OUT);

exit(0);
