#!/usr/bin/env perl

use warnings;
use strict;

use DBI;

use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::TreeUtils;

use Bio::Greg::EslrUtils;

my $url = 'mysql://ensadmin:ensembl@ens-research/gj1_gor_58';
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new( -url => $url );

my @gene_names = qw(
PSKH1
SUPT16H
DNAH1
IL6R
CHD9
RC3H1
GABRR2
APC2
CHD1L
NEURL1B
DEPDC1B
RER1
ANTXR2
SLC6A6
ATAD2B
ZNF536
);

#@gene_names = @gene_names[1];

foreach my $gene_name (@gene_names) {
  plot_sub_tree_for_gene($gene_name);
}

sub plot_sub_tree_for_gene {
  my $gene_name = shift;
  print "  plotting $gene_name\n";

  my $gene_sth = $dba->dbc->prepare("SELECT * from lnl_m_Hsap where name='${gene_name}';");
  $gene_sth->execute;
  my $gene_row = $gene_sth->fetchrow_hashref;
  $gene_sth->finish;

  my $node_id = $gene_row->{node_id};

  my $sth = $dba->prepare("SELECT * FROM lnl_m_Hsap_subs where node_id=$node_id");
  
  my $pta = $dba->get_ProteinTreeAdaptor;
  my $tree = $pta->fetch_node_by_node_id($node_id);
  
  my @species_list = (9606, 
                      9593, 
#                      9598, 
                      9600, 
                      9544, 
                      9483
  );
  my @node_ids;
  foreach my $taxid (@species_list) {
    my ($leaf) = grep {$_->taxon_id == $taxid} $tree->leaves;
    die("No leaf!") unless ($leaf);
    push @node_ids, $leaf->node_id;
  }
  
  my $sub_tree = Bio::EnsEMBL::Compara::TreeUtils->extract_subtree_from_leaves($tree, \@node_ids);
  
  my $map;
  map {
    $map->{$_->name} = $_->taxon_id;
  } $sub_tree->leaves;
  $sub_tree = Bio::EnsEMBL::Compara::TreeUtils->translate_ids($sub_tree, $map);
  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($sub_tree);
  
  map {
    $_->branch_length(0.01);
  } $treeI->nodes;
  
  $sth->execute;
  my $muts;
  my @rows = @{ $sth->fetchall_arrayref({}) }; # Fetch all rows into an arrayref of hashrefs.
  $sth->finish;

  foreach my $row (@rows) {
    my $taxid = $row->{taxon_id};
    $muts->{$taxid} = [] if (!$muts->{$taxid});
    push @{$muts->{$taxid}}, $row;
  }

  foreach my $taxid (keys %$muts) {
    # Take all mutations for a given branch, split the syn ones in half and
    # 'engulf' the nsyn's w/ the syn's.
    my @cur_branch_muts = @{$muts->{$taxid}};
    my @syn_rows = grep {$_->{mut_nsyn} == 0} @cur_branch_muts;
    my @nsyn_rows = grep {$_->{mut_nsyn} == 1} @cur_branch_muts;

    foreach my $row (@syn_rows, @nsyn_rows) {
      # Using the 'leaves_beneath' string from the table,
      # parse out the taxon IDs of the leaves beneath,
      # and match up to the node from our TreeI object
      my $leaves = $row->{leaves_beneath};
      $leaves =~ s/ens_//gi;
      my @taxids = split(",", $leaves);
      my @taxid_nodes = map {$treeI->find($_)} @taxids;
      my $tree_node;
      if (scalar(@taxid_nodes) > 1) {
        $tree_node = $treeI->lca(@taxid_nodes);
      } else {
        $tree_node = $taxid_nodes[0];
      }
      
      # Create a new internal node representing this mutation event in the tree.
      my $new_node = new $tree_node;
      $new_node->id('m_'.$row->{aln_pos});
      $new_node->set_tag_value('nsyn', $row->{mut_nsyn});
      $new_node->set_tag_value('from', $row->{nuc_from});
      $new_node->set_tag_value('to', $row->{nuc_to});
      
      $tree_node->split_branch_with_node($new_node, 1);
      $new_node->branch_length(1);
      $tree_node->branch_length(0);
    }
  }

  my $taxid_to_species;
  map {
    $taxid_to_species->{$_->taxon_id} = $_->taxon->ensembl_alias_name;
  } $sub_tree->leaves;
  
  map {
    $_->id($taxid_to_species->{$_->id});
  } $treeI->leaves;
  
  my $str = $treeI->root->as_text('nhx');
  my $tree_subs_file = "${gene_name}_subs.nhx";
  my $tree_file = "${gene_name}.nh";
  my $pdf_file = "${gene_name}.pdf";

  open(OUT, ">$tree_subs_file");
  print OUT $str."\n";
  close(OUT);

  Bio::EnsEMBL::Compara::TreeUtils->to_file($sub_tree, $tree_file);  
  
  my $rcmd = qq^
library(phylosim)
source("~/src/greg-ensembl/projects/phylosim/Plots.R")

sim <- PhyloSim()
phylo <- read.nhx("$str")
sim\$.phylo <- phylo

pdf(file="${pdf_file}")
obj <- plotTree(sim, tree.do.plot=F)
p <- obj[['grob']]
p <- p + opts(title="${gene_name}")
p <- p + scale_x_continuous(name="Number of inferred substitutions since LCA")
print(p)
dev.off()
^;
  
  Bio::Greg::EslrUtils->run_r($rcmd);
}
