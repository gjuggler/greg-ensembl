source("aln.tools.R")

#base = "http://www.ebi.ac.uk/~greg/eslr/root_nodes"
base = "http://www.ebi.ac.uk/~greg/eslr/eslr_pages/html/node_id"
node_id = 685811
alnF = sprintf("%s/%s/ptm_cmcoffee.fasta",base,node_id)
aln2F = sprintf("%s/%s/protein_tree_member.fasta",base,node_id)
treeF = sprintf("%s/%s/%s.nh",base,node_id,node_id)

library(ape)
aln = read.aln(alnF)
aln2 = read.aln(aln2F)
tree = read.tree(treeF)

aln$tree = tree
aln2$tree = tree

plot.aln(aln,square=FALSE,tree.width=10)

plot.aln.comparison(alns=list(aln,aln2),grid_lines=FALSE)