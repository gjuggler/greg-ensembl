setwd("C:/Users/Greg/Documents/Work/EBI/SVN/aln_tools/trunk/")
source("aln.tools.R")

### GJ 2009-01-23 : SIMPLE EXAMPLE
# Read in a Fasta-formatted alignment.
aln = read.aln("http://www.ebi.ac.uk/~greg/eslr/root_nodes/31/origAln.mfa")
# Plot the alignment, filling the screen entirely.
plot.aln(aln,square=FALSE)
# Plot the alignment, keeping residues square.
plot.aln(aln,square=TRUE)
# Plot the alignment all as one color.
plot.aln(aln,square=TRUE,colors="black")
### GJ 2009-01-23


### GJ 2009-01-28 : Trying to get multiple rows of alignments showing...

source("aln.tools.R")
source("plot.phylo.greg.R")
source("phylo.tools.R")
library(ape)
par(new=FALSE)
plot.new()
par(fig=c(0,1,0.666,1),new=FALSE)
aln1 = read.aln("http://www.ebi.ac.uk/~greg/eslr/dup_nodes/node_id/207514/tree_a.fasta")
aln1$tree = read.tree("http://www.ebi.ac.uk/~greg/eslr/dup_nodes/node_id/207514/tree_a.nh")
tree_len = tree_length(aln1$tree)
tree_width = 50 * tree_len;
plot.aln(aln1,square=FALSE,overlay=FALSE,tree.labels=FALSE,tree.width=tree_width,tree.space=50)

par(fig=c(0,1,0,0.333),new=TRUE)
aln2 = read.aln("http://www.ebi.ac.uk/~greg/eslr/dup_nodes/node_id/207514/tree_b.fasta")
aln2$tree = "http://www.ebi.ac.uk/~greg/eslr/dup_nodes/node_id/207514/tree_b.nh"
tree_len = tree_length(aln1$tree)
tree_width = 50 * tree_len;
plot.aln(aln2,square=FALSE,overlay=FALSE,tree.labels=FALSE,tree.width=tree_width,tree.space=50)

par(fig=c(0,1,0.333,0.666),new=TRUE)

source("slr.tools.R")
source("plot.phylo.greg.R")
source("phylo.tools.R")

# Hello!

plot.node = function(node_id,index) {
# node_id = 34
base = sprintf("http://www.ebi.ac.uk/~greg/eslr/root_nodes/%s",node_id)
f1 = sprintf("%s/%s",base,"origAln.mfa")
f2 = sprintf("%s/%s",base,"prankAln.mfa")
f3 = sprintf("%s/%s",base,"probconsAln.mfa")
treeF = sprintf("%s/%s",base,"origTree.nh")
slrF = sprintf("%s/%s",base,"origSlr.txt")
slr = read.table(slrF,header=TRUE)
aln1 = read.aln(f1)
aln2 = read.aln(f2)
aln3 = read.aln(f3)
aln1$tree = aln2$tree = aln3$tree = treeF

height = max(2.5,aln1$num_seqs / 50)
width = min(16,aln1$length/10)
#pdf(file=sprintf("%s.%s.test.pdf",format(Sys.time(),"%H.%M %S"),index),
#  width=width*5,height=height*5, family="Courier")

rowHeight = height/aln1$num_seqs
cex = rowHeight*4
aln1$colors = NULL
#plot.aln(aln1,cex=0.1,square=FALSE)

# With SLR colors.
#slr.colors = color.slr(slr)
#aln1$colors = slr.colors
aln1$tree = treeF
plot.aln(aln1,square=FALSE,tree.labels=TRUE,max.pointsize=16,tree.width=aln1$length/5)

# Alignment comparison.
alns = list(aln1,aln2,aln3)
source("aln.tools.R")
source("plot.phylo.greg.R")
plot.aln.comparison(alns,tree.labels=TRUE,draw.chars=TRUE,bar.height=.3)
dev.off()
}

for (node in c(21)) {
  plot.node(node,node)
}

# two nodes: 23, 3
# three nodes: 35
# medium: 24
# big tree: 34

x = read.tree(treeF)
max(lengths.to.root(x))

colorRow = 