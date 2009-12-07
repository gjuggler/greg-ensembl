source("plot.phylo.greg.R")

tree = "http://www.ebi.ac.uk/~greg/eslr/root_nodes/49/origTree.nh"
phylo = read.tree(file=tree)
plot(NULL,xlim=c(0,100),ylim=c(0,100))

plot.phylo.greg(phylo,
  draw.within.rect=c(0,50,0,100),
  cex=0.7)
