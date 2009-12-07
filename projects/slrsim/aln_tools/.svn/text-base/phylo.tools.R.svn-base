### phylo.tools.R -- Phylo tools for R   ###
###  written by Greg Jordan, fall 2008   ###
###  email: greg@ebi.ac.uk               ###
library(ape)

# Plots a phylogram within the specified rectangular area.
plot.phylo.within.rect = function(phylo,xleft,ybottom,xright,ytop,lwd=1,col='black',use.edge.length=TRUE) 
{  
  plot.phylo.greg(phylo,type="phylogram",
    show.tip.label=TRUE, use.edge.length=use.edge.length,
    no.margin=FALSE,
    plot=TRUE,  # Doesn't actually do any plotting.
    edge.color=col,
    edge.width=lwd
  )

  # Grab the output from the last plot command.
  p = .PlotPhyloEnv$last_plot.phylo
  xx = p$xx
  yy = p$yy
  edge = p$edge
  Nnode = p$Nnode
  Ntip = p$Ntip
  color = p$color
  width = p$width

  # Re-scale the x and y coordinates
  #xx = map(x=xx,xleft,xright)
  #yy = map(x=yy,ybottom,ytop)
  
  # Do the plotting using the real plot function.
  #phylogram.plot(edge,Ntip,Nnode,xx,yy,horizontal=TRUE,
  #  edge.color=color,edge.width=width)
}

# Get the vector of lengths-to-root for the tree.
lengths.to.root = function(x) {
  Ntip = length(x$tip.label)
  Nedge = dim(x$edge)[1]
  Nnode = x$Nnode
  dist_to_root <- .C("node_depth_edgelength", as.integer(Ntip),
                       as.integer(Nnode), as.integer(x$edge[, 1]),
                       as.integer(x$edge[, 2]), as.integer(Nedge),
                       as.double(x$edge.length), double(Ntip + Nnode),
                       DUP = FALSE, PACKAGE = "ape")[[7]]
  return(dist_to_root)
}


# Linearly maps a vector from one range to another.
map = function(x, newLo, newHi, origLo=NULL, origHi=NULL) {
  # If the user didn't specify the original range, then infer from the data.
  if (is.null(origLo)) origLo = min(x)
  if (is.null(origHi)) origHi = max(x)
  x = x - origLo
  x = x * (newHi-newLo)/(origHi-origLo)
  x = x + newLo
  return(x)
}

# Returns the maximum length from root to tip for a phylogeny.
# GJ 2009-01-28
max_length_to_root = function(phylo) {
  #max_dist = max(lengths.to.root(phylo))
  max_dist = max(dist.nodes(phylo)[root.node(phylo),])
  #if (is.rooted(phylo)) return(max_dist - phylo$root.edge)
  #else 
  return(max_dist)
}

# Alias of the max_length_to_root function.
# GJ 2009-01-28
tree_length = max_length_to_root;


# Returns a color vector where the indicated branches are given the indicated color.
# GJ 2009-01-28
color.edges = function(phylo,edges,colors=NULL,color="red") {
  if (is.null(colors)) {
    colors = rep("black",nrow(phylo$edge))
  }
  colors[edges] = color
  return(colors)
}


# Returns a list of the indices of edges within 
# the subtree defined by a given list of leaves.
# GJ 2009-01-28
extract_subtree_from_leaves = function(phylo,leaves) {
  # TODO!!
}


axis.phylo = function(tree,side=1) {
  max.len = tree_length(tree)
  if (side == 1 || side == 3) {
    axis(side=side,at=c(0,max.len),labels=c("0",max.len))
  } 
}

root.node = function(tree) {
  return(length(tree$tip.label)+1)
}