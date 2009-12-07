setwd("C:/Users/Greg/Documents/Work/EBI/SVN/aln_tools/trunk/tests/slr/2")
source("../../../aln.tools.R")
source("../../../phylo.tools.R")
source("../../../slr.tools.R")
source("../../../plot.phylo.greg.R")


dev.off()
dpar = par(no.readonly=T)

# Use layout to arrange the regions for SLR values, tracks, tree, and alignment.
num_tracks = 3
track_heights = c(.5,.1,.1)
aln_height = 1
heights=c(track_heights,aln_height)
tree_width=0.5
aln_width=1
pad_width = aln_width/10
widths=c(tree_width,aln_width,pad_width)

# Plotting order: aln, tree, track 1,2,3...n
num_rows = num_tracks+1
num_cols = 3
numbers = c()
for (i in 1:num_tracks) {
  numbers = c(numbers,3,i+3,0)
}
numbers = c(numbers,2,1,0)

layout(matrix(numbers,nrow=num_rows,ncol=num_cols,byrow=T),heights=heights,widths=widths)
layout.show(num_tracks+3)
par(xaxs='i',mai=rep(0,4))
setwd("C:/Users/Greg/Documents/Work/EBI/SVN/aln_tools/trunk/tests/slr/2")
source("../../../aln.tools.R")
source("../../../slr.tools.R")
source("../../../plot.phylo.greg.R")
source("../../../phylo.tools.R")
tree = read.tree("tree.nh")
aln = read.aln("aln.fasta")
slr = read.slr("sitewise_aln.1.txt")
a = plot.aln(aln,overlay=F,axes=F,draw.chars=F)
xl = tree_length(tree)
plot.phylo.greg2(tree,y.lim=a$ylim,x.lim=c(-1,xl*1.1),show.tip.label=T)
axis.phylo(tree,side=3)
#hist(slr$omega,xlim=c(0,2),n=nrow(slr)/5)
plot(c(1),type='n',axes=F)
# Now plot the tracks.
plot.slr.blocks(slr,xlim=a$xlim,log=T)
axis(side=2)
abline(h=mean(slr$omega))
plot.slr.types(slr,xlim=a$xlim)
plot.aln.bars(aln)



par(dpar)