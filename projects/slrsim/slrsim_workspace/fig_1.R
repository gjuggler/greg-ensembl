library(ggplot2)
library(fields)
library(RColorBrewer)

plot.f = function(data,palette="Spectral",field='fdr',limits=c(0,1),x.lim=c(0,4),y.lim=c(0,0.1),label.species=T) {
  data$z_vals = data[,field]
  p <- ggplot(data,aes(x=tree_mean_path,y=phylosim_insertrate,z=z_vals))

  if (!is.null(x.lim)) {
    p <- p + coord_cartesian(xlim = x.lim,ylim=y.lim)
  }

  p <- p + geom_tile(aes_string(fill=field))

  if (is.null(limits)) {
    limits <- c(min(data$z_vals),max(data$z_vals))
  }
  n <- 10
  breaks <- seq(from=limits[1],to=limits[2],length.out=n+1)
  p <- p + scale_fill_gradientn(colours=brewer.pal(n=n,palette),limits=limits)
  p <- p + coord_equal(ratio=20)

  if (label.species) {
    lbls <- data.frame(x=NULL,y=NULL,str=NULL)
    lbls <- rbind(lbls,data.frame(x=0.13,y=0.04,str="Primates"))
    lbls <- rbind(lbls,data.frame(x=0.25,y=0.03,str="Mammals"))
    lbls <- rbind(lbls,data.frame(x=0.75,y=0.023,str="Vertebrates"))
    lbls <- rbind(lbls,data.frame(x=1.0,y=0.015,str="Drosophilids"))
    lbls <- rbind(lbls,data.frame(x=1.25,y=0.008,str="Yeasts"))
    lbls$z_vals <- 0
    lbls$slope <- 8
    arw <- arrow(length=unit(0.05,"inches"))
    p <- p + geom_segment(data=lbls,mapping=aes(x=x+slope*y,xend=x,y=y,yend=0),colour='white',size=0.2,arrow=arw)
    p <- p + geom_text(data=lbls,mapping=aes(x=x+slope*y,y=y,label=str),hjust=0,vjust=0,size=3,colour='white')

    orig.lengths = c(0.5)
    tree.file <- data[1,]$slrsim_tree_file
    tree.file = 'anisimova_01_bglobin.nh' 
    if (tree.file == 'anisimova_01_bglobin.nh') {
      orig.lengths <- c(0.11, 0.64, 2.54)
    } else if (tree.file == 'anisimova_01_artificial.nh') {
      orig.lengths <- c(0.025, 0.25, 2.5)
    } else if (tree.file == '44mammals.nh') {
      orig.lengths <- c(1.569)
    }
    pts <- data.frame(x=orig.lengths,y=0,z_vals=0,shape=1)
    p <- p + geom_point(data=pts,aes(x=x,y=y),shape=4,size=1,colour='white')
 
  }

  print(p)
}

series = 'NO.backup/2010-10-04/1_44mammals_clustalw'
source("slrsim.functions.R")
source("slrsim.plots.R")

load(paste(series,".Rdata",sep=""))
sub.data = data
slr.threshold = 3.84
a = summarize.by.labels(sub.data,fn=fig.1.summary,thresh=slr.threshold)

pdf(file=paste(series,"_a.pdf",sep=""),width=6,height=6)
plot.f(a,field='fdr')
dev.off()

pdf(file=paste(series,"_b.pdf",sep=""),width=6,height=6)
plot.f(a,field='sens')
dev.off()

pdf(file=paste(series,"_c.pdf",sep=""),width=6,height=6)
plot.f(a,field='cor')
dev.off()

sub.data = subset(sub.data,slrsim_label == '2|0.1')


# TODO: Plot some alignments from here.
#f <- function(df) {
#  lbl <- df[1,]$slrsim_label
#  p <- plot.indel.distribution(lbl)
#  
#}
#do.by.labels(data.sub,fn=f)
