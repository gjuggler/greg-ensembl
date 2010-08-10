library(ggplot2)
library(fields)

plot.f = function(data,low='blue',high='red',field='fdr',limits=c(0,1),x.lim=c(0,20),y.lim=c(0,1)) {
  data$z_vals = data[,field]
  p <- ggplot(data,aes(x=tree_mean_path,y=phylosim_insertrate,z=z_vals))

  if (!is.null(x.lim)) {
    p <- p + coord_cartesian(xlim = x.lim,ylim=y.lim)
  }
#  if (!is.null(y.lim)) {
#    p <- p + coord_cartesian(ylim = y.lim)
#  }

  p <- p + geom_tile(aes_string(fill=field))
  if (field %in% c('fdr','fpr','sens')) {
    p <- p + stat_contour(breaks=c(0.1,0.5,0.9))
  } else if (field %in% c('cor')) {
    p <- p + stat_contour(breaks=c(0.1,0.5,0.9))
  }


  if (!is.null(limits)) {
    p <- p + scale_fill_continuous(low=low,high=high,limits=limits)
  } else {
    p <- p + scale_fill_continuous(low=low,high=high)
  }

  print(p)
}

series = '_artificial'
source("slrsim.functions.R")
source("slrsim.plots.R")

load(paste("fig_1",series,".Rdata",sep=""))
sub.data = data
slr.threshold = 0
a = summarize.by.labels(sub.data,fn=fig.1.summary,thresh=slr.threshold)

pdf(file=paste("fig_1a",series,".pdf",sep=""),width=6,height=6)
plot.f(a,low='blue',high='red')
dev.off()

pdf(file=paste("fig_1b",series,".pdf",sep=""),width=6,height=6)
plot.f(a,low='white',high='blue',field='sens')
dev.off()

pdf(file=paste("fig_1c",series,".pdf",sep=""),width=6,height=6)
plot.f(a,low='red',high='yellow',field='cor')
dev.off()

sub.data = subset(sub.data,slrsim_label == '2|0.1')


# TODO: Plot some alignments from here.
f <- function(df) {
  lbl <- df[1,]$slrsim_label
  p <- plot.indel.distribution(lbl)
  
}
do.by.labels(data.sub,fn=f)