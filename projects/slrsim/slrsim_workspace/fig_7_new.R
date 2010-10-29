library(ggplot2)
library(xtable)

folder = '~/lib/greg-ensembl/projects/slrsim/NO.backup/output/2010-10-28/2010-10-28_01'
series = paste(folder,'/fig_two',sep='')

na.rm <- TRUE
plot.x <- 'fpr'
plot.y <- 'tpr'
zoom.fpr <- 0.1

if(typeof(data) != 'list') {
  load(paste(series,".Rdata",sep=""))
}
source("slrsim.functions.R")
source("slrsim.plots.R")

print(unique(data$slrsim_label))

for (aln in unique(data$alignment_name)) {
  if (aln == 'True Alignment') {next;}
  d <- subset(data,alignment_name==aln | alignment_name=='True Alignment')
  d <- subset(d,alignment_score_threshold == 9)

  d$slrsim_label <- paste(d$slrsim_tree_file,d$alignment_name,d$filtering_name,d$alignment_score_threshold)

  # ROC plot of alignments.
  f = function(df,thresh) {
    return(slr.roc(df,na.rm=na.rm))
  }
  comb.roc <- summarize.by.labels(d,f)

  f2 = function(df,thresh) {
    roc <- slr.roc(df,na.rm=na.rm)
    auc.full <- area.under.curve(roc,x.lim=1)
    return(auc.full)
  }
  aucs <- summarize.by.labels(d,f2)
  print("Areas:");
  print(aucs)

  max.x <- max(comb.roc[,plot.y]) * zoom.fpr

  sub.roc <- subset(comb.roc,alignment_score_threshold == 9)

  sub.roc$slrsim_label <- paste(sub.roc$alignment_name,sub.roc$filtering_name,sub.roc$alignment_score_threshold)

  p <- plot.roc(sub.roc,plot=F,plot.x=plot.x,plot.y=plot.y)
  p <- p + scale_colour_brewer("Filtering Scheme",palette="Set1")

  fdr.line <- 0.2
  p <- p + geom_abline(slope=2/fdr.line,colour='gray')
  fdr.line <- 0.1
  p <- p + geom_abline(slope=2/fdr.line,colour='black')

  p <- p + facet_grid(. ~ slrsim_tree_file )

  pdf(file=paste(series,aln,"roc.pdf",sep="_"),width=20,height=10)
  print(p)
  dev.off()
  pdf(file=paste(series,aln,"roc_zoom.pdf",sep="_"),width=20,height=10)

  p <- p + scale_x_continuous(limits=c(0,max.x))
  print(p)
  dev.off()

#  pdf(file=paste(series,aln,"roc_zoom_thresholds2.pdf",sep="_"),width=50,height=10)

#  thresholds <- unique(comb.roc$alignment_score_threshold)
#  vplayout(length(thresholds),1)
#  i <- 1
#  for (filter_thresh in thresholds) {
#    sub.roc <- subset(comb.roc,alignment_score_threshold==filter_thresh)
#    #p <- p + facet_grid(. ~ alignment_score_threshold)
#    p <- plot.roc(sub.roc,plot=F,plot.x=plot.x,plot.y=plot.y)
#    fdr.line <- 0.2
#    p <- p + geom_abline(slope=2/fdr.line,colour='gray')
#    p <- p + scale_colour_brewer("Filtering Scheme",palette="Set1")
#    p <- p + scale_x_continuous(limits=c(0,max.x))
#    print(p,vp=subplot(i,1))
#    i <- i + 1
#  }
#  dev.off()

}
