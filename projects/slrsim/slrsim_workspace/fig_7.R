library(ggplot2)
library(xtable)

series = ''
fig_num = 'slrsim_seven'

load(paste(fig_num,series,".Rdata",sep=""))
source("slrsim.functions.R")
source("slrsim.plots.R")

d <- NULL
for (rate in c(0, 0.05)) {
  data.sub = subset(data,phylosim_insertrate==rate)
  if (nrow(data.sub) == 0) {
    next
  }

  # Summary table of alignments.
  d <- summarize.by.labels(data.sub,fn=fig.4.summary,thresh=3.8)
  tbl = subset(d,select=c('slrsim_label','pos_true','pos_aln','sens','spec','fpr','fdr','cor'))
  xtbl = xtable(tbl,digits=4)
  print(xtbl,type='html',file=paste(fig_num,".html",sep=""),html.table.attributes='style="border:0px;"')
  
  # ROC plot of alignments.
  pdf(file=paste(fig_num,"_roc.pdf",sep=""),width=10,height=10)
  f = function(df,thresh) {
    return(slr.roc(df,na.rm=T))
  }
  comb.roc = summarize.by.labels(data.sub,f)
  p <- plot.roc(comb.roc,plot=F,plot.x='fpr',plot.y='tpr')
  p <- p + scale_colour_discrete('Alignment program')
  p <- p + scale_x_continuous(limits=c(0,0.2))
  print(p)
  dev.off()

  # Protein plots of alignments.
  f <- function(df) {
    df.sum <- fig.4.summary(df,thresh=3.8)
    dbname='slrsim_1'

    # Important -- order by node ID so we get the first rep from each set of replicates.
    df = orderBy(~node_id,data=df)
    node_id <- df[2,]$node_id
    print(node_id)
    dump.protein(node_id=node_id,base=paste(fig_num,"_plots",sep=""),dbname=dbname)

    p <- plot.aln.overview(node_id=node_id,base=paste(fig_num,"_plots",sep=""))
    p <- p + opts(
      legend.position='none',
      axis.text.y = theme_text(size=3),
      axis.title.x = theme_blank(),
      axis.title.y = theme_blank()
    )
    p <- p + opts(title=df[1,'slrsim_label'])
    p <- p + coord_cartesian(xlim=c(0,800))
    p <- p + coord_equal(ratio=3)

    print(p,vp=subplot(c(1),c(i)))
    print(i)

    # Tick up the counter.
    assign('i',i+1,envir=.GlobalEnv)
  }

  n <- nrow(d)
  n.species = max(data.sub$true_ncod,na.rm=T)
  n.cols = max(data.sub$seq_position,na.rm=T)
  desired.aspect = (n * n.species) / (n.cols)
  print(paste("Aspect ratio:",desired.aspect,n,n.species,n.cols))

  width = 10
  height = desired.aspect * width

  print(paste("width x height:",width,height))
  pdf(file=paste(fig_num,"_alns.pdf",sep=""),width=20,height=15)
  n.cols = 1
  i <- 1
  vplayout(n.cols,n)
  do.by.labels(data.sub,fn=f)
  dev.off()
}
