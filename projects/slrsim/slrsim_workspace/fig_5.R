library(ggplot2)
library(xtable)

series = ''
fig_num = 'slrsim_five'

load(paste(fig_num,series,".Rdata",sep=""))
source("slrsim.functions.R")
source("slrsim.plots.R")

d <- NULL
for (rate in c(0, 0.05)) {
  data.sub = subset(data,phylosim_insertrate==rate)
#  data.sub = subset(data,slrsim_label=='uniform_low')
  if (nrow(data.sub) == 0) {
    next
  }

  d <- summarize.by.labels(data.sub,fn=fig.4.summary,thresh=3.8)
  tbl = subset(d,select=c('slrsim_label','pos_true','pos_aln','sens','spec','fpr','fdr','cor'))
  xtbl = xtable(tbl,digits=4)
  print(xtbl,type='html',file='fig_5.html',html.table.attributes='style="border:0px;"')
  
  pdf(file="fig_5_roc.pdf",width=10,height=10)

  f = function(df,thresh) {
    return(slr.roc(df))
  }
  comb.roc = summarize.by.labels(data.sub,f)
  p <- plot.roc(comb.roc,plot=F)
  p <- p + scale_colour_discrete('Reference species')
  print(p)
  dev.off()
}