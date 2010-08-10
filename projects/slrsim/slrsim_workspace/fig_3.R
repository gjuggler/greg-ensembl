library(ggplot2)

series = ''
fig_num = 'slrsim_three'

load(paste(fig_num,series,".Rdata",sep=""))
source("slrsim.functions.R")
source("slrsim.plots.R")

d <- NULL
for (rate in c(0, 0.05)) {
  data.sub = subset(data,phylosim_insertrate==rate)
#  data.sub = subset(data.sub,slrsim_label %in% c('NB .5 2','NB .5 3'))

  if (nrow(data.sub) == 0) {
	  next
  }
  d <- summarize.by.labels(data.sub,fn=fig.3.summary,thresh=3.8)
  pdf(file="fig_3_alns.pdf",width=15,height=20)
    
  n <- nrow(d)
  n.cols = 8
  i <- 1
  vplayout(n.cols,n*2)
  f <- function(df) {
    df.sum <- fig.3.summary(df,thresh=3.8)
    dbname='slrsim_1'
    node_id <- df[1,]$node_id
    dump.protein(node_id=node_id,base='fig_3_plots',dbname=dbname)
    
    p <- plot.aln.overview(node_id=node_id,base='fig_3_plots')
    p <- p + opts(
      legend.position='none',
      axis.text.y = theme_text(size=4),
      axis.title.x = theme_blank(),
      axis.title.y = theme_blank()
    )
    p <- p + coord_cartesian(xlim=c(0,700))
    print(p,vp=subplot(c(1:(n.cols-2)),c(i,i+1)))

    # Plot the indel distribution.
    lbl <- df[1,]$slrsim_label
    p <- plot.indel.distribution(lbl)
    p <- no.margins(p)
    p <- p + scale_x_continuous()
    print(p,vp=subplot(n.cols-1,i))

    # Indel distribution params.
    p <- ggplot(df.sum,aes(x=1,y=1))
    p <- p + geom_text(size=8,aes(label=slrsim_label))
    p <- no.margins(p)
    print(p,vp=subplot(n.cols,i))

    # FDR as color.
    p <- ggplot(df.sum,aes(x=1,y=1,fill=fdr)) + geom_tile()
    p <- no.margins(p)
    p <- p + scale_fill_gradient(name='fdr',limits=c(0,1))
    print(p,vp=subplot(n.cols-1,i+1))
    
    # FDR as text.
    p <- ggplot(df.sum,aes(x=1,y=1))
    p <- p + geom_text(size=8,aes(label=sprintf("%.2f",fdr)))
    p <- no.margins(p)
    print(p,vp=subplot(n.cols,i+1))

    print(i)

    # Tick up the counter.
    assign('i',i+2,envir=.GlobalEnv)
  }
  do.by.labels(data.sub,fn=f)

  dev.off()
}