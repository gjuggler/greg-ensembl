library(ggplot2)
library(xtable)

series = ''
fig_num = 'slrsim_four'

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
  tbl = subset(d,select=c('slrsim_label','pos_true','pos_aln','sens','spec','fpr','fdr'))
  xtbl = xtable(tbl,digits=4)
  print(xtbl,type='html',file='fig_4.html',html.table.attributes='style="border:0px;"')

  pdf(file="fig_4_omegas.pdf",width=10,height=20)

  # Create a cell-based layout.
  n.rows <- nrow(d)
  n.cols = 3
  i <- 1
  vplayout(n.cols,n.rows)

  if (nrow(data.sub) == 0) {
	  next
  }

  f <- function(df) {
    ds <- fig.4.summary(df,thresh=3.8)
    
    # Load up the distribution params (copied from LoadTrees.pm)
    lbl = ds$slrsim_label
    dist.map = list(
      'lognormal_wide_2'      = c(-1.25,0.45),
      'lognormal_wide'        = c(-1.5,0.84),
      'lognormal'             = c(-1.864,1.2007),
      'lognormal_narrow'      = c(-2.25,1.49),
      'lognormal_narrow_2'    = c(-3,1.925),
      'lognormal_low'         = c(-2.803,1.0),
      'lognormal_low_2'         = c(-4,1.84),
      'lognormal_high'        = c(-1.819,1.5),
      'lognormal_high_2'        = c(-1.0,0.784),
      'uniform_neutral'       = c(1),
      'uniform_low'           = c(0.3),
      'uniform_low2'          = c(0.1),
      'uniform_pos'           = c(1.5),
      'uniform_pos2'           = c(2)      
    )
    params <- dist.map[[lbl]]
    print(params)
    if (length(grep('lognormal',lbl))) {
      meanlog <- params[1]
      sdlog <- params[2]
      dist.f <- function(xs) dlnorm(xs,meanlog,sdlog)
    } else {
      w = params[1]
      dist.f <- function(xs) dunif(xs,min=w-.01,max=w+.01)
    }

    # Plot continuous distribution.
    bin.width = 0.02
    xs <- seq(from=0,to=2,by=bin.width)
    ys <- dist.f(xs)
    print(xs)
    print(ys)
    xy.df <- data.frame(x=xs,y=ys)
    p <- ggplot(xy.df,aes(x=x,y=y)) + geom_line()
    p <- p + scale_x_continuous('dnds')
    p <- p + scale_y_continuous('density')
    p <- p + opts(title = lbl)

    print(p,vp=subplot(c(1),i))
    
    # Plot true dnds.
    true.mean <- mean(df$true_dnds)
    p <- ggplot(df,aes(x=true_dnds))
    p <- p + geom_histogram(breaks=xs,width=bin.width,fill='black',aes(y=..count..))
    p <- p + scale_x_continuous(limits=c(0,2))
    p <- p + geom_vline(xintercept=mean(true.mean))
    print(p,vp=subplot(c(2),i))

    # Plot infered dnds for non-constant/synonymous sites.
    df.sub = df
#    df.sub = subset(df,!(aln_note %in% c('constant')))
    aln.mean <- mean(df.sub$aln_dnds)
    p <- ggplot(df.sub,aes(x=aln_dnds))
    p <- p + geom_histogram(breaks=xs,width=bin.width,fill='gray',aes(y=..count..))
#    p <- p + geom_density(colour='red',adjust=3,aes(y=..scaled..),na.rm=T)
    p <- p + scale_x_continuous(limits=c(0,2))
    p <- p + geom_vline(xintercept=aln.mean)
#    p <- p + coord_trans(y="sqrt")
    print(p,vp=subplot(c(3),i))

    # Tick up the counter.
    print(i)
    assign('i',i+1,envir=.GlobalEnv)    
  }
  do.by.labels(data.sub,fn=f)

  dev.off()
}