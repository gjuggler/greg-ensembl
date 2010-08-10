library(ggplot2)

series = 'b'
fig_num = 'fig_2'

load(paste(fig_num,series,".Rdata",sep=""))
source("slrsim.functions.R")
str(data)

for (rate in c(0, 0.025, 0.05, 0.1)) {
  data.sub = subset(data,phylosim_insertrate==rate)
  if (nrow(data.sub) == 0) {
	  next
  }
  d = summarize.by.labels(data.sub,fn=fig.2.summary,thresh=0)

  # Split the labels up by commas.
  a.counter = ddply(d,'slrsim_label',function(df) {z = unlist(strsplit(df$slrsim_label,','));z[2]})
  d$counter = a.counter[,'V1']

  a.species = ddply(d,'slrsim_label',function(df) {z = unlist(strsplit(df$slrsim_label,','));z[3]})
  #d$species = substr(a.species[,'V1'],1,5)
  d$species = a.species[,'V1']

  x.field = 'sens'
  y.field = 'spec'	
  z.field = 'min.fdr'
  d$x_vals = d[,x.field]
  d$y_vals = d[,y.field]
  d$z_vals = d[,z.field]

  pdf(file=paste(fig_num,series,rate,".pdf",sep=""))
  p <- ggplot(d,aes(x=x_vals,y=y_vals))
#  p <- p + coord_cartesian(xlim = c(0,.8),ylim=c(0.5,1))
  p <- p + geom_text(size=3,aes(label=counter),hjust=1)
  p <- p + geom_text(size=2,aes(label=species),hjust=0)
#  p <- p + scale_colour_gradient(name="fdr",limits=c(0,1))
  p <- p + scale_x_continuous(x.field) + scale_y_continuous(y.field) 
  print(p)
  dev.off()
}
