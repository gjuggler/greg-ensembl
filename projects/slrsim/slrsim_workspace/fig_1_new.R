library(ggplot2)
library(fields)
library(RColorBrewer)


plot.f = function(data,palette="Spectral",field='fdr',rev.colors=F,limits=c(0,1),x.lim=c(0,4),y.lim=c(0,0.1),label.species=T,do.plot=T) {
  data$z_vals = data[,field]
  p <- ggplot(data,aes(x=tree_mean_path,y=phylosim_insertrate*2,z=z_vals))
  p <- p + scale_x_continuous('Mean Path Length')
  p <- p + scale_y_continuous('Indel Rate')

  if (!is.null(x.lim)) {
    p <- p + coord_cartesian(xlim = x.lim,ylim=y.lim)
  }

  p <- p + geom_tile(aes_string(fill=field))

  if (is.null(limits)) {
    limits <- c(min(data$z_vals),max(data$z_vals))
  }
  n <- 10
  breaks <- seq(from=limits[1],to=limits[2],length.out=n+1)
  if (rev.colors) {
    p <- p + scale_fill_gradientn(colours=rev(brewer.pal(n=n,palette)),limits=limits)
  } else {
    p <- p + scale_fill_gradientn(colours=brewer.pal(n=n,palette),limits=limits)
  }

  p <- p + coord_equal(ratio=10)

  if (label.species) {
    lbls <- data.frame(x=NULL,y=NULL,str=NULL)
    lbls <- rbind(lbls,data.frame(x=0.13,y=0.04,str="Primates"))
    lbls <- rbind(lbls,data.frame(x=0.25,y=0.03,str="Mammals"))
    lbls <- rbind(lbls,data.frame(x=0.75,y=0.023,str="Vertebrates"))
    lbls <- rbind(lbls,data.frame(x=1.0,y=0.015,str="Drosophilids"))
    lbls <- rbind(lbls,data.frame(x=1.25,y=0.008,str="Yeasts"))
    lbls$z_vals <- 0
    lbls$slope <- 4
    arw <- arrow(length=unit(0.05,"inches"))
    p <- p + geom_segment(data=lbls,mapping=aes(x=x+slope*y,xend=x,y=y/2,yend=0),colour='white',size=0.2,arrow=arw)
    p <- p + geom_text(data=lbls,mapping=aes(x=x+slope*y,y=y/2,label=str),hjust=0,vjust=0,size=2,colour='white')

    orig.lengths = c(0.5)
    tree.file <- data[1,]$slrsim_tree_file
    if (tree.file == 'anisimova_01_bglobin.nh') {
      orig.lengths <- c(0.11, 0.64)
    } else if (tree.file == 'anisimova_01_artificial.nh') {
      orig.lengths <- c(0.025, 0.25)
    } else if (tree.file == '44mammals.nh') {
      orig.lengths <- c(1.569)
    }
    pts <- data.frame(x=orig.lengths,y=-0.002,z_vals=0,shape=1)
    p <- p + geom_point(data=pts,aes(x=x,y=y),shape=5,size=2,colour='red')
  }
  
  if(do.plot) {
    print(p)
  } else {
    return(p)
  }
}

source("slrsim.functions.R")
source("slrsim.plots.R")


folder <- 'NO.backup/2010-10-06/'
#folder <- '~/lib/greg-ensembl/projects/slrsim/NO.backup/output/2010-10-21/2010-10-21_03/'
series <- paste(folder,'/slrsim_one_new',sep='')
all.df.file <- paste(folder,'/df.list.Rdata',sep='')
if(!exists("df.list")) {
  load(all.df.file)
}

test.df.list <- list(df.list[[1]],df.list[[2]],df.list[[3]])

if(!exists("summary.df")) {
  df.f <- function(x) {
    slr.threshold <- 3.84
    a <- summarize.by.labels(x,fn=fig.1.summary,thresh=slr.threshold)
    return(a)
  }
  summary.df <- ldply(df.list,df.f)
}

bg <- 'anisimova_01_bglobin.nh'
art <- 'anisimova_01_artificial.nh'
mam <- '44mammals.nh'
stupid.df <- data.frame(
  slrsim_tree_file = c(bg,bg,art,art,mam),
  length = c(0.11,0.64,0.025,0.25,1.569),
  z_vals = ''
)

sub.df <- summary.df

tree.files <- unique(sub.df$slrsim_tree_file)
tree.labels <- c('6 artificial','17 bglobin', '44 vertebrate')
sub.df[sub.df$slrsim_tree_file == art,]$slrsim_tree_file = 0
sub.df[sub.df$slrsim_tree_file == bg,]$slrsim_tree_file = 1
sub.df[sub.df$slrsim_tree_file == mam,]$slrsim_tree_file = 2
sub.df$slrsim_tree_file <- factor(sub.df$slrsim_tree_file,labels=tree.labels)

str(sub.df)

pdf(file=paste(series,"_auc.pdf",sep=""),width=10,height=5)
p <- plot.f(sub.df,field='auc',do.plot=F,label.species=F,rev.colors=F,limits=c(0.4,1))
p <- p + facet_grid(slrsim_tree_file ~ alignment_name)
print(p)
dev.off()

pdf(file=paste(series,"_fdr.pdf",sep=""),width=10,height=5)
p <- plot.f(sub.df,field='fdr',do.plot=F,label.species=F,rev.colors=T)
p <- p + facet_grid(slrsim_tree_file ~ alignment_name)
print(p)
dev.off()

pdf(file=paste(series,"_sens.pdf",sep=""),width=10,height=5)
p <- plot.f(sub.df,field='sens',do.plot=F,label.species=F,rev.colors=F)
p <- p + facet_grid(slrsim_tree_file ~ alignment_name)
print(p)
dev.off()

#pdf(file=paste(series,"_labels.pdf"),width=12,height=4)
#vplayout(3,1)
#i <- 1
#for (aln in c('44mammals.nh','anisimova_01_artificial.nh','anisimova_01_bglobin.nh')) {
#  sub.df <- subset(summary.df,slrsim_tree_file==aln & alignment_name=='prank_codon')
#  p <- plot.f(sub.df,field='fdr',do.plot=F,label.species=T)
#  print(p,vp=subplot(i,1))
#  i <- i + 1
#}
#dev.off()