library(R.oo)
library(ape)
library(ggplot2)
source("~/src/greg-ensembl/projects/phylosim/PhyloSimSource.R")
source("~/src/greg-ensembl/projects/phylosim/examples/greg_scripts.R")

plot.sweeps <- function(prefix='true_0.01_') {

  indels <- c(0,0.01,0.02,0.1)
  lengths <- seq(from=0.1,to=2,by=0.1)

  for (i in c(1:length(lengths))) {
    file.base <- paste(prefix,lengths[i],sep='')
    print(file.base)
    aln.f <- paste(file.base,'.fasta',sep='')
    tree.f <- paste(file.base,'.nh',sep='')
    
    png(paste(file.base,".png",sep=''),width=1400,height=200)
    ggplot.aln(aln.f,tree.f,
      plot.tree=T,tree.xlim=c(0,2.5),aln.xlim=c(0,200),
      num.pages=1,plot.chars=F
    )
    dev.off()
  }
  
  make.mov(prefix)
}

make.mov <- function(prefix){
  unlink("movie.gif")
  system(paste("convert -delay 0.5 ",prefix,"*.png movie.gif",sep=''))
}

plot.all.in.dir <- function(
  aln.dir,
  ...
) {

  files <- dir(aln.dir)  
  files <- files[grep("\\.(fa|fasta)$",files)]
  for (file in files) {
    file <- paste(aln.dir,file,sep="/")
    file.base <- gsub("\\.[^\\.]+$","",file)
    tree.f <- paste(file.base,".nh",sep="")
    if(!file.exists(tree.f)) {
      tree.f=NULL
    }
    out.f <- paste(file.base,".png",sep="")    
    print(paste("tree:",tree.f))
    print(paste("aln:",file))
    print(paste("out:",out.f))

    png(out.f,width=1400,height=200)
    ggplot.aln(file,tree.f,
      plot.tree=T,tree.xlim=c(0,1.5),aln.xlim=c(0,120),
      num.pages=1,plot.chars=F,axis.text.size=3
    )
    dev.off()
  }
}
