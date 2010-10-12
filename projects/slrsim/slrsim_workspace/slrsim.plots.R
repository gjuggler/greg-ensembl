source("aln-tools/aln.tools.R")
source("aln-tools/phylo.tools.R")
source("aln-tools/plot.phylo.greg.R")
library(ape)

subplot <- function(x, y) viewport(layout.pos.col=x, layout.pos.row=y)
vplayout <- function(x, y) {
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(y,x)))
}

dump.protein = function(node_id,base=node_id,tree_file=NA,aln_file=NA,sw_file=NA,params_file=NA,sw_table=NA,dbname=NA) {
  if (is.na(tree_file)) {
    tree_file = paste(base,'/',node_id,".nh",sep="")
  }
  if (is.na(aln_file)) {
    aln_file = paste(base,'/',node_id,".fa",sep="")
  }
  if (is.na(sw_file)) {
    sw_file = paste(base,'/',node_id,".csv",sep="")
  }
  if (is.na(params_file)) {
    params_file = paste(base,'/',node_id,".txt",sep="")
  }
  if (is.na(sw_table)) {
    sw_table = "stats_sites"
  }

  url = paste("mysql://slrsim:slrsim@mysql-greg.ebi.ac.uk:4134/",
    dbname,sep="")
  print(paste("URL:",url))
  system(paste("perl ~/lib/greg-ensembl/scripts/tree_dump.pl",
    " --url=",url,
    " --id=",node_id,
    " --tree=",tree_file,
    " --aln=",aln_file,
    " --sw=",sw_file,
    " --sw_table=",sw_table,
    " --params=",params_file,
    sep="")
  )  
}

plot.indel.distribution <- function(indel_label) {
  library(VGAM)
  toks <- unlist(strsplit(indel_label,' '))
  if (toks[1] == 'NB') {
    # for NB .5 1, a=1 and b=1 - .5
    distr <- dnbinom
    a <- toks[3]
    b <- 1-as.numeric(toks[2])
  } else {
    # for POW 4 30, a=30 and b=4
    distr <- dzipf
    a <- toks[3]
    b <- toks[2]
  }
  a <- as.numeric(a)
  b <- as.numeric(b)  

  xs <- seq(from=1,to=20)
  ys <- distr(xs,a,b)
  df <- data.frame(x=xs,y=ys)
  p <- ggplot(df,aes(x=x,y=y)) + geom_bar(stat='identity')
  return(p)
}

no.margins <- function(p) {
  p <- p + opts(
      plot.margin=unit(c(0,0,0,0),'lines'),
      panel.grid.minor=theme_blank(),
      panel.grid.major=theme_blank(),
      legend.position='none',
      axis.title.x=theme_blank(),
      axis.title.y=theme_blank()
    )
  p <- p + scale_x_continuous(breaks=NA)
  p <- p + scale_y_continuous(breaks=NA)
  return(p)
}

plot.aln.overview = function(node_id,base,len=1000) {
  t_f = paste(base,'/',node_id,".nh",sep="")
  a_f = paste(base,'/',node_id,".fa",sep="")
  csv_f = paste(base,'/',node_id,".csv",sep="")

  tree = read.tree(t_f)
  aln = read.aln(a_f,seqType='protein')
  aln = sort.aln.by.tree(aln,tree)

  return(plot.aln.gg(aln,do.plot=F))
}

aln.df <- function(aln) {
  arr <- matrix(unlist(aln$seqs),nrow=aln$num_seqs,ncol=aln$length,byrow=TRUE) 
  df <- data.frame()
  for (i in 1:aln$num_seqs)
  {
    list <- arr[i,]
    name <- aln$names[i]

    # Store the indices where the gaps are.
    ids <- list == '-'
    seq.pos <- seq(1,aln$length)

    pos.nogaps <- seq.pos[ids==FALSE]
    char.nogaps <- list[ids==FALSE]
    df <- rbind(df,data.frame(
      id=rep(name,length(pos.nogaps)),
      seq_index=rep(i,length(pos.nogaps)),
      pos=pos.nogaps,
      char=char.nogaps
    ))
  }
  df$id <- factor(df$id,levels=aln$name)
  return(df)
}

plot.aln.gg <- function(aln,do.plot=T) {
  df <- aln.df(aln)

  p <- ggplot(df,aes(x=pos,y=id))
  p <- p + geom_tile(aes(fill=char))
  p <- p + scale_fill_manual(values=protein.colors())
#  p <- p + coord_equal(ratio=5)
  p <- p + opts(legend.position='none')
  if(do.plot == TRUE) {
    print(p)
  } else {
    return(p)
  }
}

plot.roc = function(data,plot.x='tn',plot.y='tp',plot.unity=F,fill.below=F,plot=T) {
  data$roc.x = data[,plot.x]
  data$roc.y = data[,plot.y]

  require(ggplot2)
  require(grid)
  print(paste("Plot.roc row count: ",nrow(data)))

  p <- ggplot(data,aes(x=roc.x,y=roc.y,colour=slrsim_label))
  p <- p + geom_line()
  p <- p + xlab(plot.x)
  p <- p + ylab(plot.y)
  if (plot) {
    print(p)
  }
  return(p)
}


protein.colors <- function(scheme='Taylor') {
  if (scheme == 'Taylor') {
    cols <- c(
       'A' = "#CCFF00",       'a' = "#CCFF00",
       'C' = "#FFFF00",       'c' = "#FFFF00",
       'D' = "#FF0000",       'd' = "#FF0000",
       'E' = "#FF0066",       'e' = "#FF0066",
       'F' = "#00FF66",       'f' = "#00FF66",
       'G' = "#FF9900",       'g' = "#FF9900",
       'H' = "#0066FF",       'h' = "#0066FF",
       'I' = "#66FF00",       'i' = "#66FF00",
       'K' = "#6600FF",       'k' = "#6600FF",
       'L' = "#33FF00",       'l' = "#33FF00",
       'M' = "#00FF00",       'm' = "#00FF00",
       'N' = "#CC00FF",       'n' = "#CC00FF",
       'P' = "#FFCC00",       'p' = "#FFCC00",
       'Q' = "#FF00CC",       'q' = "#FF00CC",
       'R' = "#0000FF",       'r' = "#0000FF",
       'S' = "#FF3300",       's' = "#FF3300",
       'T' = "#FF6600",       't' = "#FF6600",
       'V' = "#99FF00",       'v' = "#99FF00",
       'W' = "#00CCFF",       'w' = "#00CCFF",
       'Y' = "#00FFCC",       'y' = "#00FFCC",
       '2' = "#888888",       '2' = "#888888",
       'O' = "#424242",       'o' = "#424242",
       'B' = "#7D7D7D",       'b' = "#7D7D7D",
       'Z' = "#EEEEEE",       'z' = "#EEEEEE",
       'X' = "#000000",       'x' = "#000000"
    )
  } else {
    cols <- c(
       'A' = "#B8B8B8",       'a' = "#B8B8B8",
       'C' = "#E6E600",       'c' = "#E6E600",
       'D' = "#E60A0A",       'd' = "#E60A0A",
       'E' = "#E60A0A",       'e' = "#E60A0A",
       'F' = "#3232AA",       'f' = "#3232AA",
       'G' = "#C8C8C8",       'g' = "#C8C8C8",
       'H' = "#8282D2",       'h' = "#8282D2",
       'I' = "#0F820F",       'i' = "#0F820F",
       'K' = "#145AFF",       'k' = "#145AFF",
       'L' = "#0F820F",       'l' = "#0F820F",
       'M' = "#E6E600",       'm' = "#E6E600",
       'N' = "#00DCDC",       'n' = "#00DCDC",
       'P' = "#DC9682",       'p' = "#DC9682",
       'Q' = "#E60A0A",       'q' = "#E60A0A",
       'R' = "#145AFF",       'r' = "#145AFF",
       'S' = "#FA9600",       's' = "#FA9600",
       'T' = "#FA9600",       't' = "#FA9600",
       'V' = "#C8C8C8",       'v' = "#C8C8C8",
       'W' = "#C8C8C8",       'w' = "#C8C8C8",
       'Y' = "#C8C8C8",       'y' = "#C8C8C8",
       '2' = "#888888",       '2' = "#888888",
       'O' = "#424242",       'o' = "#424242",
       'B' = "#7D7D7D",       'b' = "#7D7D7D",
       'Z' = "#EEEEEE",       'z' = "#EEEEEE"
    )
  }
  return(cols)
}
