library(ggplot2)
library(phylosim)
library(plyr)
if (Sys.getenv('USER') == 'gj1') {
  source("~/src/greg-ensembl/projects/phylosim/PhyloSimPlots.R")
  source("~/src/greg-ensembl/scripts/mysql_functions.R")
} else {
  source("~/lib/greg-ensembl/projects/phylosim/PhyloSimPlots.R")
  source("~/lib/greg-ensembl/scripts/mysql_functions.R")
}

get.data <- function() {
  dbname <- 'gj1_alnscore'
  con <- connect(dbname)

  aln <- mysqlReadTable(con, 'aln')
  save(aln, file=scratch.f('alnscore.Rdata'))
  aln
}

scores.scatter <- function() {
  aln <- get.data()

  cur.df <- aln
#  cur.df$n_seqs <- factor(cur.df$n_seqs)

  stats <- c('sub_tcs', 'sub_sps', 'sub_pbl', 'sub_mmbl')
  comb.df <- data.frame()
  for (stat in stats) {
    new.df <- cur.df
    new.df$y.val <- new.df[, stat]
    new.df$stat <- stat
    comb.df <- rbind(comb.df, new.df)
  }

#  pdf(file="mpl_pbl.pdf", width=10, height=5)
#  p <- ggplot(comb.df, aes(x=mpl, y=y.val, colour=n_seqs))
#  p <- p + geom_point(size=1, alpha=1)
#  p <- p + scale_colour_brewer()
#  p <- p + facet_grid(stat ~ tree)
#  print(p)
#  dev.off()

#  pdf(file="total_bl_pbl.pdf", width=10, height=5)
#  p <- ggplot(comb.df, aes(x=total_bl, y=y.val, colour=n_seqs))
#  p <- p + geom_point(size=1, alpha=1)
#  p <- p + facet_grid(stat ~ tree)
#  print(p)
#  dev.off()

  cur.df <- subset(comb.df, aligner=='prank' & mpl > 0.9 & mpl < 1.1)
  pdf(file="n_seqs_scatter.pdf", width=10, height=5)
  p <- ggplot(cur.df, aes(x=n_seqs, y=y.val, colour=stat))
  p <- p + geom_point(size=1, alpha=1)
  p <- p + stat_smooth(aes(group=stat, colour=stat), method='lm')
  p <- p + facet_grid(. ~ tree)
  print(p)
  dev.off()

  cur.df <- subset(comb.df, n_seqs %in% c(2, 32))
  cur.df$aligner <- paste(cur.df$aligner, cur.df$n_seqs, sep=' ')

  pdf(file="aln_scatter.pdf", width=10, height=10)
  p <- ggplot(cur.df, aes(x=mpl, y=y.val, colour=aligner))
  p <- p + geom_point(size=1, alpha=1)
  p <- p + stat_smooth(aes(group=aligner, colour=aligner), method='loess')
  p <- p + facet_grid(stat ~ tree)
  print(p)
  dev.off()

  pdf(file="aln_scatter_total_bl.pdf", width=10, height=10)
  p <- ggplot(cur.df, aes(x=log(total_bl/(2*n_seqs-2)), y=y.val, colour=aligner))
  p <- p + geom_point(size=1, alpha=1)
  p <- p + stat_smooth(aes(group=aligner, colour=aligner), method='loess')
  p <- p + facet_grid(stat ~ tree)
  print(p)
  dev.off()
  

}

plot.bars <- function() {
  aln <- get.data()

  cur.df <- aln
  cur.df <- subset(cur.df, n_seqs %in% c(2, 16))
  cur.df$n_seqs <- factor(cur.df$n_seqs)
  cur.df$aligner <- factor(cur.df$aligner)

  #print(unique(cur.df$mpl))  
  cur.df$mpl <- round(cur.df$mpl * 10) / 10
  print(unique(cur.df$mpl))  
  cur.df$mpl <- factor(cur.df$mpl)

  cur.df <- ddply(cur.df, c('tree', 'aligner', 'mpl', 'n_seqs'), function(x) {
    avg.diff <- mean(x$sub_sps - x$aln_sps)
    data.frame(
      y.val=avg.diff
    )
  })

  pdf(file="aln_sub.pdf", width=10, height=5)
  p <- ggplot(cur.df, aes(x=aligner, y=y.val, fill=mpl))
  p <- p + geom_bar(position="dodge", stat='identity')
  p <- p + facet_grid(n_seqs ~ tree)
  print(p)
  dev.off()

}

scratch.f <- function(file) {
  return(paste("~/scratch/gj1_alnscore/2011-06-22_03/", file, sep=''))
}
