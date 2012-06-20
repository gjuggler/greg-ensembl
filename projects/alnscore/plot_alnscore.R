library(ggplot2)
library(phylosim)
library(plyr)

if (Sys.getenv('USER') == 'gj1') {
  source("~/src/greg-ensembl/scripts/mysql_functions.R")
} else {
  source("~/lib/greg-ensembl/scripts/mysql_functions.R")
}

get.data <- function() {
  dbname <- 'gj1_alnscore'
  con <- connect(dbname)

  ff <- scratch.f('alnscore.Rdata')
  if (!file.exists(ff)) {
    print("Getting data...")
    aln <- mysqlReadTable(con, 'aln')
    save(aln, file=ff)
    print("Done!")
  } else {
    load(ff)
  }

  aln
}

scores.scatter <- function() {
  aln <- get.data()

  cur.df <- aln
#  cur.df$n_seqs <- factor(cur.df$n_seqs)

  stats <- c('sub_tcs', 'sub_sps', 'sub_pbl', 'sub_cbl')
  comb.df <- data.frame()
  for (stat in stats) {
    new.df <- cur.df
    new.df$y.val <- new.df[, stat]
    new.df$stat <- stat
    comb.df <- rbind(comb.df, new.df)
  }

  cur.df <- subset(comb.df, aligner=='prank' & mpl > 0.9 & mpl < 1.1)
  pdf(file="n_seqs_scatter.pdf", width=10, height=5)
  p <- ggplot(cur.df, aes(x=n_seqs, y=y.val, colour=stat))
  p <- p + geom_point(size=1, alpha=1)
  p <- p + stat_smooth(aes(group=stat, colour=stat), method='lm')
  p <- p + facet_grid(. ~ tree)
  print(p)
  dev.off()

  cur.df <- subset(comb.df, n_seqs %in% c(2))
  cur.df$aligner <- paste(cur.df$aligner, cur.df$n_seqs, sep=' ')

  pdf(file="aln_scatter.pdf", width=10, height=10)
  p <- ggplot(cur.df, aes(x=mpl, y=y.val, colour=aligner))
  p <- p + geom_point(size=1, alpha=1)
  p <- p + stat_smooth(aes(group=aligner, colour=aligner), method='loess')
  p <- p + facet_grid(stat ~ tree)
  print(p)
  dev.off()

}

plot.bars <- function() {
  aln <- get.data()

  print(unique(aln$tree))
  print(unique(aln$aligner))
  print(unique(aln$n_domains))

  aln <- subset(aln, n_seqs == 2)
  aln$mpl <- round(aln$mpl * 10) / 10
  aln$mpl <- ifelse(aln$mpl > 4, 6, aln$mpl)
  aln$mpl <- factor(aln$mpl)

  print(nrow(aln))

  cur.df <- aln
  cur.df$tree <- factor(cur.df$tree)
  cur.df$n_domains <- factor(cur.df$n_domains)
  cur.df$n_seqs <- factor(cur.df$n_seqs)
  cur.df$aligner <- factor(cur.df$aligner)

  cur.df$lbl <- factor(paste(cur.df$tree, cur.df$n_domains, sep=' / '))

  relative.df <- ddply(cur.df, c('lbl', 'tree', 'aligner', 'mpl', 'n_seqs', 'n_domains'), function(x) {
    avg.diff <- mean(x$sub_sps - x$aln_sps)
    data.frame(
      y.val=avg.diff
    )
  })
  print(str(relative.df))

  relative.df$n_seqs <- factor(relative.df$n_seqs)

  pdf(file=scratch.f("aln_relative_bars.pdf"), width=10, height=5)
  p <- ggplot(relative.df, aes(x=aligner, y=y.val, fill=mpl))
  p <- p + theme_bw()
  p <- p + geom_bar(position="dodge", stat='identity')
  p <- p + scale_x_discrete("Aligner")
  p <- p + scale_y_continuous("Mean SPS score, (subsetted - realigned)")

  clr.f <- colorRampPalette(colors=c('blue', 'red'))
  n.colors <- length(unique(relative.df$mpl))
  clrs <- clr.f(n.colors)
  p <- p + scale_fill_manual(values=clrs)
  p <- p + facet_grid(lbl ~ .)
  p <- p + opts(
    strip.text.y = theme_text(angle=0, hjust=0)
  )
  print(p)

  dev.off()
}

plot.abs <- function() {

  abs.df <- ddply(cur.df, c('tree', 'aligner', 'mpl', 'n_seqs'), function(x) {
    avg.value <- mean(x$aln_sps)
    avg.hi <- quantile(x$aln_sps, probs=c(0.75))
    avg.lo <- quantile(x$aln_sps, probs=c(0.25))
    mean.bl <- mean(x$mean_bl)
    data.frame(
      mean.bl=mean.bl,
      y.val=avg.value,
      y.hi=avg.hi,
      y.lo=avg.lo
    )
  })

  abs.df <- subset(abs.df, tree == 'balanced.nh')  

  pdf(file="aln_abs_mean_bl.pdf", width=5, height=10)
  p <- ggplot(abs.df, aes(x=mean.bl, y=y.val, ymin=y.lo, ymax=y.hi, colour=n_seqs, group=n_seqs))
  p <- p + geom_point(size=1)
  p <- p + geom_errorbar(width=0.02, alpha=0.5)
  p <- p + geom_line()
  p <- p + facet_grid(aligner ~ tree)
  print(p)
  dev.off()

  pdf(file="aln_abs_mpl.pdf", width=5, height=10)
  p <- ggplot(abs.df, aes(x=mpl, y=y.val, ymin=y.lo, ymax=y.hi, colour=n_seqs, group=n_seqs))
#  p <- p + geom_point(size=1)
  p <- p + geom_errorbar(position='dodge')
#  p <- p + geom_line(alpha=0.7)
  p <- p + facet_grid(aligner ~ tree)
  print(p)
  dev.off()

  cur.abs.df <- subset(abs.df, n_seqs %in% c(2, 32))
  cur.abs.df$n_seqs <- factor(as.character(cur.abs.df$n_seqs))

  cur.abs.df$aligner <- reorder(cur.abs.df$aligner, cur.abs.df$y.val, mean)

  pdf(file="aln_abs_bars.pdf", width=5, height=10)
  p <- ggplot(cur.abs.df, aes(x=aligner, y=y.val, ymin=y.lo, ymax=y.hi, fill=n_seqs, colour=n_seqs))
  p <- p + geom_bar(stat='identity', position='dodge')
  p <- p + geom_errorbar(position='dodge')
  p <- p + facet_grid(mpl ~ tree)
  print(p)
  dev.off()


}

scratch.f <- function(file) {
  return(paste("~/scratch/gj1_alnscore/current/", file, sep=''))
}
