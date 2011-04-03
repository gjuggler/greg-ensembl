library(phylosim)
library(RColorBrewer)
library(xtable)
library(plyr)
source("~/src/greg-ensembl/projects/phylosim/PhyloSimPlots.R")
source("~/src/greg-ensembl/projects/slrsim/slrsim.functions.R")
source("~/src/greg-ensembl/projects/slrsim/slrsim.plots.R")

subplot <- function(x, y) viewport(layout.pos.col=x, layout.pos.row=y)
vplayout <- function(x, y) {
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(y,x)))
}

#
# Figure 1 - Power of three sitewise methods under zero indel rate.
#
multi.methods <- function() {
  tree.a <- read.tree('trees/artificial.nh')
  tree.b <- read.tree('trees/bglobin.nh')
  tree.c <- read.tree('trees/encode.nh') 
  tbl <- read.csv('~/scratch/gj1_fig_zero/current/table.csv')
  
  f <- function(sim) {
    p <- plotTree(sim, 
      tree.do.plot=F,
      line.width=0.75,
      axis.text.size=1.8,
      tree.xlim.expand=c(0.05, 0.4)
    )
    p <- p$grob
    p <- p + theme_bw()
    p <- p + opts(
      axis.title.x = theme_blank(),
      axis.title.y = theme_blank()
    )
    return(p)
  }

  plots <- list(
    a = list(tree.a),
    b = list(tree.b),
    c = list(tree.c)
  )

  pdf(file="fig_zero.pdf", width=20, height=8)
  vplayout(8, 3)

  for (i in 1:length(plots)) {
    cur.stuff <- plots[[i]]
    cur.tree <- cur.stuff[[1]]
    # Plot tree to the left.
    sim <- PhyloSim()
    setPhylo(sim, cur.tree)
    p <- f(sim)
    print(p, vp=subplot(1, i))
  }

  p.field <- function(field, lbl, keep.legend=F) {
    tbl[, 'z_value'] <- tbl[, field]
    p <- ggplot(tbl, aes(x=length, y=z_value, colour=analysis))
    p <- p + theme_bw()
    p  <- p + geom_line(alpha=0.8)
    p <- p + scale_colour_brewer(name="Analysis Method", palette="Set1")
    p <- p + scale_x_continuous(name="Tree Length", limits=c(0.1, 2), breaks=c(0.1, 1.0, 2.0))
    p <- p + scale_y_continuous(name=lbl, limits=c(0, 1), breaks = c(0, 0.5, 1))
    p <- p + facet_grid(tree ~ .)
    p <- p + opts(
      axis.title.x = theme_blank(),
      axis.title.y = theme_blank(),
      strip.text.x = theme_blank(),
      strip.text.y = theme_blank(),
      strip.background = theme_blank(),
      title = paste(lbl,sep='')
    )
    if (!keep.legend) {
      p <- p + opts(legend.position = 'none')
    }
    return(p)
  }

  p <- p.field('cor', 'Sitewise Correlation')
  print(p, vp=subplot(2, 1:3))
  
  p <- p.field('auc', 'AUC (FPR<0.1)')
  print(p, vp=subplot(3, 1:3))

  p <- p.field('tpr_at_fpr', 'TPR at FPR<0.05')
  print(p, vp=subplot(4, 1:3))

  p <- p.field('tpr_at_fdr', 'TPR at FDR<0.1')
  print(p, vp=subplot(5, 1:3))

  p <- p.field('tpr_at_thresh', 'TPR at Default Threshold')
  print(p, vp=subplot(6, 1:3))

  p <- p.field('tpr_at_bh', 'TPR at Adjusted Threshold')
  print(p, vp=subplot(7, 1:3))

  p <- p.field('tpr_at_fdr', '', keep.legend=T)
  print(p, vp=subplot(8, 1:3))

  dev.off()  
}

#
# Figure 2 - Tree shapes, TPR_FPR<0.05, and ROC plots at 1.0/0.1
#
multi.plot <- function() {
  tree.a <- read.tree('trees/artificial.nh')
  tree.b <- read.tree('trees/bglobin.nh')
  tree.c <- read.tree('trees/encode.nh')

  tbl.a <- read.csv('~/scratch/gj1_fig_one_a/current/table.csv')
  tbl.b <- read.csv('~/scratch/gj1_fig_one_b/current/table.csv')
  tbl.c <- read.csv('~/scratch/gj1_fig_one_c/current/table.csv')

  f <- function(sim) {
    p <- plotTree(sim, 
      tree.do.plot=F,
      line.width=0.75,
      axis.text.size=1.8,
      tree.xlim.expand=c(0.05, 0.4)
    )
    p <- p$grob
    p <- p + theme_bw()
    p <- p + opts(
      axis.title.x = theme_blank(),
      axis.title.y = theme_blank()
    )
    return(p)
  }

  plots <- list(
    a = list(tree.a, tbl.a),
    b = list(tree.b, tbl.b),
    c = list(tree.c, tbl.c)
  )

  pdf(file='fig_multi.pdf', width=24, height=10)
  vplayout(6, 3)

  all.tbls <- data.frame()
  for (i in 1:length(plots)) {
    cur.stuff <- plots[[i]]
    cur.tree <- cur.stuff[[1]]
    cur.tbl <- cur.stuff[[2]]

    cur.tbl <- subset(cur.tbl, analysis == 'SLR Sitewise')
    all.tbls <- rbind(all.tbls, cur.tbl)

    # Plot tree to the left.
    sim <- PhyloSim()
    setPhylo(sim, cur.tree)
    p <- f(sim)
    print(p, vp=subplot(1, i))
  }

  # Write all-table to a file.
  write.csv(all.tbls, file="tabl_multi.csv", row.names=FALSE)

  # Length / indel sweep of TPR at FPR.
  z_val <- 'tpr_at_fpr'
  all.tbls[,'z_val'] <- all.tbls[, z_val]
  rect.df <- expand.grid(unique(all.tbls[,'aligner']), unique(all.tbls[,'tree']))
  colnames(rect.df) <- c('aligner', 'tree')
  rect.df[,'xmin'] <- 0.9
  rect.df[,'xmax'] <- 1.1
  rect.df[,'ymin'] <- 0.045
  rect.df[,'ymax'] <- 0.055
  rect.df[,'ins_rate'] <- 0
  rect.df[,'length'] <- 0
  rect.df[,'z_val'] <- 0
  p <- ggplot(all.tbls, aes(x=length, y=ins_rate, z=z_val))
  p <- p + theme_bw()
  x.brks <- c(0.2, 1.0, 2.0)
  y.brks <- c(0, 0.05, 0.10)
  p <- p + scale_x_continuous('Mean Path Length', breaks=x.brks)
  p <- p + scale_y_continuous('Indel Rate', breaks=y.brks)
  p <- p + coord_cartesian(xlim=c(0.1, 2.1), ylim=c(-0.005, 0.105))
  p <- p + geom_tile(aes(fill=z_val))
  p <- p + geom_rect(data=rect.df, fill=rgb(0,0,0,alpha=0), aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, colour=aligner))
  p <- p + scale_colour_brewer(name="Alignment Method", palette="Set1")
  n <- 5
  clr.lim <- c(0, 1)
  breaks <- seq(from=clr.lim[1], to=clr.lim[2], length.out=n+1)
  p <- p + scale_fill_gradientn(z_val, colours=brewer.pal(n=n,'Spectral'), limits=clr.lim, breaks=breaks)
  p <- p + facet_grid(tree ~ aligner)
  p <- p + opts(
    panel.grid.major = theme_blank(),
    panel.grid.minor = theme_blank(),
    strip.text.x = theme_blank(),
    strip.text.y = theme_blank(),
    strip.background = theme_blank(),
    axis.title.x = theme_blank(),
    axis.title.y = theme_blank()
  )
  print(p, vp=subplot(2:5, 1:3))

  # ROC plot at default length / indel.
  a <- '~/scratch/gj1_fig_one_a/current'
  b <- '~/scratch/gj1_fig_one_b/current'
  c <- '~/scratch/gj1_fig_one_c/current'
  sets <- c(a,b,c)
  all.sites <- data.frame()
  for (i in 1:length(sets)) {
    cur.dir <- sets[i]
    print(cur.dir)
    all.sites.f <- paste(cur.dir, '/', 'sites.Rdata', sep='')
    default.sites.f <- paste(cur.dir, '/', 'sites_at_default.Rdata', sep='')
    if (!file.exists(default.sites.f)) {
      all.sites <- load(all.sites.f)
      sub.sites <- subset(sites, tree_length == 1 & ins_rate == 0.05)
      save(sub.sites, file=default.sites.f)
    }
    load(default.sites.f)
    all.sites <- rbind(all.sites, sub.sites)
  }
  all.sites[,'tree'] <- factor(as.character(all.sites[,'tree']))
  f <- function(df, thresh) {return(slr.roc(df, na.rm=T))}
  roc <- summarize.by.labels(all.sites, f)
  p <- plot.roc(roc, plot=F, plot.x='fpr', plot.y='tpr', color.by='aligner')
  p <- p + theme_bw()
  p <- p + scale_colour_brewer(name="Alignment Method", palette="Set1")
  p <- p + geom_vline(xintercept=0.05,colour='black')
  p <- p + geom_vline(xintercept=0.01,colour='gray')
  max.x <- max(roc[, 'fpr']) * 0.1
  p <- p + scale_x_continuous(name="False Positive Rate", limits=c(0,max.x))
  p <- p + scale_y_continuous(name="True Positive Rate", limits=c(0,1))
  p <- p + facet_grid(tree ~ .)
  p <- p + opts(
    legend.position='none',
    panel.grid.major = theme_blank(),
    panel.grid.minor = theme_blank(),
    strip.text.x = theme_blank(),
    strip.text.y = theme_blank(),
    strip.background = theme_blank(),
    axis.title.x = theme_blank(),
    axis.title.y = theme_blank()
  )

  print(p, vp=subplot(6, 1:3))

  dev.off()

}



#
# Figure 4 - FPR vs TPR ROC plot for filters in three trees (A/B/C)
#
figure.four <- function() {
  .figure.four('figure_four.pdf', '~/scratch/gj1_fig_two_b/current',
    'TPR vs FPR for 17-Taxon Tree, MAFFT Alignments'
  )
}

figure.six <- function() {
  .figure.four('figure_six.pdf','~/scratch/gj1_fig_two_e/current', 
    'TPR vs FPR for 17-Taxon Tree, PRANK(codon) Alignments'
  )
}

filter.rocs <- function() {
  roc.a <- load.roc('~/scratch/gj1_fig_two_a/current')
  roc.b <- load.roc('~/scratch/gj1_fig_two_b/current')
  roc.c <- load.roc('~/scratch/gj1_fig_two_c/current')
  roc.d <- load.roc('~/scratch/gj1_fig_two_d/current')  

}

figure.four.new <- function() {

  roc <- load.roc('gj1_fig_two_a', return.merged=F)
  roc <- rbind(roc, load.roc('gj1_fig_two_b', return.merged=F))
  roc <- rbind(roc, load.roc('gj1_fig_two_c', return.merged=F))
  roc <- rbind(roc, load.roc('gj1_fig_two_clustalw', return.merged=F))

  keep.filters <- c('none', 'branchlength_hi', 'guidance_hi', 'tcoffee_hi', 'optimal_hi', 'true')
  roc <- subset(roc, filter %in% keep.filters)

  roc <- ins.mpl.factors(roc)
  roc <- aln.factors(roc)
  roc <- filter.factors(roc)
  roc <- roc[ order(roc$filter, roc$aligner), ]
  roc$analysis <- paste(roc$aligner, roc$filter, sep=" / ")
  roc$analysis <- factor(roc$analysis, levels=unique(roc$analysis))
  roc$simulation <- paste(roc$tree_length, roc$ins_rate, sep="\n")

  roc <- subset(roc, fpr <= 0.1)

  p <- plot.roc(roc,plot=F,plot.x='fpr',plot.y='tpr', color.by='filter')
  p <- p + theme_bw()
  p <- p + geom_vline(x=0.05, colour='gray', linetype='dashed')
  p <- p + geom_hline(y=0.5, colour='gray', linetype='dashed')
  p <- p + scale_colour_brewer(name="Filtering Method", palette="Set1")
  p <- p + facet_grid(aligner ~ simulation)
  p <- p + scale_y_continuous(name="True Positive Rate", limits=c(0,1.0), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))
  p <- p + scale_x_continuous(name="False Positive Rate", limits=c(0,0.1), breaks=c(0, 0.02, 0.04, 0.06, 0.08, 0.1))
  p <- p + opts(
    panel.margin=unit(rep(0.5, times=4), "lines")
  )

  pdf(file="figure_four_new.pdf",width=16,height=12)
    print(p)
  dev.off()
}

.figure.four <- function(file, dir, title) {
  sites.file <- paste(dir, '/', 'sites.Rdata', sep='')

  load.roc <- function(d) {
    figure.roc.file <- paste(d, '/', 'figure_roc.Rdata', sep='')
    if (!file.exists(figure.roc.file)) {
      load(sites.file)
    
      f <- function(df, thresh) {return(slr.roc(df, na.rm=T))}
      roc <- summarize.by.labels(sites, f)
      save(roc, file=figure.roc.file)
    }
    load(figure.roc.file)
    return(roc)
  }
  
  roc <- load.roc(dir)

  roc <- subset(roc, filter != 'columns')
  roc <- filter.factors(roc)
  roc <- ins.mpl.factors(roc)

  p <- plot.roc(roc,plot=F,plot.x='fpr',plot.y='tpr', color.by='filter')
  p <- p + geom_vline(x=0.01, colour='gray', linetype='dashed')
  p <- p + geom_vline(x=0.05, colour='black', linetype='dashed')
  p <- p + scale_colour_brewer(name="Filtering Method", palette="Set1")
  p <- p + facet_grid(ins_rate ~ tree_length)
  p <- p + opts(
    panel.margin=unit(rep(0.5, times=4), "lines"),
    title=title
  )
  p <- p + scale_y_continuous(name="True Positive Rate", limits=c(0,1.0), breaks=c(0,0.5,1))
  p <- p + scale_x_continuous(name="False Positive Rate", limits=c(0,0.1), breaks=c(0, 0.05, 0.1))

  pdf(file=file,width=10,height=10)
  print(p)
  dev.off()
}

supp.fdr.figures <- function() {
#  .supp.fdr.figure('~/scratch/gj1_fig_two_a/current', 'supp_fdr_a.pdf', '6-taxon Tree')
  .supp.fdr.figure('~/scratch/gj1_fig_two_b/current', 'supp_fdr_b.pdf', '17-taxon Tree')
#  .supp.fdr.figure('~/scratch/gj1_fig_two_c/current', 'supp_fdr_c.pdf', '44-taxon Tree')
}

.supp.fdr.figure <- function(dir, file, tree_label) {
  print(dir)

  sites.file <- paste(dir, '/', 'sites.Rdata', sep='')
  figure.roc.file <- paste(dir, '/', 'figure_roc.Rdata', sep='')
  if (!file.exists(figure.roc.file)) {
    load(sites.file)

    f <- function(df, thresh) {return(slr.roc(df, na.rm=T))}
    roc <- summarize.by.labels(sites, f)
    save(roc, file=figure.roc.file)
  }
  load(figure.roc.file)

  roc <- filter.factors(roc)

  print("  plotting")
  p <- plot.roc(roc,plot=F,plot.x='p',plot.y='fp', color.by='filter', alpha=0.8)
  fdr <- 0.5
  p <- p + geom_abline(intercept=0, slope=fdr, colour='gray', linetype='dashed')
  fdr <- 0.1
  p <- p + geom_abline(intercept=0, slope=fdr, colour='black', linetype='dashed')
  p <- p + scale_colour_brewer(name="Filtering Method", palette="Set1")
  p <- p + facet_grid(ins_rate ~ tree_length)
  p <- p + opts(
    panel.margin=unit(rep(0.5, times=4), "lines"),
    title=paste("False Positives vs Positives for ",tree_label,sep="")
  )
  max.x <- max(roc[,'p']) * 0.02
  sub.below.max <- subset(roc, p < max.x)
  max.y <- max(sub.below.max[,'fp']) / 3
  p <- p + scale_y_continuous(name="Number of False Positives", limits=c(0,max.y))
  p <- p + scale_x_continuous(name="Number of Positives Identified", limits=c(0,max.x))

  pdf(file=file,width=10,height=10)
  print(p)
  dev.off()
}

#
# Figure 5 - TPR_FDR<0.1 sensitivity for filters in three trees (A/B/C)
#

figure.five <- function() {
  b <- '~/scratch/gj1_fig_two_b/current'
  c <- '~/scratch/gj1_fig_two_e/current'

  pdf(file="figure_five.pdf", width=8, height=6)
  p <- .figure.five('tp_at_fdr', 'TP at FDR<0.1 (Relative to True Alignment)', FALSE, c(b,c), show.filtered.fraction=T)
  print(p)
  dev.off()
}

supp.figure.two.through.four <- function() {
  b <- '~/scratch/gj1_fig_two_b/current'
  e <- '~/scratch/gj1_fig_two_e/current'

  pdf(file="supp_figure_two.pdf", width=10, height=10)
  p <- .figure.five('fdr_at_thresh', 'FDR at Threshold', FALSE, c(b,e), show.filtered.fraction=TRUE)
  print(p)
  dev.off()

  pdf(file="supp_figure_three.pdf", width=10, height=10)
  p <- .figure.five('tp_at_thresh', 'TP at Threshold', FALSE, c(b,e), show.filtered.fraction=TRUE)
  print(p)
  dev.off()

  pdf(file="supp_figure_four.pdf", width=10, height=10)
  p <- .figure.five('tp_at_thresh', 'FP at Threshold', FALSE, c(b,e), show.filtered.fraction=TRUE)
  print(p)
  dev.off()

}

.figure.five <- function(field, lbl, relative.to.true, dirs, show.filtered.fraction=FALSE) {

  tbl <- data.frame()
  for (i in 1:length(dirs)) {
    cur.dir <- dirs[i]
    cur.tbl <- read.csv(paste(cur.dir, '/', 'table.csv', sep=''))
    tbl <- rbind(tbl, cur.tbl)
  }

  print(lbl)

  tbl <- subset(tbl, filter != 'columns')
  
  tbl[,'z_value'] <- tbl[, field]
  tbl[is.nan(tbl$z_value), 'z_value'] <- 1
  tbl[tbl$z_value < 0, 'z_value'] <- 0

  # If relative to the true value, go through and normalize.
  if (relative.to.true) {
    for (i in 1:nrow(tbl)) {
      cur.row <- tbl[i,]
      true.row <- subset(tbl,
        ins_rate == cur.row$ins_rate &
        length == cur.row$length &
        tree == cur.row$tree &
        aligner == cur.row$aligner &
        filter == 'True Alignment'
      )
      ratio <- cur.row$z_value / true.row$z_value
      tbl[i, 'z_value'] <- ratio
    }
  }

  tbl <- ins.mpl.factors(tbl)
  tbl <- aln.factors(tbl)
  tbl <- filter.factors(tbl)

  filter.lbl <- levels(tbl$filter)
  tbl[, 'filter_x'] <- as.numeric(tbl[, 'filter'])

  tbl <- tree.factors(tbl)
  tbl <- aln.factors(tbl)

  split.by <- 'aligner'
  values <- unique(tbl[, split.by])
  layout.width <- 0.8
  n <- length(values)
  for (i in 1:length(values)) {
    ind <- tbl[, split.by] == values[i]
    tbl[ind, 'filter_x'] <- tbl[ind, 'filter_x'] - (1/2)*layout.width/n + (i-1)*layout.width/n
  }

  flank <- 0.3 # flank is now half of the width of our bars, for flush layout.
  flank <- flank / 2 # this adds more space between.
  tbl[,'flank'] <- flank

  main.y.min <- 0
  filter.fraction.height <- 0.2
  if (field %in% c('tp_at_fpr', 'tp_at_fdr', 'tp_at_thresh')) {
    filter.fraction.height <- max(tbl[,'z_value']) * 0.2
  }  

  ymax <- max(tbl[, 'z_value'])
  y.lim <- c(0,ymax)
  y.breaks <- round(c(0, ymax/4, ymax/4*2, ymax/4*3, ymax), digits=2)

  y.labels <- as.character(c(0, 0.5, 1.0))
  if (show.filtered.fraction) {
    main.y.min <- filter.fraction.height + filter.fraction.height * .8
  }

  tbl[,'main.y.min'] <- main.y.min
  tbl[, 'filter.fraction.height'] <- filter.fraction.height

  p <- ggplot(tbl, aes(x=filter_x, y=z_value, fill=filter))
  p <- p + theme_bw()
  p <- p + geom_rect(aes(xmin=filter_x-flank, xmax=filter_x+flank, ymin=main.y.min, ymax=main.y.min+z_value))
  p <- p + scale_fill_brewer(name="Filtering Method", palette="Set1")

  if (show.filtered.fraction) {
    clr <- 'gray'
    p <- p + geom_rect(colour=clr, fill=clr, aes(xmin=filter_x-flank, xmax=filter_x+flank, ymin=0, ymax=filtered.fraction * filter.fraction.height))
  }  

  if (field == 'tp_at_fpr' | field == 'tp_at_fdr' | field == 'tp_at_thresh') {
    y.lim <- c(0, max(tbl[, 'z_value']))
    y.breaks <- c(0, as.integer(y.lim[2]/2), y.lim[2])
    y.labels <- as.character(y.breaks)
  }
  if (relative.to.true == TRUE) {
    y.lim <- c(0, max(tbl[, 'z_value']))
    y.breaks <- unique(c(0, 1, as.integer(y.lim[2]/2), as.integer(y.lim[2])))
    y.labels <- as.character(y.breaks)
  }

  if (show.filtered.fraction) {
    y.lim[2] <- y.lim[2] + main.y.min
    y.labels <- as.character(y.breaks)
    y.breaks <- y.breaks + main.y.min
    y.breaks <- c(0, filter.fraction.height, y.breaks)
    y.labels <- as.character(c(0, 1, y.labels))
  }

  print(y.breaks)
  print(y.labels)

  p <- p + scale_y_continuous(name=lbl, limits=y.lim, breaks=y.breaks, labels=y.labels)
  brk <- c(0, 1:(length(filter.lbl)+1))
  lbls <- c('', filter.lbl, '')
  p <- p + scale_x_continuous(name='Filtering Method', breaks=c(0, 1:(length(filter.lbl)+1)), labels=c('a',filter.lbl,'b'))
  p <- p + facet_grid(ins_rate ~ length)

  tree.string = "for 17-Taxon Tree"
  if (length(dirs) == 3) {
    tree.string = "for 6-, 17-, and 44-taxon Trees"
  }
  p <- p + opts(
    panel.margin=unit(rep(0.5, times=4), "lines"),
    title=paste(lbl, tree.string, sep=' '),
    axis.text.x = theme_text(angle=90, hjust=1),
    legend.position = 'none'
  )
  return(p)
}


false.pos.stats <- function() {

  plot.distributions <- function(roc, field, xlim=c(0,15), plot.diff=F) {
    print(paste("Plotting FP distributions for",field))
    df <- roc

    # Keep only sites that are FPs in the inferred alignment and TNs in the true alignment.
    # IMPORTANT DETAIL: I use the FPR<0.05 threshold from the TRUE alignment to define
    # our threshold for FPs and TPs...

    roc$aligner <- roc$aligner.aln
    roc$filter <- roc$filter.aln
    roc <- ins.mpl.factors(roc)
    roc <- aln.factors(roc)
    roc <- filter.factors(roc)
    roc <- roc[ order(roc$filter, roc$aligner), ]
    roc$analysis <- paste(roc$aligner, roc$filter, sep=" / ")
    roc$analysis <- factor(roc$analysis, levels=unique(roc$analysis))
    roc$simulation <- paste(roc$tree_length, roc$ins_rate, sep="\n")

    for (cur.set in unique(roc$analysis)) {
      print(paste("  ",cur.set, ":", nrow(subset(roc, analysis == cur.set))))
    }
    for (cur.set in unique(roc$simulation)) {
      print(paste("  ",cur.set, ":", nrow(subset(roc, simulation == cur.set))))
    }
 
    roc[, 'bl_mismatch.aln'] <- roc$bl_aligned.aln - roc$bl_match.aln
    roc[, 'bl_mismatch.true'] <- 0

    cur.roc <- roc
    cur.roc$stat_aln <- cur.roc[, paste(field,'.aln',sep='')]
    cur.roc$stat_true <- cur.roc[, paste(field,'.true',sep='')]

    if (plot.diff) {
      cur.roc$stat_aln <- cur.roc$stat_aln - cur.roc$stat_true
    }

    # Color by the aln - true score difference.
    cur.roc$score_diff <- cur.roc$score.aln - cur.roc$score.true;
    cur.roc$score_diff <- round(cur.roc$score_diff)
    cur.roc[cur.roc$score_diff > 4, 'score_diff'] <- 4
    cur.roc$score_diff <- factor(cur.roc$score_diff, levels=rev(c(0:4)))
    cur.roc <- cur.roc[order(cur.roc$score_diff),]

    bin.width <- diff(xlim) / 20
    ramp <- colorRampPalette(c('blue','red'))

    p <- ggplot(cur.roc, aes(x=stat_aln))
    p <- p + theme_bw()
    p <- p + geom_bar(binwidth=bin.width)
    #p <- p + scale_fill_manual(values=rev(ramp(5)))
    #p <- p + geom_vline(xintercept=0)
    x.lbl <- paste(field, "(aln)", sep=' ')
    if (plot.diff) {
      x.lbl <- paste(field, "(aln - true)", sep=' ')
    }
    p <- p + scale_x_continuous(x.lbl)
    p <- p + facet_grid( analysis ~ simulation)
    p <- p + opts(strip.text.y = theme_text())
    p <- p + opts(
     legend.position = 'none'
    )
    return(p)
  }

  combine.fp.rocs.from.dirs <- function(dbs) {
    all.stats <- data.frame()
    for (db in dbs) {
      stats <- load.roc(db, return.fps=T)
      all.stats <- rbind(all.stats, stats)  
    }
    return(all.stats)
  }

  db.a <- 'gj1_fig_two_a'
  db.b <- 'gj1_fig_two_b'
  db.c <- 'gj1_fig_two_c'
  db.clustalw <- 'gj1_fig_two_clustalw'

  comb.stats <- combine.fp.rocs.from.dirs(c(db.a, db.b, db.c, db.clustalw))
  print(str(comb.stats))

  comb.stats$filter <- comb.stats$filter.aln
  comb.stats$aligner <- comb.stats$aligner.aln
  
  comb.stats <- subset(comb.stats, 
    filter == 'none'
  )

  pdf("fp.stats.pdf", width=10, height=14)
    vplayout(1, 3)
    p <- plot.distributions(comb.stats, 'bl_aligned')
    print(p, vp=subplot(1,1))
    p <- plot.distributions(comb.stats, 'bl_aligned', plot.diff=T)
    print(p, vp=subplot(1,2))
    p <- plot.distributions(comb.stats, 'bl_mismatch')
    print(p, vp=subplot(1,3))
  dev.off()

}


false.positives <- function(fp_index=1, file="test.pdf") {
  dbs <- c('gj1_fig_two_a', 'gj1_fig_two_b', 
           'gj1_fig_two_c',  'gj1_fig_two_clustalw')

  losers.f <- 'losers.Rdata'
  if (!file.exists(losers.f)) {
    print("  loading genes")
    genes <- load.all.genes(dbs)
    print("  loading rocs")

    keep.fields <- c('truth.aln', 'score.aln', 'score.true',
      'tree_length', 'ins_rate',
      'slrsim_rep', 'seq_position',
      'aligner.aln', 'filter.aln',
      'db'
    )
    all.sites <- load.all.rocs(dbs, keep.fields=keep.fields)
    print(nrow(all.sites))

    fp.db <- 'gj1_fig_two_c'
    ins <- 0.05
    mpl <- 1

    all.sites <- subset(all.sites, tree_length == mpl & ins_rate == ins)
    print(nrow(all.sites))

    cur.sub <- subset(all.sites, 
      filter.aln == 'none' & db == fp.db
    )
    print(nrow(cur.sub))
    losers <- subset(cur.sub, truth.aln == 0 & score.aln > 2.71 & score.true < 2.71)
    print(nrow(losers))

    loser.others <- data.frame()
    for (i in 1:nrow(losers)) {
      cur.loser <- losers[i,]
      print(paste(cur.loser$slrsim_rep, cur.loser$seq_position))

      cur.others <- subset(all.sites, 
        slrsim_rep == cur.loser$slrsim_rep & 
        seq_position == cur.loser$seq_position
      )
      print(paste(i, nrow(cur.others)))
      print(unique(cur.others$filter.aln))
      loser.others <- rbind(loser.others, cur.others)
    }

    genes<- subset(genes, tree_length == mpl & ins_rate == ins)
    save(losers, loser.others, genes, ins, mpl, file=losers.f)
  }
  load(losers.f)

  print(nrow(losers))
  cur.fp <- losers[fp_index,]

  load.window <- function(x, seq_id, seq_position, flank.width) {
    db <- x$db
    data.dir <- paste('~/scratch/',db,'/current/data', sep='')
    tree.f <- paste(data.dir, '/', x$tree_file, sep='')
    aln.f <- paste(data.dir, '/', x$masked_aln_file, sep='')
    print(aln.f)
    aln <- read.aln(aln.f)
    pep.aln <- aln.tx(aln)

    tree <- NULL
    if (file.exists(tree.f)) {
      print(tree.f)
      tree <- read.tree(tree.f)
    }

    seq_id = 'human'
    aln.col <- aln.col.position(pep.aln, seq_id, seq_position)

    col.lo <- max(1, aln.col - flank.width)
    col.hi <- min(aln.length(pep.aln), aln.col + flank.width)

    below <- aln.slice(pep.aln, col.lo, aln.col-1)
    at <- aln.slice(pep.aln, aln.col, aln.col)
    above <- aln.slice(pep.aln, aln.col+1, col.hi)
        
    spliced <- aln.splice(below, at, empty.space=0)
    spliced <- aln.splice(spliced, above, empty.space=0)
    return(list(aln=spliced, tree=tree))
  }

  window.width <- 19

  filters.alns <- list(
    c('none', 'mafft'),
    c('none', 'clustalw'),
    c('none', 'prank'),
    c('none', 'prank_codon'),
    c('branchlength_hi', 'prank_codon'),
    c('tcoffee_hi', 'prank_codon'),
    c('true', 'prank_codon')
  )

  pdf.width <- 1.5 * (window.width*2+3)
  pdf.height <- 17 * length(filters.alns)
  pdf(file=file, width=pdf.width/5, height=pdf.height/5)
  vplayout(1, length(filters.alns))

  for (i in 1:length(filters.alns)) {
    filt.aln <- filters.alns[[i]]
    cur.filter <- filt.aln[1]
    cur.aln <- filt.aln[2]
    print(paste(cur.filter, cur.aln))

    cur.site <- cur.fp
    cur_seq_position <- cur.site$seq_position
    cur_slrsim_rep <- cur.site$slrsim_rep

    cur.gene <- subset(genes,
      slrsim_rep == cur_slrsim_rep & 
      filter == cur.filter &
      aligner == cur.aln
    )

    cur.site <- subset(loser.others,
      slrsim_rep == cur_slrsim_rep & 
      filter.aln == cur.filter &
      aligner.aln == cur.aln &
      seq_position == cur_seq_position
    )

    if (nrow(cur.gene) != 1) {
      if (nrow(cur.gene) == 0) {
         stop("No gene found!")
      }
      stop("Need just one gene!")
    }

    if (cur.filter == 'true') {
      cur.score <- cur.fp$score.true
    } else if (nrow(cur.site) != 1 && cur.filter != 'none') {
      cur.score <- -1
    } else if (nrow(cur.site) == 1) {
      cur.score <- cur.site$score.aln
    } else {
      if (nrow(cur.site) == 0) {
        stop("No site found!")
      } else {
        stop("More than one site found!")
      }
    }

    lbl <- sprintf("Aligner=%s / Filter=%s / score=%.3f", cur.aln, cur.filter, cur.score)

    alntree <- load.window(cur.gene, 'human', cur_seq_position, window.width)
    aln <- alntree$aln
    tree <- alntree$tree

    n.seqs <- aln.num.seqs(aln)
    aln.modifier <- function(x) {
        x = x + geom_rect(xmin=window.width+.5, xmax=window.width+1.5, ymin=0.5, ymax=n.seqs+0.5, fill=NA, colour='black', border=1, alpha=0.5)
        return(x)
    }

    p <- aln.plot(aln,
      aln.plot.xlabel='', 
      tree=tree,
      aln.mod=aln.modifier,
      tree.labels.in.tree=F,
      aln.to.tree.size = 5,
      plot.title = lbl
    )
    pushViewport(subplot(1,i))
    grid.draw(p$grob)
    upViewport()
  }
  dev.off()
  
}

lots.of.fps <- function() {
  for (i in 1:20) {
    false.positives(i, file=paste("fps/",i,".pdf",sep=''))
  }
}

merged.datasets <- function(tbl) {

  merged.data.f <- 'merged_data.csv'
  if (!file.exists(merged.data.f)) {
    print(paste("Creating merged datasets..."))

    load.s <- function(db) {return(subset(load.sites(db), filter=='none'))} 
    sites <- load.s('gj1_fig_two_c')
    mafft <- load.s('gj1_fig_two_a')
    prank <- load.s('gj1_fig_two_b')
    probcons <- load.s('gj1_fig_two_d')
    pagan <- load.s('gj1_fig_two_e')
    clustalw <- load.s('gj1_fig_two_clustalw')
    tc_lo <- subset(load.sites('gj1_fig_two_c'), filter=='tcoffee_lo')
    bl_lo <- subset(load.sites('gj1_fig_two_c'), filter=='branchlength_lo')
    tc_hi <- subset(load.sites('gj1_fig_two_c'), filter=='tcoffee_hi')
    bl_hi <- subset(load.sites('gj1_fig_two_c'), filter=='branchlength_hi')

    true <- subset(load.sites('gj1_fig_two_c'), filter=='true')
    true$aligner <- 'True Alignment'
    aln <- subset(load.sites('gj1_fig_two_c'), filter=='none')

    print(nrow(sites))
    print(nrow(tc_lo))
    print(nrow(bl_lo))

    print("  Combining pvals...")
  
    add.rows <- function(input, a, b, lbl) {
      return(rbind(input, combine.pvals(a, b, lbl)))
    }

    rows <- data.frame()
    
    rows <- rbind(rows, true)
    rows <- rbind(rows, aln)    

    rows <- add.rows(rows, sites, clustalw, 'prank(codon)+clustalw')
    rows <- add.rows(rows, sites, mafft,    'prank(codon)+mafft')
    rows <- add.rows(rows, sites, prank,    'prank(codon)+prank')
    rows <- add.rows(rows, sites, pagan,    'prank(codon)+pagan')
    rows <- add.rows(rows, sites, probcons, 'prank(codon)+probcons')

    codon.pagan <- combine.pvals(sites, pagan, 'prank(codon)+pagan')
    codon.pagan.prank <- combine.pvals(codon.pagan, prank, 'prank(codon)+pagan+prank')
    codon.pagan.tc_lo <- combine.pvals(codon.pagan, tc_lo, 'prank(codon)+pagan+tc(low)')
    codon.pagan.prank.tc_lo <- combine.pvals(codon.pagan.prank, tc_lo, 'prank(codon)+pagan+prank+tc(low)')

    rows <- add.rows(rows, codon.pagan, prank,             'prank(codon)+pagan+prank')
    rows <- add.rows(rows, codon.pagan, tc_lo,             'prank(codon)+pagan+tc(low)')
    rows <- add.rows(rows, codon.pagan, bl_lo,             'prank(codon)+pagan+bl(low)')
    rows <- add.rows(rows, codon.pagan.tc_lo, bl_lo,       'prank(codon)+pagan+tc(low)+bl(low)')
    rows <- add.rows(rows, codon.pagan.prank, tc_lo,       'prank(codon)+pagan+prank+tc(low)')
    rows <- add.rows(rows, codon.pagan.prank, bl_lo,       'prank(codon)+pagan+prank+bl(low)')
    rows <- add.rows(rows, codon.pagan.prank.tc_lo, bl_lo, 'prank(codon)+pagan+prank+tc(low)+bl(low)')

    codon.prank <- combine.pvals(sites, prank, 'prank(codon)+prank')
    codon.prank.tc_lo <- combine.pvals(codon.prank, tc_lo, 'prank(codon)+prank+tc(low)')

    rows <- add.rows(rows, codon.prank, mafft,       'prank(codon)+prank+mafft')
    rows <- add.rows(rows, codon.prank, tc_lo,       'prank(codon)+prank+tc(low)')
    rows <- add.rows(rows, codon.prank, bl_lo,       'prank(codon)+prank+bl(low)')
    rows <- add.rows(rows, codon.prank.tc_lo, bl_lo, 'prank(codon)+prank+tc(low)+bl(low)')

    rows <- add.rows(rows, sites, bl_lo, 'prank(codon)+bl(low)')
    rows <- add.rows(rows, sites, bl_hi, 'prank(codon)+bl(high)')
    rows <- add.rows(rows, sites, tc_lo, 'prank(codon)+tc(low)')
    rows <- add.rows(rows, sites, tc_hi, 'prank(codon)+tc(high)')
    rows <- add.rows(rows, tc_lo, bl_lo, 'prank(codon)+tc(low)+bl(low)')
    rows <- add.rows(rows, tc_lo, bl_hi, 'prank(codon)+tc(low)+bl(high)')
    rows <- add.rows(rows, tc_hi, bl_lo, 'prank(codon)+tc(high)+bl(low)')
    rows <- add.rows(rows, tc_hi, bl_hi, 'prank(codon)+tc(high)+bl(high)')

    print("  Gathering data from combined sites...")
    pt <- summarize.by.labels(rows, paper.table)
    write.csv(pt, file=merged.data.f, row.names=F)
  }
  tbl <- read.csv(merged.data.f)
  tbl$filter <- ''

#  tbl <- aln.factors(tbl)
#  tbl <- filter.factors(tbl)
  #tbl <- ins.mpl.factors(tbl)
  output.tables.long(tbl, "merged")

  return(tbl)
}

filter.tables <- function() {
  tbl.a <- load.tbl('gj1_fig_two_a')
  tbl.b <- load.tbl('gj1_fig_two_b')
  tbl.c <- load.tbl('gj1_fig_two_c')
  tbl.clustalw <- load.tbl('gj1_fig_two_clustalw')
  tbl <- rbind(tbl.a, tbl.b, tbl.c, tbl.clustalw)

  tbl <- aln.factors(tbl)
  tbl <- filter.factors(tbl)
  #tbl <- ins.mpl.factors(tbl)

  output.tables.long(tbl, "filters")
}

output.tables.long <- function(tbl, file) {
  tbl$tree_length <- tbl$length
  tbl <- tbl[with(tbl, order(tree_length, ins_rate, aligner, filter)),]
  tbl$conditions <- paste(tbl$tree_length, tbl$ins_rate, sep=' / ')

  id.fields <- c(
    'aligner', 'filter'
  )
  variable.fields <- c(
    'n_sites', 'tpr_at_fpr', 'tp_at_fdr', 'fdr_at_thresh', 'tp_at_thresh'
  )

  df <- unique(tbl[, id.fields])
  tables <- list()
  for (cond in unique(tbl$conditions)) {
    print(cond)
    cur.tbl <- subset(tbl, conditions == cond)
    cur.tbl <- cur.tbl[, c(id.fields, variable.fields)]
    names(cur.tbl) <- c(
      names(cur.tbl[, id.fields]),
      paste(cond, "\n", names(cur.tbl[,variable.fields]), sep=' ')
    )
    df <- merge(df, cur.tbl, by=id.fields)
  }

  df <- df[with(df, order(aligner, filter)),]
  print.tables(list(df), file=paste(file,'.html',sep=''))

  write.csv(df, file=paste(file,'.csv', sep=''), row.names=F)
  
}

print.tables <- function(tables, file) {
  all.strings <- ''
  for (i in 1:length(tables)) {
    cur.tbl <- tables[[i]]
    print(str(cur.tbl))
    x.tbl <- xtable(cur.tbl, digits=3)
    tbl.string <- print(x.tbl, type='html')
    all.strings <- paste(all.strings, tbl.string)
  }
  cat(all.strings, file=file)
}


combined.tests <- function() {

  load.sites <- function(db) {
    dir = paste('~/scratch/',db,'/current', sep='')
    sites.file <- paste(dir, '/', 'sites.Rdata', sep='')
    load(sites.file)
    sites <- subset(sites, ins_rate == 0.025 & tree_length == 2 & tree == 'bglobin.nh')
    return(sites)
  }

  sites.b <- load.sites('gj1_fig_two_b')
  sites.e <- load.sites('gj1_fig_two_e')
  sites.h <- load.sites('gj1_fig_two_h')

  sites.clustalw <- subset(sites.h, filter == 'No filter')
  sites.mafft <- subset(sites.b, filter == 'No filter')
  sites.prank <- subset(sites.e, filter == 'No filter')
  sites.tcoffee <- subset(sites.e, filter == 'tcoffee')
  sites.bl <- subset(sites.e, filter == 'branchlength')
  sites.true <- subset(sites.e, filter == 'True Alignment')

  print.vals <- c('slrsim_label', 'auc', 'tpr_at_fpr', 'tp_at_fdr', 'fdr_at_thresh')

  merged <- combine.pvals(sites.clustalw, sites.mafft)
  pt <- paper.table(merged)
  pt$slrsim_label <- 'clustalw_mafft'
  print(pt[, print.vals])

  merged <- combine.pvals(sites.clustalw, sites.prank)
  pt <- paper.table(merged)
  pt$slrsim_label <- 'clustalw_prank'
  print(pt[, print.vals])

  merged <- combine.pvals(sites.mafft, sites.prank)
  pt <- paper.table(merged)
  pt$slrsim_label <- 'mafft_prank'
  print(pt[, print.vals])

  merged <- combine.pvals(sites.tcoffee, sites.bl) 
  pt <- paper.table(merged)
  pt$slrsim_label <- 'tcoffee_bl'
  print(pt[, print.vals])

  pt <- paper.table(sites.prank)
  pt$slrsim_label <- 'prank'
  print(pt[, print.vals])

  pt <- paper.table(sites.true)
  pt$slrsim_label <- 'true'
  print(pt[, print.vals])

}



aln.factors <- function(data) {
  aligners <- unique(data$aligner)
  filters <- unique(data$filter)

  ordered.f <- c('clustalw', 'mafft', 'probcons', 'fsa', 'fsa_careful', 'prank', 'pagan', 'prank_codon', 'true')
  aligners <- unique(c(ordered.f[ordered.f %in% aligners], as.character(aligners)))

  levels <- c()
  labels <- c()
  for (aligner in aligners) {
    levels <- c(levels, aligner)
    labels <- c(labels, 
      switch(aligner,
        clustalw = 'ClustalW',
        mafft = 'MAFFT',
        fsa = 'FSA',
        fsa_careful = 'FSA (high specificity)',
        probcons = 'Probcons',
        prank = 'PRANK',
        prank_codon = 'PRANK(codon)',
        pagan = 'PAGAN',
        pagan_groups = 'PAGAN(groups)',
        pagan_codon = 'PAGAN(codon)',
        True_Alignment = 'True Alignment',
        'True Alignment' = 'True Alignment',
        true = 'True Alignment',
        none = 'True Alignment',
        aligner
      )
    )
  }
  data$aligner <- factor(data$aligner, levels=levels, labels=labels)
  return(data)
}

ins.mpl.factors <- function(data, short=F) {
  if (is.null(data$tree_length)) {
    data$tree_length <- data$length
  }
  # First, take the insertion rate and multiply by two to get the indel rate
  data$ins_rate <- data$ins_rate * 2
  if (short) {
    data$ins_rate <- factor(paste("Indel =",data$ins_rate))
    data$tree_length <- paste("MPL = ",data$tree_length)
  } else {
    data$ins_rate <- factor(paste("Indel Rate = ",data$ins_rate))
    data$tree_length <- paste("Divergence = ",data$tree_length)
  }
  data$ins_rate <- factor(data$ins_rate, levels=rev(levels(data$ins_rate)))
  return(data)
}

tree.factors <- function(data) {
  trees <- unique(data$tree)

  ordered.f <- c('artificial.nh', 'bglobin.nh', 'encode.nh')
  trees <- ordered.f[ordered.f %in% trees]
  
  levels <- c()
  labels <- c()
  for (f in trees) {
    levels <- c(levels, f)
    labels <- c(labels, 
      switch(f,
        artificial.nh = '6-Taxon',
        bglobin.nh = '17-Taxon',
        encode.nh = '44-Taxon'
      )
    )
  }
  data$tree <- factor(data$tree, levels=levels, labels=labels)
  return(data)
} 

filter.factors <- function(data) {
  filters <- unique(data$filter)

  ordered.f <- c('No filter', 'none',
    'gblocks', 
    'guidance', 'guidance_lo', 'guidance_hi',
    'columns', 
    'branchlength', 'branchlength_lo', 'branchlength_hi',
    'tcoffee', 'tcoffee_lo', 'tcoffee_hi',
    'optimal', 'optimal_lo', 'optimal_hi',
    'True Alignment', 'True_Alignment', 'true')
  filters <- unique(c(ordered.f[ordered.f %in% filters], filters))
  
  levels <- c()
  labels <- c()
  for (f in filters) {
    levels <- c(levels, f)
    labels <- c(labels, 
      switch(f,
        gblocks = 'Gblocks',
        optimal = 'Optimal',
        'optimal_lo' = 'Optimal(low)',
        'optimal_hi' = 'Optimal(high)',
        tcoffee = 'T-Coffee',
        'tcoffee_lo' = 'T-Coffee(low)',
        'tcoffee_hi' = 'T-Coffee(high)',
        guidance = 'GUIDANCE',
        'guidance_lo' = 'GUIDANCE(low)',
        'guidance_hi' = 'GUIDANCE(high)',
        columns = 'Non-gap Sequences',
        branchlength = 'Branch Length',
        'branchlength_lo' = 'Branch Length(low)',
        'branchlength_hi' = 'Branch Length(high)',
        'No filter' = 'None',
        'none' = 'None',
        'True Alignment' = 'True Alignment',
        True_Alignment = 'True Alignment',
        true = 'True Alignment'
      )
    )
  }
  data$filter <- factor(data$filter, levels=levels, labels=labels)
  return(data)
}

load.tbl <- function(db) {
  dir = paste('~/scratch/',db,'/current', sep='')
  tbl.file <- paste(dir, '/', 'table.csv', sep='')
  tbl <- read.csv(file=tbl.file, stringsAsFactors=F)
  return(tbl)
}

load.sites <- function(db) {
  dir = paste('~/scratch/',db,'/current', sep='')
  sites.file <- paste(dir, '/', 'sites.Rdata', sep='')
  load(sites.file)
  return(sites)
}

load.genes <- function(db) {
  dir = paste('~/scratch/',db,'/current', sep='')
  genes.file <- paste(dir, '/', 'genes.Rdata', sep='')
  load(genes.file)
  return(genes)
}

load.all.genes <- function(dbs) {
  all.genes <- data.frame()
  for (i in 1:length(dbs)) {
    db <- dbs[i]
    dir = paste('~/scratch/',db,'/current', sep='')
    genes.file <- paste(dir, '/', 'genes.Rdata', sep='')
    load(genes.file)
    genes$db <- db
    all.genes <- rbind(all.genes, genes)
  }
  return(all.genes)
}

load.all.rocs <- function(dbs, return.merged=T, return.fps=F, keep.fields=NULL) {

  df <- data.frame()
  for (i in 1:length(dbs)) {
    cur.df <- load.roc(dbs[i], return.merged=return.merged, return.fps=return.fps)
    cur.df$db <- dbs[i]
    if (!is.null(keep.fields)) {
      cur.df <- subset(cur.df, select=keep.fields)
    }
    df <- rbind(df, cur.df)
  }
  return(df)
}

load.roc <- function(db, return.merged=T, return.fps=F) {
  dir <- paste('~/scratch/', db, '/current', sep='')
  sites.file <- paste(dir, '/sites.Rdata', sep='')
  rocs.file <- paste(dir, '/rocs.Rdata', sep='')
  merged.file <- paste(dir, '/rocs_merged.Rdata', sep='')
  fps.file <- paste(dir, '/rocs_fps.Rdata', sep='')

  if (!file.exists(rocs.file)) {
    print(paste("Making roc for",db,"..."))
    load(sites.file)
    roc <- summarize.by.labels(sites, slr.roc)
    save(roc, file=rocs.file)
  }
  print(paste("Loading roc for",db,"..."))

  if (!return.merged && !return.fps) {
    load(rocs.file)
    return(roc)
  }

  if (!file.exists(merged.file)) {
    print(paste("Making merged roc for",db,"..."))
    load(rocs.file)
    roc$alnfilter <- paste(roc$align, roc$filter, sep='')

    roc.aln <- subset(roc, filter != 'true')
    roc.true <- subset(roc, filter == 'true')
    
    merge.fields <- c('tree_length', 'ins_rate', 'slrsim_rep', 'seq_position')
  
    # For some reason we get duplicates for a given rep and seq_position. Remove them for both sets.
    remove.dups <- function(df) {
      return(df[!duplicated(df[, merge.fields]),])
    }
    roc.true <- remove.dups(roc.true)
  
    all.merged <- data.frame()

    filters <- unique(roc.aln$alnfilter)
    for (i in 1:length(filters)) {
      cur.f <- filters[i]
      sub.roc.aln <- subset(roc.aln, alnfilter == cur.f)
      sub.roc.aln <- remove.dups(sub.roc.aln)
      roc.merged <- merge(
        sub.roc.aln, 
        roc.true,
        by=merge.fields,
        suffixes=c('.aln', '.true')
      )
      print(paste("  ",cur.f, nrow(roc.merged)))
      all.merged <- rbind(all.merged, roc.merged)
    }
    roc <- all.merged
    rm(all.merged)
    save(roc, file=merged.file)
  }
  if (!return.fps) {
    load(merged.file)
    return(roc)
  }

  if (!file.exists(fps.file)) {
    print(paste("Making FPs roc for",db,"..."))
    load(merged.file)
    roc <- subset(roc, truth.aln == 0 & score.aln > thresh_at_fpr.true & score.true < thresh_at_fpr.true)
    save(roc, file=fps.file)
  }    
  load(fps.file)

  return(roc)
}
