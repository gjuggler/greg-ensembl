library(phylosim)
library(RColorBrewer)
library(xtable)
library(plyr)
if (Sys.getenv('USER') == 'gj1') {
  source("~/src/greg-ensembl/projects/phylosim/PhyloSimPlots.R")
  source("~/src/greg-ensembl/projects/slrsim/slrsim.functions.R")
  source("~/src/greg-ensembl/projects/slrsim/slrsim.plots.R")
} else {
  source("~/lib/greg-ensembl/projects/phylosim/PhyloSimPlots.R")
  source("~/lib/greg-ensembl/projects/slrsim/slrsim.functions.R")
  source("~/lib/greg-ensembl/projects/slrsim/slrsim.plots.R")
}


subplot <- function(x, y) viewport(layout.pos.col=x, layout.pos.row=y)
vplayout <- function(x, y) {
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(y,x)))
}

#
# Figure 1 - Power of three sitewise methods under zero indel rate.
#
multi.methods <- function() {
  tbl <- read.csv('~/scratch/gj1_fig_zero/current/table.csv')
  tbl[, 'tdr_at_thresh'] <- 1 - tbl[, 'fdr_at_thresh']

  comb.df <- data.frame()
  y.fields <- c('tpr_at_thresh', 'fpr_at_thresh', 'tpr_at_fpr2', 'cor_pearson')
  for (field in y.fields) {
    cur.df <- tbl
    cur.df$stat <- field
    cur.df$y.val <- cur.df[, field]
    if (field == 'fpr_at_thresh') {
      cur.df$y.val <- cur.df[, field] * 100
    }
    if (field == 'fdr_at_thresh') {
      cur.df$y.val <- cur.df[, field] * 3
    }
    
    comb.df <- rbind(comb.df, cur.df)
  }

  comb.df$stat <- factor(comb.df$stat, levels=y.fields)

  pdf(file="fig_zero.pdf", width=6, height=5)
  p <- ggplot(comb.df, aes(x=length, y=y.val, colour=analysis, linetype=analysis))
  p <- p + theme_bw()
  p  <- p + geom_line(alpha=1, size=0.5)
  clrs <- c('gray10', 'gray10', 'gray10')
  line.types <- c('dashed', 'dotted', 'solid')
  p <- p + scale_colour_manual(name="Analysis Method", values=clrs)
  p <- p + scale_linetype_manual(name="Analysis Method", values=line.types)
  p <- p + scale_y_continuous(limits=c(0, 1))
  p <- p + facet_grid(stat ~ tree)
  p <- p + coord_equal(ratio=2)
  p <- p + opts(
    axis.text.x = theme_text(size=5, angle=90, hjust=1.2),
    axis.text.y = theme_text(size=5, angle=0, hjust=1.2),
    axis.title.x = theme_blank(),
    axis.title.y = theme_blank(),
    strip.text.x = theme_blank(),
    strip.text.y = theme_blank(),
    strip.background = theme_blank(),
#    panel.grid.minor = theme_line(colour='gray60', size=0.2),
    panel.grid.major = theme_line(colour='gray80', size=0.2)
  )
  print(p)
  dev.off()
}

multi.aligners <- function() {
  tbl.a <- read.csv('~/scratch/gj1_fig_one_a/current/table.csv', stringsAsFactors=F)
  tbl.b <- read.csv('~/scratch/gj1_fig_one_b/current/table.csv', stringsAsFactors=F)
  tbl.c <- read.csv('~/scratch/gj1_fig_one_c/current/table.csv', stringsAsFactors=F)
  
  tbl.a$tree <- 'artificial'
  tbl.b$tree <- 'bglobin'
  tbl.c$tree <- 'encode'

  tbl.a <- subset(tbl.a, ins_rate == 0.05 | (ins_rate == 0 & aligner == 'none'))
  tbl.b <- subset(tbl.b, ins_rate == 0.05 | (ins_rate == 0 & aligner == 'none'))
  tbl.c <- subset(tbl.c, ins_rate == 0.05 | (ins_rate == 0 & aligner == 'none'))
  
  tbl <- rbind(tbl.a, tbl.b, tbl.c)
  tbl <- subset(tbl, aligner %in% c('clustalw', 'mafft', 'prank_codon', 'none'))
  tbl[tbl$aligner == 'none' & tbl$ins_rate == 0, 'aligner'] <- 'no_indels'
  tbl[, 'tdr_at_thresh'] <- 1 - tbl[, 'fdr_at_thresh']

  tbl$indel_rate <- tbl$ins_rate * 2
  tbl$indel_rate <- factor(tbl$indel_rate)

  tbl <- aln.factors(tbl)

  comb.df <- data.frame()
  y.fields <- c('tpr_at_thresh', 'fpr_at_thresh', 'tpr_at_fpr2', 'cor_pearson')
  for (field in y.fields) {
    cur.df <- tbl
    cur.df$stat <- field
    cur.df$y.val <- cur.df[, field]
    if (field == 'fpr_at_thresh') {
      cur.df$y.val <- cur.df[, field] * 100
    }
    comb.df <- rbind(comb.df, cur.df)
  }

  comb.df$stat <- factor(comb.df$stat, levels=y.fields)

  pdf(file="fig_aligners.pdf", width=6, height=5)
  p <- ggplot(comb.df, aes(x=length, y=y.val, linetype=aligner, colour=indel_rate, size=indel_rate))
  p <- p + theme_bw()
  p  <- p + geom_line(alpha=1)
  line.sizes <- c(1.5, 0.5)
  p <- p + scale_size_manual(values=line.sizes)
  line.clrs <- c('gray80', 'black', '#FF3333')
  p <- p + scale_colour_manual(name="Indel Rate", values=line.clrs)
  line.types <- c(1, 1, 2, 3, 4)
  p <- p + scale_linetype_manual(name="Aligner", values=line.types)
  p <- p + scale_y_continuous(limits=c(0, 1))
  p <- p + scale_x_continuous(limits=c(0, 2))
  p <- p + facet_grid(stat ~ tree)
  p <- p + coord_equal(ratio=2)
  p <- p + opts(
    axis.text.x = theme_text(size=5, angle=90, hjust=1.2),
    axis.text.y = theme_text(size=5, angle=0, hjust=1.2),
    axis.title.x = theme_blank(),
    axis.title.y = theme_blank(),
    strip.text.x = theme_blank(),
    strip.text.y = theme_blank(),
    strip.background = theme_blank(),
#    panel.grid.minor = theme_line(colour='gray60', size=0.2),
    panel.grid.major = theme_line(colour='gray80', size=0.2)
  )
  print(p)
  dev.off()
}

multi.filt <- function() {
  tbl.a <- read.csv('~/scratch/gj1_fig_three_a/current/table.csv')
  tbl.b <- read.csv('~/scratch/gj1_fig_three_b/current/table.csv')
  
  tbl.a$aligner <- 'prank_c'
  tbl.b$aligner <- 'clustalw'

  tbl.a <- subset(tbl.a, ins_rate %in% c(0, 0.04))
  tbl.b <- subset(tbl.b, ins_rate %in% c(0, 0.04))
  
  tbl <- rbind(tbl.a, tbl.b)
  tbl <- subset(tbl, !(filter %in% c('branchlength', 'optimal_a', 'optimal_b', 'meta')))
  tbl <- subset(tbl, filter %in% c('true', 'none', 'optimal_c', 'tcoffee'))
  tbl[, 'tdr_at_thresh'] <- 1 - tbl[, 'fdr_at_thresh']

  tbl$indel_rate <- tbl$ins_rate * 2
  tbl$indel_rate <- factor(tbl$indel_rate)
#  tbl <- aln.factors(tbl)

  comb.df <- data.frame()
  y.fields <- c('tpr_at_thresh', 'fpr_at_thresh', 'tpr_at_fpr', 'cor')
  for (field in y.fields) {
    print(field)
    cur.df <- tbl
    cur.df$stat <- field
    cur.df$y.val <- cur.df[, field]
    if (field == 'fpr_at_thresh') {
      cur.df$y.val <- cur.df[, field] * 100
    }
    comb.df <- rbind(comb.df, cur.df)
  }
  comb.df$stat <- factor(comb.df$stat, levels=y.fields)

  comb.df <- filter.factors(comb.df)

  pdf(file="fig_filters.pdf", width=6, height=5)
  p <- ggplot(comb.df, aes(x=length, y=y.val, linetype=filter, colour=indel_rate))
  p <- p + theme_bw()
  p  <- p + geom_line(alpha=0.7)
  line.clrs <- c('black', '#FF3333', '#FF3333')
  p <- p + scale_colour_manual(name="Indel Rate", values=line.clrs)
  line.types <- c(1, 2, 3, 4, 5, 6, 7, 8)
  p <- p + scale_linetype_manual(name="Filter", values=line.types)
  p <- p + scale_y_continuous(limits=c(0, 1))
  p <- p + scale_x_continuous(limits=c(0, 2))
  p <- p + facet_grid(stat ~ aligner)
  p <- p + coord_equal(ratio=2)
  p <- p + opts(
    axis.text.x = theme_text(angle=90, hjust=1),
    axis.title.x = theme_blank(),
    axis.title.y = theme_blank(),
#    strip.text.x = theme_blank(),
#    strip.text.y = theme_blank(),
    strip.background = theme_blank(),
#    panel.grid.minor = theme_line(colour='gray60', size=0.2),
    panel.grid.major = theme_line(colour='gray80', size=0.2)
  )
  print(p)
  dev.off()
}


multi.trees <- function() {
  tree.a <- read.tree('trees/artificial.nh')
  tree.b <- read.tree('trees/bglobin.nh')
  tree.c <- read.tree('trees/encode.nh')
  pdf(file='fig_multi_trees.pdf', width=15, height=5)
  vplayout(3, 1)

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
  for (i in 1:length(plots)) {
    cur.stuff <- plots[[i]]
    cur.tree <- cur.stuff[[1]]

    # Plot tree to the left.
    sim <- PhyloSim()
    setPhylo(sim, cur.tree)
    p <- f(sim)
    print(p, vp=subplot(i, 1))
  }
  dev.off()
}

multi.filters <- function() {

  _multi.filters(include.mafft=FALSE)
  _multi.filters(include.mafft=TRUE, file='filt_fig_multi_supp.pdf')
}

_multi.filters <- function(include.mafft=FALSE, file='filt_fig_multi.pdf') {
  tbl.a <<- read.csv('~/scratch/gj1_fig_three_a/current/table.csv')
  tbl.b <<- read.csv('~/scratch/gj1_fig_three_b/current/table.csv')

  tbl.tca <<- read.csv('~/scratch/gj1_fig_three_a/new_tcoffee/table.csv', stringsAsFactors=F)
  tbl.tcb <<- read.csv('~/scratch/gj1_fig_three_b/new_tcoffee/table.csv', stringsAsFactors=F)

  if (include.mafft) {
    tbl.c <<- read.csv('~/scratch101/gj1_fig_three_e/current/table.csv')
    tbl.c$aligner <- 'mafft'
  }


  reorder.f <- function(tbl) {
    print(levels(tbl$filter))
    unq <- unique(c('gblocks', 'guidance', 'tcoffee', 'optimal_c', 'true', as.character(tbl$filter)))
    
    tbl$filter <- factor(as.character(tbl$filter), levels=unq, labels=unq, ordered=T)
    tbl <- tbl[order(tbl$filter),]
    print(levels(tbl$filter))
    return(tbl)
  }
  tbl.a <- reorder.f(tbl.a)
  tbl.b <- reorder.f(tbl.b)
  if (include.mafft) {
    tbl.c <- reorder.f(tbl.c)
    tbl.c$aln.filt <- factor(tbl.c$filter, labels=paste(tbl.c[1,]$aligner, levels(tbl.c$filter), sep=' '))
  }

  tbl.a$aln.filt <- factor(tbl.a$filter, labels=paste(tbl.a[1,]$aligner, levels(tbl.a$filter), sep=' '))
  tbl.b$aln.filt <- factor(tbl.b$filter, labels=paste(tbl.b[1,]$aligner, levels(tbl.b$filter), sep=' '))
  #tbl.b$aln.filt <- paste(tbl.b$aligner, tbl.b$filter, sep=' ')
  
  tbl.a <- rbind(tbl.a, tbl.b)
  
  #print(head(tbl.a))

  pdf(file="filt_fig_multi.pdf", width=width+3, height=height)
  vplayout(2, 1)
  print(plot.filter.f(tbl.a, 'a'), vp=subplot(1,1))
  print(plot.filter.f(tbl.b, 'b'), vp=subplot(2,1))
  dev.off()

  pdf(file="filt_fig_multi_bw.pdf", width=width+3, height=height)
  vplayout(2, 1)
  print(plot.filter.f(tbl.a, 'a', no.color=T), vp=subplot(1,1))
  print(plot.filter.f(tbl.b, 'b', no.color=T), vp=subplot(2,1))
  dev.off()

}

plot.filter.f <- function(df, lbl, plot=F, no.color=F) {
  print(lbl)
  df <- subset(df, ins_rate > 0)
  df$remaining_fraction <- 1 - df$filtered_fraction
  df$tdr_at_thresh <- 1 - df$fdr_at_thresh

  none.aln <- subset(df, filter=='none')
  df <- subset(df, !(filter %in% c('none', 'branchlength')))

  m.by <- c('length', 'ins_rate')
  filter.fields <- c('remaining_fraction', 'fpr_at_thresh', 'tpr_at_thresh', 'tpr_at_fpr2')

  tbl.m <- df
  comb.df <- data.frame()
  for (fld in filter.fields) {
    n.fld <- paste(fld, '.none', sep='')
    none.aln[, n.fld] <- none.aln[, fld]
    tbl.m <- merge(tbl.m, none.aln[, c(n.fld, m.by)], by=m.by)

    tbl.m[, fld] <- tbl.m[, fld] / tbl.m[, n.fld]

    tbl.m[, 'cur.z'] <- tbl.m[, fld]
    tbl.m[, n.fld] <- NULL
    tbl.m[, 'field'] <- fld
    comb.df <- rbind(comb.df, tbl.m)
  }

  comb.df$tree <- factor(comb.df$field, levels=filter.fields, ordered=T)
  comb.df <- filter.factors(comb.df)

  comb.df$aligner <- comb.df$filter

  p <- multi.plot(comb.df, color.by='cur.z', prefix='', relative=T, plot=F, no.color=no.color)
  p <- p + coord_equal(ratio=10)  
  return(p)
}

multi.plots <- function() {
  tbl.a <<- read.csv('~/scratch/gj1_fig_one_a/current/table.csv')
  tbl.b <<- read.csv('~/scratch/gj1_fig_one_b/current/table.csv')
  tbl.c <<- read.csv('~/scratch/gj1_fig_one_c/current/table.csv')

  tbl.a$tree <- 'artificial'
  tbl.b$tree <- 'bglobin'
  tbl.c$tree <- 'encode'
  plots <- rbind(tbl.a, tbl.b, tbl.c)

  plots <- subset(plots, !(length %in% c(0.05, 0.1)))

  multi.plot(plots, color.by='tpr_at_fpr2', keep=c('clustalw', 'mafft', 'prank', 'prank_codon', 'none'))
  multi.plot(plots, color.by='tpr_at_thresh', keep=c('clustalw', 'prank_codon', 'none'))
  multi.plot(plots, color.by='fpr_at_thresh', keep=c('clustalw', 'prank_codon', 'none'))
  multi.plot(plots, color.by='tdr_at_thresh', keep=c('clustalw', 'prank_codon', 'none'))
}

#
# Figure 2 - Tree shapes, TPR_FPR<0.05, and ROC plots at 1.0/0.1
#
multi.plot <- function(all.tbls, color.by='tpr_at_fdr',
  facet.y = 'aligner',
  facet.x = 'tree',
  keep = NULL,
  prefix = '',
  plot = T,
  relative=F,
  no.color=F
) {

  print(color.by)
  print(nrow(all.tbls))

   if (!any(is.null(keep))) {
     all.tbls <- subset(all.tbls, aligner %in% keep)
     all.tbls <- aln.factors(all.tbls)
   }

  all.tbls[, 'tdr_at_thresh'] <- 1 - all.tbls[, 'fdr_at_thresh']

  all.tbls$fx <- all.tbls[, facet.x]
  all.tbls$fy <- all.tbls[, facet.y]

  width <- length(unique(all.tbls$fx))
  height <- length(unique(all.tbls$fy))

  if (plot) {
    pdf(file=paste(prefix, 'fig_multi_', color.by, '.pdf', sep=''), width=width+3, height=height)
  }

  # Write all-table to a file.
  write.csv(all.tbls, file=paste(prefix, 'fig_multi_.csv', sep=''), row.names=FALSE)

  z.val <- color.by
  all.tbls[, 'z_val'] <- all.tbls[, z.val]

  n <- 10
  if (z.val == 'fpr_at_thresh') {
    #all.tbls[, 'z_val'] <- log10(all.tbls[, 'z_val'])
    #clr.lim <- c(-3.2, -2)
    clr.lim <- c(0, 0.005)
  } else if (z.val == 'tpr_at_thresh') {
    clr.lim <- c(0, 1.0)
  } else if (z.val == 'tpr_at_fpr' || z.val == 'tpr_at_fdr') {
    clr.lim <- c(0, 1.0)
  } else {
    clr.lim <- c(0, 1.0)
  }

  if (relative) {
    clr.lim <- c(0, max(all.tbls$z_val))
    clr.val <- c(
      rgb(1, 0, 0),
      rgb(1, 0.5, 0.5),
      rgb(1, 0.75, 0.75),
      rgb(1, 0.9, 0.9),
      rgb(1, 1, 1),
      rgb(0.9, 0.9, 1),
      rgb(0.75, 0.75, 1),
      rgb(0.5, 0.5, 1),
      rgb(0, 0, 1)
    )  
    if (no.color) {
      clr.val <- c(
      gray(0.2),
      gray(0.4),
      gray(0.6),
      gray(0.8),
      gray(1),
      gray(0.8),
      gray(0.6),
      gray(0.4),
      gray(0.2)
      )
    }

    clr.brk <- c(0, 0.3, 0.65, 0.85, 0.95, 1.05, 1.15, 1.35, 1.7, 10)
    clr.lbl <- sprintf("%.2f - %.2f", clr.brk[-length(clr.brk)], clr.brk[-1])
    #clr.lbl <- paste(clr.brk[-length(clr.brk)], clr.brk[-1], sep=' - ')
  } else {
    clr.brk <- seq(from=clr.lim[1], to=clr.lim[2], length.out=n+1)
    if (z.val == 'fpr_at_thresh') {
      clr.lbl <- list()
      for (i in 1:(length(clr.brk)-1)) {
        #clr.lbl[[i]] <- bquote(10^.(sprintf("%.2f",clr.brk[i])) - 10^.(sprintf("%.2f",clr.brk[i+1])))
        #clr.lbl[[i]] <- sprintf("%.4f - %.4f", 10^(clr.brk[i]), 10^(clr.brk[i+1]))
        clr.lbl[[i]] <- sprintf("%.4f - %.4f", clr.brk[i], clr.brk[i+1])
      }
    } else {
      clr.lbl <- paste(sprintf("%.1f",clr.brk[-length(clr.brk)]), sprintf("%.1f", clr.brk[-1]), sep=' - ')
    }

    clrs <- brewer.pal(n=5, name='Spectral')
    clr.rmp <- colorRampPalette(clrs, bias=1.5)
    clr.val <- clr.rmp(n)
  }

  if(any(is.nan(all.tbls$z_val))) {
    all.tbls[is.nan(all.tbls[, 'z_val']), 'z_val'] <- 1
  }
  if(any(is.infinite(all.tbls$z_val))) {
    all.tbls[is.infinite(all.tbls[, 'z_val']), 'z_val'] <- 1
  }
  all.tbls[(all.tbls[, 'z_val'] > clr.lim[2]), 'z_val'] <- clr.lim[2] - abs(clr.lim[2])*0.01
  all.tbls[(all.tbls[, 'z_val'] < clr.lim[1]), 'z_val'] <- clr.lim[1] + abs(clr.lim[1])*0.01

  for (i in 1:length(clr.val)) {
    within.rng <- all.tbls$z_val >= clr.brk[i] & all.tbls$z_val < clr.brk[i+1]
    if (!any(is.na(within.rng))) {
      print(paste(i,sum(within.rng)))
      all.tbls[within.rng, 'new.z'] <- i
    }
  }
  all.tbls$z_val <- factor(all.tbls$new.z, levels=1:length(clr.val))
  clr.scale <- scale_fill_manual(color.by, values=clr.val, breaks=1:length(clr.lbl), labels=clr.lbl)

  # When doing relative plots w/ ClustalW alignments, we seem to get some NA z-vals...
  all.tbls[is.na(all.tbls$z_val), 'z_val'] <- 9

  all.tbls[as.numeric(all.tbls$z_val) < 5, 'dir'] <- 1
  all.tbls[as.numeric(all.tbls$z_val) >= 5, 'dir'] <- 2
  all.tbls$dir <- as.factor(all.tbls$dir)
  changed.tbls <- all.tbls[all.tbls$z_val != 5,]

  p <- ggplot(all.tbls, aes(x=length, y=ins_rate*2))
  p <- p + theme_bw()
  if (relative && no.color) {
    #p <- p + geom_tile(aes(fill=z_val))
    p <- p + geom_point(data=changed.tbls, aes(colour=z_val, shape=dir), size=3)
    clr.scale <- scale_colour_manual(color.by, values=clr.val, breaks=1:length(clr.lbl), labels=clr.lbl)
    p <- p + clr.scale
    p <- p + scale_shape_manual(breaks=c(1,2), values=c(18, 15))
  } else {
    p <- p + geom_tile(aes(fill=z_val))
    p <- p + clr.scale
  }

  if (relative) {
    x.brks <- c(0.4, 0.8, 1.2, 1.6, 2)
    y.brks <- c(0.04, 0.08, 0.12, 0.16, 0.2)
    p <- p + scale_y_continuous(breaks=y.brks, labels=sprintf("%.2f",y.brks))
    p <- p + scale_x_continuous(breaks=x.brks, labels=sprintf("%.1f",x.brks))
  }
  p <- p + facet_grid(fy ~ fx)
  p <- p + coord_equal(ratio=10)
  p <- p + opts(
#    legend.position='none',
    panel.grid.major = theme_blank(),
    panel.grid.minor = theme_blank(),
#    strip.text.x = theme_blank(),
#    strip.text.y = theme_blank(),
    strip.background = theme_blank(),
    axis.text.x = theme_text(size=6, angle=90, hjust=1.2),
    axis.text.y = theme_text(size=6, hjust=1.2),
    axis.title.x = theme_blank(),
    axis.title.y = theme_blank()
  )
  if (plot) {
    print(p)
    dev.off()
  } else {
    return(p)
  }
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

  keep.filters <- c('none', 'gblocks', 'branchlength_hi', 'guidance_hi', 'tcoffee_hi', 'optimal_hi', 'true')
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
      'thresh_at_fpr.true',
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

    #losers <- subset(cur.sub, truth.aln == 0 & score.aln > 2.71 & score.true < 2.71)
    losers <- subset(cur.sub, truth.aln == 0 & score.aln > thresh_at_fpr.true & score.true < thresh_at_fpr.true)

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

  write.csv(losers, file="losers.csv", row.names=F)

  prank.losers <- 0
  for (i in 1:nrow(losers)) {
    cur.loser <- losers[i,]
    prank.other <- subset(loser.others, 
      aligner.aln=='prank' &
      seq_position==cur.loser$seq_position &
      slrsim_rep==cur.loser$slrsim_rep &
      filter.aln=='none'
    )
    if (prank.other$score.aln > prank.other$thresh_at_fpr.true) {
      prank.losers <- prank.losers + 1
    }
  }
  print(prank.losers)


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

    #below <- aln.slice(pep.aln, col.lo, aln.col-1)
    #at <- aln.slice(pep.aln, aln.col, aln.col)
    #above <- aln.slice(pep.aln, aln.col+1, col.hi)        
    #spliced <- aln.splice(below, at, empty.space=0)
    #spliced <- aln.splice(spliced, above, empty.space=0)
    
    sub.aln <- aln.slice(pep.aln, col.lo, col.hi)
    return(list(aln=sub.aln, tree=tree))
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

    #print(cur_seq_position)
    #print(window.width)
    alntree <- load.window(cur.gene, 'human', cur_seq_position, window.width)
    aln <- alntree$aln
    tree <- alntree$tree

    aln <- sort.aln.by.tree(aln, tree)

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
      plot.title = lbl,
      aln.plot.chars=T,
      aln.char.text.size=2
    )
    pushViewport(subplot(1,i))
    grid.draw(p$grob)
    upViewport()
  }
  dev.off()
  
}

lots.of.fps <- function() {
  for (i in 2:2) {
    false.positives(i, file=paste("fps/",i,".pdf",sep=''))
  }
}

example.roc <- function() {
  if (!exists('example.roc.df')) {
    roc <- load.roc('gj1_fig_two_a', return.merged=F, na.rm=F)
    example.roc.df <<- roc
  }
  roc <- example.roc.df

  keep.filters <- c('tcoffee_hi')
  roc <- subset(roc, filter %in% keep.filters)
  roc <- subset(roc, tree_length == 1 & ins_rate == 0.05)

  print(nrow(roc))

  roc <- ins.mpl.factors(roc)
  roc <- aln.factors(roc)
  roc <- filter.factors(roc)
  roc <- roc[ order(roc$filter, roc$aligner), ]

  roc$filter <- as.character(roc$filter)
  roc$aligner <- as.character(roc$aligner)
#  roc[roc$filter == 'None', 'aligner'] <- 'MAFFT'
#  roc[roc$filter == 'True Alignment', 'aligner'] <- 'True Alignment'
  roc[roc$filter == 'Branch Length(high)', 'aligner'] <- 'MAFFT + Branch Length filter'

  roc$analysis <- roc$aligner
#  roc$analysis <- factor(roc$aligner, levels=unique(roc$aligner))

  n.fp <- max(roc$fp)
  n.tp <- max(roc$tp)

  add.fdrs <- function(df, p) {
    fdr <- 0.1
    p <- p + geom_abline(intercept=0, slope=(1 - fdr)/fdr * n.fp/n.tp, colour='black', size=0.5, linetype=1)
#    fdr <- 0.5
#    p <- p + geom_abline(intercept=0, slope=(1 - fdr)/fdr * n.fp/n.tp, colour='black', linetype=2)
  }
  add.fps <- function(df, p) {
    p <- p + geom_vline(x=0.05, colour='black', linetype=3)
  }

  roc$is_filtered <- FALSE
  roc[roc$score == -9999, 'is_filtered'] <- TRUE
  roc[roc$score == -9999, 'score'] <- 0

  p <- ggplot(roc, aes(x=fpr, y=tpr))
  p <- p + theme_bw()

  p <- p + geom_abline(intercept=0, slope=1, colour='gray', linetype=1)

  add.areas <- function(roc, p) {
    # AUC polygons.
    below.fpr <- subset(roc, is_filtered == FALSE)
    first.row <- below.fpr[1, ]
    last.row <- below.fpr[nrow(below.fpr),]
    last.row$tpr <- 0
    first.row$tpr <- 0
    below.fpr <- rbind(first.row, below.fpr, last.row)

    auc.tpr <- max(below.fpr$tpr)
    whole.rectangle <- below.fpr
    whole.rectangle[whole.rectangle$tpr > 0, 'tpr'] <- auc.tpr
    p <- p + geom_polygon(data=whole.rectangle, alpha=0.3, fill=rgb(.8, .8, 1), aes(colour=NULL))
    p <- p + geom_polygon(data=below.fpr, alpha=0.5, fill=rgb(.5, .5, 1), aes(colour=NULL))
  }

  add.points <- function(roc, p) {
    # Threshold point marker.
    t <- qchisq(0.975, df=1)
    for (lbl in unique(roc$analysis)) {
      df.sub <- subset(roc, analysis==lbl)
      row.index <- min(which(df.sub$score <= t))
      row <- df.sub[row.index, ]
      p <- p + geom_segment(data=row, xend=-0.1, yend=row$tpr, colour=I('black'), size=0.25, linetype=2)
  
      row.index <- max(which(df.sub$fdr <= 0.1))
      row <- df.sub[row.index, ]
      p <- p + geom_segment(data=row, xend=-0.1, yend=row$tpr, colour=I('black'), size=0.25, linetype=2)

      row.index <- max(which(df.sub$fpr <= 0.05))
      row <- df.sub[row.index, ]
      p <- p + geom_segment(data=row, xend=-0.1, yend=row$tpr, colour=I('black'), size=0.25, linetype=2)
    }
    p
  }

#  p <- p + scale_colour_brewer(name="Alignment", palette="Set1")
#   p <- p + scale_colour_manual(values=c('blue', 'red'))
  p <- p + scale_y_continuous(name="True Positive Rate", limits=c(0,1.0), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))
  p <- p + opts(legend.position = 'none') 

  p.full <- p + scale_x_continuous(name="False Positive Rate", 
    limits=c(-0.1,1.0), 
    breaks=c(-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1),
    labels=c('', 0, 0.2, 0.4, 0.6, 0.8, 1),
    expand = c(0, 0)
  )
  p.zoom <- p + scale_x_continuous(name="False Positive Rate", 
    limits=c(-0.1,1.1) * 0.1, 
    breaks=c(-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1) * 0.1,
    labels=c('', c(0, 0.2, 0.4, 0.6, 0.8, 1) * 0.1),
    expand = c(0, 0)
  )

  p.full <- add.areas(roc, p.full)
#  p.zoom <- add.areas(roc, p.zoom)

#  p.full <- add.fps(roc, p.full)
  p.zoom <- add.fps(roc, p.zoom)

  p.full <- p.full + geom_point(size=1, aes(colour=is_filtered))
  p.zoom <- p.zoom + geom_point(size=1, aes(colour=is_filtered))

  p.full <- p.full + scale_colour_manual(values=c('blue', 'black'))
  p.zoom <- p.zoom + scale_colour_manual(values=c('blue', 'black'))

  p.zoom <- add.fdrs(roc, p.zoom)
  p.zoom <- add.points(roc, p.zoom)

  pdf(file="roc.example.pdf",width=8,height=4)
    vplayout(2, 1)
    print(p.full, vp=subplot(1, 1))
    print(p.zoom, vp=subplot(2, 1))
  dev.off()
}

meta.analysis <- function() {
    asdf <- load.sites('gj1_fig_three_a')
    asdf$label <- T
    asdf$tree <- T
    asdf$analysis <- T
    asdf$true_type <- as.factor(asdf$true_type)

    print(object.sizes())

    comb.df <- data.frame()
    for (ln in unique(asdf$tree_length)) {
      for (ins in unique(asdf$ins_rate)) {
        print(paste(ln, ins, '...'))
        cur.asdf <- subset(asdf, tree_length == ln & ins_rate == ins)
        cur.df <- meta.sub(cur.asdf)
        comb.df <- rbind(comb.df, cur.df)
      }
    }   
    write.csv(comb.df, file="meta.analysis.csv")
}

meta.sub <- function(asdf) {
    prank_c <- subset(asdf, filter=='none')
    bl <- subset(asdf, filter=='branchlength')
    tc <- subset(asdf, filter=='tcoffee')
    gb <- subset(asdf, filter=='gblocks')
    true <- subset(asdf, filter=='true')

    merge.by <- c('tree_length', 'ins_rate', 'slrsim_rep', 'seq_position')

    prep.df <- function(df, lbl) {
      df <- df[!duplicated(df[, merge.by]),]
      df <- df[, c('lrt_stat', merge.by)]
      df[, paste('lrt.', lbl, sep='')] <- df$lrt_stat
      df$lrt_stat <- NULL
      return(df)
    }

    collect.scores <- function(df.list) {
      df.names <- names(df.list)
      all.sites <- prep.df(df.list[[1]], df.names[1])
  
      for (i in 2:length(df.list)) {
        cur.df <- prep.df(df.list[[i]], df.names[i])
        #print(str(cur.df))
        all.sites <- merge(all.sites, cur.df, by=merge.by, all.x=T, all.y=T)
      }
      return(all.sites)
    }

    combine.median <- function(df) {
      combine.fields <- grep("lrt\\.", colnames(df), value=T)
      #print(combine.fields)
  
      meta.matrix <- as.matrix(df[, combine.fields])
      lrts <- apply(meta.matrix, 1, median, na.rm=T)
      df[, 'lrt_stat'] <- lrts
        
      for (i in 1:length(combine.fields)) {
        df[, combine.fields[i]] <- NULL
      }
      return(df)
    }
  
    meta <- list(
      none = prank_c,
      bl = bl,
      tc = tc,
      gb = gb
    )
    #print("  collecting scores")
    meta <- collect.scores(meta)
    #print("  combining median")
    meta <- combine.median(meta)

    #print(colnames(meta))
    #print(colnames(true))
    #print("  merging true")
    true$lrt_stat <- NULL
    test.df <- merge(true, meta, by=merge.by, all.x=T)
    return(test.df)
}

add.len.ins.jobs <- function() {
  source("~/src/greg-ensembl/scripts/mysql_functions.R")

  dbname <- 'gj1_fig_three_a'
  con <- connect(dbname)

  genes <- mysqlReadTable(con, 'genes')

  len <- unique(genes$slrsim_tree_mean_path)
  ins <- unique(genes$phylosim_insertrate)

  for (l in len) {
    for (i in ins) {
      sql <- sprintf("INSERT INTO analysis_job (analysis_id, input_id, status) VALUES (4, \"{phylosim_insertrate => %.3f, slrsim_tree_mean_path => %.3f}\", 'READY')", i, l)
      print(sql)
      mysqlExecStatement(con, sql)
    }    
  }

  print(len)
  print(ins)


}

object.sizes <- function() {
  return(rev(sort(sapply(ls(envir=.GlobalEnv), function (object.name)
        object.size(get(object.name))))))
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

  ordered.f <- rev(c('clustalw', 'mafft', 'probcons', 'fsa', 'fsa_careful', 'prank', 'pagan', 'prank_codon', 'true', 'none', 'no_indels'))
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
        prank = 'PRANK_A',
        prank_codon = 'PRANK_C',
        pagan = 'PAGAN',
        pagan_groups = 'PAGAN(groups)',
        pagan_codon = 'PAGAN(codon)',
        True_Alignment = 'True Alignment',
        'True Alignment' = 'True Alignment',
        true = 'True Alignment',
        none = 'True Alignment',
        'no_indels' = 'No Indels',
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
  filters <- unique(as.character(data$filter))

  ordered.f <- c('No filter', 'none',
    'True Alignment', 'True_Alignment', 'true',
    'optimal', 'optimal_lo', 'optimal_hi', 'optimal_a', 'optimal_b', 'optimal_c',
    'tcoffee', 'tcoffee_lo', 'tcoffee_hi',
    'guidance', 'guidance_lo', 'guidance_hi',
    'prank_mean',
    'columns', 
    'branchlength', 'branchlength_lo', 'branchlength_hi',
    'gblocks'
    )
  filters <- unique(c(ordered.f[ordered.f %in% filters], filters))
  
  levels <- c()
  labels <- c()
  for (f in filters) {
    levels <- c(levels, f)
    labels <- c(labels, 
      switch(f,
        gblocks = 'Gblocks',
        guidance = 'GUIDANCE',
        'guidance_lo' = 'GUIDANCE(low)',
        'guidance_hi' = 'GUIDANCE(high)',
        prank_mean = 'PRANK',
        columns = 'Non-gap Sequences',
        branchlength = 'Branch Length',
        'branchlength_lo' = 'Branch Length(low)',
        'branchlength_hi' = 'Branch Length(high)',
        tcoffee = 'T-Coffee',
        'tcoffee_lo' = 'T-Coffee(low)',
        'tcoffee_hi' = 'T-Coffee(high)',
        'No filter' = 'None',
        'none' = 'None',
        optimal = 'Optimal',
        'optimal_c' = 'Optimal',
        'optimal_lo' = 'Optimal(low)',
        'optimal_hi' = 'Optimal(high)',
        'True Alignment' = 'True Alignment',
        True_Alignment = 'True Alignment',
        true = 'True Alignment'
      )
    )
  }
  data$filter <- factor(data$filter, levels=levels, labels=labels, ordered=T)
  return(data)
}

load.tbl <- function(db) {
  dir = paste('~/scratch/',db,'/current', sep='')
  tbl.file <- paste(dir, '/', 'table.csv', sep='')
  tbl <- read.csv(file=tbl.file, stringsAsFactors=F)
  return(tbl)
}

load.sites <- function(db) {
  if (exists(db)) {
    print(paste("Loading cached sites from",db))
    merged <- get(db)
    return(merged)
  }
  dir = paste('~/scratch/',db,'/current', sep='')
  sites.file <- paste(dir, '/', 'merged.Rdata', sep='')
  load(sites.file)
  assign(db, merged, envir=.GlobalEnv)
  return(merged)
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

load.roc <- function(db, return.merged=T, return.fps=F, na.rm=F) {
  dir <- paste('~/scratch/', db, '/current', sep='')
  sites.file <- paste(dir, '/sites.Rdata', sep='')
  rocs.file <- paste(dir, '/rocs.Rdata', sep='')
  merged.file <- paste(dir, '/rocs_merged.Rdata', sep='')
  fps.file <- paste(dir, '/rocs_fps.Rdata', sep='')

  if (!file.exists(rocs.file)) {
    print(paste("Making roc for",db,"..."))
    load(sites.file)
    roc.f <- function(df) {slr.roc(df, na.rm=na.rm)}
    roc <- summarize.by.labels(sites, roc.f)
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
      print(paste("  ", cur.f))
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
