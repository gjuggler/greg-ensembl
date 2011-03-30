library(phylosim)
library(RColorBrewer)
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

figure.x <- function() {
  tbl.a <- read.csv('~/scratch/gj1_fig_one_a/current/table.csv')
  tbl.b <- read.csv('~/scratch/gj1_fig_one_b/current/table.csv')
  tbl.c <- read.csv('~/scratch/gj1_fig_one_c/current/table.csv')
  all.tbls <- data.frame()
  for (cur.tbl in list(tbl.a, tbl.b, tbl.c)) {
    cur.tbl <- subset(cur.tbl, analysis == 'SLR Sitewise')
    all.tbls <- rbind(all.tbls, cur.tbl)
  }

  power.field <- 'tpr_at_fpr'

  tbl <- all.tbls
  for (i in 1:nrow(tbl)) {
    cur.row <- tbl[i,]
    true.row <- subset(tbl, 
      ins_rate==cur.row$ins_rate &
      length==cur.row$length &
      aligner=='True_Alignment' &
      tree==cur.row$tree
    )
    true.row[, power.field] <- pmax(0.001, true.row[, power.field])
    tbl[i, 'rel_power'] <- cur.row[, power.field] / true.row[, power.field]
  }

  p.f <- function(aln, cmp) {

    tbl <- subset(tbl, length == 1)
    tbl <- subset(tbl, aligner != 'True_Alignment')
    tbl <- subset(tbl, aligner == aln)

    if (cmp == '6-17') {
      tbl <- subset(tbl, tree != 'encode.nh')
    } else {    
      tbl <- subset(tbl, tree != 'artificial.nh')
    }
    tbl <- aln.factors(tbl) 
    tbl <- tree.factors(tbl) 

    print(str(tbl))

    combined.tbl <- data.frame()
    lines <- data.frame()
    
    tbl[, 'match_fraction'] <- (tbl$mean_bl_match) / tbl$mean_bl_aligned
    tbl[, 'bl_mismatch'] <- (tbl$mean_bl_mismatch)
    for (field in c('total_column_score', 'sum_of_pairs_score', 'match_fraction', 'bl_mismatch')) {
      cur.tbl <- tbl
      cur.tbl$score <- cur.tbl[,field]
      cur.tbl$score_field <- field
  
      sm <- subset(cur.tbl, tree == '6-Taxon')
      md <- subset(cur.tbl, tree == '17-Taxon')
      lg <- subset(cur.tbl, tree == '44-Taxon')

      if (cmp == '6-17') {  
      line.df <- merge(sm, md, by=c('ins_rate', 'length', 'aligner', 'score_field'))
      line.df$tree <- line.df$tree.x
      lines <- rbind(lines, line.df)
      } else {
        line.df <- merge(md, lg, by=c('ins_rate', 'length', 'aligner', 'score_field'))
        line.df$tree <- line.df$tree.x
        lines <- rbind(lines, line.df)
      }
  
      combined.tbl <- rbind(combined.tbl, cur.tbl)
    }

    p <- ggplot(combined.tbl, aes(x=score, y=rel_power, colour=tree, group=tree))    
    p <- p + theme_bw()
#    p <- p + geom_point(size=0.5, alpha=0.8)

    p <- p + stat_smooth(method='lm', alpha=0.3, fullrange=FALSE, se=TRUE, linewidth=0.5)
    p <- p + geom_segment(data=lines,
      aes(x=score.x, y=rel_power.x, xend=score.y, yend=rel_power.y),
      linewidth=0.1, alpha=0.6, colour='black',
      arrow=arrow(length=unit(0.15,"cm"))
      )
    p <- p + scale_colour_hue(name="Tree", l=50)
  
    p <- p + scale_x_continuous("Alignment Score")
    p <- p + scale_y_continuous("Relative Sitewise Power")
    p <- p + facet_grid(aligner ~ score_field, scales="free")
    p <- p + opts(
      title=paste(combined.tbl[1,'aligner'])
    )
    return(p)
  }

  pdf("figure_x.pdf", height=10, width=10)
  vplayout(1, 3)
  p <- p.f('clustalw', '6-17')
  print(p, vp=subplot(1,1))
  p <- p.f('mafft', '6-17')
  print(p, vp=subplot(1,2))
  p <- p.f('prank_codon', '6-17')
  print(p, vp=subplot(1,3))
  dev.off()

  pdf("figure_y.pdf", height=10, width=10)
  vplayout(1, 3)
  p <- p.f('clustalw', '17-44')
  print(p, vp=subplot(1,1))
  p <- p.f('mafft', '17-44')
  print(p, vp=subplot(1,2))
  p <- p.f('prank_codon', '17-44')
  print(p, vp=subplot(1,3))
  dev.off()

}

try.test <- function() {

  f <- 3

  try({
    f <- 10 + 3
  })
  
  print(f)
}


#
# Figure 3 - Sitewise power vs alignment accuracy, and equivalent alignment
# accuracy in the 17-taxon vs 44-taxon trees.
#
figure.three <- function() {

  tbl.a <- read.csv('~/scratch/gj1_fig_one_a/current/table.csv')
  tbl.b <- read.csv('~/scratch/gj1_fig_one_b/current/table.csv')
  tbl.c <- read.csv('~/scratch/gj1_fig_one_c/current/table.csv')
  all.tbls <- data.frame()
  for (cur.tbl in list(tbl.a, tbl.b, tbl.c)) {
    cur.tbl <- subset(cur.tbl, analysis == 'SLR Sitewise')
    all.tbls <- rbind(all.tbls, cur.tbl)
  }

  tbl <- all.tbls
  #tbl <- subset(tbl, aligner != 'clustalw')

  for (i in 1:nrow(tbl)) {
    cur.row <- tbl[i,]

    true.row <- subset(tbl, 
      ins_rate==cur.row$ins_rate &
      length==cur.row$length &
      aligner=='True_Alignment' &
      tree==cur.row$tree
    )
    tbl[i, 'tpr_reduction'] <- cur.row$tpr_at_fpr / (true.row$tpr_at_fpr+0.001)
  }
  tbl[, 'tree'] <- factor(tbl[, 'tree'], 
    levels=c('artificial.nh', 'bglobin.nh', 'encode.nh'), 
    labels=c('6-taxon', '17-taxon', '44-taxon')
  )

  plot.f <- function(df, aln_score_field, lbl) {
    df[,'score'] <- df[,aln_score_field]
    p <- ggplot(df, aes(
      x=score,
      y=tpr_reduction,
      group=tree,
      colour=tree
    ))
    p <- p + stat_smooth(method='loess', alpha=0.4)
    p <- p + geom_point(size=1, alpha=0.9)
    p <- p + scale_colour_hue(name="Tree", l=50)
    p <- p + scale_x_continuous(lbl)
    p <- p + scale_y_continuous("Power (Relative to True Alignment)", limits=c(0,1))
    p <- p + opts(
      title=paste("Loss in Power vs ", lbl, sep='')
    )
    p <- p + facet_grid( . ~ aligner)
    return(p)
  }

  plot.g <- function(df, tree_a, tree_b, field, lbl) {
    # Find the equivalent rows for tree A and tree B
    a <- subset(df, tree==tree_a)
    a <- subset(a, aligner!='True_Alignment')

    for (i in 1:nrow(a)) {
      row.b <- subset(df, 
        tree==tree_b &
        ins_rate==a[i,]$ins_rate &
        length==a[i,]$length &
        aligner==a[i,]$aligner
      )
      a[i, 'tree_b_score'] <- row.b[,field]
    }
    
    a[,'value'] <- a[,field]

    a <- aln.factors(a)

    p <- ggplot(a, aes(x=value, y=tree_b_score, group=aligner, colour=aligner))
    p <- p + geom_abline(intercept=0, slope=1, colour='black', linetype='dashed')
    p <- p + scale_colour_hue("Alignment Method", l=50)
    p <- p + stat_smooth(method='loess', alpha=0.4)
    p <- p + geom_point(size=1, alpha=0.9)
    p <- p + scale_x_continuous(paste(lbl, " in ", tree_a, " tree", sep=''))
    p <- p + scale_y_continuous(paste(lbl, " in ", tree_b, " tree", sep=''))
    p <- p + opts(
      title=paste(lbl, " in ",tree_b," vs ", tree_a, " Tree", sep='')
    )

    return(p)
  }

  plot.h <- function(df, tree_a, tree_b, field, lbl) {
    # Find the equivalent rows for tree A and tree B
    a <- subset(df, tree==tree_a)
    a <- subset(a, aligner!='True_Alignment')

    a[, 'tree_a_score'] <- a[, field]
    for (i in 1:nrow(a)) {
      row.b <- subset(df, 
        tree==tree_b &
        ins_rate==a[i,]$ins_rate &
        length==a[i,]$length &
        aligner==a[i,]$aligner
      )
      a[i, 'tree_b_score'] <- row.b[,field]
    }    
    a[, 'tree_score_ratio'] <- a$tree_b_score - a$tree_a_score
    a <- aln.factors(a)

    p <- ggplot(a, aes(x=tree_score_ratio, fill=aligner))
    p <- p + geom_histogram()
    p <- p + geom_vline(xintercept=0)
    p <- p + scale_fill_hue("Aligner", l=50)
    p <- p + scale_x_continuous(substitute(paste(Delta, lbl))) #, limits=c(-0.1, 0.1))
    p <- p + facet_grid(aligner ~ .)
    p <- p + opts(
      title=substitute(paste(Delta, lbl, " (",tree_b," vs ", tree_a, ")", sep=''))
    )
    return(p)
  }

  plot.i <- function(df, low.tree, field, lbl) {
    power.field <- 'tpr_at_fpr'

    # Find the equivalent rows for tree A and tree B
    get.tree.comparison <- function(df, tree_a, tree_b) {
      a <- subset(df, tree==tree_a)
      a <- subset(a, aligner!='True_Alignment')
      a[, 'tree_a_score'] <- a[, field]
      a[, 'tree_a_power'] <- a[, power.field]
      for (i in 1:nrow(a)) {
        row.b <- subset(df, 
          tree==tree_b &
          ins_rate==a[i,]$ins_rate &
          length==a[i,]$length &
          aligner==a[i,]$aligner
        )
        a[i, 'tree_b_score'] <- row.b[,field]
        a[i, 'tree_b_power'] <- row.b[, power.field]
      }    
      a[, 'tree_score_ratio'] <- a$tree_b_score - a$tree_a_score
      a <- aln.factors(a)
    a[, 'tree_power_ratio'] <- a$tree_b_power - a$tree_a_power
    return(a)
  }

  tree_a <- '6-taxon'
  tree_b <- '17-taxon'
  if (low.tree == '17-taxon') {
    tree_a <- '17-taxon'
    tree_b <- '44-taxon'
  }
  all.df <- get.tree.comparison(df, tree_a, tree_b)

    p <- ggplot(all.df, aes(x=tree_score_ratio, y=tree_power_ratio, colour=aligner))
    p <- p + geom_hline(yintercept=0, colour='black')
    p <- p + geom_vline(yintercept=0, colour='black')
    p <- p + scale_colour_hue("Aligner", l=50)
    p <- p + stat_smooth(method='lm', alpha=0.4, fullrange=FALSE, se=FALSE)
    p <- p + geom_point(size=1.3, alpha=0.8)
    p <- p + scale_x_continuous(substitute(paste(Delta, lbl))) #, limits=c(-0.1, 0.1))
    p <- p + scale_y_continuous(substitute(paste(Delta, power.field)))
    p <- p + opts(
      title=substitute(paste(Delta, power.field, " vs ", Delta, lbl, " (", tree_b, " vs ",tree_a,")",sep=""))
    )
    return(p)
  }

  tbl <- subset(tbl, aligner != 'True_Alignment')

  print(str(tbl))
  tbl[, 'new_score'] <- tbl$mean_bl_match * tbl$aln_length - tbl$mean_bl_mismatch * tbl$aln_length

  pdf(file="figure_three.pdf", width=20, height=16)
  vplayout(3, 4)

#  sub.tbl <- subset(tbl, aligner != 'clustalw')
  sub.tbl <- tbl
  p <- plot.i(sub.tbl, '17-taxon', 'total_column_score', 'TCS')
  print(p, vp=subplot(1,4))
  p <- plot.i(sub.tbl, '17-taxon', 'sum_of_pairs_score', 'SPS')
  print(p, vp=subplot(2,4))
  p <- plot.i(sub.tbl, '17-taxon', 'new_score', 'new_score')
  print(p, vp=subplot(3,4))

  sub.tbl <- tbl
  p <- plot.i(sub.tbl, '6-taxon', 'total_column_score', 'TCS')
  print(p, vp=subplot(1,2))
  p <- plot.i(sub.tbl, '6-taxon', 'sum_of_pairs_score', 'SPS')
  print(p, vp=subplot(2,2))
  p <- plot.i(sub.tbl, '6-taxon', 'new_score', 'new_score')
  print(p, vp=subplot(3,2))

  p <- plot.h(tbl, '6-taxon', '17-taxon', 'total_column_score', 'TCS')
  print(p, vp=subplot(1,1))
  p <- plot.h(tbl, '17-taxon', '44-taxon', 'total_column_score', 'TCS')
  print(p, vp=subplot(1,3))

  p <- plot.h(tbl, '6-taxon', '17-taxon', 'sum_of_pairs_score', 'SPS')
  print(p, vp=subplot(2,1))
  p <- plot.h(tbl, '17-taxon', '44-taxon', 'sum_of_pairs_score', 'SPS')
  print(p, vp=subplot(2,3))

  p <- plot.h(tbl, '6-taxon', '17-taxon', 'new_score', 'new_score')
  print(p, vp=subplot(3,1))
  p <- plot.h(tbl, '17-taxon', '44-taxon', 'new_score', 'new_score')
  print(p, vp=subplot(3,3))

  dev.off()

  return()

  pdf(file="supp_figure_one.pdf", width=24, height=24)
  vplayout(3, 3)

  p <- plot.g(tbl, '6-taxon', '17-taxon', 'total_column_score', "Total Column Score")
  print(p, vp=subplot(1,1))
  p <- plot.g(tbl, '6-taxon', '17-taxon', 'sum_of_pairs_score', "Sum of Pairs Score")
  print(p, vp=subplot(2,1))
  p <- plot.g(tbl, '6-taxon', '17-taxon', 'match_bl_score', "Match BL Score")
  print(p, vp=subplot(3,1))

  p <- plot.g(tbl, '17-taxon', '44-taxon', 'total_column_score', "Total Column Score")
  print(p, vp=subplot(1,2))
  p <- plot.g(tbl, '17-taxon', '44-taxon', 'sum_of_pairs_score', "Sum of Pairs Score")
  print(p, vp=subplot(2,2))
  p <- plot.g(tbl, '17-taxon', '44-taxon', 'match_bl_score', "Match BL Score")
  print(p, vp=subplot(3,2))

  p <- plot.f(tbl, 'total_column_score', 'Total Column Score')
  print(p, vp=subplot(1, 3))
  p <- plot.f(tbl, 'sum_of_pairs_score', 'Sum of Pairs Score')
  print(p, vp=subplot(2, 3))
  p <- plot.f(tbl, 'match_bl_score', 'Match BL Score')
  print(p, vp=subplot(3, 3))

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

figure.four.new <- function() {

  load.roc <- function(d) {
    sites.file <- paste(d, '/', 'sites.Rdata', sep='')
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
  
  roc.b <- load.roc('~/scratch/gj1_fig_two_b/current')
  print(str(roc.b))
  roc.e <- load.roc('~/scratch/gj1_fig_two_e/current')

  roc <- rbind(roc.b, roc.e)

  roc <- subset(roc, 
    (ins_rate == 0.025 & tree_length == 1) |
    (ins_rate == 0.05 & tree_length == 1) |
    (ins_rate == 0.05 & tree_length == 2) |
    (ins_rate == 0.1 & tree_length == 2)
  )
  roc <- subset(roc, filter != 'columns')
  roc <- filter.factors(roc)
  roc <- ins.mpl.factors(roc, short=T)
  roc <- aln.factors(roc)

  roc$grp <- paste(roc$tree_length, roc$ins_rate, sep=" ")
  print(unique(roc$grp))
  print(unique(roc$aligner))

  p <- plot.roc(roc,plot=F,plot.x='fpr',plot.y='tpr', color.by='filter')
  p <- p + geom_vline(x=0.01, colour='gray', linetype='dashed')
  p <- p + geom_vline(x=0.05, colour='black', linetype='dashed')
  p <- p + scale_colour_brewer(name="Filtering Method", palette="Set1")
  p <- p + facet_grid(aligner ~ grp)
  p <- p + opts(
    panel.margin=unit(rep(0.5, times=4), "lines"),
    title='TPR vs FPR for 17-Taxon Tree'
  )
  p <- p + scale_y_continuous(name="True Positive Rate", limits=c(0,1.0), breaks=c(0,0.5,1))
  p <- p + scale_x_continuous(name="False Positive Rate", limits=c(0,0.1), breaks=c(0, 0.05, 0.1))

  pdf(file="figure_four_new.pdf",width=16,height=8)
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

figure.ari <- function() {
  dir <- '~/scratch/gj1_ari_indels/current'
  sites.file <- paste(dir, '/', 'sites.Rdata', sep='')
  roc.file <- paste(dir, '/', 'roc.Rdata', sep='')
  tbl.file <- paste(dir, '/', 'table.csv', sep='')

  print("  scatter")
  pdf(file="ari_bars.pdf")  
  tbl <- read.csv(file=tbl.file)

  tbl <- aln.factors(tbl)
  tbl <- ins.mpl.factors(tbl)

  tbl$x <- as.integer(tbl$aligner)
  tbl$dx <- 0.3
    p <- ggplot(tbl, aes(fill=aligner))
    p <- p + geom_rect(aes(xmin=x-dx, xmax=x+dx, ymin=0, ymax=tpr_at_fpr))
    x.lbl <- aln.labels
    p <- p + scale_x_continuous("Aligner", breaks=c(1:length(x.lbl)), labels=x.lbl)
    p <- p + scale_y_continuous("Sitewise TPR at FPR<0.05")
    p <- p + facet_grid(ins_rate ~ length)
    p <- p + opts(
    axis.text.x = theme_text(angle=90, hjust=1),
    legend.position = 'none'
    )
    print(p)
  dev.off()

  print("  ROC")
  if (!file.exists(roc.file)) {
    load(sites.file)

    f <- function(df, thresh) {return(slr.roc(df, na.rm=T))}
    roc <- summarize.by.labels(sites, f)
    save(roc, file=roc.file)
  }
  load(roc.file)

  roc <- subset(roc, aligner != 'clustalw')
  roc <- subset(roc, aligner != 'mafft')

  roc <- aln.factors(roc)
  roc <- ins.mpl.factors(roc)

  p <- plot.roc(roc,plot=F,plot.x='fpr',plot.y='tpr', color.by='aligner')
  p <- p + geom_vline(x=0.01, colour='gray', linetype='dashed')
  p <- p + geom_vline(x=0.05, colour='black', linetype='dashed')
  p <- p + facet_grid(ins_rate ~ tree_length)
  p <- p + opts(
    panel.margin=unit(rep(0.5, times=4), "lines"),
    title="True Positive Rate vs False Positive Rate for 17-Taxon Tree"
  )
  p <- p + scale_y_continuous(name="True Positive Rate", limits=c(0,1.0), breaks=c(0,0.5,1))
  p <- p + scale_x_continuous(name="False Positive Rate", limits=c(0,0.1), breaks=c(0, 0.05, 0.1))

  pdf(file="ari_roc.pdf",width=10,height=10)
  print(p)
  dev.off()
}

false.pos.stats <- function() {

  p.s <- function(dir, filter, field, xlim, plot.diff=F) {
    roc.file <- paste(dir, '/', 'false.pos.stats.Rdata', sep='')
    sites.file <- paste(dir, '/', 'sites.Rdata', sep='')
    if (!file.exists(roc.file)) {
      load(sites.file)
  
      sites <- subset(sites, tree_length == 1 & ins_rate == 0.025)
      print(nrow(sites))
      roc <- summarize.by.labels(sites, slr.roc)
  
      roc.aln <- subset(roc, filter != 'True Alignment')
      print(nrow(roc.aln))
      roc.true <- subset(roc, filter == 'True Alignment')
      print(nrow(roc.true))
  
      # For some reason we get duplicates for a given rep and seq_position. Remove them for both sets.
      remove.dups <- function(df) {
        return(df[!duplicated(df[, c('slrsim_rep', 'seq_position')]),])
      }
      roc.true <- remove.dups(roc.true)

      all.merged <- data.frame()
      filters <- unique(roc.aln$filter)
      for (i in 1:length(filters)) {
        cur.f <- filters[i]
        sub.roc.aln <- subset(roc.aln, filter == cur.f)
        sub.roc.aln <- remove.dups(sub.roc.aln)
        roc.merged <- merge(
          sub.roc.aln, 
          roc.true,
          by=c('tree_length', 'ins_rate', 'slrsim_rep', 'seq_position'),
          suffixes=c('.aln', '.true')
        )
        print(paste(cur.f, nrow(roc.merged)))
        all.merged <- rbind(all.merged, roc.merged)
      }
      roc <- all.merged
      print(nrow(roc))
      save(roc, file=roc.file)
    }
    load(roc.file)

    df <- roc

    print(paste(filter, field))
    roc <- subset(df, filter.aln == filter)
    #print(nrow(roc))

    score.thresh <- roc[1, 'thresh_at_fpr.true']
    roc <- subset(roc, truth.aln == 0 & score.aln > score.thresh & score.true < score.thresh)

    print(nrow(roc))
 
    roc[, 'bl_mismatch.aln'] <- roc$bl_aligned.aln - roc$bl_match.aln
    roc[, 'bl_mismatch.true'] <- 0

    cur.roc <- roc
    cur.roc[, 'stat_aln'] <- cur.roc[, paste(field,'.aln',sep='')]
    cur.roc[, 'stat_true'] <- cur.roc[, paste(field,'.true',sep='')]

    cur.roc$score_diff <- cur.roc$score.aln - cur.roc$score.true;

    cur.roc$score_diff <- round(cur.roc$score_diff)
    cur.roc[cur.roc$score_diff > 4, 'score_diff'] <- 4
    cur.roc$score_diff <- factor(cur.roc$score_diff, levels=rev(c(0:4)))
    cur.roc <- cur.roc[order(cur.roc$score_diff),]

    bin.width <- diff(xlim) / 20

    ramp <- colorRampPalette(c('blue','red'))

    p <- ggplot(cur.roc, aes(x=stat_aln, fill=score_diff))
    p <- p + theme_bw()
    p <- p + geom_bar(binwidth=bin.width)
    p <- p + scale_fill_manual(values=rev(ramp(5)))
    p <- p + geom_vline(xintercept=0)
    p <- p + scale_x_continuous(paste(field, "(aln)",sep=' '), limits=xlim)
    p <- p + opts(
     title=paste(filter, field)
    )
    return(p)
  }

  p.f <- function(dir, filter, row) {
    p <- p.s(dir, filter, 'bl_aligned', c(-0.5, 4.5))
    print(p, vp=subplot(1,row))
    p <- p.s(dir, filter, 'bl_mismatch', c(-0.5,4.5))
    print(p, vp=subplot(2,row))
  }

  dir.b <- '~/scratch/gj1_fig_two_b/current'
  dir.e <- '~/scratch/gj1_fig_two_e/current'

  pdf("figure_seven_mafft.pdf")
  vplayout(2, 4)
  p.f(dir.b, 'No filter', 1)
  p.f(dir.b, 'branchlength', 2)
  p.f(dir.b, 'tcoffee', 3)
  p.f(dir.b, 'optimal', 4)
  dev.off()

  pdf("figure_seven_prank.pdf")
  vplayout(2, 4)
  p.f(dir.e, 'No filter', 1)
  p.f(dir.e, 'branchlength', 2)
  p.f(dir.e, 'tcoffee', 3)
  p.f(dir.e, 'optimal', 4)
  dev.off()
}


false.positives <- function(cur_aligner) {
  dir <- '~/scratch/gj1_ari_indels/current'
  data.dir <- paste(dir, '/', 'data', sep='')
  sites.file <- paste(dir, '/', 'sites.Rdata', sep='')
  genes.file <- paste(dir, '/', 'genes.Rdata', sep='')
  roc.file <- paste(dir, '/', 'roc.Rdata', sep='')
  tbl.file <- paste(dir, '/', 'table.csv', sep='')

  if (!file.exists(roc.file)) {
    load(sites.file)
    roc <- summarize.by.labels(sites, slr.roc)
    save(roc, file=roc.file)
  }
  load(roc.file)
  load(genes.file)

  for (cur_ins_rate in unique(roc$ins_rate)) {
    for (cur_tree_length in unique(roc$tree_length)) {
      if (cur_tree_length != 1) {
        next
      }
      if (cur_ins_rate != 0.025) {
        next
      }

      cur.sub <- subset(roc, tree_length == cur_tree_length & ins_rate == cur_ins_rate & aligner == cur_aligner)
      true.sub <- subset(roc, tree_length == cur_tree_length & ins_rate == cur_ins_rate & aligner == "True_Alignment")
      test <- merge(cur.sub, true.sub[, c('slrsim_rep', 'seq_position', 'score')], by=c('slrsim_rep', 'seq_position'), suffixes=c('.aln', '.true'))
      losers <- subset(test, truth == 0 & score.aln > 2.71 & score.true < 2.71)

      top.fps <- losers
      if (nrow(top.fps) > 5) {
        top.fps <- top.fps[1:5,]
      }

      load.window <- function(x, seq_id, seq_position, flank.width) {
        tree.f <- paste(data.dir, '/', x$tree_file, sep='')
        aln.f <- paste(data.dir, '/', x$aln_file, sep='')
#        if (file.exists(tree.f)) {
#          tree <- read.tree(tree.f)
#        }
        aln <- read.aln(aln.f)
        pep.aln <- aln.tx(aln)

        seq_id = 'human'
        aln.col <- aln.col.position(pep.aln, seq_id, seq_position)

        col.lo <- max(1, aln.col - flank.width)
        col.hi <- min(aln.length(pep.aln), aln.col + flank.width)

        below <- aln.slice(pep.aln, col.lo, aln.col-1)
        at <- aln.slice(pep.aln, aln.col, aln.col)
        above <- aln.slice(pep.aln, aln.col+1, col.hi)
        
        spliced <- aln.splice(below, at, empty.space=0)
        spliced <- aln.splice(spliced, above, empty.space=0)
        return(spliced)
      }

      window.width <- 20

      all.alignments <- list()
      aligners <- unique(c(cur_aligner, 'mafft', 'prank', 'prank_codon', 'True_Alignment'))
      for (i in 1:length(aligners)) {
        a <- aligners[i]
        cur.list <- list()

        for (j in 1:nrow(top.fps)) {
          cur.site <- top.fps[j,]
          cur_seq_position <- cur.site$seq_position
          cur_slrsim_rep = cur.site$slrsim_rep

          cur.gene <- subset(genes,
            ins_rate == cur_ins_rate &
            tree_length == cur_tree_length &
            slrsim_rep == cur_slrsim_rep & 
            aligner == a
          )
          aln <- load.window(cur.gene, 'human', cur_seq_position, window.width)

          cur.list[[j]] <- aln
        }
        all.alignments[[i]] <- cur.list
      }
      names(all.alignments) <- aligners

      pdf.width <- nrow(top.fps) * (window.width*2+3)
      pdf.height <- 17 * length(all.alignments)
      pdf(file="test.pdf", width=pdf.width/5, height=pdf.height/5)
      vplayout(nrow(top.fps), length(all.alignments))
      for (i in 1:length(all.alignments)) {
        cur.alns <- all.alignments[[i]]
        aln.name <- names(all.alignments)[i]
        print(aln.name)
        for (j in 1:nrow(top.fps)) {
          print(j)
          cur.aln <- cur.alns[[j]]
          n.seqs <- aln.num.seqs(cur.aln)
          p <- aln.plot(cur.aln)
          p <- p + geom_rect(xmin=window.width+.5, xmax=window.width+1.5, ymin=0.5, ymax=n.seqs+0.5, fill=NA, colour='black', border=1, alpha=0.5)
          p <- p + opts(title=aln.name)
          print(p, vp=subplot(j,i))
        }
      }
      dev.off()
      return()
    }
  }
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

  if (is.factor(aligners)) {
    #print(levels(aligners))
    #return(data)
  }

  levels <- c()
  labels <- c()
  for (aligner in aligners) {
    levels <- c(levels, aligner)
    labels <- c(labels, 
      switch(aligner,
        clustalw = 'ClustalW',
        mafft = 'MAFFT',
        prank = 'PRANK',
        prank_codon = 'PRANK(codon)',
        pagan = 'PAGAN',
        pagan_groups = 'PAGAN(groups)',
        pagan_codon = 'PAGAN(codon)',
        True_Alignment = 'True Alignment',
        'True Alignment' = 'True Alignment',
        true = 'True Alignment'
      )
    )
  }
  data$aligner <- factor(data$aligner, levels=levels, labels=labels)
  return(data)
}

ins.mpl.factors <- function(data, short=F) {
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

  ordered.f <- c('No filter', 'gblocks', 'guidance', 
    'columns', 'branchlength', 'tcoffee',
    'optimal', 'True Alignment', 'True_Alignment', 'true')
  filters <- ordered.f[ordered.f %in% filters]
  
  levels <- c()
  labels <- c()
  for (f in filters) {
    levels <- c(levels, f)
    labels <- c(labels, 
      switch(f,
        gblocks = 'Gblocks',
        optimal = 'Optimal',
        tcoffee = 'T-Coffee',
        guidance = 'GUIDANCE',
        'No filter' = 'None',
        columns = 'Non-gap Sequences',
        branchlength = 'Branch Length',
        'True Alignment' = 'True Alignment',
        True_Alignment = 'True Alignment',
        true = 'True Alignment'
      )
    )
  }
  data$filter <- factor(data$filter, levels=levels, labels=labels)
  return(data)
}