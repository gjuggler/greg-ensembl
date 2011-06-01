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
