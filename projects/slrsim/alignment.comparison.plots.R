source("combined.plot.R")

figure.ari <- function() {
  dir <- '~/scratch/gj1_ari_indels/current'
  sites.file <- paste(dir, '/', 'sites.Rdata', sep='')
  tbl.file <- paste(dir, '/', 'table.csv', sep='')

  print("  table")
  tbl <- read.csv(file=tbl.file)
  tbl <- tbl[with(ret.df, order(length,ins_rate, aligner)),]
  tbl <- tbl[, c(
    'sum_of_pairs_score', 'total_column_score', 'mean_bl_mismatch', 'aln_length', 'lambda', 
    'auc', 'tpr_at_fpr', 'tpr_at_fdr', 'fdr_at_thresh'
  )]
  tbl.table <- xtable(tbl)
  print(tbl.table, file="ari_indels_table.html", type='html')
  return()

  print("  bars")
  pdf(file="ari_indels_bars.pdf")  
  tbl <- read.csv(file=tbl.file)
  tbl$lbl <- paste(tbl$aligner, tbl$ins_rate, tbl$length, sep='--')

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
  roc.file <- paste(dir, '/', 'fig_roc.Rdata', sep='')
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

  pdf(file="ari_indels_roc.pdf",width=10,height=10)
  print(p)
  dev.off()
}
