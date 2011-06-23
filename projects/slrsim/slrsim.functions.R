library(plyr)
library(doBy)
library(qvalue)

adj.paper.table <- function(df) {
  paper.table(df, adjust.pvals=T)
}

paper.table <- function(df, adjust.pvals=F) {

  print("  cor")
  print(head(df))
  cor = cor(df$true_dnds, df$aln_dnds, method='spearman', use='complete.obs')
  cor_pearson = cor(df$true_dnds, df$aln_dnds, method='pearson', use='complete.obs')

  print("  roc")
  roc <- slr.roc(df)

  # Use the technique where the p-value can be doubled because we're just looking
  # at positive selection. This gives thresholds of 2.71 and 5.41 for 5% and 1% FPR.
  print("  stats")
  stats.05 <- df.stats(df, thresh=2.71, paml_thresh=0.95, adjust.pvals=adjust.pvals)

  rrow <- roc[1,]
  row <- df[1,]

  n.genes <- length(unique(df$slrsim_rep))

  ret.df <- data.frame(
    label = row$label,
    tree     = row$tree,
    analysis = row$analysis,
    length   = row$tree_length,
    ins_rate    = row$ins_rate,
    aligner  = row$aligner,
    filter   = row$filter,

    `tpr_at_fpr`    = rrow$tpr_at_fpr,
    `tp_at_fpr`     = rrow$tp_at_fpr,
    `fp_at_fpr`     = rrow$fp_at_fpr,
    `tpr_at_fpr2`    = rrow$tpr_at_fpr2,
    `tp_at_fpr2`     = rrow$tp_at_fpr2,
    `fp_at_fpr2`     = rrow$fp_at_fpr2,
    `tpr_at_thresh` = stats.05$tpr,
    `fpr_at_thresh` = stats.05$fpr,
    `fdr_at_thresh`    = stats.05$fdr,
    `tp_at_thresh` = stats.05$tp,
    `fp_at_thresh` = stats.05$fp,
    `thresh_at_fpr`    = rrow$thresh_at_fpr,
    `thresh_at_fpr2`    = rrow$thresh_at_fpr2,
    `cor`   = cor,
    `cor_pearson` = cor_pearson,
    `n_sites` = nrow(roc),
    `n_genes` = n.genes
  )

  #aln.acc.df <- df.alignment.accuracy(df)
  #ret.df <- cbind(ret.df, aln.acc.df)

  print(ret.df[, c('tpr_at_fpr', 'tpr_at_thresh', 'fpr_at_thresh')])
  return(ret.df)
}


slr.roc = function(df, na.value=-9999) {
  df$score = df$lrt_stat

  # Fix NA rows to a very low score.
  df[is.na(df$lrt_stat), 'score'] <- na.value

  df$truth = as.integer( df$true_dnds > 1 )
  df <- df[order(-df$score), ]

  df$tp = cumsum(df$truth)
  df$tn = cumsum(1-df$truth)
  df$count = cumsum(rep(1,nrow(df)))
  df$p <- 1:nrow(df) # Count of positive calls.
  df$fp = df$tn # Count of false positive calls.

  df$fpr = df$tn / max(df$tn)
  df$tpr = df$tp / max(df$tp)

  df$fdr = df$fp/(df$count)

  n <- nrow(df)

  df$auc_full <- 0
  df$auc <- 0
  get.auc <- FALSE
  if (n > 1 && get.auc) {
    auc.df <- subset(df, score > na.value)
    auc <- area.under.curve(auc.df, x.lim=0.1)
    df[df$score > na.value, ]$auc <- auc
    auc_full <- area.under.curve(auc.df, x.lim=1)
    df[df$score > na.value, ]$auc_full <- auc_full
  }

  df$tpr_at_fpr <- 0
  df$fdr_at_fpr <- 1
  df$tp_at_fpr <- 0
  df$fp_at_fpr <- 0
  df$thresh_at_fpr <- max(df$lrt_stat)
  if (n > 1) {
    fpr.sub <- subset(df, fpr < 0.05)
    if (nrow(fpr.sub) > 0) {
      # Take the last (rightmost) row with fpr < 0.05
      row <- fpr.sub[nrow(fpr.sub), ]
      df$tpr_at_fpr <- row$tpr
      df$fdr_at_fpr <- row$fdr
      df$tp_at_fpr <- row$tp
      df$fp_at_fpr <- row$fp
      df$thresh_at_fpr <- row$lrt_stat
    }
  }

  df$tpr_at_fpr2 <- 0
  df$fdr_at_fpr2 <- 1
  df$tp_at_fpr2 <- 0
  df$fp_at_fpr2 <- 0
  df$thresh_at_fpr2 <- max(df$lrt_stat)
  if (n > 1) {
    fpr.sub <- subset(df, fpr < 0.01)
    if (nrow(fpr.sub) > 0) {
      # Take the last (rightmost) row with fpr < 0.01
      row <- fpr.sub[nrow(fpr.sub), ]
      df$tpr_at_fpr2 <- row$tpr
      df$fdr_at_fpr2 <- row$fdr
      df$tp_at_fpr2 <- row$tp
      df$fp_at_fpr2 <- row$fp
      df$thresh_at_fpr2 <- row$lrt_stat
    }
  }

  return(df)
}

area.under.curve <- function(roc,x.field='tn',y.field='tp',x.lim=1) {
  roc$tmp <- roc[,x.field]
  max.tmp <- max(roc$tmp)
  roc <- subset(roc,(tmp/max.tmp) <= x.lim)
  widths <- c(0,diff(roc[,x.field]))
  heights <- roc[,y.field]
  
  areas <- heights * widths
  filled.area <- sum(areas)

  total.area <- max(roc[,x.field]) * max(roc[,y.field])
  total.rows <- nrow(roc)

  roc$tmp <- NULL
  return(filled.area/total.area)
}

summarize.by.labels = function(data,fn,thresh=3.8) {
  library(plyr)
  a = ddply(data,.(slrsim_label),fn)
  return(a)
}

do.by.labels = function(data,fn) {
  library(plyr)
  d_ply(data,.(slrsim_label),fn)
}

is.paml = function(df) {
  lrt.range <- range(df$lrt_stat, na.rm=T)
  if (lrt.range[2] <= 1) {
    return(TRUE)
  }
  if (length(grep("paml", df[1, 'analysis'], ignore.case=T)) > 0) {
    return(TRUE)
  }
  return(FALSE)
}

combine.pvals <- function(df.a, df.b, label) {  
  merge.by <- c('tree_length', 'ins_rate', 'slrsim_rep', 'seq_position')

  remove.dups <- function(df) {
    return(df[!duplicated(df[, merge.by]),])
  }
  df.a <- remove.dups(df.a)
  df.b <- remove.dups(df.b)

  merged <- merge(df.a, df.b[, c('lrt_stat', merge.by)],
    by=merge.by,
    suffixes=c('.a', '.b')
  )

  print(paste("before removing NAs",nrow(merged)))
  merged <- merged[!is.na(merged$lrt_stat.a),]
  merged <- merged[!is.na(merged$lrt_stat.b),]
  print(paste("after removing NAs",nrow(merged)))
  merged$aligner <- label  
  merged$slrsim_label <- paste(merged$slrsim_label, merged$aligner, sep=' / ')

  merged$lrt_stat = pmin(merged$lrt_stat.a, merged$lrt_stat.b)
  merged$lrt_stat.a <- NULL
  merged$lrt_stat.b <- NULL

  return(merged)
}

df.stats = function(df,
  thresh=3.8,
  paml_thresh=0.95,
  adjust.pvals=F) {

  if (adjust.pvals) {
    aa <- df[1, ]$analysis_action
    if (aa == 'paml_m8' || aa == 'paml_m2') {
      print("PAML!!!")

      lnl_alt <- df$lnl_alt
      lnl_null <- df$lnl_null
      lnl_diff <- 2* (lnl_alt - lnl_null)
      chisq.t <- qchisq(0.95, df=2)

      n.nonsig <- sum(lnl_diff < chisq.t)
      print(paste(n.nonsig, "non-significant LRT sites!"))

      # set all posterior probs to 0 if the LRT isn't significant at 5%
      df[lnl_diff < chisq.t, 'lrt_stat'] <- 0
    } else {
      # Do the gene-by-gene Hochberg '88 correction to p-values.
      df[is.na(df$lrt_stat), 'lrt_stat'] <- 0      
      df$pval.tmp <- pchisq(abs(df$lrt_stat), df=1, lower.tail=F)
      df <- ddply(df, c('label', 'slrsim_rep'), function(x) {
        x$pval.tmp <- p.adjust(x$pval.tmp, method='hochberg')
        x
      })
      # Set all LRT values to 0 if the adjusted p-value isn't significant at 5%
      df[df$pval.tmp > 0.05, 'lrt_stat'] <- 0
      df$pval.tmp <- NULL
    }
  }

  pos_pos = nrow(subset(df,true_type=="positive1" & lrt_stat>thresh))
  neg_pos = nrow(subset(df,true_type!="positive1" & lrt_stat>thresh))
  neg_neg = nrow(subset(df,true_type!="positive1" & !(lrt_stat>thresh)))
  pos_neg = nrow(subset(df,true_type=="positive1" & !(lrt_stat>thresh)))

  pos_all = nrow(subset(df,true_type=="positive1"))
  neg_all = nrow(subset(df,true_type!="positive1"))
  all_pos = nrow(subset(df,lrt_stat>thresh))
  all_neg = nrow(subset(df,!(lrt_stat>thresh)))

  all = nrow(df)

  # see http://en.wikipedia.org/wiki/Receiver_operating_characteristic
  sens = pos_pos/pos_all  # also: tpr, power, hit rate, recall (anisimova 2002 as power)
  spec = neg_neg/neg_all  # also: tnr
  ppv = pos_pos/all_pos   # also: precision, (anisimova 2002 has it as accuracy)
  npv = neg_neg/all_neg
  fdr = neg_pos/all_pos
  fpr = neg_pos/neg_all   # also: false alarm rate, fallout
  dlr = sens/(1-spec)     # diagnostic likelihood ratio, see http://www.rapid-diagnostics.org/accuracy.htm#posdlr

  pos_true = pos_all/all # proportion of true positives
  pos_inf = all_pos/all # proportion of inferred positives

  return(list(
    sens=sens,
    tpr=sens,
    spec=spec,
#    ppv=ppv,
#    npv=npv,
    fdr=fdr,
    fpr=fpr,
#    dlr=dlr,
    pos_true=pos_true,
    pos_aln=pos_inf,
    pos_inf=pos_inf,
    true_pos_count = pos_all,
    true_neg_count = neg_all,
    aln_pos_count = all_pos,
    aln_neg_count = all_neg,
    tp = pos_pos,
    tn = neg_neg,
    fp = neg_pos,
    fn = pos_neg
  ))
}

format.numeric.df <- function(x,digits=3) {  
  nums <- unlist(lapply(x,is.double))
  #print(nums)
  for (col in names(x)) {
    if (nums[col] == TRUE) {
      x[,col] <- formatC(x[,col],digits=digits,format='fg')
    }
  }

  return(x)
}

df.gene.detection.rate = function(df) {
  # Note: the detection rate is based on PAML's LRT test or SLR's multiple testing-corrected estimates.
  # (This means that a different SLR threshold will *not* change these calculated detection rates.)
  df.detected.positive = function(df) {
    if (is.paml(df)) {
      return(nrow(subset(df, paml_lrt > 5.99)) >= 1)
    } else {
      return(nrow(subset(df, aln_type == 'positive3' | aln_type == 'positive4')) >= 1)
    }
  }

  num.detected = by(df,df$node_id,df.detected.positive)
  return(list(detected=sum(num.detected)/length(num.detected)))
}

df.cor = function(df) {
  return(list(cor=cor(df$true_dnds,df$aln_dnds,method='spearman',use='complete.obs')))
}

df.alignment.accuracy = function(df) {
  first.r <- function(x,field) {
    return(x[1, field])
  }                      
  mean.f <- function(field) {
    return(mean(by(df, df$node_id, first.r, field)))
  }
  sum.f <- function(field) {
    return(sum(by(df, df$node_id, first.r, field)))
  }

  out.df <- data.frame(
    sum_of_pairs_score=mean.f('sum_of_pairs_score'),
    total_column_score=mean.f('total_column_score'),
    mean_bl_aligned=mean.f('mean_bl_aligned'),
    mean_bl_match=mean.f('mean_bl_match'),
    mean_bl_mismatch=mean.f('mean_bl_mismatch'),
    match_bl_score=mean.f('match_bl_score'),
    mismatch_bl_score=mean.f('mismatch_bl_score'),
    entropy=mean.f('entropy'),
    lambda=mean.f('lambda'),
    aln_length=mean.f('aln_length'),
    filtered_fraction=1 - (sum.f('residue_count') / sum.f('unmasked_residue_count'))
  )
  return(out.df)
}

meta.sub <- function(asdf) {
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
