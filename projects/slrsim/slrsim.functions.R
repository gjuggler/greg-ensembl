library(plyr)
library(doBy)

paper.table <- function(df) {

  # Sum up the total and unfiltered site counts.
  one.index.per.rep <- which(!duplicated(df$slrsim_rep))
  one.row.per.rep <- df[one.index.per.rep,]
  all.residues <- sum(one.row.per.rep$site_count)
  left.residues <- sum(one.row.per.rep$unfiltered_site_count)   

  row <- df[1,]

  # Use the technique where the p-value can be doubled because we're just looking
  # at positive selection. This gives thresholds of 2.71 and 5.41 for 5% and 1% FPR.
  stats.05 <- df.stats(df,thresh=2.71)
  stats.01 <- df.stats(df,thresh=5.41)
  
  stats.05 <- data.frame(stats.05)
  stats.01 <- data.frame(stats.01)

  cor = cor(df$true_dnds,df$aln_dnds,method='spearman',use='complete.obs')

  roc <- slr.roc(df,na.rm=TRUE)
  n.sites <- nrow(roc)
  rrow <- roc[1,]

  indel = row$phylosim_insertrate * 2

  ret.df <- data.frame(
    tree     = row$slrsim_tree_file,
    aligner  = row$alignment_name,
    length   = row$tree_mean_path,
    indel    = indel,
    filter   = row$filtering_name,

    auc      = rrow$auc,
    `fpr_tpr`    = rrow$tpr_at_fpr,
    `fpr_fdr`    = rrow$tpr_at_fpr,
    `fpr_tp`    = rrow$tp_at_fpr,
    `fpr_fp`    = rrow$fp_at_fpr,

    `fdr_tpr`    = rrow$tpr_at_fdr,
    `fdr_fpr`    = rrow$fpr_at_fdr,
    `fdr_tp`    = rrow$tp_at_fdr,
    `fdr_fp`    = rrow$fp_at_fdr,

    `cor`   = cor
  )

  ret.df <- data.frame(ret.df,
    sites = all.residues,
    filtered.fraction = 1 - (left.residues/all.residues)
  )

  ret.df <- ret.df[with(ret.df, order(tree,aligner,length,indel,filter)),]
  return(ret.df)
}

fig.1.summary <- function(df,thresh) {
  
  stats <- df.stats(df,thresh=thresh)
  roc <- slr.roc(df,na.rm=TRUE)

  # Need to copy these ROC-derived summary stats over to the 'stats' object.
  stats$auc <- roc[1,]$auc
  stats$tpr_at_fpr <- roc[1,]$tpr_at_fpr

  print(paste("summarizing",df[1,]$experiment_name,df[1,]$slrsim_label))
  return(cbind(df[1,],stats))
}

min.fdr.score <- function(data) {
    roc = slr.roc(data,na.rm=F)
    roc = subset(roc, tp > 0)
		min.fdr = min(roc$fdr)
    which.min = which(roc$fdr <= min.fdr)
    best.index = max(which.min)
    min.score = roc[best.index,]$score
#    print(paste(roc[1,]$slrsim_label,min.fdr,best.index,min.score))
    return(data.frame(min.score=min.score,min.fdr=min.fdr))
}

slr.roc = function(df,na.rm=F) {
  library(doBy)
  library(plyr)

  if (na.rm) {
    n.na = nrow(subset(df,is.na(aln_lrt)))
    df = subset(df,!is.na(aln_lrt))
    n.left <- nrow(df)
    #print(sprintf("Removed %d NA rows (%d remaining)",n.na,n.left))
    if(n.left == 0) {
      # Think of something to do here. Return empty df?
    }
  }

  if (!is.paml(df)) {
    # Create a signed LRT if it's SLR-based data.
    # We'll replace NAs with an extremely low score.
    df[is.na(df$aln_lrt),'aln_lrt'] = 1000
    df[is.na(df$aln_dnds),'aln_dnds'] = 0
    df$score = sign(df$aln_dnds-1)*df$aln_lrt
  } else {
    df$score = df$aln_lrt
  }
  df$truth = as.integer( df$true_dnds > 1 )
  df <- orderBy(~-score,data=df)

  df$tp = cumsum(df$truth)
  df$tn = cumsum(1-df$truth)
  df$count = cumsum(rep(1,nrow(df)))
  df$fp = df$tn

  df$fpr = df$tn / max(df$tn)
  df$tpr = df$tp / max(df$tp)

  df$fdr = df$fp/(df$count)

  if (!na.rm) {
    df <- subset(df,aln_lrt != 1000)
  }

  df$auc_full <- -1
  df$auc <- -1
  if (nrow(df) > 1) {
    auc <- area.under.curve(df,x.lim=0.2)
    df$auc <- auc
    auc_full <- area.under.curve(df,x.lim=1)
    df$auc_full <- auc_full
  }

  df$tpr_at_fpr <- -1
  df$fdr_at_fpr <- -1
  df$tp_at_fpr <- -1
  df$fp_at_fpr <- -1
  if (nrow(df) > 1) {
    fpr.sub <- subset(df,fpr < 0.05)
    if (nrow(fpr.sub) > 0) {
      # Take the last (rightmost) row with fpr < 0.05
      row <- fpr.sub[nrow(fpr.sub),]
      df$tpr_at_fpr <- row$tpr
      df$fdr_at_fpr <- row$fdr
      df$tp_at_fpr <- row$tp
      df$fp_at_fpr <- row$fp
    }
  }

  df$tpr_at_fdr <- -1
  df$fpr_at_fdr <- -1
  df$tp_at_fdr <- -1
  df$fp_at_fdr <- -1
  if (nrow(df) > 1) {
    fpr.sub <- subset(df,fdr < 0.1)
    if (nrow(fpr.sub) > 0) {
      # Take the last (rightmost) row with fdr < 0.1
      row <- fpr.sub[nrow(fpr.sub),]
      df$tpr_at_fdr <- row$tpr
      df$fpr_at_fdr <- row$fpr
      df$tp_at_fdr <- row$tp
      df$fp_at_fdr <- row$fp
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
  a = ddply(data,.(slrsim_label),fn,thresh=thresh)
  return(a)
}

do.by.labels = function(data,fn) {
  library(plyr)
  d_ply(data,.(slrsim_label),fn)
}

is.paml = function(df) {
  if (!any(colnames(df) %in% c('sitewise_action'))) {
    return(FALSE)
  }
  if (grepl("paml",df[1,'sitewise_action'],ignore.case=T)) {
    return(TRUE)
  } else {
    return(FALSE)
  } 
}

df.stats = function(df,
  thresh=3.8,
  paml_thresh=0.95,type='all') {

  aln_thresh = 1
  if (is.paml(df)) {
    thresh = paml_thresh
    aln_thresh = -1
  }

  # Collect stats for SLR-type runs.
  pos_pos = nrow(subset(df,true_type=="positive1" & aln_lrt>thresh & aln_dnds>aln_thresh))
  neg_pos = nrow(subset(df,true_type!="positive1" & aln_lrt>thresh & aln_dnds>aln_thresh))
  neg_neg = nrow(subset(df,true_type!="positive1" & !(aln_lrt>thresh & aln_dnds>aln_thresh)))
  pos_neg = nrow(subset(df,true_type=="positive1" & !(aln_lrt>thresh & aln_dnds>aln_thresh)))

  pos_all = nrow(subset(df,true_type=="positive1"))
  neg_all = nrow(subset(df,true_type!="positive1"))
  all_pos = nrow(subset(df,aln_lrt>thresh & aln_dnds>aln_thresh))
  all_neg = nrow(subset(df,!(aln_lrt>thresh & aln_dnds>aln_thresh)))

  all = nrow(df)

  if (nrow(subset(df,!is.na(true_dnds) & !is.na(aln_dnds))) > 0) {
    cor = cor(df$true_dnds,df$aln_dnds,method='spearman',use='complete.obs')
  } else {
    cor = 0
  }

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
    cor = cor
  ))
}
