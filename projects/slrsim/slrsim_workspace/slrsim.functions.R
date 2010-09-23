fig.1.summary <- function(df,thresh) {
  
  #min.score = min.fdr.score(df)
  stats = df.stats(df,thresh=thresh)
  return(cbind(df[1,],stats))
}

fig.2.summary <- function(df,thresh) {
  min.score = min.fdr.score(df)
  stats = df.stats(df,thresh=thresh)
  return(cbind(df[1,],stats,min.score))
}

fig.3.summary <- function(df,thresh) {
  stats = df.stats(df,thresh=thresh)
  return(cbind(df[1,],stats))
}

fig.4.summary <- function(df,thresh) {
  stats = df.stats(df,thresh=thresh)
  return(cbind(df[1,],stats))
}

min.fdr.score <- function(data) {
    roc = slr.roc(data,na.rm=F)
    roc = subset(roc, tp > 0)
		min.fdr = min(roc$fdr)
    which.min = which(roc$fdr <= min.fdr)
    best.index = max(which.min)
    min.score = roc[best.index,]$score
    print(paste(roc[1,]$slrsim_label,min.fdr,best.index,min.score))
    return(data.frame(min.score=min.score,min.fdr=min.fdr))
}

slr.roc = function(df,na.rm=F) {
  library(doBy)
  library(plyr)

  if (na.rm) {
    n.na = nrow(subset(df,is.na(aln_lrt)))
    df = subset(df,!is.na(aln_lrt))
    print(sprintf("Removed %d NA rows",n.na))
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

	# Last row before fpr is 0.1
  last.row.index <- max(which(df$fpr <= 0.1))
  last.row <- df[last.row.index,]
  df$`tpr.at.0.1` <- last.row$tpr

  return(df)
}

summarize.by.labels = function(data,fn,thresh=3.8) {
  library(plyr)
  a = ddply(data,.(slrsim_label),fn,thresh=thresh)
  if (!is.null(a$`tpr.at.0.1`)) {
    #levels(a$slrsim_label) <- 
#    print("Ordering!")
#    print(unique(a$`tpr.at.0.1`))
#    a <- orderBy(~-`tpr.at.0.1`,data=a)
#    print(unique(a$`tpr.at.0.1`))
  }
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

  cor = cor(df$true_dnds,df$aln_dnds,method='spearman',use='complete.obs')

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

  true_pos_count = pos_pos
  true_neg_count = neg_neg


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
    true_pos_count = pos_all,
    true_neg_count = neg_all,
    aln_pos_count = all_pos,
    aln_neg_count = all_neg,
    cor = cor
  ))
}
