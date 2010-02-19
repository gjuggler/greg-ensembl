if (exists('drv')) {
  lapply(dbListConnections(drv),dbDisconnect)
} else {
  library(RMySQL)
  drv <- dbDriver('MySQL')
}

con <- dbConnect(drv, host='mysql-greg.ebi.ac.uk', port=4134, user='slrsim', password='slrsim', dbname='slrsim_anisimova')

get.vector = function(con,query,columns=1) {
  res = dbSendQuery(con,query)
  if (columns == 'all') {
    data = fetch(res,n=-1)  # Important: the n=-1 causes the entire resultset to be fetched, not just the first 500.
  } else {
    data = fetch(res,n=-1)[[columns]]
  }
  dbClearResult(res)
  return(data)
}

# Grabs from the database the true and inferred omegas for the given reference sequence.
get.all.data = function() {
  query = sprintf("SELECT * FROM stats_slrsim");
  data = get.vector(con,query,columns='all')
  return(data)
}

get.test.data = function() {
  query = sprintf("SELECT * FROM stats_slrsim LIMIT 500000");
  data = get.vector(con,query,columns='all')
  return(data)
}

is.paml = function(df) {
  if (grepl("paml",df[1,]$sitewise_name,ignore.case=T)) {
    return(TRUE)
  } else {
    return(FALSE)
  } 
}

df.stats = function(df,thresh=3.8,paml_thresh=0.95,type='all') {

  aln_thresh = 1
  if (is.paml(df)) {
    thresh = paml_thresh
    aln_thresh = -1
  }

    # Collect stats for SLR-type runs.
    pos_pos = nrow(subset(df,true_type=="positive1" & lrt>thresh & aln_dnds>aln_thresh))
    neg_pos = nrow(subset(df,true_type!="positive1" & lrt>thresh & aln_dnds>aln_thresh))
    neg_neg = nrow(subset(df,true_type!="positive1" & !(lrt>thresh & aln_dnds>aln_thresh)))
    pos_neg = nrow(subset(df,true_type=="positive1" & !(lrt>thresh & aln_dnds>aln_thresh)))

    pos_all = nrow(subset(df,true_type=="positive1"))
    neg_all = nrow(subset(df,true_type!="positive1"))
    all_pos = nrow(subset(df,lrt>thresh & aln_dnds>aln_thresh))
    all_neg = nrow(subset(df,!(lrt>thresh & aln_dnds>aln_thresh)))

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
    spec=spec,
#    ppv=ppv,
#    npv=npv,
    fdr=fdr,
    fpr=fpr,
#    dlr=dlr,
    pos_true=pos_true,
    pos_inf=pos_inf
  ))
}

df.swfwer = function(df,thresh=3.8,paml_thresh=0.95) {
  aln_thresh = 1
  if (is.paml(df)) {
    thresh = paml_thresh
    aln_thresh = -1
  }

  df.fp = function(df,thresh) {
    fps = nrow(subset(df,true_type!="positive1" & lrt>thresh & aln_dnds>aln_thresh))
    return(fps)
  }

  error_by_protein = by(df,df$node_id,df.fp,thresh=thresh)
  max_errors = max(error_by_protein)
  min_errors = min(error_by_protein)
  return(list(swfwer=sum(error_by_protein >= 1)/length(error_by_protein),
              max=max_errors,min=min_errors))
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
  return(list(cor=cor(df$true_dnds,df$aln_dnds,method='spearman')))
}

df.alignment.accuracy = function(df) {
  df.sps = function(df) {
    return(df[1,]$sum_of_pairs_score)
  }
  df.tcs = function(df) {
    return(df[1,]$total_column_score)
  }
  sps.by.gene = by(df,df$node_id,df.sps)
  tcs.by.gene = by(df,df$node_id,df.tcs)

  return(list(
    sum_of_pairs=mean(sps.by.gene),
    total_column_score=mean(tcs.by.gene),
    aln_entropy=mean(df$aln_entropy),
    true_entropy=mean(df$true_entropy)
  ))
}

# Go through each experiment, calculate summary stats, and build a data frame.
summarize.results = function(data,thresh=2,paml_thresh=0.95) {

  # Paste together some metadata so that we have one ID per experiment.
  attrs = c('slrsim_file','alignment_name','filtering_name','species_name','sitewise_name','phylosim_ins_rate','slrsim_tree_length')

  ids = rep("",nrow(data))
  for (attr in attrs) {
    ids = paste(ids,data[[attr]],sep=" ")
  }
  data$id = ids

  split.sets = split(data,data$id)
  for (i in 1:length(split.sets)) {
    df = split.sets[[i]]
    print(df[1,]$id)
    #df = data[ids==my_id,]

    # Calculate the summaries.
    stats = df.stats(df,thresh=thresh,paml_thresh=paml_thresh)
    gene.detection.rate = df.gene.detection.rate(df)
    swfwer = df.swfwer(df,thresh=thresh,paml_thresh=paml_thresh)
    cor = df.cor(df)
    acc = df.alignment.accuracy(df)

    # Put them all together into a data frame.
    attr.list = list()
    for (attr in attrs) {
      attr.list[[attr]] = df[[attr]][1]
    }   

    if(!exists('res.df')) {
      res.df = data.frame(attr.list,stats,gene.detection.rate,swfwer,cor,acc)
    } else {
      res.df = rbind(res.df,data.frame(attr.list,stats,gene.detection.rate,swfwer,cor,acc))
    }
  } 
  return(res.df)
}

plot.roc = function(df,color='black',lty='solid',lwd=1,overlay=T) {
  require(RColorBrewer)       

  if (!overlay) {
    plot(c(),xlim=c(0,0.1),ylim=c(0,1))
  }

  if (!is.paml(df)) {
    # Create a signed LRT if it's SLR-based data.
    df$lrt = sign(df$aln_dnds-1)*df$lrt
  }
  a = get.perf(df)
  plot(a,col=color,lty=lty,lwd=lwd,add=overlay)
}

get.perf = function(df,...) {
  truth = as.integer( df$true_dnds > 1 )
  
  require(ROCR)
  pred = prediction(df$lrt,truth)
  perf = performance(pred,measure="tpr",x.measure="fpr")
  return(perf)
}

empty.plot = function() {
  plot(c(),xlim=c(0,0.1),ylim=c(0,1),xlab='',ylab='')
}

do.some.plots = function(data) {

  slr.lty = 'solid'
  paml.lty = 'dashed'
  leg.txt = c('SLR','PAML M3')
  leg.lty = c(slr.lty,paml.lty)
  # Massingham 2005 simulation A
  empty.plot()
  title(main="Massingham 2005 A",xlab="Type I Error",ylab="True Positives Rate")
  legend('bottomright',legend=leg.txt,lty=leg.lty)
  d = subset(data,sim_name=="mass_05_A" & parameter_set_name=="SLR")
  plot.roc(d,lty=slr.lty)
  d = subset(data,sim_name=="mass_05_A" & parameter_set_name=="PAML M3")
  plot.roc(d,lty=paml.lty)
  # Massingham 2005 simulation B
  empty.plot()
  title(main="Massingham 2005 B",xlab="Type I Error",ylab="True Positives Rate")
  legend('bottomright',legend=leg.txt,lty=leg.lty)
  d = subset(data,sim_name=="mass_05_B" & parameter_set_name=="SLR")
  plot.roc(d,lty=slr.lty)
  d = subset(data,sim_name=="mass_05_B" & parameter_set_name=="PAML M3")
  plot.roc(d,lty=paml.lty)


  # Xenopus vs Human bglobin simulations
  empty.plot()
  rate.cols = 'Set1'
  hum.lty = 'solid'
  xen.lty = 'dashed'
  leg.txt = c('Human','Xenopus')
  leg.lty = c(hum.lty,xen.lty)

  bglobin.data = subset(data,sim_name=='bglobin_ref_hum' | sim_name=='bglobin_ref_xen')
  bglobin.data = subset(bglobin.data,parameter_set_name=="SLR")
  rate.list = unique(bglobin.data$ins_rate)
  cols = brewer.pal(length(rate.list),rate.cols)
  for (i in 1:length(rate.list)) {
    rate = rate.list[i]
    data.sub = subset(bglobin.data,ins_rate == rate)

    d = subset(data.sub,sim_name=='bglobin_ref_hum')
    plot.roc(d,lty=hum.lty,col=cols[i])

    d = subset(data.sub,sim_name=='bglobin_ref_xen')
    plot.roc(d,lty=xen.lty,col=cols[i])
  }

  title(main='Anisimova-M3 with B-globain tree (n=17)',xlab='Type I Error',ylab='Sensitivity')
  legend('bottomright',title='Reference species',bg='white',
    legend=leg.txt,lty=leg.lty)
  legend('bottomright',title='Indel rate',inset=c(0,0.2),bg='white',
    legend=unlist(rate.list),lty='solid',col=cols)
  # Try plotting a tree.
  library(ape)
  library(Hmisc)
  plot.tree = function() {
    tree = read.tree('trees/anisimova_01_bglobin.nh')
    plot.phylo(tree,cex=0.3)
  }
  subplot(plot.tree,x=0.087,y=0.65,size=c(1.5,1.5))

}