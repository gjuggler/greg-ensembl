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
  query = sprintf("SELECT * FROM stats_sites")
  sites = get.vector(con,query,columns='all')
  
  query = sprintf("SELECT * FROM stats_genes")
  genes = get.vector(con,query,columns='all')

  all = merge(sites,genes,by=c('data_id','node_id','parameter_set_id'))
  
  return(all)
}

get.all.nodes = function() {
  query = sprintf("SELECT * FROM stats_sites GROUP BY node_id")
  data = get.vector(con,query,columns='all')
  return(data)
}

get.all.first.reps = function() {
  nodes = get.all.nodes()
  first.reps = subset(nodes,slrsim_rep == 1);
  return(first.reps)
}

dump.protein = function(node_id,base=node_id,tree_file=NA,aln_file=NA,sw_file=NA,params_file=NA) {
  if (is.na(tree_file)) {
    tree_file = paste(base,".nh",sep="")
  }
  if (is.na(aln_file)) {
    aln_file = paste(base,".fa",sep="")
  }
  if (is.na(sw_file)) {
    sw_file = paste(base,".csv",sep="")
  }
  if (is.na(params_file)) {
    params_file = paste(base,".txt",sep="")
  }

  url = "mysql://slrsim:slrsim@mysql-greg.ebi.ac.uk:4134/slrsim_anisimova"
  system(paste("perl ~/lib/greg-ensembl/scripts/tree_dump.pl",
    " --url=",url,
    " --id=",node_id,
    " --tree=",tree_file,
    " --aln=",aln_file,
    " --sw=",sw_file,
    " --sw_table=","stats_slrsim",
    " --params=",params_file,
    sep="")
  )  
}

dump.nodes = function(nodes,base) {
  for (i in 1:nrow(nodes)) {
    node = nodes[i,]
    dump.protein(node$node_id,base=paste(base,node$node_id,sep=""))
  }
}

dump.all.nodes = function(base) {
  nodes = get.all.nodes()
  dump.nodes(nodes,base)
}

get.test.data = function() {
  query = sprintf("SELECT * FROM stats_slrsim LIMIT 500000");
  data = get.vector(con,query,columns='all')
  return(data)
}

is.paml = function(df) {
  if (grepl("paml",df[1,]$sitewise_action,ignore.case=T)) {
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

df.swfwer = function(df,
    thresh=3.8,
    paml_thresh=0.95) {

  aln_thresh = 1
  if (is.paml(df)) {
    thresh = paml_thresh
    aln_thresh = -1
  }

  df.fp = function(df,thresh) {
    fps = nrow(subset(df,true_type!="positive1" & aln_lrt>thresh & aln_dnds>aln_thresh))
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
  return(list(cor=cor(df$true_dnds,df$aln_dnds,method='spearman',use='complete.obs')))
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

  df.unfiltered = function(df){return(df[1,]$unfiltered_site_fraction)}
  unfiltered.by.gene = by(df,df$node_id,df.unfiltered)

  return(list(
    sum_of_pairs=mean(sps.by.gene),
    total_column_score=mean(tcs.by.gene),
    aln_entropy=mean(df$aln_entropy),
    true_entropy=mean(df$true_entropy),
    unfiltered_fraction=mean(unfiltered.by.gene)
  ))
}

# Go through each experiment, calculate summary stats, and build a data frame.
summarize.results = function(data,thresh=3.8,paml_thresh=0.95) {

  # Paste together some metadata so that we have one ID per experiment.
  attrs = c('slrsim_file','alignment_name','filtering_name','alignment_score_threshold','slrsim_ref','sitewise_action','phylosim_insertrate','slrsim_tree_length')

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

  a = get.perf(df)
  plot(a,col=color,lty=lty,lwd=lwd,add=overlay,colorize=F)
}

get.perf = function(df,...) {
  if (!is.paml(df)) {
    #print("Signed LRT...")
    # Create a signed LRT if it's SLR-based data.
    df$aln_lrt = sign(df$aln_dnds-1)*df$aln_lrt
  }
  truth = as.integer( df$true_dnds > 1 )
  
  require(ROCR)
  pred = prediction(df$aln_lrt,truth)
  perf = performance(pred,measure="tpr",x.measure="fpr")
  return(perf)
}

get.perf.df = function(label="unlabeled",df) {
  perf = get.perf(df)
  df = data.frame(label=label,fpr=unlist(perf@x.values),tpr=unlist(perf@y.values))
  return(df)
}

empty.plot = function() {
  plot(c(),xlim=c(0,0.1),ylim=c(0,1),xlab='',ylab='')
}

filtering.roc = function(data) {
  require(ggplot2)
  png("~/public_html/blank.png")

  comb = data.frame()
#  for (aln in c('mcoffee','prank_f','True Alignment')) {
#    for (filt in c('None','tcoffee')) {
  for (aln in c('prank_f','True Alignment')) {
    for (filt in c('None','prank_mean','prank_treewise','indelign','tcoffee')) {
      if (aln == 'True Alignment' && filt != 'None') {
        next;
      }
      sub = subset(data, alignment_name==aln & filtering_name==filt)
      sub.roc = slr.roc(sub)
      sub.roc$filter = as.factor(filt)
      sub.roc$aln = as.factor(aln)
      sub.roc$label = paste(aln,filt,sep="/")
      comb = rbind(sub.roc,comb)
    }
  }
  print(nrow(comb))
  p <- ggplot(comb,aes(x=fp,y=tp,group=label,colour=label)) + geom_line()
  ggsave("~/public_html/slrsim_filt.png",width=8,height=6)
  dev.off()
}

filter.sweep.roc = function(data) {
  comb = data.frame()
  for (aln in sort(unique(data$alignment_name))) {
    for (filt in sort(unique(data$filtering_name))) {
      for (thresh in sort(unique(data$alignment_score_threshold))) {
#        if (aln == 'True Alignment' && (filt != 'None' || thresh != 5)) {
#          next;
#        }

        sub = subset(data, alignment_name==aln & filtering_name==filt & alignment_score_threshold==thresh)
        print(paste(aln,filt,thresh,nrow(sub),sep="/"))
        if (nrow(sub)==0) {next}
        sub.roc = slr.roc(sub)
        sub.roc$filter = as.factor(filt)
        sub.roc$aln = as.factor(aln)
        sub.roc$thresh = as.factor(thresh)
        sub.roc$label = paste(aln,filt,thresh,sep="/")
        comb = rbind(sub.roc,comb)      
      }
    }
  }
  return(comb)
}

test.plot.roc = function(data) {
  png(file="~/public_html/roc.png",width=600,height=600)
  plot.roc(data)
  dev.off()
}

plot.roc = function(data) {
  require(ggplot2)
  require(grid)
  print(nrow(data))

  # See http://learnr.wordpress.com/2009/05/26/ggplot2-two-or-more-plots-sharing-the-same-legend/
  lay = grid.layout(nrow=2,ncol=2,
    widths=unit(c(2,0.5),c("null","null")),
    heights=unit(c(1,1),c("null","null")))
  vplayout = function(...) {grid.newpage();pushViewport(viewport(layout=lay))}
  subplot = function(x,y){viewport(layout.pos.row=x,layout.pos.col=y)}
    
  vplayout()

  p <- ggplot(data,aes(x=tn,y=tp,group=label,colour=label)) + geom_line()
  no.leg = opts(legend.position="none")

  p1 <- p + no.leg
  p2 <- p + xlim(0,1000) + no.leg
  leg <- p + opts(keep="legend_box")

  print(p1,vp=subplot(1,1))
  print(p2,vp=subplot(2,1))
  print(leg,vp=subplot(1:2,2))
}

slr.roc = function(df) {
  library(doBy)
  library(plyr)

  if (!is.paml(df)) {
    # Create a signed LRT if it's SLR-based data.
    df$score = sign(df$aln_dnds-1)*df$aln_lrt
  } else {
    df$score = df$aln_lrt
  }
  df$truth = as.integer( df$true_dnds > 1 )
  df <- orderBy(~-score,data=df)
  #print(df[1:10,])
  df$tp = cumsum(df$truth)
  df$tn = cumsum(1-df$truth)
  
  return(df)
}

my.own.roc = function(data, score.field='score', truth.fun=NULL) {
  library(doBy)
  library(plyr)

  #png("~/public_html/test.png")
  #n <- 2000
  #truth <- rbinom(n,1,0.8)
  #score <- rnorm(n)

  truth <- adply(df, 1, truth.fun);
  score <- data[[score.field]]

  data$truth = truth
  data$score = score
  data <- orderBy(~+score,data=data)

  data$tp = cumsum(data$truth)
  data$fp = cumsum(1-data$truth)

  return(data)
}


gg.df.cor = function(df) {
  return(cor(df$true_dnds,df$aln_dnds,method='spearman',use='complete.obs'))
}

gg.df.fpr = function(df) {
  df.res = df.stats(df)
  return(df.res[['fpr']])
}

correlation.boxplot = function(data,variable,n=5,title="") {
  png("~/public_html/blank.png")  

  dx = (max(data[[variable]]) - min(data[[variable]])) / (n-1)
  data$var_orig = data[[variable]]
  data$var = round_any(data[[variable]],dx)

  cor = ddply(data,.(var),'gg.df.cor')

  p <- qplot(x=var,y=gg.df.cor,data=cor) + xlab(variable) + ylab("Spearman dN/dS Correlation") + ylim(0,1)
  gr <- geom_rug(data=data,aes(x=var_orig,alpha=0.05,size=0.1))
  p + gr + opts(legend.position="none") + opts(title=title)

  #p + geom_histogram(data=data,aes(x=var, y=..density..))
  #p + stat_summary(fun.data=fun,colour="red",geom="point",size=3)
  #p + geom_boxplot() + xlab(variable) + ylab("Error")
  ggsave("~/public_html/slrsim_cor.png",width=8,height=6)
  dev.off()
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