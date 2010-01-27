if (exists('drv')) {
  lapply(dbListConnections(drv),dbDisconnect)
} else {
  library(RMySQL)
  drv <- dbDriver('MySQL')
}

con <- dbConnect(drv, host='mysql-greg.ebi.ac.uk', port=4134, user='slrsim', password='slrsim', dbname='slrsim_anisimova')

get_vector = function(con,query) {
  res = dbSendQuery(con,query)
  data = fetch(res,n=-1)[[1]]  # Important: the n=-1 causes the entire resultset to be fetched, not just the first 500.
  dbClearResult(res)
  return(data)
}

# Get the parameter sets.
query = 'SELECT parameter_set_id AS id,parameter_value AS name FROM parameter_set where parameter_name="name";'
param_sets = get_vector(con,query)

# Grabs from the database the true and inferred omegas for the given reference sequence.
get_data_alt = function() {
  header_q = "SELECT meta_value FROM meta WHERE meta_key='slrsim_stats_header'";
  header = get_vector(con,header_q)
  query = sprintf("SELECT stats FROM sitewise_stats");
  data = get_vector(con,query)
  data = c(header,data)
  toks = paste(data,collapse="\n")
  cat(toks,file="tmp.txt")
  tbl = read.table(file="tmp.txt",sep="\t",header=T)
  return(tbl)
}

if (!exists('all_data')) {
  all_data = get_data_alt()

  # Touch-ups.
  #all_data[is.na(all_data$ins_rate),]$ins_rate = 0
  #all_data[is.na(all_data$del_rate),]$del_rate = 0
}

is.paml = function(df) {
  if (grepl("paml",df[1,]$parameter_set_name,ignore.case=T)) {
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
    pos_pos = nrow(subset(df,true_type=="positive1" & lrt>thresh & aln>aln_thresh))
    neg_pos = nrow(subset(df,true_type!="positive1" & lrt>thresh & aln>aln_thresh))
    neg_neg = nrow(subset(df,true_type!="positive1" & !(lrt>thresh & aln>aln_thresh)))
    pos_neg = nrow(subset(df,true_type=="positive1" & !(lrt>thresh & aln>aln_thresh)))

    pos_all = nrow(subset(df,true_type=="positive1"))
    neg_all = nrow(subset(df,true_type!="positive1"))
    all_pos = nrow(subset(df,lrt>thresh & aln>aln_thresh))
    all_neg = nrow(subset(df,!(lrt>thresh & aln>aln_thresh)))

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

  if(all_pos==0) {
#    ppv = NA
#    fdr = NA
  }
  if (pos_all==0) {
#    sens = NA
  } 

  return(list(sens=sens,spec=spec,ppv=ppv,npv=npv,fdr=fdr,fpr=fpr,dlr=dlr,pos_true=pos_true,pos_inf=pos_inf))
}

df.fwer = function(df,thresh=3.8,paml_thresh=0.95) {
  aln_thresh = 1
  if (is.paml(df)) {
    thresh = paml_thresh
    aln_thresh = -1
  } 

  df.fp = function(df,thresh) {
    fps = nrow(subset(df,true_type!="positive1" & lrt>thresh & aln>aln_thresh))
    return(fps)
  }

  error_by_protein = by(df,df$node_id,df.fp,thresh=thresh)
  max_errors = max(error_by_protein)
  min_errors = min(error_by_protein)
  return(list(fwer=sum(error_by_protein >= 1)/length(error_by_protein),
              max=max_errors,min=min_errors))
}

df.cor = function(df) {
  return(list(cor=cor(df$true,df$aln,method='spearman')))
}


# Go through each experiment, calculate summary stats, and build a data frame.
summarize.results = function(data,thresh=3.8,paml_thresh=0.95) {

  # Paste together some metadata so that we have one ID per experiment.
  attrs = c('slrsim_scheme_name','alignment_name','filtering_name','species_name','sitewise_name','phylosim_ins_rate','slrsim_tree_length')

  ids = rep("",nrow(data))
  for (attr in attrs) {
    ids = paste(ids,data[[attr]],sep=" ")
  } 
  unique_ids = unique(ids)

  for (my_id in unique_ids) {
    print(my_id)
    df = data[ids==my_id,]

    # Calculate the summaries.
    stats = df.stats(df,thresh=thresh,paml_thresh=paml_thresh)
    fwer = df.fwer(df,thresh=thresh,paml_thresh=paml_thresh)
    cor = df.cor(df)

    # Put them all together into a data frame.
    attr.list = list()
    for (attr in attrs) {
      attr.list[[attr]] = df[[attr]][1]
    }   

    if(!exists('res.df')) {
      res.df = data.frame(attr.list,stats,fwer,cor)
    } else {
      res.df = rbind(res.df,data.frame(attr.list,stats,fwer,cor))
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
    df$lrt = sign(df$aln-1)*df$lrt
  }
  a = get.perf(df)
  plot(a,col=color,lty=lty,lwd=lwd,add=overlay)
}

get.perf = function(df,...) {
  truth = as.integer( df$true > 1 )
  
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