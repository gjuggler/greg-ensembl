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

# Get the sim_sets.
sim_sets = get_vector(con,paste(
  'SELECT distinct value FROM protein_tree_tag where tag="sim_name"'
))

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
  all_data$program = rep('slr',nrow(all_data))
  all_data[all_data$parameter_set==3,]$program = 'paml'
  all_data$program = as.factor(all_data$program)
  all_data[is.na(all_data$ins_rate),]$ins_rate = 0
  all_data[is.na(all_data$del_rate),]$del_rate = 0
}

df.stats = function(df,thresh=3.8,paml_thresh=0.95,type='all') {

  aln_thresh = 1
  if (df$program[1] == 'paml') {
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
  if (df$program[1] == 'paml') {
    thresh = paml_thresh
    aln_thresh = -1
  } 

  df.fp = function(df,thresh) {
    fps = nrow(subset(df,true_type!="positive1" & lrt>thresh & aln>aln_thresh))
    if(fps>0)
      return(1)
    return(0)
  }

  error_by_protein = by(df,df$node_id,df.fp,thresh=thresh)
  return(list(fwer=sum(error_by_protein)/length(error_by_protein)))
}

df.cor = function(df) {
  return(list(cor=cor(df$true,df$aln,method='spearman')))
}

# Go through each experiment, calculate summary stats, and build a data frame.
summarize.results = function(data,thresh=3.8,paml_thresh=0.95) {

  # Paste together some metadata so that we have one ID per experiment.
  attrs = c('sim_name','program','ins_rate','sim_length')
  ids = paste(data$sim_name,data$parameter_set,data$ins_rate,data$sim_length)
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

slr.roc = function(df) {
  require(RColorBrewer)       

  lengths = unique(df$sim_length)
  indels = unique(df$ins_rate)
  psets = unique(df$parameter_set)

  length.col = brewer.pal(length(lengths),"Greens")
  indel.col = brewer.pal(length(indels),"Reds")
  pset.lty = c('solid','dashed')
  plot(c(),xlim=c(0,0.1),ylim=c(0,1))

  # Hold tree size const, vary indel rate.
  for (i in 1:length(psets)) {
    pset = psets[i]
    for (j in 1:length(indels)) {
      x = indels[i]
      slr = subset(df,parameter_set==pset & sim_name=='art' & ins_rate==x & sim_length==1.1)
      if (pset == 2) {
        slr$lrt = sign(slr$aln-1)*slr$lrt
      }
      a = get.perf(slr)
      plot(a,col=indel.col[i],lty=pset.lty[i],add=T)
    }
  }


  # Hold indel rate const, vary tree size.
  for (i in 1:length(lengths)) {
    x = lengths[i]
    slr = subset(df,parameter_set==2 & sim_name=='art' & ins_rate==0 & sim_length==x)
    slr$lrt = sign(slr$aln-1)*slr$lrt
    a = get.perf(slr)
    #plot(a,col=length.col[i],add=T)
  }

}

get.perf = function(df,...) {
  truth = as.integer( df$true > 1 )
  
  require(ROCR)
  pred = prediction(df$lrt,truth)
  perf = performance(pred,measure="tpr",x.measure="fpr")
  return(perf)
}