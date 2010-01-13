library(RMySQL)
drv <- dbDriver('MySQL')
con <- dbConnect(drv, host='mysql-greg.ebi.ac.uk', port=4134, user='slrsim', password='slrsim', dbname='slrsim_anisimova')

get_vector = function(con,query) {
  res = dbSendQuery(con,query)
  data = fetch(res)[[1]]
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

# Executes a Perl script which extracts the true and inferred omegas for a reference sequence.
get_data = function(node_id,param_set)  {
  cmd = sprintf("perl ./slrsim_stats_for_node.pl %s %s",node_id,param_set)
  tryCatch(
    {
      tbl = read.table(pipe(cmd),header=T,sep="\t")
      return(tbl)
    },
    error=function(err){return()}
  )
}

if (!exists('all_data')) {
  all_data = data.frame();
  for (parameter_set_id in c(2,3)) {
    node_ids = get_vector(con,'SELECT distinct node_id FROM protein_tree_tag WHERE tag="sim_name"')
    for (node_id in node_ids) {
      print(sprintf("%s %s",parameter_set_id,node_id))

      # Grab the site-wise information into a data frame.
      data = get_data(node_id,parameter_set_id)

      if (!is.data.frame(data) || nrow(data)==0)
        next

      # Add some useful info to the data frame.
      n = nrow(data)

      # Concatenate to the end of the data frame.
      all_data = rbind(all_data,data)
    } 
  }
  dbDisconnect(con)
}

df.stats = function(df,thresh=6,type='all') {
  pos_pos = nrow(subset(df,true_type=="positive1" & lrt>thresh & aln>1))
  neg_pos = nrow(subset(df,true_type!="positive1" & lrt>thresh & aln>1))
  neg_neg = nrow(subset(df,true_type!="positive1" & !(lrt>thresh & aln>1)))
  pos_neg = nrow(subset(df,true_type=="positive1" & !(lrt>thresh & aln>1)))

  pos_all = nrow(subset(df,true_type=="positive1"))
  neg_all = nrow(subset(df,true_type!="positive1"))
  all_pos = nrow(subset(df,lrt>thresh & aln>1))
  all_neg = nrow(subset(df,!(lrt>thresh & aln>1)))
  all = nrow(df)

  # see http://en.wikipedia.org/wiki/Receiver_operating_characteristic
  sens = pos_pos/pos_all  # also: tpr, power, hit rate, recall
  spec = neg_neg/neg_all  # also: tnr
  ppv = pos_pos/all_pos   # also: precision
  npv = neg_neg/all_neg
  fdr = neg_pos/all_pos
  fpr = neg_pos/neg_all   # also: false alarm rate, fallout
  dlr = sens/(1-spec)     # diagnostic likelihood ratio, see http://www.rapid-diagnostics.org/accuracy.htm#posdlr

  pr = pos_all/all # just the proportion of positives

  if(all_pos==0) {
    ppv = NA
    fdr = NA
  }
  if (pos_all==0) {
    sens = NA
  } 

  return(list(sens=sens,spec=spec,ppv=ppv,npv=npv,fdr=fdr,fpr=fpr,pr=pr))
}

df.fwer = function(df,thresh=6) {
  df.fp = function(df,thresh) {
    fps = nrow(subset(df,true_type!="positive1" & lrt>thresh & aln>1))
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
summarize.results = function(data) {
  ids = paste(data$sim_name,data$sim_length,data$parameter_set)

  for (my_id in unique(ids)) {
    df = data[ids==my_id,]

    # These values should be the same for each subset.
    sim_name = df$sim_name[1]
    sim_length = df$sim_length[1]
    parameter_set = df$parameter_set[1]
    sim_file = df$sim_file[1]
    sim_ref = df$sim_ref[1]
    
    # Calculate the summaries.
    stats = df.stats(df,thresh=3.8)
    fwer = df.fwer(df,thresh=3.8)
    cor = df.cor(df)

    # Put them all together into a data frame.
    attrs = list(sim_name=sim_name,sim_length=sim_length,parameter_set=parameter_set,sim_file=sim_file,sim_ref=sim_ref)
    res.list = c(stats,fwer,cor,attrs)
    if(!exists("res.df")) {
      res.df = data.frame(attrs,stats,fwer,cor)
    }
    res.df = rbind(res.df,data.frame(attrs,stats,fwer,cor))
  } 
  return(res.df)
}