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
  'SELECT distinct value FROM protein_tree_tag where tag="sim_set"'
))

# Get the parameter sets.
res = dbSendQuery(con,'SELECT parameter_set_id AS id,parameter_value AS name FROM parameter_set where parameter_name="name";')
param_sets = as.vector(fetch(res))
dbClearResult(res)

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

all_data = data.frame();
for (parameter_set_id in c(1,2)) {
  for (sim_set in sim_sets) {
    node_ids = get_vector(con,sprintf('SELECT distinct node_id FROM protein_tree_tag WHERE tag="sim_set" AND value="%s";',sim_set))
    for (node_id in node_ids) {
    print(sprintf("%s %s %s",parameter_set_id,sim_set,node_id))
      # Grab the site-wise information into a data frame.
      data = get_data(node_id,parameter_set_id)

      if (!is.data.frame(data) || nrow(data)==0)
        next

      # Add some useful info to the data frame.
      data$sim_set = rep(sim_set,n=nrow(data))
      data$parameter_set = rep(parameter_set_id,n=nrow(data))

      # Concatenate to the end of the data frame.
      all_data = rbind(all_data,data)
    } 
  }
}
dbDisconnect(con)

df.accuracy = function(df) {
  pos_vec = c('positive1','positive2','positive3','positive4')
  pos_pos = nrow(subset(df,true_type=="positive1" & aln_type %in% pos_vec))
  all_pos = nrow(subset(df,aln_type %in% pos_vec))
  all_pos = max(all_pos,1)
  pos_pos/all_pos
}

