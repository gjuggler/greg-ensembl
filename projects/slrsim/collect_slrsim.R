library(RMySQL)
drv <- dbDriver('MySQL')
con <- dbConnect(drv, host='mysql-greg.ebi.ac.uk', port=4134, user='greg', password='TMOqp3now', dbname='gj1_slrsim_1')

get_vector = function(con,query) {
  res = dbSendQuery(con,query)
  data = as.vector(fetch(res))
  dbClearResult(res)
  return(data)
}

# Get the sim_sets.
res = dbSendQuery(con,paste(
  'SELECT distinct value FROM protein_tree_tag where tag="sim_set"'
))
sim_sets = as.vector(fetch(res))
dbClearResult(res)

# Get the parameter sets.
res = dbSendQuery(con,'SELECT parameter_set_id AS id,parameter_value AS name FROM parameter_set where parameter_name="name";')
param_sets = as.vector(fetch(res))
dbClearResult(res)

get_data = function(node_id,param_set)  {
  cmd = sprintf("perl ./slrsim_stats_for_node.pl %s %s",node_id,param_set)
  tbl = read.table(pipe(cmd),header=T,sep="\t")
  return(tbl)
}

for (sim_set in sim_sets) {
  sim_results = 
    
  node_ids = get_vector(con,sprintf('SELECT distinct node_id FROM protein_tree_tag WHERE tag="sim_set" AND value="%s";',sim_set))

  for (node_id in node_ids) {
    data = get_data(node_id,2)
    
  } 
  
}

dbDisconnect(con)


