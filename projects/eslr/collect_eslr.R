if (exists('drv')) {
  lapply(dbListConnections(drv),dbDisconnect)
} else {
  library(RMySQL)
  drv <- dbDriver('MySQL')
}

con <- dbConnect(drv, host='ens-research', port=3306, user='ensro', password='', dbname='gj1_eslr')

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

# Get the parameter sets.
query = 'SELECT parameter_set_id AS id,parameter_value AS name FROM parameter_set where parameter_name="name";'
param.sets = get.vector(con,query,columns='all')

get.genes = function(parameter.set.id=1) {
  query = sprintf("SELECT * FROM stats_genes where parameter_set_id=%s",parameter.set.id);
  data = get.vector(con,query,columns='all')

  # Temoporary fix.
  psc = data$num_pscs
  psc_weak = data$num_pscs_weak
  data$num_pscs = psc_weak
  data$num_pscs_weak = psc

  return(data)
}

if (!exists('all.genes')) {
  all.genes = get.genes(1)
}
