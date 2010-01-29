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
#print(param.sets)

get.data.alt = function() {
  query = sprintf("SELECT * FROM stats_genes where parameter_set_id=1");
  data = get.vector(con,query,columns='all')
  return(data)
}

if (!exists('all.data')) {
  all.data = get.data.alt()
}



ize = function(izer) {
  function(data, columns=names(data)) {
    data[columns] = lapply(data[columns], izer)
    data
  }
}
logicalize = ize(as.logical)
characterize = ize(as.character)
factorize = ize(as.factor)
