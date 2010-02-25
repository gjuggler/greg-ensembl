if (exists('drv')) {
  lapply(dbListConnections(drv),dbDisconnect)
} else {
  library(RMySQL)
  drv <- dbDriver('MySQL')
}
con <- dbConnect(drv, host='ens-research', port=3306, user='ensro', password='', dbname='gj1_eslr_57')

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
get.psets = function() {
  query = paste('SELECT n.parameter_set_id AS id,sn.parameter_value AS shortname, n.parameter_value AS name FROM parameter_set n, parameter_set sn',
    'WHERE n.parameter_set_id=sn.parameter_set_id AND n.parameter_name="name" AND sn.parameter_name="shortname";')
  return(get.vector(con,query,columns='all'))
}

get.genes = function(parameter.set.id=1) {
  query = sprintf("SELECT * FROM stats_genes where parameter_set_id=%s",parameter.set.id);
  data = get.vector(con,query,columns='all')

  return(data)
}

get.all.merged = function() {
  all.node.ids = get.vector(con,'SELECT DISTINCT(node_id) FROM stats_genes')
  
  param.sets = get.psets()
  for (pset in param.sets$id) {
    genes = get.genes(pset)
    col.dnds = paste(param.sets[pset,]$shortname,'_dnds',sep="")
    col.psc = paste(param.sets[pset,]$shortname,'_psc',sep="")
    col.psc.weak = paste(param.sets[pset,]$shortname,'_psc_weak',sep="")
    genes.subset = data.frame(a=genes$node_id,b=genes$omega_mean,c=genes$psc_count,d=genes$weak_psc_count)
    colnames(genes.subset) = c('node_id',col.dnds,col.psc,col.psc.weak)
    
    if (!exists('merged.df')) {
      merged.df = data.frame(node_id=all.node.ids)
    }
    print(merged.df[1,])
    merged.df = merge(merged.df,genes.subset,all.x=T)
  }
  return(merged.df);
}
