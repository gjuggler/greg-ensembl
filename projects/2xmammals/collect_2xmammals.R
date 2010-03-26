if (exists('drv')) {
  lapply(dbListConnections(drv),dbDisconnect)
} else {
  library(RMySQL)
  drv <- dbDriver('MySQL')
}
connect.db = function(dbName) {
  con <- dbConnect(drv, host='ens-research', port=3306, user='ensro', password='', dbname='gj1_eslr_57')
  return(con)
}
con = connect.db("gj1_2x_57")

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

get.genes = function(parameter.set.id=1,db="gj1_eslr_57") {
  query = sprintf("SELECT * FROM %s.stats_genes where parameter_set_id=%s",db,parameter.set.id);
  data = get.vector(con,query,columns='all')

  return(data)
}

get.genes.merged = function(db="gj1_eslr_57") {
  all.node.ids = get.vector(con,sprintf('SELECT DISTINCT(node_id) FROM %s.stats_genes',db))
  
  param.sets = get.psets()

  genes = get.genes(1)
  
  for (pset in param.sets$id) {
    cur.genes = get.genes(pset,db=db)

    create.name = function(ext) {paste(param.sets[pset,]$shortname,ext,sep='')}
    
    col.dnds = create.name('.dnds')
    col.dnds.m0 = create.name('.dnds.m0')
    col.psc.count = create.name('.psc.count')
    col.weak.psc.count = create.name('.weak.psc.count')

    genes.subset = data.frame(
      a=cur.genes$node_id,
      b=cur.genes$omega_mean,
      c=cur.genes$omega_m0,
      d=cur.genes$psc_count,
      e=cur.genes$weak_psc_count
      )
    colnames(genes.subset) = c(
              'node_id',
              col.dnds,
              col.dnds.m0,
              col.psc.count,
              col.weak.psc.count
              )
    
    genes = merge(genes,genes.subset,all.x=T)
  }
  return(genes);
}

factorize = function(data,columns=names(data)) {
  data[columns] = lapply(data[columns],as.factor)
  return(data)
}

get.all.sites = function(parameter.set.id=1) {
  query = sprintf("SELECT * FROM stats_sites where parameter_set_id=%s",parameter.set.id)
  data = get.vector(con,query,columns='all')
  return(data)
}
