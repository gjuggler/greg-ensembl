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

get.vector = function(con,query,columns='all') {
  res = dbSendQuery(con,query)
  if (columns == 'all') {
    data = fetch(res,n=-1)  # Important: the n=-1 causes the entire resultset to be fetched, not just the first 500.
  } else {
    data = fetch(res,n=-1)[[columns]]
  }
  dbClearResult(res)
  return(data)
}

factorize = function(data,columns=names(data)) {
  data[columns] = lapply(data[columns],as.factor)
  return(data)
}

# Get the parameter sets.
get.psets = function(db='gj1_eslr_57') {
  query = sprintf(
    paste('SELECT n.parameter_set_id AS id,sn.parameter_value AS shortname, n.parameter_value AS name FROM %s.parameter_set n, %s.parameter_set sn',
      'WHERE n.parameter_set_id=sn.parameter_set_id AND n.parameter_name="name" AND sn.parameter_name="shortname";'),
    db,db)
  return(get.vector(con,query,columns='all'))
}

get.genes = function(parameter.set.id=1,db="gj1_eslr_57") {
  query = sprintf("SELECT * FROM %s.stats_genes where parameter_set_id=%s",db,parameter.set.id);
  data = get.vector(con,query,columns='all')

  return(data)
}

get.genes.list = function(db="gj1_eslr_57") {
  query = sprintf("SELECT * FROM %s.stats_genes",db);
  data = get.vector(con,query,columns='all')
  return(data)
}

get.genes.merged = function(db="gj1_eslr_57") {
  all.node.ids = get.vector(con,sprintf('SELECT DISTINCT(node_id) FROM %s.stats_genes',db))
  
  param.sets = get.psets(db=db)
  print(param.sets)
  
  genes = get.genes(1,db=db)
  
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

get.sites = function(parameter.set.id=1,db="gj1_eslr_57",limit="") {
  query = sprintf("SELECT * from %s.stats_sites WHERE parameter_set_id=%s %s",db,parameter.set.id,limit)
  sites = get.vector(con,query,columns='all')
  return(sites)
}

get.sites.slim = function(parameter.set.id=1,db="gj1_eslr_57",limit="") {
  fields = c('node_id','parameter_set_id','chr_name','chr_start','lrt_stat','ncod','omega','note','type')
  fields.string = paste(fields,collapse=",")
  query = sprintf("SELECT %s from %s.stats_sites WHERE parameter_set_id=%s %s",fields.string,db,parameter.set.id,limit)
  sites = get.vector(con,query,columns='all')
  return(sites)
}