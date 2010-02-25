if (exists('drv')) {
  lapply(dbListConnections(drv),dbDisconnect)
} else {
  library(RMySQL)
  drv <- dbDriver('MySQL')
}

con <- dbConnect(drv, host='ens-research', port=3306, user='ensadmin', password='ensembl',dbname='gj1_eslr_nodesets')

# Get the sim_sets.
get.nodesets = function() {

  res = dbSendQuery(con,paste(
    'select name,ns.node_set_id,node_id',
    ', leaf_count(node_id) AS size',
    ', num_human_genes(node_id) AS num_human',
    ', human_gene_list(node_id) AS human_genes',
    ' FROM node_set_member nsm JOIN node_set ns USING(node_set_id)'
  ))
  node.sets = fetch(res,n=-1)
  dbClearResult(res)

  return(node.sets)
}

export.lists = function(node.sets) {
  split.list = split(node.sets,node.sets["name"])

  for (i in 1:length(split.list)) {
    df = split.list[[i]]
    write.csv(df,file=paste(df[1,]$name,".csv",sep=""))
  }
}

summarize.results = function(node.sets) {
  l = list(id=node.sets$node_set_id)
  s = node.sets$size
  h = node.sets$num_human
  genes = node.sets$human_genes

  counts = aggregate(s,l,length)$x
  medians = aggregate(s,l,median)$x
  means = aggregate(s,l,mean)$x
  mins = aggregate(s,l,min)$x
  maxes = aggregate(s,l,max)$x
  hum_count = aggregate(h,l,sum)$x
  hum_mean = aggregate(h,l,mean)$x
  
  unq = subset(node.sets,!duplicated(node_set_id))
  all = data.frame(id=unq$node_set_id,name=unq$name)
  data = data.frame(id=1:length(medians),count=counts,median=medians,
    mean=means,min=mins,max=maxes,hum_count=hum_count,hum_mean=hum_mean)
  all = merge(all,data,by='id')
  return(all)
}
