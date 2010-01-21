library(RMySQL)
drv <- dbDriver('MySQL')
if (!exists('con')) {
  con <- dbConnect(drv, host='compara2', port=3306, user='ensadmin', password='ensembl',dbname='gj1_57')
}

# Get the sim_sets.
if (!exists('node_sets')) {
  res = dbSendQuery(con,paste(
    'select name,ns.node_set_id,node_id',
    ', leaf_count(node_id) AS size',
    ', num_human_genes(node_id) AS hums',
    'FROM node_set_member nsm JOIN node_set ns USING(node_set_id)'
  ))
  node_sets = fetch(res,n=-1)
  dbClearResult(res)
}

l = list(id=node_sets$node_set_id)
s = node_sets$size
h = node_sets$hums
counts = aggregate(s,l,length)$x
medians = aggregate(s,l,median)$x
means = aggregate(s,l,mean)$x
mins = aggregate(s,l,min)$x
maxes = aggregate(s,l,max)$x
hum_count = aggregate(h,l,sum)$x
hum_mean = aggregate(h,l,mean)$x

unq = subset(node_sets,!duplicated(node_set_id))
all = data.frame(id=unq$node_set_id,name=unq$name)
data = data.frame(id=1:length(medians),count=counts,median=medians,
  mean=means,min=mins,max=maxes,hum_count=hum_count,hum_mean=hum_mean)
all = merge(all,data,by='id')

dbDisconnect(con)
