library(RMySQL)

drv <- dbDriver('MySQL')
con <- dbConnect(drv, host='compara2', user='ensro', dbname='gj1_57')

factorize = function(data,columns=names(data)) {
  data[columns] = lapply(data[columns],as.factor)
  return(data)
}

### Get the site-wise info.
getData = function() {
query = paste(
  'SELECT straight_join np.node_id AS node_id, nc1.node_id AS c1_id, nc2.node_id AS c2_id, sp.aln_position AS aln_position',
    ', sp.omega AS p_omega, sc1.omega AS c1_omega, sc2.omega AS c2_omega, sp.type AS p_type, sc1.type AS c1_type, sc2.type AS c2_type',
    ', sp.ncod AS p_ncod, sc1.ncod AS c1_ncod, sc2.ncod AS c2_ncod',

  'FROM protein_tree_node np JOIN (sitewise_aln sp, protein_tree_node nc1, protein_tree_node nc2)',
    ' ON (np.node_id=sp.node_id AND nc1.parent_id=np.node_id AND nc2.parent_id=np.node_id)',
    ' JOIN (sitewise_aln sc1, sitewise_aln sc2) ON (sc1.node_id=nc1.node_id AND sc2.node_id=nc2.node_id)',

  'WHERE nc2.node_id > nc1.node_id',
    'AND sp.aln_position=sc1.aln_position AND sc2.aln_position=sc1.aln_position',
    'AND sp.parameter_set_id=sc1.parameter_set_id AND sc2.parameter_set_id=sc1.parameter_set_id',
    'AND sc1.ncod >= 6 AND sc2.ncod >= 6'
#  'ORDER BY sp.node_id,sp.aln_position'
#  'LIMIT 5000'
)

res = dbSendQuery(con,query)
data = fetch(res,n=-1)

dup_data = factorize(data,columns=c('node_id','c1_id','c2_id','p_type','c1_type','c2_type'))
save(dup_data,file='~/scratch/dup_data.Rdata')
dbClearResult(res);
} # End getData


### Get the dupldiv tag values for the applicable node ids.
getTags = function() {
query = paste(
  'SELECT t.node_id AS node_id, t.tag AS tag, t.value AS value'
  ,'FROM protein_tree_tag t'
  ,'WHERE t.tag LIKE "dupldiv%"'
#  ,'LIMIT 100'
)
res = dbSendQuery(con,query)
data = fetch(res,n=-1)
data = as.data.frame(data)

print(nrow(data))
dup_tags = factorize(data,columns=c('node_id','tag'))
save(dup_tags,file='~/scratch/dup_tags.Rdata')

} # end getTags

#getData()
getTags()

dbDisconnect(con);
