if (exists('drv')) {
  lapply(dbListConnections(drv),dbDisconnect)
} else {
  library(RMySQL)
  drv <- dbDriver('MySQL')
}

if (!exists('dbname')) {
  dbname = 'slrsim_anisimova'
  host = 'mysql-greg.ebi.ac.uk'
  port=4134
  user='slrsim'
  password='slrsim'
  userpass='slrsim:slrsim'
} else if (Sys.getenv('USER') == 'gj1') {
  host = 'ens-research'
  port=3306
  user='ensro'
  password=''
  userpass='ensro'
} else {
  host = 'mysql-greg.ebi.ac.uk'
  port=4134
  user='slrsim'
  password='slrsim'
  userpass='slrsim:slrsim'
}
con <- dbConnect(drv, host=host, port=port, user=user, password=password, dbname=dbname)
dbURL = paste("mysql://",userpass,"@",host,":",port,"/",dbname,sep="")
print(paste("Connected to:",user,"@",host,":",port,"/",dbname))
print(paste("[",dbURL,"]"))
#url = paste("mysql://slrsim:slrsim@mysql-greg.ebi.ac.uk:4134/",dbname,sep="")

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

# Grabs from the database the true and inferred omegas for the given reference sequence.
get.all.data = function(sites.cols=NULL,genes.cols=NULL,where=NULL) {
  if (is.null(sites.cols)) {
    sites.cols = "s.aln_dnds,s.lrt_stat,s.aln_position,s.seq_position,s.true_dnds,s.true_type";
  }
  if (is.vector(sites.cols)) {
    sites.cols = paste(sites.cols,collapse=",")
  }
  if (is.null(genes.cols)) {
    genes.cols = "g.node_id,g.slrsim_label";
  }
  if (is.vector(genes.cols)) {
    genes.cols = paste(genes.cols,collapse=",")
  }

  if (is.null(where)) {
    where.statement = ''
  } else {
    where.statement = where
  }

  for (necessary in c(
    'g.slrsim_label',
    'g.node_id',
    'g.phylosim_insertrate',
    'g.phylosim_deleterate',
    'g.alignment_name',
#    'g.aligner',
    'g.alignment_score_threshold',
    'g.filtering_name',
    'g.filter',
    'g.slrsim_analysis_name',
    'g.slrsim_tree_file',
    'g.slrsim_rep',
    'g.tree_total_length',
    'g.tree_mean_path',
    'g.unmasked_residue_count',
    'g.residue_count',
    'g.sum_of_pairs_score',
    'g.total_column_score',
    'g.match_bl_score',
    'g.mismatch_bl_score',
    'g.mean_bl_aligned',
    'g.mean_bl_match',
    'g.mean_bl_mismatch',
    'g.mean_entropy',
    'g.lambda',
    'g.aln_length'
    )) {
    if (length(grep(necessary, genes.cols)) == 0) {
      genes.cols = paste(genes.cols,necessary,sep=',')
    }
  }

  for (necessary in c(
    's.seq_position',
    's.aln_position',
    's.aln_dnds',
    's.lrt_stat',
    's.true_dnds',
    's.true_type',
    's.bl_aligned',
    's.bl_match',
    's.entropy'
    )) {
    if (length(grep(necessary, sites.cols)) == 0) {
      sites.cols = paste(sites.cols,necessary,sep=',')
    }
  }

  query <- sprintf("select %s,%s from genes g JOIN sites s ON g.node_id=s.node_id %s",genes.cols,sites.cols,where.statement)
  all <- get.vector(con,query,columns='all')
 
  new.col.names <- list(
    c('slrsim_tree_file', 'tree'),
    c('slrsim_analysis_name', 'analysis'),
    c('alignment_name', 'aligner'),
    c('tree_mean_path', 'tree_length'),
    c('phylosim_insertrate', 'ins_rate'),
    c('phylosim_deleterate', 'del_rate')
  )
  all <- rename.cols(all, new.col.names)

  analysis.nm <- all[1, 'analysis']
  if (grepl('SLR', analysis.nm)) {
    all[, 'lrt_stat'] <- all[, 'lrt_stat'] * sign(all[, 'aln_dnds'] - .9999)
  }

  print(head(all$lrt_stat, n=200))

  return(all)
}

get.all.genes <- function() {
  query <- sprintf("select * from genes;")
  all <- get.vector(con,query,columns='all')
  new.col.names <- list(
    c('slrsim_tree_file', 'tree'),
    c('slrsim_analysis_name', 'analysis'),
    c('tree_mean_path', 'tree_length'),
    c('phylosim_insertrate', 'ins_rate'),
    c('phylosim_deleterate', 'del_rate')
  )
  all <- rename.cols(all, new.col.names)
  return(all)                          
}

rename.cols <- function(data, old.new.list) {
  all.names <- names(data)
  for (old.new in old.new.list) {
    old <- old.new[1]
    new <- old.new[2]
    all.names[all.names == old] <- new
  }
  names(data) <- all.names
  return(data)
}
