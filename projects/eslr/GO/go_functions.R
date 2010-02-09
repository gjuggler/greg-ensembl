top.terms.table = function(df.list) {
  query = "SELECT go_id AS id,term AS desc,ontology FROM go_term WHERE ontology in ('BP','CC','MF');"
  retVal = dbGetQuery(GO_dbconn(),query)
  retVal$id = as.character(retVal$id)

  rv = data.frame()
  for (i in 1:length(df.list)) {
    df = df.list[[i]]
    df = merge(retVal,df,by=c('id'))
    ordered = df[order(df$pval),]
    top_ten = ordered[1:10,]
    top_ten$rank = 1:10
    rv = rbind(rv,top_ten)
  }
  rv$ont <- NULL
  rv$go <- NULL
  return(rv)
}

create.go.table = function (df.list,names=NULL,threshold=0.05,corr.pvals=F) {
  # Get all GO terms from the DB.
  query = "SELECT go_id AS id,term AS desc,ontology FROM go_term WHERE ontology in ('BP','CC','MF');"
  retVal = dbGetQuery(GO_dbconn(),query)
  retVal$id = as.character(retVal$id)

  # Gather a list of all IDs.
  all_terms = c()
  for (i in 1:length(df.list)) {
    df = df.list[[i]]
    all_terms = c(all_terms,as.character(df$id))
  }
  unique_terms = unique(all_terms)

  indices = retVal$id %in% unique_terms
  retVal = retVal[indices,]
  print(paste("Number of GO terms:",nrow(retVal)))

  # Now, construct the new table.
  table.df = data.frame(id=retVal$id,desc=retVal$desc,ont=retVal$ont)
  pval_cols = c()
  for (i in 1:length(df.list)) {
    df = df.list[[i]]

    ids = df$id
    pvals = df$pval
    if (corr.pvals) {
      pvals = df$pval_corr
    }
    if (!is.null(names)) {
      col.name = names[i]
    } else if (!is.null(df[['type']])) {
      col.name = df[1,]$type
    } else {

    }
    pval_cols = c(pval_cols,col.name)
    new.df = data.frame(id=ids)
    new.df[[col.name]] = pvals
    table.df = merge(table.df,new.df,by=c('id'),all.x=T)
  }

  # Store a column which is the # of enriched (p<=0.05) terms:
  patterns = c()
  sums = rep(0,nrow(table.df))
  for (col.name in pval_cols) {
    col = table.df[,col.name]
    col[is.na(col)] = 1
    table.df[,col.name] = col

    col = as.integer(col <= threshold)
    patterns = paste(patterns,col,sep="")
  }
  table.df$patterns = as.character(patterns)
  only_ones = gsub("0","",table.df$patterns)
  n_ones = sapply(only_ones,nchar)
  table.df$pattern_sums = n_ones

  # Store a column which is the sum of all the p-values of constituent datasets.
  x = table.df[,pval_cols]
  #x[is.na(x)] = 1
  sums = apply(x,1,sum)
  table.df$sums = sums
  return(table.df)
}
