dbname <- 'gj1_subs'
source("~/src/greg-ensembl/scripts/collect_sitewise.R")

query <- 'SELECT other_taxon_id AS id,length FROM stats_subs'
lengths = get.vector(con,query=query,columns='all')

library(plyr)
ddply(lengths,c("id"),function(df) table(df$length)[c(1:20)])