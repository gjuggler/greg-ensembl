source("collect_2xmammals.R")

get.merged.from.db = function(limit="") {
  query = sprintf(
    paste(
          "SELECT a.node_id AS node_id,a.chr_name AS chr_name,a.chr_start AS chr_start,",
          "a.ncod AS a_ncod,a.lrt_stat AS a_lrt_stat,a.omega_lower AS a_omega_lower,a.omega AS a_omega,a.omega_upper AS a_omega_upper,a.note AS a_note, a.type AS a_type,",
          "b.ncod AS b_ncod,b.lrt_stat AS b_lrt_stat,b.omega_lower AS b_omega_lower,b.omega AS b_omega,b.omega_upper AS b_omega_upper,b.note AS b_note, b.type AS b_type",
          "FROM %s a, %s b",
          "WHERE a.parameter_set_id=1 AND b.parameter_set_id=1 AND",
          "a.node_id=b.node_id AND a.chr_name=b.chr_name AND a.chr_start=b.chr_start %s",
          sep=" "),
    'gj1_2x_57.stats_sites',
    'gj1_2xmc_57.stats_sites',
    limit)

  print("Getting merged...");
  merged = get.vector(con,query,columns='all')
  #print("Checking names...");
  #merged = data.frame(merged,check.names=T,stringsAsFactors=T)
  #merged = factorize(merged,c("chr_name","a_note","b_note"))
  print(merged[1,])
  print("Saving file...");
  save(merged,file="~/scratch/2x_57_merged.Rdata")  
}

factorize.merged = function() {
  load(file="~/scratch/2x_57_merged.Rdata")
  gcinfo(TRUE)
  for (col in c('chr_name','a_note','a_type','b_note','b_type')) {
    print(paste("Factorizing ",col,sep=""))
    merged[[col]] = as.factor(merged[[col]])
    gc()
  }
  save(merged,file="~/scratch/2x_57_merged_factors.Rdata")
}

find.overlap = function() {
  #load(file="~/scratch/2x_57_merged_factors.Rdata",.GlobalEnv)

  pos = c('positive1','positive2','positive3','positive4')
  pos.hi = c('positive3','positive4')

  conf = subset(merged,a_type %in% pos & b_type %in% pos)
  conf.hi = subset(merged,a_type %in% pos.hi & b_type %in% pos.hi)

  save(conf,conf.hi,file="~/scratch/2x_57_merged_confident.Rdata")
}

get.from.db = function(limit="") {
  fields = paste('node_id','chr_name','chr_start','ncod','lrt_stat','omega_lower','omega','omega_upper','note','type',sep=',')
  
  print("Getting prank aln...")
  query = sprintf(
    paste(
          "SELECT %s FROM %s s",
          "WHERE s.parameter_set_id=1 AND s.chr_name IS NOT NULL %s",
          sep=" "),
    fields,
    'gj1_2x_57.stats_sites',
    limit)
  pr = get.vector(con,query,columns='all')

  print("Getting mc aln...")
  query = sprintf(
    paste(
          "SELECT %s FROM %s s",
          "WHERE s.parameter_set_id=1 AND s.chr_name IS NOT NULL %s",
          sep=" "),
    fields,
    'gj1_2xmc_57.stats_sites',
    limit)
  mc = get.vector(con,query,columns='all')

  print("Saving files...")
  save(pr,mc,file="~/scratch/2x_57_sites.Rdata")
}

load.from.file = function() {
  load(file="~/scratch/2x_57_sites.Rdata",.GlobalEnv)
}
