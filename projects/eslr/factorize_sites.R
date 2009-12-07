ize = function(izer) {
  function(data, columns=names(data)) {
    data[columns] = lapply(data[columns], izer)
    data
  }
}

logicalize = ize(as.logical)
characterize = ize(as.character)
factorize = ize(as.factor)

load('~/scratch/ensembl_sites.tsv.Rdata')
#data = data[1:1000,]
data = factorize(data,columns=c(
  'chr','node_id',
  'g_note','l_note','v_note','p_note','no2_note','nof_note',
  'noseq_note','only2_note'
))

data = subset(data,select=-c(
  stable_id,node_id,chr,chr_start,chr_end,
  g_lrt,l_lrt,p_lrt,v_lrt,
  no2_lrt,nof_lrt,noseq_lrt,only2_lrt,
  g_note,l_note,p_note,v_note,no2_note,nof_note,noseq_note,only2_note
))

save(data,file="~/scratch/sites_test.Rdata")