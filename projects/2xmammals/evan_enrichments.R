load("evan_enrichments.Rdata")
source("go_enrichments.R")

genes$stable_id = genes$human_protein
for (list in c('list1','list2','list3','list4','list5','list6')) {
  print(list)

  data = read.csv(paste("evan_",list,".txt",sep=""),header=F)
  names(data) = 'name'
  data$name = gsub("(ENSG.*)_(ENST.*)","\\1",data$name)

  merge.names = merge(genes,data,by='name')
  merge.ensg = merge(genes,data,by.x='human_gene',by.y='name')
  merged = rbind(merge.names,merge.ensg)

  print(paste("before",nrow(data),"after",nrow(merged)))

#  tbl = get.enrich.by.subset(
#    subset=merged$stable_id,
#    all=genes$stable_id,
#    go.df=go.hs.iea,
#    go.field.name='stable_id'
#  )
#  head(tbl)
#  write.csv(tbl,file=paste("evan_enrich_",list,".csv",sep=""),row.names=F)
}