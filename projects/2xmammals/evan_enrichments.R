load("evan_enrichments_new.Rdata")
source("go_enrichments.R")

genes$stable_id = genes$human_protein
for (iea in c('incl_iea')) {
  for (list in c('list7','list8','list9')) {
    print(list)

    if (iea == 'incl_iea') {
      go.data = get("go.hs.iea")
    } else {
      go.data = get("go.hs")
    }   

    data = read.csv(paste("evan_",list,".txt",sep=""),header=F)
    names(data) = 'name'
    data$name = gsub("(ENSG.*)_(ENST.*)","\\1",data$name)

    merge.names = merge(genes,data,by='name')
    merge.ensg = merge(genes,data,by.x='human_gene',by.y='name')
    merged = rbind(merge.names,merge.ensg)

    print(paste("before",nrow(data),"after",nrow(merged)))

    tbl = get.enrich.by.subset(
      subset=merged$stable_id,
      all=genes$stable_id,
      go.df=go.data,
      go.field.name='stable_id'
    )
    write.csv(tbl,file=paste("evan_enrich_",list,"_",iea,".csv",sep=""),row.names=F)
  }
}