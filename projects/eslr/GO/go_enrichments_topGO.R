# To create the all.genes dataset, source the ../collect_eslr.R file.
if (!exists('all.genes')) {
  source("../collect_eslr.R",echo=T)
}
ensp.genes = all.genes[grep("ENSP0",all.genes$human_gene),]
nrow(ensp.genes)
ensp.genes[1:10,]$human_gene

hi.genes = subset(ensp.genes,avg_omega > 0.4)
psc.genes = subset(ensp.genes,num_pscs_weak >= 1)

go.cmd = 'select node_id, stable_id, group_concat(DISTINCT go_term) AS go from go_terms WHERE source_taxon=9606 AND evidence_code != "IEA" group by stable_id';
go.hs = get.vector(con,go.cmd,columns='all')

go.cmd = 'select node_id, stable_id, group_concat(DISTINCT go_term) AS go from go_terms WHERE source_taxon=10090 AND evidence_code != "IEA" group by stable_id';
go.mm = get.vector(con,go.cmd,columns='all')

library(biomaRt)
library(topGO)

get.go.table = function(int.symbols,all.symbols,symbols.to.go,ontology) {
  myGeneList <- factor(as.integer(all.symbols %in% int.symbols))
  names(myGeneList) <- all.symbols
  GOdata <- new("topGOdata",ontology=ontology,annot=annFUN.gene2GO,gene2GO=symbols.to.go,allGenes=myGeneList)
  
  test.stat <- new("elimCount",testStatistic=GOFisherTest,name="Fisher count")
  res.Fis <- getSigGroups(GOdata,test.stat)
  
  ts <- termStat(GOdata)
  df.terms = data.frame(
    id=row.names(ts),
    ann=ts$Annotated,
    sig=ts$Significant,
    exp=ts$Expected
    )
  df.scores = data.frame(
    id=names(res.Fis@score),
    pval=res.Fis@score
    )
  df.extra = data.frame(ont=rep(ontology,nrow(df.scores)))
  my.df = merge(df.terms,df.scores)
  my.df = cbind(my.df,df.extra)
  return(my.df)
}

for (species in c('mm')) {
  go.df = get(paste("go.",species,sep=""))
  go.vec = strsplit(go.df$go,split=",",fixed=T)
  names(go.vec) = go.df$node_id
  
  for (subset in c('psc')) {
    gene.subset = get(paste(subset,".genes",sep=""))
    
    bp = get.go.table(gene.subset$node_id,ensp.genes$node_id,go.vec,"BP")
    mf = get.go.table(gene.subset$node_id,ensp.genes$node_id,go.vec,"MF")
#    cc = get.go.table(gene.subset$human_gene,ensp.genes$human_gene,go.vec,"CC")
    
    all = rbind(bp,mf)
    
    write.csv(all,file=paste("go",species,subset,"csv",sep="."))
  }
}
