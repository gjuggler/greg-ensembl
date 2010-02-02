# To create the all.genes dataset, source the ../collect_eslr.R file.
source("../collect_eslr.R",echo=T)

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

do.enrichments = function() {
  objs = c()
  for (param.set in 1:8) {
    pset.name = param.sets[param.set,]$name
    all.genes = get.genes(param.set)
                                        #  ensp.genes = all.genes[grep("ENSP0",all.genes$human_gene),]
    
    hi.genes = subset(all.genes,avg_omega > 0.4)
    psc.genes = subset(all.genes,num_pscs >= 1)
    
    print(paste("Param set:",pset.name,
                "Num genes:",nrow(all.genes),
                "Hi:",nrow(hi.genes),
                "PSC:",nrow(psc.genes)
                ))
    
    for (species in c('hs','mm')) {
      go.df = get(paste("go.",species,sep=""))
      go.vec = strsplit(go.df$go,split=",",fixed=T)
      names(go.vec) = go.df$node_id
      
      for (subset in c('psc','hi')) {
        gene.subset = get(paste(subset,".genes",sep=""))
        
        bp = get.go.table(gene.subset$node_id,all.genes$node_id,go.vec,"BP")
        mf = get.go.table(gene.subset$node_id,all.genes$node_id,go.vec,"MF")
        cc = get.go.table(gene.subset$node_id,all.genes$node_id,go.vec,"CC")
        
        all = rbind(bp,mf)
        
        object.name = paste("go",param.set,species,subset,sep=".")
        assign(object.name,all)
        objs = c(objs,object.name)
                                        #write.csv(all,file=paste("go",species,subset,"csv",sep="."))
      }
    }
  }
}

save.objs = function(objs) {
  # Save all objects to a Rdata file.
  save(list=objs,file="go_objs_eslr.Rdata")

  df.list = list(go.1.hs.psc,go.1.hs.hi,go.1.mm.psc,go.1.mm.hi)
  df.names = c('hs_psc','hs_hi','mm_psc','mm_hi')
  tbl = create.go.table(df.list=df.list,names=df.names)
  write.csv(tbl,file="go_table_eslr.csv",row.names=F)
}
