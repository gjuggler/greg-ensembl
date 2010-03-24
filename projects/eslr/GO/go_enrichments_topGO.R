# To create the all.genes dataset, source the ../collect_eslr.R file.
source("../collect_eslr.R",echo=T)
source("go_functions.R",echo=T)

go.cmd = 'select node_id, stable_id, group_concat(DISTINCT go_term) AS go from go_terms WHERE source_taxon=%s AND evidence_code != "%s" group by stable_id';

go.hs = get.vector(con,sprintf(go.cmd,9606,'IEA'),columns='all')
go.hs.iea = get.vector(con,sprintf(go.cmd,9606,''),columns='all')

go.mm = get.vector(con,sprintf(go.cmd,'10090','IEA'),columns='all')
go.mm.iea = get.vector(con,sprintf(go.cmd,'10090',''),columns='all')

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


get.enrich.df = function(subset,all,go.df) {
  go.vec = strsplit(go.df$go,split=",",fixed=T)
  names(go.vec) = go.df$node_id
  
  bp = get.go.table(subset,all,go.vec,"BP")
  mf = get.go.table(subset,all,go.vec,"MF")
  cc = get.go.table(subset,all,go.vec,"CC")
  return(rbind(bp,mf,cc))
}

do.vizbi.enrichments = function() {
  all.genes = get.all.merged()

  hi.thresh = 0.3
  psc.thresh = 1

  p = subset(all.genes,p.psc.count >= psc.thresh)
  g = subset(all.genes,g.psc.count >= psc.thresh)
  l = subset(all.genes,l.psc.count >= psc.thresh)

  print(nrow(p))
  print(nrow(g))
  print(nrow(l))

  go.df = go.hs.iea
  p.enriched = get.enrich.df(p$node_id,all.genes$node_id,go.df)
  g.enriched = get.enrich.df(g$node_id,all.genes$node_id,go.df)
  l.enriched = get.enrich.df(l$node_id,all.genes$node_id,go.df)

  write.csv(p.enriched,file="p_psc.csv",row.names=F)
  write.csv(g.enriched,file="g_psc.csv",row.names=F)
  write.csv(l.enriched,file="l_psc.csv",row.names=F)

  p = subset(all.genes,p.dnds.m0 >= hi.thresh)
  g = subset(all.genes,g.dnds.m0 >= hi.thresh)
  l = subset(all.genes,l.dnds.m0 >= hi.thresh)

  print(nrow(p))
  print(nrow(g))
  print(nrow(l))
  
  go.df = go.hs.iea
  p.enriched = get.enrich.df(p$node_id,all.genes$node_id,go.df)
  g.enriched = get.enrich.df(g$node_id,all.genes$node_id,go.df)
  l.enriched = get.enrich.df(l$node_id,all.genes$node_id,go.df)

  write.csv(p.enriched,file="p_hi.csv",row.names=F)
  write.csv(g.enriched,file="g_hi.csv",row.names=F)
  write.csv(l.enriched,file="l_hi.csv",row.names=F)

}

do.enrichments = function() {
  all.genes = get.all.merged()
  objs = c()

  # Add one for primate-greater and glires-greater genes.
  primates.higher = subset(all.genes, p.dnds > np.dnds)
  primates.lower = subset(all.genes, p.dnds < np.dnds)
  primates.morepsc = subset(all.genes, p.psc.count > np.psc.count & np.psc.count == 0)
  glires.higher = subset(all.genes, g.dnds > ng.dnds)
  glires.lower = subset(all.genes, g.dnds < ng.dnds)

  assign('go.primates.higher',get.enrich.df(primates.higher$node_id,all.genes$node_id,go.hs),pos=.GlobalEnv)
  assign('go.primates.lower',get.enrich.df(primates.lower$node_id,all.genes$node_id,go.hs),pos=.GlobalEnv)
  assign('go.primates.morepsc',get.enrich.df(primates.morepsc$node_id,all.genes$node_id,go.hs),pos=.GlobalEnv)
  assign('go.glires.higher',get.enrich.df(glires.higher$node_id,all.genes$node_id,go.mm),pos=.GlobalEnv)
  assign('go.glires.lower',get.enrich.df(glires.lower$node_id,all.genes$node_id,go.mm),pos=.GlobalEnv)
  objs = c(objs,'go.primates.higher','go.primates.lower','go.primates.morepsc','go.glires.higher','go.glires.lower')

  param.sets = get.psets()
  
  for (param.set in param.sets$id) {
    pset.name = param.sets[param.set,]$name
    pset.shortname = param.sets[param.set,]$shortname
    
    col.dnds = paste(pset.shortname,'.dnds',sep="")
    col.psc.count = paste(pset.shortname,'.psc.count',sep="")
    hi.genes = all.genes[all.genes[[col.dnds]] > 0.4,]
    psc.genes = all.genes[all.genes[[col.psc.count]] >= 1,]
    
    print(paste("Param set:",pset.name,
                "Num genes:",nrow(all.genes),
                "Hi:",nrow(hi.genes),
                "PSC:",nrow(psc.genes)
                ))
    
    for (species in c('hs.iea','mm.iea','hs','mm')) {
      go.df = get(paste("go.",species,sep=""))
      
      for (subset in c('psc','hi')) {
        gene.subset = get(paste(subset,".genes",sep=""))

        all = get.enrich.df(gene.subset$node_id,all.genes$node_id,go.df)
        
        object.name = paste("go",param.set,species,subset,sep=".")
        print(object.name)
        assign(object.name,all,pos=.GlobalEnv)
        objs = c(objs,object.name)
      }
    }
  }      
}

save.enrichments = function() {
  ids = 1:8
  names = c('Mammals','Primates','Glires','Laurasiatheria','No 2x','2x Only','No Primates','No Glires')
  param.sets = data.frame(id=ids,name=names)
  
  df.list = list()
  df.names = c()
  i = 1
  for (go in c('iea','')) {
    for (species in c('hs','mm')) {
      for (ps in param.sets$id) {
        ps.name = param.sets[ps,]$name
        if (nchar(go) > 0) {
          species.full = paste(species,go,sep=".")
        } else {
          species.full = species
        }
        for (type in c('psc','hi')) {
          obj.name = paste("go",ps,species.full,type,sep=".")
          obj = get(obj.name)
                                        # P-value correction:
          df.list[[i]] = obj
          display.name = sprintf("%s %s (GO:%s %s)",
            ps.name, type, species, go)
          df.names = c(df.names,display.name)
          
          i = i + 1
        }
      }
    }
  }

  tbl = create.go.table(df.list=df.list,names=df.names)
  write.csv(tbl,file="eslr_go.csv")
}

save.special = function() {
  df.list = list(
    go.primates.higher, go.primates.lower, go.primates.morepsc,
    go.glires.higher, go.glires.lower
    )
  df.names = c(
    'primates_higher','primates_lower','primates_morepsc',
    'glires_higher','glires_lower'
    )
  tbl = create.go.table(df.list=df.list,names=df.names)
  write.csv(tbl,file="go_table_special.csv",row.names=F)
}

save.objs = function(objs) {
  # Save all objects to a Rdata file.
  save(list=objs,file="go_objs_eslr.Rdata")

  df.list = list(
    go.1.hs.psc,go.1.hs.hi,go.1.mm.psc,go.1.mm.hi,
    go.1.hs.iea.psc,go.1.hs.iea.hi,go.1.mm.iea.psc,go.1.mm.iea.hi)
  df.names = c(
    'hs_psc','hs_hi','mm_psc','mm_hi',
    'hs_iea_psc','hs_iea_hi','mm_iea_psc','mm_iea_hi')
  tbl = create.go.table(df.list=df.list,names=df.names)
  write.csv(tbl,file="go_table_eslr.csv",row.names=F)
}


main.done = FALSE
main = function() {
  do.enrichments()
  save.enrichments()
  save.special()
  
  assign("main.done",TRUE,pos=.GlobalEnv)
}

#if (!main.done) {main()}
