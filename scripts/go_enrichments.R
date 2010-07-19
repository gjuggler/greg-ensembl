go.cmd = 'select data_id, protein_id AS stable_id, group_concat(DISTINCT go_term) AS go from go_terms 
  WHERE taxon_id=%s
  AND (evidence_code is null || evidence_code != "%s")
  AND subset = "%s"
  AND namespace = "%s"
  AND ancestral_mapping = 0
  GROUP BY protein_id';

if (!exists('go.hs')) {
  go.hs = get.vector(con,sprintf(go.cmd,9606,'IEA','go','biological_process'),columns='all')
  go.hs.iea = get.vector(con,sprintf(go.cmd,9606,'','go','biological_process'),columns='all')

  go.mm = get.vector(con,sprintf(go.cmd,'10090','IEA','go','biological_process'),columns='all')
  go.mm.iea = get.vector(con,sprintf(go.cmd,'10090','','go','biological_process'),columns='all')
}

library(biomaRt)
library(topGO)

get.go.table.subset = function(int.symbols,all.symbols,symbols.to.go,ontology,nodeSize=5,description='') {
 scores <- factor(as.integer(all.symbols %in% int.symbols))
 names(scores) <- all.symbols
 return(get.go.table(scores,symbols.to.go,ontology,nodeSize,description))
}

get.go.table = function(scores,symbols.to.go,ontology,nodeSize=5,description='',scoresAreBinary=FALSE) {
  if (scoresAreBinary) {
    geneSelectionFun = function(score){return(score >= 1)}
  } else {
    geneSelectionFun = function(score){return(score < 0.1)}
  }

  GOdata <- new("topGOdata",
    ontology = ontology,
    allGenes = scores,
    annot = annFUN.gene2GO,
    gene2GO = symbols.to.go,
    geneSelectionFun = geneSelectionFun,
    nodeSize = nodeSize,
    description = description
  )

  go.terms = usedGO(GOdata)

  if (scoresAreBinary) {
    res.Fis <- runTest(GOdata,algorithm="classic",statistic="Fisher")
    res.Fis.weight <- runTest(GOdata,algorithm="weight",statistic="Fisher")
    res.Fis.elim <- runTest(GOdata,algorithm="elim",statistic="Fisher")
    res.Fis.parentChild <- runTest(GOdata,algorithm="parentChild",statistic="Fisher")
    test.df = GenTable(GOdata, 
      pval.fis = res.Fis,
      pval.fis.weight = res.Fis.weight,
      pval.fis.elim = res.Fis.elim,
      pval.fis.parentchild = res.Fis.parentChild,
      topNodes=length(go.terms)
    )
  } else {
    res.KS <- runTest(GOdata,algorithm="classic",statistic="KS")
    res.KS.elim <- runTest(GOdata,algorithm="elim",statistic="KS")
    res.T <- runTest(GOdata,algorithm="classic",statistic="t")
    res.Fis <- runTest(GOdata,algorithm="classic",statistic="Fisher")
    test.df = GenTable(GOdata,
      pval.ks.elim = res.KS.elim,
      pval.ks = res.KS,
      pval.t = res.T,
      pval.fis = res.Fis,
      topNodes=length(go.terms)
    )
  }


  return(test.df)
}

get.enrich.by.score = function(named.scores,go.df,go.field.name='node_id',nodeSize=5) {
  go.vec = strsplit(go.df$go,split=",",fixed=T)
  names(go.vec) = go.df[,go.field.name]

  bp = get.go.table(
    scores = named.scores,
    symbols.to.go = go.vec,
    ontology = 'BP',
    nodeSize = nodeSize,
    description = '',
    scoresAreBinary = FALSE
  )
  return(bp)
}

get.enrich.by.subset = function(subset,all,go.df,go.field.name='node_id',nodeSize=5) {
  go.vec = strsplit(go.df$go,split=",",fixed=T)
  names(go.vec) = go.df[,go.field.name]

  scores = as.integer(all %in% subset)
  names(scores) = all

  bp = get.go.table(
    scores = scores,
    symbols.to.go = go.vec,
    ontology = 'BP',
    nodeSize = nodeSize,
    description = '',
    scoresAreBinary = TRUE
  )
  return(bp)
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
