# Load packages

library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(Category)
library(GOstats)

# Create the different sets of genes.
genes = subset(all.genes,!is.na(avg_omega))
genes = genes[grep("ENSP0",genes$human_gene),]

hi.threshold = 0.3;
psc.threshold = 1;

hi.omega.genes = subset(genes,avg_omega > hi.threshold)
psc.genes = subset(genes,num_pscs >= psc.threshold)

get.entrez = function(ensps) {
  egs <- unlist(mget(ensps,revmap(org.Hs.egENSEMBLPROT),ifnotfound=NA),use.names=F)
  # Remove duplicate Entrez IDs.
  egs <- egs[!duplicated(egs)]
  return(egs)
}

get.go = function(genes) {
  source("collect_eslr.R")
  
  query = sprintf("SELECT * from go_terms where source_taxon=10090 and node_id=?");

}

category.go.enrichments = function() {
                                        # Entrez IDs for the universe.
  all.egs	<- get.entrez(genes$human_gene)
  
                                        # Entrez IDs for the subsets
  hi.egs = get.entrez(hi.omega.genes$human_gene)
  psc.egs = get.entrez(psc.genes$human_gene)

  for (g in c('hi','psc')) {
    eg	<- get(paste(g, ".egs", sep=""))
    for (go in c("BP","MF","CC")) {
      print(paste("Running ",g,go,"..."))
      
      para <- new("GOHyperGParams", geneIds=eg, universeGeneIds=all.egs, annotation="org.Hs.eg.db", ontology=go, pvalueCutoff = 1, conditional=F)
      hgTest = hyperGTest(para)
      summary = summary(hgTest)
      print(summary[1:10,])
      df = data.frame(
        GO=summary[[paste("GO",go,"ID",sep="")]],
        Pvalue=summary$Pvalue,
        OddsRatio=summary$OddsRatio,
        ExpCount=summary$ExpCount,
        Count=summary$Count,
        Size=summary$Size,
        Term=summary$Term,
        ontology=go
        )
      if (!exists("all.summaries")) {
        all.summaries = df
      } else {
        all.summaries = rbind(all.summaries,df)
      }
    } # for go
    write.csv(all.summaries,file=paste("eslr_enrichments_",g,".csv",sep=""),row.names=F)
    rm(all.summaries)
  } # for g in groups
  
}
