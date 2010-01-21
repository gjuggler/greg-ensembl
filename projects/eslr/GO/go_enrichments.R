# Load packages

library(org.Hs.eg.db)
library(Category)
library(GOstats)

if (!exists("all_genes")) {
  all_genes= read.table("~/src/ensembl_svn/ensembl-greg/scripts/eslr/data_dump/data/NO.backup/genes/ensembl_genes.tsv",sep="\t",header=T,stringsAsFactors=F)
}

# Create the different sets of genes.
genes = subset(all_genes,!is.na(v_omega))
genes = genes[grep("ENSP0",genes$stable_id),]

hi.threshold = 0.43;
lo.threshold = 0.12;
psc.threshold = 1;

hi.omega.genes = subset(genes,v_omega > hi.threshold)
lo.omega.genes = subset(genes,v_omega <= lo.threshold)
psc.genes = subset(genes,v_pscs >= psc.threshold)

hi.egs = mget(hi.omega.genes$stable_id,org.Hs.egENSEMBLPROT2EG,ifnotfound=NA)
lo.egs = mget(lo.omega.genes$stable_id,org.Hs.egENSEMBLPROT2EG,ifnotfound=NA)
psc.egs = mget(psc.genes$stable_id,org.Hs.egENSEMBLPROT2EG,ifnotfound=NA)

# Defining the universe IDs = all genes in this case
egUNI	<- mget(genes$stable_id,org.Hs.egENSEMBLPROT2EG,ifnotfound=NA)

                                        # Run through all groups
for (g in c('hi','lo','psc')) {
                                        # Get the data
  eg	<- get(paste(g, ".egs", sep=""))
                                        # File name and other parameters
  add	<- list(categorySize=3)
                                        # Running though the GO categories
  for (go in c("BP","MF","CC")) {
    para <- new("GOHyperGParams", geneIds=eg, universeGeneIds=egUNI, annotation="org.Hs.eg.db", ontology=go, pvalueCutoff = Inf, conditional=T)
    hgTest = hyperGTest(para)
    summary = summary(hgTest)
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
  write.csv(all.summaries,file=paste("go_enrich_",g,".csv",sep=""))
  rm(all.summaries)
} # for g in groups
