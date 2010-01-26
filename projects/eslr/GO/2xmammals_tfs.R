# Load packages

library(org.Hs.eg.db)
library(Category)
library(GOstats)

for (num in c(1,2,3)) {
  print(paste(num,"..."))
                                        # Auto-generate the file name by pasting together the strings and numbers.
  file = paste("2xmammals_data","/","list",num,".txt",sep="")
  
                                        # Read in the list of genes.
  all.genes = read.table(file,sep="\t",header=F,stringsAsFactors=F)
                                        # Rename the column label (for clarity)
  names(all.genes) = "id"
  
                                        # Some of the genes are ENSG___ names... separate those out.
  ens.genes = subset(all.genes,grepl("ENSG",id))
  ens.genes = unlist(strsplit(ens.genes$id,"_"))
  ens.genes = ens.genes[seq(from=1,to=length(ens.genes),by=2)]
                                        # Non-ensembl genes.
  rest.genes = subset(all.genes,!grepl("ENSG",id))
  
                                        # Get the Entrez Gene identifiers for Ensembl gene stable IDs.
  ens.egs = mget(ens.genes,org.Hs.egENSEMBL2EG,ifnotfound=NA)
                                        # Get the Entrez Gene identifiers for non-Ensembl gene names.
  name.egs = mget(rest.genes$id,revmap(org.Hs.egSYMBOL),ifnotfound=NA)  
  subset.egs = c(ens.egs,name.egs)
                                        # Create the gene "universe" by taking the L keys of the Ensembl 2 EG mapping.
  all.egs = Lkeys(org.Hs.egENSEMBL2EG)
                                        # Require at least 3 GO annotations for a term to be enriched.
  add	<- list(categorySize=3)
  
  for (go in c("BP","MF","CC")) { # For each GO ontology...
    print(paste(go,"..."))
    para <- new("GOHyperGParams", geneIds=subset.egs, universeGeneIds=all.egs, annotation="org.Hs.eg.db", ontology=go, pvalueCutoff = 0.1, conditional=T)
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
  write.csv(all.summaries,file=paste("enrich_list",num,".csv",sep=""))
  rm(all.summaries)  
}
