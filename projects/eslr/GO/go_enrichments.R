# Load packages

library(org.Mm.eg.db)
library(Category)
library(GOstats)

# NB: All mappings have to be via Entrez gene ID!
# If you have Ensembl ID, use org.Mm.egENSEMBL2EG to map between IDs

# Defining the universe IDs = all genes in this case
egUNI	<- Lkeys(org.Mm.egGENENAME)

# Random groups of genes to test here
group1	<- sample(egUNI, 200)
group2	<- sample(egUNI, 300)
group3	<- egUNI[1:200]

# Run through all groups
for (g in 1:3) {
	# Get the data
	eg	<- get(paste("group", g, sep=""))

	# File name and other parameters
	file	<- paste("enrichment_group", g, ".html", sep="")
	caption	<- paste("<br><br>Genes from group", g, "<br>Number of genes (Entrez ID):", length(eg), "<br>")
	add	<- list(categorySize=3, "htmlLinks"=TRUE)

	# Run KEGG:
	# Define parameters
	para	<- new("KEGGHyperGParams", geneIds=eg, universeGeneIds=egUNI, annotation="org.Mm.eg.db", pvalueCutoff = 0.05)
	# Run test
	HGtest	<- hyperGTest(para)
	# Write to output file
	htmlReport(HGtest, file=file, append=TRUE, label=caption, summary.args=add)

	# Running though the GO categories
	for (go in c("BP", "MF", "CC")) {
		para	<- new("GOHyperGParams", geneIds=eg, universeGeneIds=egUNI, annotation="org.Mm.eg.db", ontology=go, pvalueCutoff = 0.05, conditional=TRUE)
		HGtest	<- hyperGTest(para)
		htmlReport(HGtest, file=file, summary.args=add, append=TRUE, label=caption)
	} # for go

	# Testing for chr location
	para	<- new("ChrMapHyperGParams", geneIds=eg, universeGeneIds=egUNI, annotation="org.Mm.eg.db", pvalueCutoff = 0.05, conditional=TRUE)
	HGtest	<- hyperGTest(para)
	htmlReport(HGtest, file=file, summary.args=add[1], append=TRUE, label=caption)

} # for g in groups

