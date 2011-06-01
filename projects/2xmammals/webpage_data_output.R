setwd("/nfs/users/nfs_g/gj1/scratch/gj1_2x_57/2011-01-21_01")

# Load data from the table.
dbname <- 'gj1_2x_57'
source('~/src/greg-ensembl/scripts/collect_sitewise.R')
query <- 'SELECT * FROM web_data WHERE sites_csv IS NOT NULL;'
web.genes <- get.vector(con,query)

query <- 'SELECT * FROM stats_genes;'
genes <- get.vector(con, query)

# Split into sub-clades
print("  splitting genes")
psets <- list(
  m = 1,
  p = 2,
  g = 3,
  l = 4
)
for (i in 1:length(psets)) {
  nm <- names(psets)[i]
  ps <- psets[[i]]
  cur.genes <- subset(genes, parameter_set_id==ps)

  cur.genes$chr_name <- NULL
  cur.genes$chr_start <- NULL
  cur.genes$chr_end <- NULL

  print(paste(nm, ps, nrow(cur.genes)))
  
  merged <- merge(cur.genes, web.genes, by=c('data_id'), suffixes=c('','.web'))

  cols <- colnames(merged)
  non.web.fields <- cols[grep("web",cols, invert=T)]
  merged <- subset(merged, select=non.web.fields)
  print(paste(nm, ps, nrow(merged)))
  assign(paste('genes_',nm,sep=''), merged, envir=.GlobalEnv)
}

source("~/src/greg-ensembl/scripts/liftOver.R")
genes_m <- lift.over(genes_m, 'hg19ToHg18')
genes_p <- lift.over(genes_p, 'hg19ToHg18')
genes_g <- lift.over(genes_g, 'hg19ToHg18')
genes_l <- lift.over(genes_l, 'hg19ToHg18')

lift.sites <- function(x) {
  print(str(x))
  if (! 'chr_start' %in% colnames(x)) {
    return(x)
  }
  x$chr_end <- x$chr_start + 1
  return(lift.over(x, 'hg19ToHg18'))
}

# Load sites from files.
print("  loading sites")
load("sites_1.Rdata")
sites_m <- lift.sites(sites)
print(head(sites_m))
print(nrow(sites_m))
load("sites_2.Rdata")
sites_p <- lift.sites(sites)
load("sites_3.Rdata")
sites_g <- lift.sites(sites)
load("sites_4.Rdata")
sites_l <- lift.sites(sites)

sites_info <- list(
node_id="the internal Ensembl database ID used to fetch the tree and alignment. Not useful for external purposes.",
parameter_set_id="the internal ID used to distinguish between different clade datasets. Not useful for external purposes.",
domain="the Pfam domain associated with a given site.",
filter_value="contains a value representing the stringency  of filter which the given site passed. In this data, only sites with the max filter_value of 3 were retained.",
omega="the maximum-likelihood estimate of dN/dS",
omega_lower="",
omega_upper="the lower and upper range of the 95% confidence interval on dN/dS",
lrt_stat="the likelihood ratio statistic calculated by SLR. see http://www.ebi.ac.uk/goldman/SLR for more info.",
note="a note indicating whether the site was constant (no changes in any sequence at the site) or synonymous (only synonymous changes in any sequence at the site)",
type="a string indicating the amount of evidence at the given site for negative or positive selection: positive1-4 indicates positive selection at increasing confidence, and negative1-4 indicates negative selection at increasing confidence",
chr_name="the chromosome where the site resides (in the ref. sequence)",
chr_start="the start position of the site's codon (in the ref. sequence)",
pval="a p-value for non-neutral evolution at this site. Calculated by comparing the lrt_stat to a chi-squared distribution with df=1"
)

genes_info <- list(
data_id="an internal ID used for data tracking.",
parameter_set_id="the internal ID used to distinguish between different clade datasets. Not useful for external purposes.",
aln_file="the location of the alignment file in the bulk dataset output.",
aln_pdf="the location of the alignment PDF plot in the bulk dataset output.",
aln_png="the location of the thumbnail alignment plot in the bulk dataset output.",
chr_name="the chromosome where the gene resides (in the ref. sequence)",
data_file="the location of a text file with miscellaneous data about the pipeline job",
f_neg="the proportion of sitewise values which showed evidence for negative selection at a nominal 1% false positive rate",
f_neut="the proportion of sitewise values which did not show significant evidence for negative or positive selection",
f_pos="the proportion of sitewise values which showed evidence for positive selection at a nominal 1% false positive rate",
gc_content_mean="the mean GC content over all sequences in the alignment",
gene_name="the name of the reference sequence gene",
leaf_count="the number of sequences in the alignment",
node_id="an internal Ensembl ID used to fetch sequences and trees",
params_file="the location of a file containing various parameters used by the pipeline job",
pval_fisher="a combined score for positive selection in the alignment. Fisher's method for combining independent p-values was employed.",
seq_length_mean="the mean length of all coding sequences in the alignment",
sites_csv="the location of a file containing the sitewise output for this gene in the bulk data output",
slr_dnds="the overall dN/dS parameter obtained from SLR, optimizing a M0 model of codon evolution on the given tree and alignment",
slr_kappa="the kappa parameter obtained from SLR",
stable_id="the stable ID of the reference protein",
tree_file="the location of a file with the sub-tree structure from Ensembl",
tree_mean_path="the mean root-to-tip path length of the tree",
chr_start="the start chromosome coordinate of the coding sequence",
chr_end="the end chromosome coordinate of the coding sequence",
chr_strand="the strand on which the gene resides",
sitewise_value_count="the number of sitewise values collected for this gene"
)


files <- c(ls(pattern="genes_"),ls(pattern="sites_"),ls(pattern="info"))
print("  saving all")
save(
  sites_m,sites_p,sites_g,sites_l,
  genes_m,genes_p,genes_g,genes_l,
  sites_info,
  genes_info,
  file="mammals_e57_sitewise_tables.Rdata"
)

print("  done!")