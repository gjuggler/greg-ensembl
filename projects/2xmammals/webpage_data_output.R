setwd("/nfs/users/nfs_g/gj1/scratch/gj1_2x_57/2011-01-21_01")

# Load data from the table.
dbname <- 'gj1_2x_57'
source('~/src/greg-ensembl/scripts/collect_sitewise.R')
query <- 'SELECT * FROM web_data WHERE sites_csv IS NOT NULL;'
genes <- get.vector(con,query)

# Split into sub-clades
print("  splitting genes")
genes_m <- subset(genes,parameter_set_id==1)
genes_p <- subset(genes,parameter_set_id==2)
genes_g <- subset(genes,parameter_set_id==3)
genes_l <- subset(genes,parameter_set_id==4)

# Load sites from files.
print("  loading sites")
load("sites_1.Rdata")
sites_m <- sites
load("sites_2.Rdata")
sites_p <- sites
load("sites_3.Rdata")
sites_g <- sites
load("sites_4.Rdata")
sites_l <- sites

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
  file="mammals_slr_tables.Rdata"
)

print("  done!")