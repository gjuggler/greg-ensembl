library(R.oo)
library(ape)
library(plyr)
source("~/src/greg-ensembl/scripts/mysql_functions.R")
source("~/src/greg-ensembl/scripts/xtable_utils.R")
source("~/src/greg-ensembl/projects/orthologs/analyze_orthologs.R")


get.optic.data <- function() {
  if (!file.exists("groups.gz")) {
    system("wget http://genserv.anat.ox.ac.uk/downloads/clades/sauropsid_clade_version1/orthologs/groups.gz")
  }
  grps <- read.table(file="groups.gz", header=T, stringsAsFactor=F)
  token.list <- strsplit(grps$gene, split="|", fixed=T)
  grps$species <- sapply(token.list,function(y) y[1])
  grps$gene <- sapply(token.list,function(y) y[2])
  df <- grps

  if (!file.exists("optic.Rdata")) {
    all.species <- unique(df$species)
    optic <- ddply(df, .(group_id), function(x) {
      tmp.df <- data.frame(
        method_id = 'Optic',
        group_id = x[1, 'group_id'],
        species_count = length(unique(x$species)),
        leaf_count = nrow(x),
        tree_mpl = 0
      )
  
      for (spc in all.species) {
        tmp.df[, spc] <- sum(x$species == spc)
      }
      tmp.df$human_count <- tmp.df[, 'mm_hsapiens5']
      tmp.df
    })
    print(head(optic))
    save(optic, file="optic.Rdata")
  }
  load("optic.Rdata")
  optic
}

get.optic.summary <- function() {
  optic <- get.optic.data()
  optic.summary <- summarize.set(optic)
  optic.summary
}