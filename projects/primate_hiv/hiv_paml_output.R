library(plyr)
library(ggplot2)
dbname <<- 'gj1_hiv_full_57_c'
source("~/src/greg-ensembl/scripts/collect_sitewise.R")

get.genes <- function() {
  genes <- get.vector(con, "select * from genes;")

  genes$lrt <- 2* (genes$m8_lnl - genes$m7_lnl)
  genes$pval <- pchisq(genes$lrt, df=2, lower.tail=F)
  genes$pval.bh <- p.adjust(genes$pval, method='BH')

  save(genes, file=wd.f("genes.Rdata"))
  assign("genes", genes, envir=.GlobalEnv)
}

interaction.dnds <- function() {
  x <- genes

  # Downloaded from Biomart on Ensembl 62.
  # http://www.ensembl.org/biomart/martview/89acb8d7f0b4b12f6f1bd121e6d490a1/89acb8d7f0b4b12f6f1bd121e6d490a1/89acb8d7f0b4b12f6f1bd121e6d490a1/89acb8d7f0b4b12f6f1bd121e6d490a1?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.refseq_peptide|hsapiens_gene_ensembl.default.feature_page.chromosome_name|hsapiens_gene_ensembl.default.feature_page.start_position|hsapiens_gene_ensembl.default.feature_page.end_position|hsapiens_gene_ensembl.default.feature_page.external_gene_id|hsapiens_gene_ensembl.default.feature_page.strand&FILTERS=hsapiens_gene_ensembl.default.filters.source."ensembl"|hsapiens_gene_ensembl.default.filters.biotype."protein_coding"&VISIBLEPANEL=resultspanel
  ens.ncbi <- read.csv(wd.f("ens_ncbi_ids.txt"))
  ens.ncbi$gene_name <- ens.ncbi[, 'Associated.Gene.Name']
  assign("ens.ncbi", ens.ncbi, envir=.GlobalEnv)

  unique.genes <- function(x) {
    x <- x[!duplicated(x$gene_name), ]    
    return(x)
  }

  summarize.genes <- function(x) {
    x <- x[!duplicated(x$gene_name), ]    
    x <- subset(x, paml_dnds <= 5)

    n.pos <- nrow(subset(x, pval.bh < 0.05))
    f.pos <- n.pos / nrow(x)
    dnds.median <- median(x$paml_dnds)
    dnds.mean <- mean(x$paml_dnds)
    dnds.var <- var(x$paml_dnds)

    top.genes <- x[order(x$pval),]
    top.genes <- head(top.genes, n=5)
    top.gene.names <- paste(top.genes$gene_name, collapse=' ')
    if (nrow(x) > 5) {
      top.gene.names <- paste(top.gene.names, '...')
    }

    key.tbl <- table(x$Keyword)
    keys <- names(key.tbl)
    top.keys <- keys[order(-key.tbl)]
    top.keys <- head(top.keys, n=3)
    top.keys <- paste(top.keys, collapse=', ')
    if (length(keys) > 3) {
      top.keys <- paste(top.keys, '...')
    }

    return(data.frame(
      n = nrow(x),
      dnds_median = dnds.median,
      dnds_mean = dnds.mean,
      dnds_var = dnds.var,
      n.pos = n.pos,
      f.pos = f.pos,
      genes = top.gene.names,
      interactions = top.keys
    ))
  }

  int.tbl <- function(interaction.file, lbl) {
    print(interaction.file)
    ints <- read.table(wd.f(interaction.file), sep="\t", header=T)
    ints$X <- NULL
    ints$ncbi_id <- ints[, 'Human_Prot_Acc']
    ints$ncbi_id <- sub("\\.\\d+", "", ints$ncbi_id)

    # merge genes with ens-ncbi mapping.
    gene.ens <- merge(x, ens.ncbi)
    gene.ens$ncbi_id <- gene.ens[, 'RefSeq.Protein.ID']
    # merge genes/ncbi with interaction data.
    gene.hiv <- merge(gene.ens, ints)
    gene.hiv$hiv_name <- gene.hiv[, 'HIV.1_Prot_Name']

    total <- summarize.genes(x)
    total$Keyword <- 'All genes'
    total$hiv_name <- ''

    present <- summarize.genes(gene.hiv)
    present$Keyword <- 'All interacting genes'
    present$hiv_name <- ''

    # Get all unique gene-pep interactions.
    cur.genes <- ddply(gene.hiv, c('Keyword', 'hiv_name'), unique.genes)
    all.ints <<- rbind(all.ints, cur.genes)

    all.merged <<- rbind(all.merged, gene.hiv)

    # Summarize by each pep-interaction combo.
    by.pep.int <- ddply(gene.hiv, c('Keyword', 'hiv_name'), summarize.genes)
    by.pep.int <- by.pep.int[order(by.pep.int$n), ]
    all.pep.ints <<- rbind(all.pep.ints, by.pep.int)

    # Summarize by all interactions for each pep.
    pep.summary <- ddply(gene.hiv, c('hiv_name'), summarize.genes)
    print(pep.summary)
    all.pep.sums <<- rbind(all.pep.sums, pep.summary)
  }

  all.pep.ints <<- data.frame()
  all.pep.sums <<- data.frame()
  all.ints <<- data.frame()
  all.merged <<- data.frame()
  for (suffix in c('env', 'gag', 'nef', 'pol', 'rev', 'tat', 'vif', 'vpr', 'vpu')) {
    int.tbl(paste('interactions_', suffix, '.txt', sep=''), suffix)
  }

  tot.present <- summarize.genes(all.merged)
  tot.genes <- summarize.genes(x)

  tot.present$Keyword <- ''
  tot.genes$Keyword <- ''
  tot.present$hiv_name <- 'All interacting genes'
  tot.genes$hiv_name <- 'All genes'
  all.pep.ints <- rbind(all.pep.ints, tot.present, tot.genes)
  tot.present$Keyword <- NULL
  tot.genes$Keyword <- NULL
  all.pep.sums <- rbind(all.pep.sums, tot.present, tot.genes)


  assign("all.pep.ints", all.pep.ints, envir=.GlobalEnv)
  write.csv(all.pep.ints, file=wd.f("interactions.keywords.csv"), row.names=F)
  assign("all.ints", all.ints, envir=.GlobalEnv)
  write.csv(all.ints, file=wd.f("interactions.all.csv"), row.names=F)
  assign("all.pep.sums", all.pep.sums, envir=.GlobalEnv)
  write.csv(all.pep.sums, file=wd.f("interactions.peps.csv"), row.names=F)

}

interaction.plots <- function() {

  tbl <- read.csv(wd.f("interactions.peps.csv"))

  tbl <- subset(tbl, n >= 10)

  tbl[tbl$f.pos > 0.5, 'f.pos'] <- 0.5
  tbl[, 'hiv_name'] <- paste(tbl[, 'hiv_name'], ' (', tbl[, 'n'], ')', sep='')
  tbl$hiv_name <- reorder(tbl$hiv_name, tbl$f.pos)
  tbl$hiv_name <- reorder(tbl$hiv_name, !grepl('All', tbl$hiv_name))

  p <- ggplot(tbl, aes(x=hiv_name, 
    y=dnds_mean, ymin=dnds_mean-dnds_var, ymax=dnds_mean+dnds_var,
    colour=f.pos))
  p <- p + geom_linerange(size=2)
  p <- p + geom_point(size=4)

  p <- p + scale_colour_gradient("Fraction with sitewise LRT p<0.05", limits=c(0.1, 0.5))
  p <- p + scale_y_continuous("Mean dN/dS of human interacting proteins")  
  p <- p + scale_x_discrete("HIV interacting protein")
  p <- p + opts(
    axis.text.x = theme_text(angle=90, hjust=1)
  )
  pdf(file=wd.f("interactions.peps.pdf"), width=8, height=5)
  print(p)
  dev.off()

  tbl <- read.csv(wd.f("interactions.keywords.csv"))
  tbl[, 'Keyword'] <- paste(tbl[, 'hiv_name'], tbl[, 'Keyword'], sep=' ')

  tbl <- subset(tbl, n >= 10)

  tbl[tbl$f.pos > 0.5, 'f.pos'] <- 0.5
  tbl[, 'Keyword'] <- paste(tbl[, 'Keyword'], ' (', tbl[, 'n'], ')', sep='')
  tbl$Keyword <- reorder(tbl$Keyword, tbl$f.pos)
  tbl$Keyword <- reorder(tbl$Keyword, !grepl('All', tbl$Keyword))

  p <- ggplot(tbl, aes(x=Keyword, 
    y=dnds_mean, ymin=dnds_mean-dnds_var, ymax=dnds_mean+dnds_var,
    colour=f.pos))
  p <- p + geom_linerange(size=1)
  p <- p + scale_colour_gradient("Fraction with sitewise LRT p<0.05", limits=c(0, 0.5))
  p <- p + scale_y_continuous("Mean dN/dS of human interacting proteins")  
  p <- p + scale_x_discrete("HIV interacting protein and interaction keyword")
  p <- p + facet_grid(. ~ hiv_name, scales="free", space="free")
  p <- p + opts(
    axis.text.x = theme_text(angle=90, hjust=1)
  )
  pdf(file=wd.f("interactions.keywords.pdf"), width=12, height=6)
  print(p)
  dev.off()

}

wd.f <- function(file) {
  return(paste("~/src/greg-ensembl/projects/primate_hiv/", file, sep=''))
}

