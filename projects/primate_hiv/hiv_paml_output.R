library(plyr)
library(ggplot2)
dbname <<- 'gj1_hiv_full_57_c'
source("~/src/greg-ensembl/scripts/collect_sitewise.R")

bsub.candidates <- function() {
  dbname <<- 'gj1_hiv_62'
  source("~/src/greg-ensembl/scripts/mysql_functions.R")
  con <- connect(dbname)
  
  rows <- read.table(wd.f("candidates_final_ids.txt"))  
  rows$stable_id <- rows[, 1]
  for (i in 1:nrow(rows)) {
    cur.id <- as.character(rows[i, 'stable_id'])
    job.id <- -1
    job.id <- dbGetQuery(con, sprintf("select analysis_job_id from analysis_job where input_id like '%%%s%%'", cur.id))
    job.id <- as.numeric(job.id)
    
    if (job.id == -1) {
      stop("Job not found!")
    }
    sql <- sprintf("update analysis_job set status='READY', job_claim=NULL, retry_count=0 where analysis_job_id=%d", job.id)
    print(sql)
    dbSendQuery(con, sql)
    #print(cur.nm)
    #cmd <- sprintf("bsub -q normal runWorker.pl -url mysql://ensadmin:ensembl@ens-research:3306/gj1_hiv_62 -debug 1 -job_id %d", job.id)
    #print(cmd)
    #system(cmd)
  }  
}

new.genes <- function() {
  dbname <<- 'gj1_hiv_62'
  source("~/src/greg-ensembl/scripts/collect_sitewise.R")
  genes <- get.vector(con, "select * from genes;")

  genes$lrt <- 2* (genes$m8_lnl - genes$m7_lnl)
  genes$pval <- pchisq(genes$lrt, df=2, lower.tail=F)
  genes$pval.bh <- p.adjust(genes$pval, method='BH')

  genes[genes$paml_dnds > 3, 'paml_dnds'] <- 3

  save(genes, file=wd.f("genes.e62.Rdata"))
  write.csv(genes[, manuscript.fields()], file=wd.f("genes.e62.csv"), row.names=F)
  new.genes <<- genes[, manuscript.fields()]
}

manuscript.fields <- function() {
  c(
    'stable_id_gene', 'gene_name', 'aln_length', 'masked_nucs', 'paml_dnds', 'slr_dnds', 'pval', 'pval.bh'
  )
}

get.genes <- function() {
  genes <- get.vector(con, "select * from genes;")

  genes$lrt <- 2* (genes$m8_lnl - genes$m7_lnl)
  genes$pval <- pchisq(genes$lrt, df=2, lower.tail=F)
  genes$pval.bh <- p.adjust(genes$pval, method='BH')

  genes[genes$paml_dnds > 3, 'paml_dnds'] <- 3

  save(genes, file=wd.f("genes.Rdata"))
  assign("genes", genes, envir=.GlobalEnv)
}

interaction.dnds <- function() {
  if (!exists('genes')) {
    get.genes()
  }  

  x <- genes

  # Downloaded from Biomart on Ensembl 62.
  # http://www.ensembl.org/biomart/martview/89acb8d7f0b4b12f6f1bd121e6d490a1/89acb8d7f0b4b12f6f1bd121e6d490a1/89acb8d7f0b4b12f6f1bd121e6d490a1/89acb8d7f0b4b12f6f1bd121e6d490a1?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.refseq_peptide|hsapiens_gene_ensembl.default.feature_page.chromosome_name|hsapiens_gene_ensembl.default.feature_page.start_position|hsapiens_gene_ensembl.default.feature_page.end_position|hsapiens_gene_ensembl.default.feature_page.external_gene_id|hsapiens_gene_ensembl.default.feature_page.strand&FILTERS=hsapiens_gene_ensembl.default.filters.source."ensembl"|hsapiens_gene_ensembl.default.filters.biotype."protein_coding"&VISIBLEPANEL=resultspanel
  ens.ncbi <- read.csv(wd.f("ens_ncbi_ids.txt"))
  ens.ncbi$gene_name <- ens.ncbi[, 'Associated.Gene.Name']
  assign("ens.ncbi", ens.ncbi, envir=.GlobalEnv)

  summarize.genes <- function(x) {
    x <- unique.genes(x)

    n.pos <- nrow(subset(x, pval.bh < 0.05))
    f.pos <- n.pos / nrow(x)
    dnds.median <- median(x$paml_dnds)
    dnds.mean <- mean(x$paml_dnds)
    dnds.var <- var(x$paml_dnds)

    qs <- quantile(x$paml_dnds, c(0.25, 0.5, 0.75))

    top.genes <- x[order(x$pval),]
    #top.genes <- head(top.genes, n=5)
    top.gene.names <- paste(top.genes$gene_name, collapse=' ')
    #if (nrow(x) > 5) {
    #  top.gene.names <- paste(top.gene.names, '...')
    #}

    key.tbl <- table(x$Keyword)
    keys <- names(key.tbl)
    top.keys <- keys[order(-key.tbl)]
    #top.keys <- head(top.keys, n=3)
    top.keys <- paste(top.keys, collapse=', ')
    #if (length(keys) > 3) {
    #  top.keys <- paste(top.keys, '...')
    #}

    # Test #1: mann-whitney U test for dnds distribution versus all genes
    all.genes <- unique.genes(genes)
    mwu.all <- wilcox.test(x$paml_dnds, all.genes$paml_dnds, conf.int=T)

    # Test #2: mann-whitney U test for dnds distribution versus interacting genes
    int.genes <- interacting.genes
    mwu.int <- wilcox.test(x$paml_dnds, int.genes$paml_dnds, conf.int=T)

    # Test #3: Fisher's exact test for contingency table vs all genes
    sig.p <- factor(all.genes$pval.bh < 0.05, levels=c(FALSE, TRUE), labels=c("Non-significant", "Significant"))
    cur.nms <- x$gene_name
    cur.int <- factor(all.genes$gene_name %in% cur.nms, levels=c(FALSE, TRUE), labels=c("Non-interaction", "Interaction"))
    if(length(unique(sig.p)) > 1 && length(unique(cur.int)) > 1) {
      all.fish <- fisher.test(sig.p, cur.int)
    } else {
      all.fish <- data.frame(p.value = 1, estimate=1)
    }
    
    # Test #4: Fisher's exact test for contingency table vs interacting genes
    sig.p <- factor(int.genes$pval.bh < 0.05, levels=c(FALSE, TRUE), labels=c("Non-significant", "Significant"))
    cur.nms <- x$gene_name
    cur.int <- factor(int.genes$gene_name %in% cur.nms, levels=c(FALSE, TRUE), labels=c("Non-interaction", "Interaction"))
    if(length(unique(sig.p)) > 1 && length(unique(cur.int)) > 1) {
      int.fish <- fisher.test(sig.p, cur.int)
    } else {
      int.fish <- data.frame(p.value = 1, estimate=1)
    }

    return(data.frame(
      n = nrow(x),
      dnds_median = dnds.median,
      dnds_mean = dnds.mean,
      dnds_var = dnds.var,
      q_lo = qs[1],
      q_mid = qs[2],
      q_hi = qs[3],
      n.pos = n.pos,
      f.pos = f.pos,
      genes = top.gene.names,
      interactions = top.keys,
      t1 = mwu.all$p.value * sign(mwu.all$estimate),
      t2 = mwu.int$p.value * sign(mwu.int$estimate),
      t3 = all.fish$p.value * sign(all.fish$estimate - 1),
      t4 = int.fish$p.value * sign(int.fish$estimate - 1)
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
    gene.hiv$hiv_name <- tx.hiv.names(gene.hiv$hiv_name)

    # Summarize by HIV protein and interaction keyword.
    pep.keyword.summary <- ddply(gene.hiv, c('Keyword', 'hiv_name'), summarize.genes)
    all.pep.keyword <<- rbind(all.pep.keyword, pep.keyword.summary)

    # Summarize by HIV protein.
    pep.summary <- ddply(gene.hiv, c('hiv_name'), summarize.genes)
    all.pep <<- rbind(all.pep, pep.summary)

    all.interactions <<- rbind(all.interactions, gene.hiv)
  }

  # Load up all interacting genes.
  all.interacting.genes()

  all.interactions <<- data.frame()
  all.pep.keyword <<- data.frame()
  all.pep <<- data.frame()
  for (suffix in c('env', 'gag', 'nef', 'pol', 'rev', 'tat', 'vif', 'vpr', 'vpu')) {
    int.tbl(paste('interactions_', suffix, '.txt', sep=''), suffix)
  }

  tot.present <- summarize.genes(interacting.genes)
  tot.genes <- summarize.genes(x)

  tot.present$Keyword <- ''
  tot.genes$Keyword <- ''
  tot.present$hiv_name <- 'All interacting genes'
  tot.genes$hiv_name <- 'All genes'
  all.pep.keyword <- rbind(all.pep.keyword, tot.present, tot.genes)
  tot.present$Keyword <- NULL
  tot.genes$Keyword <- NULL
  all.pep <- rbind(all.pep, tot.present, tot.genes)

  keep.fields <- c(
    'HIV.1_Prot_Name', 'HIV.1_Prot_Acc', 'HIV.1_GeneID',
    'Keyword',
    'Human_Prot_Name', 'Human_Prot_Acc', 'Human_GeneId',
    'gene_name', 'Ensembl.Gene.ID', 'RefSeq.Protein.ID',
    'paml_dnds', 'pval', 'pval.bh'
  )
  write.csv(all.interactions[, keep.fields], file=wd.f("interactions.all.csv"), row.names=F)

  assign("all.pep.keyword", all.pep.keyword, envir=.GlobalEnv)
  write.csv(all.pep.keyword, file=wd.f("interactions.keywords.csv"), row.names=F)

  assign("all.pep", all.pep, envir=.GlobalEnv)
  write.csv(all.pep, file=wd.f("interactions.peps.csv"), row.names=F)
}

output.all <- function() {
  keep.fields <- c(
    'HIV.1_Prot_Name', 'HIV.1_Prot_Acc', 'HIV.1_GeneID',
    'Keyword',
    'Human_Prot_Name', 'Human_Prot_Acc', 'Human_GeneId',
    'gene_name', 'Ensembl.Gene.ID', 'RefSeq.Protein.ID',
    'paml_dnds', 'pval', 'pval.bh'
  )
  write.csv(all.interactions[, keep.fields], file=wd.f("interactions.all.csv"), row.names=F)
}

unique.genes <- function(x) {
  x <- x[!duplicated(x$gene_name), ]    
  x
}

all.interacting.genes <- function() {
  if (!exists('genes')) {
    get.genes()
  }  
  x <- genes

  if (exists('interacting.genes')) {
    return(interacting.genes)
  }

  print("  loading all interacting genes...")

  ens.ncbi <- read.csv(wd.f("ens_ncbi_ids.txt"))
  ens.ncbi$gene_name <- ens.ncbi[, 'Associated.Gene.Name']
  assign("ens.ncbi", ens.ncbi, envir=.GlobalEnv)

  comb.df <- data.frame()
  for (suffix in c('env', 'gag', 'nef', 'pol', 'rev', 'tat', 'vif', 'vpr', 'vpu')) {
    interaction.file <- paste('interactions_', suffix, '.txt', sep='')

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
    gene.hiv$hiv_name <- tx.hiv.names(gene.hiv$hiv_name)

    comb.df <- rbind(comb.df, gene.hiv)
  }

  comb.df <- unique.genes(comb.df)
  assign("interacting.genes", comb.df, envir=.GlobalEnv)
  comb.df
}

interaction.plots <- function() {

  tbl <- read.csv(wd.f("interactions.peps.csv"))

  print(nrow(tbl))
  tbl <- subset(tbl, n > 1)
  print(nrow(tbl))

  # Benjamini-Hochberg correction for p-values.
  for (test in 1:4) {
    cur.t <- paste('t', test, sep='')
    cur.t.adj <- paste('t', test, '.adj', sep='')
    tbl[, cur.t.adj] <- p.adjust(tbl[, cur.t], method='BH')
  }

  tbl[, 'hiv_name'] <- paste(tbl[, 'hiv_name'], ' (', tbl[, 'n'], ')', sep='')
  tbl$hiv_name <- reorder(tbl$hiv_name, tbl$f.pos)
  tbl$hiv_name <- reorder(tbl$hiv_name, !grepl('All', tbl$hiv_name))

  tbl <- tbl[order(tbl$hiv_name), ]
  write.csv(tbl, file=wd.f("interactions.peps.tbl.csv"), row.names=F)  

  tbl$category <- tbl$hiv_name
  tbl[tbl$f.pos > 0.5, 'f.pos'] <- 0.5

  p <- ggplot(tbl, aes(x=category, 
     y=q_mid, ymin=q_lo, ymax=q_hi,
    colour=f.pos))
  p <- p + theme_bw()

  all.tbl <- tbl[grepl('All genes', tbl$category), ]
  p <- p + geom_hline(yintercept = all.tbl$q_mid, colour='black', linetype='dashed')

  p <- p + geom_linerange(size=2)
  p <- p + geom_point(size=4)

  p <- p + scale_colour_gradient("Fraction with site-specific positive selection", limits=c(0.1, 0.5))

  add.test.results <- function(tbl, p) {
    y.min <- min(tbl$q_lo)
    y.max <- max(tbl$q_hi)
    sig.space <- 0.5
    sig.ylim <- c(0 - sig.space * .5, 0 - sig.space*1.5)
    sig.lims <- c()
    sig.labs <- c(
      expression(MWU[all]),
      expression(MWU[int]),
      expression(Fisher[all]),
      expression(Fisher[int])
    )
    for (i in 1:4) {
      cur.t <- paste('t', i, sep='')
      cur.t.adj <- paste('t', i, '.adj', sep='')
  
      sig.sub <- tbl[abs(tbl[, cur.t]) < 0.05, ]
      sig.sub$sig.level <- 'gray70'
      sig.sub[abs(sig.sub[, cur.t.adj]) < 0.05, 'sig.level'] <- 'black'

      sig.sub[sig.sub[, cur.t] > 0, 'sig.shape'] <- 24
      sig.sub[sig.sub[, cur.t] < 0, 'sig.shape'] <- 25
    
      cur.y <- sig.ylim[1] + diff(sig.ylim) * (i/4)
      sig.lims[i] <- cur.y
      sig.sub[, 'dnds_mean'] <- cur.y
      p <- p + geom_point(data=sig.sub, aes(x=category, y=dnds_mean, fill=sig.level, shape=sig.shape), colour=NA)
      p <- p + scale_fill_identity()
      p <- p + scale_shape_identity()
    }

    y.lim <- c(sig.ylim[1], y.max)
    y.labs <- seq(from=0, to=1, by=0.2)
    y.brks <- y.labs
    y.labs <- c(sig.labs, y.labs)
    y.brks <- c(sig.lims, y.brks)

    p <- p + scale_y_continuous("Mean dN/dS of human interacting proteins", limits=range(y.brks), breaks=y.brks, labels=y.labs)
    p
  }

  p <- add.test.results(tbl, p)
  p <- p + scale_x_discrete("HIV interacting protein")
  p <- p + opts(
    axis.text.x = theme_text(angle=90, hjust=1)
  )
  pdf(file=wd.f("interactions.peps.pdf"), width=8, height=5)
  print(p)
  dev.off()

  tbl <- read.csv(wd.f("interactions.keywords.csv"))

  # Benjamini-Hochberg correction for p-values.
  for (test in 1:4) {
    cur.t <- paste('t', test, sep='')
    cur.t.adj <- paste('t', test, '.adj', sep='')
    tbl[, cur.t.adj] <- p.adjust(tbl[, cur.t], method='BH')
  }

  tbl[, 'Keyword'] <- paste(tbl[, 'hiv_name'], tbl[, 'Keyword'], sep=' ')

  tbl[, 'Keyword'] <- paste(tbl[, 'Keyword'], ' (', tbl[, 'n'], ')', sep='')
  tbl$Keyword <- reorder(tbl$Keyword, tbl$f.pos)
  tbl$Keyword <- reorder(tbl$Keyword, !grepl('All', tbl$Keyword))

  tbl <- tbl[order(tbl$Keyword), ]
  write.csv(tbl, file=wd.f("interactions.keywords.tbl.csv"), row.names=F)  

  tbl$category <- tbl$Keyword
  tbl[tbl$f.pos > 0.5, 'f.pos'] <- 0.5

  # Only show categories with more than 1 interacting gene.
  print(nrow(tbl))
  tbl <- subset(tbl, n > 1)
  print(nrow(tbl))

  p <- ggplot(tbl, aes(x=category, 
     y=q_mid, ymin=q_lo, ymax=q_hi,
    colour=f.pos))
  p <- p + theme_bw()

  all.tbl <- tbl[grepl('All genes', tbl$category), ]
  p <- p + geom_hline(yintercept = all.tbl$q_mid, colour='black', linetype='dashed')

  p <- p + geom_linerange(size=1.5)
  p <- p + geom_point(size=3)
  p <- p + scale_colour_gradient("Fraction with site-specific positive selection", limits=c(0, 0.5))

  p <- add.test.results(tbl, p)

  p <- p + scale_x_discrete("HIV interacting protein and interaction keyword")
  p <- p + facet_grid(. ~ hiv_name, scales="free", space="free")
  p <- p + opts(
    axis.text.x = theme_text(angle=90, hjust=1)
  )
  pdf(file=wd.f("interactions.keywords.pdf"), width=40, height=6)
  print(p)
  dev.off()

}

tx.hiv.names <- function(x) {
  tx.list <- list(
    'gp41' = 'gp41',
    'gp120' = 'gp120',
    'gp160' = 'gp160',
    'retropepsin' = 'protease'
  )
  nms <- names(tx.list)
  for (i in 1:length(tx.list)) {
    tx.grep <- nms[i]
    tx.sub <- tx.list[[i]]
    x <- sub(paste('.*', tx.grep, '.*', sep=''), tx.sub, x)
  }
  x
}


wd.f <- function(file) {
  return(paste("~/src/greg-ensembl/projects/primate_hiv/", file, sep=''))
}

