library(plyr)
library(ggplot2)
library(boot)
#library(Hmisc)
#library(goseq)
#source("~/src/greg-ensembl/scripts/go_enrichments.R")

dbname <- 'gj1_gorilla'

pst <- function(...) { return(paste(..., sep=''))}

gene.sub.counts <- function(x) {
  x <- subset(x, !is.na(mut_nsyn))
  ns <- subset(x, mut_nsyn == 1)
  s <- subset(x, mut_nsyn == 0)
  return(data.frame(
    g.ns = nrow(subset(ns, taxon_id == 9593)),
    g.s = nrow(subset(s, taxon_id == 9593)),
    h.ns = nrow(subset(ns,  taxon_id == 9606)),
    h.s = nrow(subset(s,taxon_id == 9606)),
    c.ns = nrow(subset(ns,  taxon_id == 9598)),
    c.s = nrow(subset(s,taxon_id == 9598)),
    hc.ns = nrow(subset(ns,  taxon_id == 1234)),
    hc.s = nrow(subset(s,taxon_id == 1234)),
    hcg.ns = nrow(subset(ns,  taxon_id == 207598)),
    hcg.s = nrow(subset(s,taxon_id == 207598)),
    o.ns = nrow(subset(ns,  taxon_id == 9600)),
    o.s = nrow(subset(s,taxon_id == 9600)),
    hcgo.ns = nrow(subset(ns,  taxon_id == 9604)),
    hcgo.s = nrow(subset(s,taxon_id == 9604)),
    m.ns = nrow(subset(ns,  taxon_id == 9544)),
    m.s = nrow(subset(s,taxon_id == 9544)),
    hcgom.ns = nrow(subset(ns,  taxon_id == 376913)),
    hcgom.s = nrow(subset(s,taxon_id == 376913)),
    r.ns = nrow(subset(ns,  taxon_id == 9483)),
    r.s = nrow(subset(s,taxon_id == 9483)),
    all.ns = nrow(ns),
    all.s = nrow(s)
  ))
}


gene.grantham.pval <- function(x, taxon_ids=c(9593), dnds.field='omega', all.neutral=FALSE, muts.per.codon=0.001) {
  x$dnds <- x[, dnds.field]
  if (all.neutral) {
    x$dnds = 1
  }
  x$dnds <- pmax(0.001, x$dnds)
  x$dnds <- pmin(3, x$dnds)

  x <- subset(x, is.na(mut_nsyn) | taxon_id %in% taxon_ids)
  mut <- subset(x, mut_nsyn == 1 & !is.na(mut_nsyn))

  if (nrow(mut) > 0) {
    gr.scores <- round(mut$grantham / 25)
    dnds.scores <- round(-log(mut$dnds))
    dnds.scores <- pmin(9, dnds.scores)
    dnds.scores <- pmax(-1, dnds.scores)

    total.scores <- (2 + gr.scores + dnds.scores ) / nrow(x)
    total.sum <- sum(total.scores)
  } else {
    total.sum <- 0
  }

  return(data.frame(
    pval=-total.sum
  ))
}


gene.sitewise.pval <- function(x, taxon_ids=c(9593), dnds.field='omega', all.neutral=FALSE, muts.per.codon=0.001) {
  x$dnds <- x[, dnds.field]
  if (all.neutral) {
    x$dnds = 1
  }
  x$dnds <- pmax(0.001, x$dnds)
  x$dnds <- pmin(3, x$dnds)
  #print(muts.per.codon*2/3*200)
  muts.per.codon <- 0.33
  x$exp <- x$dnds * muts.per.codon

  x$exp <- pmax(1e-12, x$exp)
  x$exp <- pmin(0.99, x$exp)

  x <- subset(x, is.na(mut_nsyn) | taxon_id %in% taxon_ids)

  lnl.without.muts <- sum(log(1 - x$exp))

  nomut <- subset(x, mut_nsyn == 0 | is.na(mut_nsyn))
  mut <- subset(x, mut_nsyn == 1 & !is.na(mut_nsyn))

  nomut.p <- log(1 - nomut$exp)
  mut.p <- log(1)
  if (nrow(mut) > 0) {
    mut.p <- log(mut$exp)
  }

  total.p <- sum(mut.p) + sum(nomut.p)

  lnl.with.muts <- total.p

  total.p <- total.p / nrow(x)
  total.p <- exp(total.p)

  #print(-2*(lnl.with.muts - lnl.without.muts) / nrow(x))

  lnl <- 2*(lnl.with.muts - lnl.without.muts) / nrow(x)

  return(data.frame(
    lnl=lnl,
    pval=lnl,
    pval.old=total.p
  ))
}

gene.binomial.pval <- function(x, taxon_ids=c(9593)) {
  aln.length <- x[1, 'aln_length']
  cur.dnds <- x[1, 'slr_dnds']

  # Extract the count of gorilla non-synonymous substitutions.
  ns.subs <- subset(x, mut_nsyn==1 & !is.na(mut_nsyn))
  g.subs <- subset(ns.subs, taxon_id %in% taxon_ids)
  n.g <- nrow(g.subs)

  #####    
  mut.rate <- 0.5e-9
  n.years <- 10e6
  mut.per.site <- mut.rate * n.years
  #####

  p.per.site <- mut.per.site * cur.dnds
  # Rough estimate: 2 non-synonymous sites per codon.
  ns.sites <- aln.length * 3 * 2/3

  res <- binom.test(n.g, ns.sites, p.per.site, "greater")
  pval <- res$p.value
  data.frame(
    pval=pval,
    mut.per.site=mut.per.site,
    p.per.site=p.per.site,
    ns.sites = ns.sites
  )
}

test.gene.prob <- function() {
  test.data <- data.frame(
    omega = 1,
    mut_nsyn = c(rep(0, times=100)),
    taxon_id=9593
  )
  print(gene.sitewise.pval(test.data))
  test.data <- data.frame(
    omega = 1,
    mut_nsyn = c(rep(0, times=200)),
    taxon_id=9593
  )
  print(gene.sitewise.pval(test.data))
  test.data <- data.frame(
    omega = 0.5,
    mut_nsyn = c(1, rep(0, times=99)),
    taxon_id=9593
  )
  print(gene.sitewise.pval(test.data))
  test.data <- data.frame(
    omega = 0.5,
    mut_nsyn = c(1, 1, rep(0, times=198)),
    taxon_id=9593
  )
  print(gene.sitewise.pval(test.data))
  test.data <- data.frame(
    omega = 0.5,
    mut_nsyn = c(1, 1, 1, rep(0, times=297)),
    taxon_id=9593
  )
  print(gene.sitewise.pval(test.data))
  print('---')
  test.data <- data.frame(
    omega = 0.1,
    mut_nsyn = c(1, 1, 1, rep(0, times=297)),
    taxon_id=9593
  )
  print(gene.sitewise.pval(test.data))
  print('---')
  test.data <- data.frame(
    omega = 1,
    mut_nsyn = c(rep(1,times=2), rep(0, times=98)),
    taxon_id=9593
  )
  print(gene.sitewise.pval(test.data))
  test.data <- data.frame(
    omega = 0.1,
    mut_nsyn = c(rep(1,times=2), rep(0, times=98)),
    taxon_id=9593
  )
  print(gene.sitewise.pval(test.data))
  test.data <- data.frame(
    omega = 1,
    mut_nsyn = c(rep(1,times=3), rep(0, times=97)),
    taxon_id=9593
  )
  print(gene.sitewise.pval(test.data))
}

add.grantham <- function(x) {

  y <- ddply(merged, .(data_id), gene.grantham.pval, taxon_ids=c(9593))
  y$pval.gr.gorilla <- y$pval
  y$pval <- NULL

  return(merge(x, y[,c('data_id', 'pval.gr.gorilla')]))
}

get.pvals <- function(nofilter=F) {
  prefix <- ''
  subs.tbl <- 'subs'
  if (nofilter) {
    prefix='nofilter_'
    subs.tbl <- 'subs_nofilters'
  }

  # Collect genes and PAML inferred substitutions from the SQL tables.
  if (!exists("ggenes")) {
    source("~/src/greg-ensembl/scripts/collect_sitewise.R")
    print("  loading genes")
    ggenes <- get.vector(con, 'select * from genes')
    ggenes$slr_dnds <- as.numeric(ggenes$slr_dnds)
    ggenes$m0_dnds <- as.numeric(ggenes$m0_dnds)
    assign("ggenes", ggenes, envir=.GlobalEnv)
    save(ggenes, file=pdf.f('ggenes.Rdata'))
  }
  if (!exists("subs")) {
    source("~/src/greg-ensembl/scripts/collect_sitewise.R")
    print("  loading subs")
    subs <- get.vector(con, sprintf('select data_id, aln_pos, mut_nsyn, taxon_id, pattern, confidence, nongap_count from %s', subs.tbl))

    subs <- subset(subs, taxon_id %in% c(9593, 9606, 9598, 9600, 9544, 9478, 9483, 1234, 207598, 9604, 376913, 9526))
    assign("subs", subs, envir=.GlobalEnv)
    save(subs, file=pdf.f(paste(prefix, 'subs.Rdata', sep='')))
  }
  if (!exists("sites")) {
    source("~/src/greg-ensembl/scripts/collect_sitewise.R")
    print("  loading sites")
    sites <- get.vector(con, 'select data_id, aln_position as aln_pos, omega, lrt_stat, nongap_count as nongap_mammal from sites')
    sites[sites$omega < 1, 'lrt_stat'] <- -sites[sites$omega < 1, 'lrt_stat'] # Create signed LRT.
    assign("sites", sites, envir=.GlobalEnv)
    save(sites, file=pdf.f(paste(prefix, 'sites.Rdata', sep='')))
  }
  if (!exists("subs.sites")) {
    print("  merging subs & sites")
    subs.sites <- merge(sites, subs, by=c('data_id', 'aln_pos'), all.x=TRUE)
    print(nrow(subs.sites))
    assign("subs.sites", subs.sites, envir=.GlobalEnv)
    save(subs.sites, file=pdf.f(paste(prefix, 'subs.sites.Rdata', sep='')))
  }
  if (!exists("merged")) {
    print("  merging genes & subs")
    merged <- merge(ggenes[, c('data_id', 'aln_length', 'm0_dnds', 'slr_dnds')], 
      subs.sites, by=c('data_id'))
    assign("merged", merged, envir=.GlobalEnv)
    save(merged, file=pdf.f(paste(prefix, 'merged.Rdata', sep='')))
  }
  g <- ggenes
  s <- subs.sites

  return()

  muts.per.codon <- nrow(subset(s, taxon_id==9593 & !is.na(mut_nsyn))) / nrow(s)

  test.subset <- FALSE
  if (test.subset) {
    data_ids <- unique(s$data_id)
    data_ids <- head(data_ids, n=100)
    print(nrow(merged))
    merged <- subset(merged, data_id %in% data_ids)
    print(nrow(merged))
  }

  # Get an estimated likelihood value for observed gorilla substitutions based on 
  # the mutation rate, divergence time, and mammalian sitewise dN/dS.
  add.pvals <- FALSE
  if (add.pvals) {
    print("  sitewise p-value testing...")  
    res <- ddply(merged, .(data_id), gene.sitewise.pval, muts.per.codon=muts.per.codon, taxon_ids=c(9593))
    res$pval.sw.gorilla <- res$pval
    res$pval <- NULL
    cur.res <- ddply(merged, .(data_id), gene.sitewise.pval, muts.per.codon=muts.per.codon, taxon_ids=c(9593), all.neutral=TRUE)
    res$pval.sw.gorilla.neutral <- cur.res$pval
    cur.res <- ddply(merged, .(data_id), gene.sitewise.pval, muts.per.codon=muts.per.codon, taxon_ids=c(9606,1234))
    res$pval.sw.human <- cur.res$pval
    cur.res <- ddply(merged, .(data_id), gene.sitewise.pval, muts.per.codon=muts.per.codon, taxon_ids=c(9598,1234))
    res$pval.sw.chimp <- cur.res$pval
    # Use dN/dS and grantham scores to score each gene's substitutions
    #
    cur.res <- ddply(merged, .(data_id), gene.grantham.pval, taxon_ids=c(9593))
    res$pval.gr.gorilla <- cur.res$pval
    cur.res <- ddply(merged, .(data_id), gene.grantham.pval, taxon_ids=c(9606, 1234))
    res$pval.gr.human <- cur.res$pval
    cur.res <- ddply(merged, .(data_id), gene.grantham.pval, taxon_ids=c(9598, 1234))
    res$pval.gr.chimp <- cur.res$pval
    # Use the mutation rate, divergence time estimate, gene length, and gene-wide dN/dS to construct
    # a binomial test for greater than the expected number of nonsynonymous substitutions in
    # the Gorilla lineage.
    print("  whole-gene binomial testing...")
    cur.res <- ddply(merged, .(data_id), gene.binomial.pval, taxon_ids=c(9593))
    res$pval.bn.gorilla <- cur.res$pval
    cur.res <- ddply(merged, .(data_id), gene.binomial.pval, taxon_ids=c(9606,1234))
    res$pval.bn.human <- cur.res$pval
    cur.res <- ddply(merged, .(data_id), gene.binomial.pval, taxon_ids=c(9598,1234))
    res$pval.bn.chimp <- cur.res$pval  
  } else {
    print("  counting ns & s subs...")
    res <- ddply(merged, .(data_id), gene.sub.counts)
    res.merged <- merge(ggenes, res)
    print("  merging res & stats")
    x <- merge(res.merged, stats.all)
    print("  getting gene ranks")
    y <- get.top.genes(x, stat='lrt.5', n=nrow(x))
    return(y)
  }
}

get.lrts <- function() {
  if (!exists("ggenes")) {
    source("~/src/greg-ensembl/scripts/collect_sitewise.R")
    print("  loading genes")
    ggenes <- get.vector(con, 'select * from genes')
    ggenes$slr_dnds <- as.numeric(ggenes$slr_dnds)
    ggenes$m0_dnds <- as.numeric(ggenes$m0_dnds)
    assign("ggenes", ggenes, envir=.GlobalEnv)
  }

  g <- ggenes

  for (cur.type in c('all')) {
    cur.stats <- g
    # Branch models.
    for (i in 1:10) {
      lrt.key <- paste('lrt', '.', i, sep='')
      pval.key <- paste('pval', '.', i, sep='')
      pval.adj.key <- paste('pval.adj','.',i,sep='')
      omega.key <- paste('omega', '.', i, sep='')
      print(paste(cur.type, lrt.key))

      null.lnl <- g[, pst(cur.type, '_', '0_lnL')]
      null.omega <- g[, pst(cur.type, '_', '0_w_0')]
      
      exp.lnl <- g[, pst(cur.type, '_', i, '_lnL')]
      w.b.key <- pst(cur.type, '_', i, '_w_1')
      exp.omega.a <- g[, pst(cur.type, '_', i, '_w_0')]
      exp.omega.b <- exp.omega.a
      if (w.b.key %in% colnames(ggenes)) {
        exp.omega.b <- g[, w.b.key]
      }
      signed.lrt <- 2 * pmax(0, exp.lnl - null.lnl) * sign(exp.omega.b - exp.omega.a)
      pval <- 1 - pchisq( 2 * (exp.lnl - null.lnl), df=1)

      has.pv <- !is.na(pval)
      pval.adj <- pval
      if (any(has.pv)) {
        pval.adj[has.pv] <- p.adjust(pval.adj[has.pv], 'BH')
      }

      cur.stats[, lrt.key] <- signed.lrt
      cur.stats[, pval.key] <- pval
      cur.stats[, pval.adj.key] <- pval.adj
    }

    # Branch-sites models.
    get.branchsites <- TRUE
    if (get.branchsites) {
      for (i in 11:18) {    
        if (i %in% c(13,14)) {next}

        lrt.key <- paste('lrt','.',i,sep='')
        pval.key <- paste('pval','.',i,sep='')
        pval.adj.key <- paste('pval.adj','.',i,sep='')
        print(paste(cur.type, lrt.key))
  
        alt.lnl <- g[, pst(cur.type, '_', i, '_', 'alt_lnL')]
        null.lnl <- g[, pst(cur.type, '_', i, '_', 'null_lnL')]
  
        alt.lnl <- as.numeric(alt.lnl)
        null.lnl <- as.numeric(null.lnl)
   
        lrt <- 2 * pmax(0, alt.lnl - null.lnl)
        pval <- 1 - pchisq( 2 * (alt.lnl - null.lnl), df=1)
  
        has.pv <- !is.na(pval)
        pval.adj <- pval
        if (any(has.pv)) {
          pval.adj[has.pv] <- p.adjust(pval.adj[has.pv], 'BH')
        }
  
        cur.stats[, lrt.key] <- lrt
        cur.stats[, pval.key] <- pval
        cur.stats[, pval.adj.key] <- pval.adj
      }
    }

    cur.stats <- cur.stats[order(-cur.stats$lrt.5),]
    assign(paste("stats.", cur.type, sep=''), cur.stats, envir=.GlobalEnv)
  }

  cols <- colnames(stats.all)
  pval.fields <- cols[grepl('pval',cols)]
  lrt.fields <- cols[grepl('lrt',cols)]

#  for (prefix in c('noC', 'noH')) {
#    cur.stats <- get(paste('stats.', prefix,sep=''))
#    for (field in c(pval.fields, lrt.fields)) {
#      prefix.lbl <- paste(prefix, '.', field, sep='')
#      stats.all[, prefix.lbl] <- cur.stats[, field]
#    }
#  }
  
  stats.all$ils.hg <- chunk.b.gene(stats.all)
  stats.all$ils.cg <- chunk.c.gene(stats.all)

  assign('stats.all', stats.all, envir=.GlobalEnv)
  save(stats.all, file=pdf.f("stats.all.Rdata"))
  stats.all
}

reload.lrts.into.results <- function() {
  stats.all <- get.lrts()
  
  y <- get('y', envir=.GlobalEnv)
  stats.cols <- colnames(stats.all)
  for (field in stats.cols) {
    if (field %in% c('data_id')) {
      next
    }
    y[, field] <- NULL
  }

  y <- merge(y, stats.all, by='data_id')
  y
}

collect.results <- function() {
  rm(ggenes, envir=.GlobalEnv); rm(subs, envir=.GlobalEnv);
  rm(sites, envir=.GlobalEnv); rm(subs.sites, envir=.GlobalEnv); rm(merged, envir=.GlobalEnv);

  get.lrts()
  x <- get.pvals(nofilter=T)

  rm(ggenes, envir=.GlobalEnv); rm(subs, envir=.GlobalEnv);
  rm(sites, envir=.GlobalEnv); rm(subs.sites, envir=.GlobalEnv); rm(merged, envir=.GlobalEnv);

  get.lrts()
  y <- get.pvals(nofilter=F)
  y <- remove.false.positives(y)

  write.csv(y, file="~/src/greg-ensembl/projects/gorilla/table.csv",row.names=F)
  save(y, file=pdf.f("y.Rdata"))
  assign("y", y, envir=.GlobalEnv)
}

remove.false.positives <- function(x) {
  bad.genes <- c('ITPK1', 'POLR2A', 'ATN1', 'GAS6')
                       
  print("  removing false positives...")
  print(paste("  before:", nrow(x)))
  x <- subset(x, !(gene_name %in% bad.genes))
  print(paste("  after:", nrow(x)))
  return(x)
}


orthologs.table <- function() {
  if (!exists('orthologs')) {
    source("~/src/greg-ensembl/scripts/collect_sitewise.R")
    print("  loading orthologs")
    orthologs <- get.vector(con, 'select gene_name, gene_id, o_9606, o_9598, o_9593, o_9600, o_9544, o_9483 from gj1_orthologs.orthologs')
    assign("orthologs", orthologs, envir=.GlobalEnv)
  }

  if (!exists('orth.df')) {
    print('  categorizing orthologs')
    orth.df <- ddply(orthologs, 'o_9606', function(x) {
      orth.fields <- c('o_9606', 'o_9598', 'o_9593', 'o_9600', 'o_9544', 'o_9483')
  
      n.deletions <- sum(is.na(x[, orth.fields]))
      n.duplications <- sum(!is.na(as.numeric(x[, orth.fields])))
      n.paralogs <- sum(!is.na(x[, orth.fields]) & x[, orth.fields] == 'human_dup')
      if (n.deletions == 0 && n.duplications == 0 && n.paralogs == 0) {
        x$category <- 'one2one'
      } else {
        ctg <- c()
        if (n.deletions > 0) ctg <- c('deletion')
        if (n.duplications > 0) ctg <- c(ctg, 'duplication')
        if (n.paralogs > 0) ctg <- c(ctg, 'paralog')
        x$category <- paste(ctg, collapse=' & ')
      }
      
      x$deletions <- n.deletions
      x$duplications <- n.duplications
      x$paralogs <- n.paralogs
      return(x)
    })
    print(head(orth.df, n=3))
    print(table(orth.df$category))
    save(orth.df, file=pdf.f("orth.df.Rdata"))
    assign('orth.df', orth.df, envir=.GlobalEnv)
  }

  orth.fields <- c('o_9598', 'o_9593', 'o_9600', 'o_9544')
  full.orth.fields <- c('o_9598', 'o_9593', 'o_9600', 'o_9544', 'o_9483')

  deletion.df <- orth.df[, orth.fields]
  deletion.df <- !is.na(deletion.df)
  deletion.df <- !deletion.df
  write.csv(deletion.df, file=pdf.f('deletions.csv'), row.names=F)
  print(table(apply(deletion.df, 1, function(x) {paste(as.numeric(x), collapse='')})))

  deletion.df <- orth.df[, full.orth.fields]
  deletion.df <- !is.na(deletion.df)
  deletion.df <- !deletion.df
  write.csv(deletion.df, file=pdf.f('deletions_all.csv'), row.names=F)
  print(table(apply(deletion.df, 1, function(x) {paste(as.numeric(x), collapse='')})))

  dup.df <- orth.df[, orth.fields]
  for (fld in orth.fields) {
    dup.df[, fld] <- !is.na(as.numeric(dup.df[, fld]))
  }   
  write.csv(dup.df, file=pdf.f('duplications.csv'), row.names=F)
  print(table(apply(dup.df, 1, function(x) {paste(as.numeric(x), collapse='')})))

  dup.df <- orth.df[, full.orth.fields]
  for (fld in full.orth.fields) {
    dup.df[, fld] <- !is.na(as.numeric(dup.df[, fld]))
  }   
  write.csv(dup.df, file=pdf.f('duplications_all.csv'), row.names=F)
  print(table(apply(dup.df, 1, function(x) {paste(as.numeric(x), collapse='')})))


  for (fld in orth.fields) {
    if (any(is.na(orth.df[, fld]))) {
      orth.df[is.na(orth.df[, fld]), fld] <- 0
    }
  }       
  write.csv(orth.df, file=pdf.f('orthologs.csv'), row.names=F)
}

diff.datasets <- function() {
  compare.top <- function(x, x2, fld) {
    mg <- merge(get(x), get(x2), by=c('data_id','gene_name'))
    mg <- merge(mg, y[, c('data_id', 'ptrn_110', 'ptrn_101', 'ptrn_011', 'c.ns', 'h.ns', 'g.ns')])

    mg$ils1 <- (mg$ptrn_011) / (mg$ptrn_110+1)
    mg$ils2 <- (mg$ptrn_101) / (mg$ptrn_110+1)

    mg <- subset(mg, !is.na(lrt.5.x) & !is.na(lrt.5.y))
    f.x <- paste(fld,'.x',sep='')
    f.y <- paste(fld,'.y',sep='')
    mg.x <- mg[order(-mg[,f.x]),]
    mg.y <- mg[order(-mg[,f.y]),]

    n.top <- floor(nrow(mg) / 20)
    top.x <- head(mg.x$gene_name, n=n.top)
    top.y <- head(mg.y$gene_name, n=n.top)

    mg$x.sig <- FALSE
    mg$y.sig <- FALSE
    mg[mg$gene_name %in% top.x, 'x.sig'] <- TRUE
    mg[mg$gene_name %in% top.y, 'y.sig'] <- TRUE
    mg$ptrn <- paste(as.numeric(mg$x.sig), as.numeric(mg$y.sig), sep='')
    print(table(mg$ptrn))

    only.in.x <- mg[mg$ptrn == '10',]
    only.in.y <- mg[mg$ptrn == '01',]

    print(mean(only.in.x$g.ns))
    print(mean(only.in.x$h.ns))
    print(mean(only.in.x$c.ns))
    print(mean(only.in.y$g.ns))
    print(mean(only.in.y$h.ns))
    print(mean(only.in.y$c.ns))

#    print(mg[mg$ptrn == '01', c('gene_name', 'data_prefix.y', 'ptrn',
#      'lrt.5.x', 'lrt.5.y', 'lrt.9.x', 'lrt.10.x')])

    #print(head(mg[, c('gene_name', 'x.sig', 'y.sig')]))

    spr.cor <- cor(mg[, f.x], mg[, f.y], method='spearman')
    print(paste("Comparison of", fld, "from", x, "and", x2))
    print(paste("  spearman rank correlation:", spr.cor))
    #print(t.test(mg[, f.x], mg[, f.y]))
  }

  compare.top('stats.all', 'stats.noC', 'lrt.5')
  compare.top('stats.all', 'stats.noH', 'lrt.5')

  #compare.top('stats.all', 'stats.noC', 'lrt.1')
  #compare.top('stats.all', 'stats.noH', 'lrt.2')
}

branchsite.checks <- function() {
  data <- get('y', envir=.GlobalEnv)

  cor.f <- function(data1) {
    min.pv <- min(data1, na.rm=T)
    print(min.pv)
  }

 y <- subset(data, rank.lrt.5 <= 50)
 print(sum(y$pval.15 < 0.05))
 y <- subset(data, rank.lrt.1 < 50)
 print(sum(!is.na(y$pval.11) & y$pval.11 < 0.05))
 y <- subset(data, rank.lrt.2 < 50)
 print(sum(y$pval.12 < 0.05))
 y <- subset(data, rank.lrt.6 < 50)
 print(sum(y$pval.16 < 0.05))
 y <- subset(data, rank.lrt.7 < 50)
 print(sum(y$pval.17 < 0.05))
}

tiny.fields <- function() {
  return(c('gene_name', 'job_id', 'data_prefix'))
}

med.fields <- function() {
  return(c('gene_name', 'job_id', 'data_prefix', 'aln_length', 'g.ns', 'g.s'))
}

distribution.fields <- function() {
    fields <- c('gene_name', 'gene_id', 'protein_id', 'slr_dnds', 'm0_dnds',
      'aln_length', 'g.ns', 'g.s', 'h.ns', 'h.s', 'c.ns', 'c.s',      
      'lrt.1', 'lrt.2', 'lrt.5', 'lrt.6', 'lrt.7', 'lrt.8', 'lrt.9', 'lrt.10',
      'lrt.11', 'lrt.12', 'lrt.15', 'lrt.16', 'lrt.17',
      'pval.1', 'pval.adj.1',
      'pval.2', 'pval.adj.2',
      'pval.5', 'pval.adj.5', 
      'pval.6', 'pval.adj.6', 
      'pval.7', 'pval.adj.7', 
      'pval.8', 'pval.adj.8', 
      'pval.9', 'pval.adj.9', 
      'pval.10', 'pval.adj.10',
      'pval.11', 'pval.adj.11', 
      'pval.12', 'pval.adj.12', 
      'pval.15', 'pval.adj.15', 
      'pval.16', 'pval.adj.16', 
      'pval.17', 'pval.adj.17', 
      'pval.18', 'pval.adj.18',
      'masked_nucs', 'masked_ils',
      'ptrn_000', 'ptrn_100', 'ptrn_010', 'ptrn_001', 'ptrn_110', 'ptrn_011', 'ptrn_101', 'ptrn_111')
  return(fields)
}

table.fields <- function() {
    fields <- c('gene_name', 'job_id', 'data_prefix', 'slr_dnds', 'm0_dnds',
      'aln_length', 'g.ns', 'g.s', 'h.ns', 'h.s', 'c.ns', 'c.s', 'hc.ns', 'hc.s',
#      'g.f.nsyn','c.f.nsyn', 'h.f.nsyn', 'hc.f.nsyn',
       'ils.hg', 'ils.cg',
      'rank.lrt.1', 'rank.lrt.2', 'rank.lrt.5', 'rank.lrt.6', 'rank.lrt.7', 'rank.lrt.8', 'rank.lrt.9', 'rank.lrt.10',
      'lrt.1', 'lrt.2', 'lrt.5', 'lrt.6', 'lrt.7', 'lrt.8', 'lrt.9', 'lrt.10',
      'pval.1', 'pval.adj.1',
      'pval.2', 'pval.adj.2',
      'pval.5', 'pval.adj.5',
      'pval.6', 'pval.adj.6',
      'pval.7', 'pval.adj.7',
      'pval.8', 'pval.adj.8',
      'pval.9', 'pval.adj.9',
      'pval.10', 'pval.adj.10',
      'pval.11', 'pval.adj.11',
      'pval.12', 'pval.adj.12',
      'pval.15', 'pval.adj.15',
      'pval.16', 'pval.adj.16',
      'pval.17', 'pval.adj.17',
#      'pval.18', 'pval.adj.18',
      'masked_nucs', 'masked_ids', 'poor_coverage', 'masked_ils')
  return(fields)
}

submit.jobs <- function(rws) {
  print(rws[, c('gene_name', 'job_id')])
  l_ply(rws$job_id,
    function(x){
      cmd <- paste('bsub -q normal runWorker.pl -url mysql://ensadmin:ensembl@ens-research:3306/gj1_gorilla -debug 1 -job_id ',x,sep='')
      print(cmd)
      system(cmd)
    }
  )
}

renew.jobs <- function(rws) {
  print(rws[, c('gene_name', 'job_id')])
  l_ply(rws$job_id,
    function(x){
      cmd <- paste("mysql -hens-research -uensadmin -pensembl -A gj1_gorilla -e 'update analysis_job set status=\"READY\", job_claim=NULL where analysis_job_id=", x, ";'", sep='')
      print(cmd)
      system(cmd)
    }
  )
}

copy.pdfs <- function(rws) {
  print(rws[, c('gene_name', 'job_id')])            
  for (i in 1:nrow(rws)) {
    row <- rws[i,]
    dir <- paste("/nfs/users/nfs_g/gj1/scratch/gj1_gorilla/2011-04-18_01/data/genomic_primates/",
    row$data_prefix, sep='')
    aln.f <- paste(dir, "/", row$gene_name, "_plot_aln.pdf", sep='')
    tree.f <- paste(dir, "/", row$gene_name, "_plot_tree.pdf", sep='')
    fasta.f <- paste(dir, "/", row$gene_name,"_aln_masked_subs.fasta", sep='')
    unmasked.f <- paste(dir, "/", row$gene_name,"_aln_unmasked.fasta", sep='')
    print(aln.f)
    file.copy(aln.f, pdf.f(paste('pdfs/', row$gene_name, "_plot_aln.pdf", sep='')))
    print(tree.f)
    file.copy(tree.f, pdf.f(paste('pdfs/', row$gene_name, "_plot_tree.pdf", sep='')))
    print(fasta.f)
    file.copy(fasta.f, pdf.f(paste('pdfs/', row$gene_name, "_aln_masked_subs.fasta", sep='')))
    print(unmasked.f)
    file.copy(unmasked.f, pdf.f(paste('pdfs/', row$gene_name, "_aln_unmasked.fasta", sep='')))
  }
}

top.genes.table <- function(x, source="rank.lrt.5", n=50, dir='up', remove.masked=F) {
  comb.df <- data.frame()

  if (remove.masked) {
    x <- subset(x, is.na(masked_ids))
  }

  if (dir == 'up') {
    #vars <- c('rank.lrt.9', 'rank.lrt.10', 'rank.lrt.5', 'rank.g.f.nsyn')
    vars <- c(source)
    #vars <- c('rank.lrt.5', 'rank.g.f.nsyn', 'rank.pval.bn.gorilla', 'rank.pval.sw.gorilla', 'rank.pval.gr.gorilla')
  } else {
    vars <- c('rank.lrt.5', 'rank.lrt.7', 'rank.lrt.6')
  }

  available.fields <- colnames(x)
  for (var in vars) {
    fields <- table.fields()
    fields <- fields[fields %in% available.fields]
    cur.df <- head(x[order(x[, var]), fields], n=n)
    if (dir == 'down') {
      cur.df <- tail(x[order(x[, var]), fields], n=n)
    }
    cur.df$source.list <- var
    comb.df <- rbind(comb.df, cur.df)
  }
  comb.df <- comb.df[!duplicated(comb.df$gene_name),]
  comb.df <- format.numeric.df(comb.df, digits=3)

  if (dir == 'up') {
    write.csv(comb.df, file=pdf.f(paste("top.genes.", source, ".csv", sep='')), row.names=F)
  } else {
    write.csv(comb.df, file=pdf.f("bot.genes.", source, "csv"), row.names=F)
  }
  return(comb.df)
}

all.genes.table <- function() {
  y <- get("y", envir=.GlobalEnv)

  fields <- table.fields()
  fields <- fields[fields %in% colnames(y)]
    
  out.df <- y[, fields]

  write.csv(out.df, file=pdf.f("all.genes.csv"), row.names=F)

  fields <- distribution.fields()
  fields <- fields[fields %in% colnames(y)]
  out.df <- y[, fields]  
  write.csv(out.df, file=pdf.f("all.distribution.csv"), row.names=F)
}

top.tens.table <- function(x) {

  top.ten.fields <- c('cur.nm', 'gene_name', 'aln_length', 'slr_dnds', 'm0_dnds', 
    'cur.ns', 'cur.s',
    'cur.lrt', 'cur.pval', 'cur.pval.adj', 
    'cur.bs.lrt', 'cur.bs.pval', 'cur.bs.pval.adj'
  )
  
  tests.lists <- list(
    Human_branch = c(1, 'h.ns', 'h.s', 1, 11),
    Chimpanzee_branch = c(2, 'c.ns', 'c.s', 2, 12),
    Gorilla_branch = c(5, 'g.ns', 'g.s', 5, 15),
    AGA_Stem_branch = c(6, '', '', 6, 16),
    AGA_Clade_branch = c(7, '', '', 7, 17),
    Human_Chimpanzee_branch = c(8, '', '', 8, -1),
    Human_Gorilla_branch = c(9, '', '', 9, -1),
    Gorilla_Chimpanzee_branch = c(10, '', '', 10, -1),
    Human_branchsite = c(11, 'h.ns', 'h.s', 1, 11),
    Chimpanzee_branchsite = c(12, 'c.ns', 'c.s', 2, 12),
    Gorilla_branchsite = c(15, 'g.ns', 'g.s', 5, 15),
    AGA_Stem_branchsite = c(16, '', '', 6, 16),
    AGA_Clade_branchsite = c(17, '', '', 7, 17)
  )

  all.df <- data.frame()
  for (i in 1:length(tests.lists)) {
    cur.nm <- names(tests.lists)[i]
    cur.list <- tests.lists[[i]]
    sort.test <- cur.list[1]
    ns.fld <- cur.list[2]
    s.fld <- cur.list[3]
    b.test <- cur.list[4]
    bs.test <- cur.list[5]

    sort.fld <- paste('lrt.', sort.test, sep='')

    test.fld <- paste('lrt.', b.test, sep='')
    pval.fld <- paste('pval.', b.test, sep='')
    pval.adj.fld <- paste('pval.adj.', b.test, sep='')

    sorted.x <- x[order(-x[,sort.fld]),]
    top.ten <- head(sorted.x, n=10)

    if (bs.test != -1) {
      pval.bs.lrt.fld <- paste('lrt.', bs.test, sep='')
      pval.bs.fld <- paste('pval.', bs.test, sep='')
      pval.bs.adj.fld <- paste('pval.adj.', bs.test, sep='')
    
      top.ten[, 'cur.bs.lrt'] <- top.ten[, pval.bs.lrt.fld]
      top.ten[, 'cur.bs.pval'] <- top.ten[, pval.bs.fld]
      top.ten[, 'cur.bs.pval.adj'] <- top.ten[, pval.bs.adj.fld]
    } else {
      top.ten[, 'cur.bs.lrt'] <- '-'
      top.ten[, 'cur.bs.pval'] <- '-'
      top.ten[, 'cur.bs.pval.adj'] <- '-'
    }

    print(top.ten)
    if (ns.fld != '') {
      top.ten[, 'cur.ns'] <- top.ten[, ns.fld]
      top.ten[, 'cur.s'] <- top.ten[, s.fld]
    } else {
      top.ten[, 'cur.ns'] <- '-'
      top.ten[, 'cur.s'] <- '-'
    }

    top.ten[, 'cur.nm'] <- cur.nm    
    top.ten[, 'cur.lrt'] <- top.ten[, test.fld]
    top.ten[, 'cur.pval'] <- top.ten[, pval.fld]
    top.ten[, 'cur.pval.adj'] <- top.ten[, pval.adj.fld]

    all.df <- rbind(all.df, top.ten[, top.ten.fields])
  }

  assign("top.ten.df", all.df, envir=.GlobalEnv)
  write.csv(all.df, file=pdf.f("top.ten.csv"), row.names=F)
}

venn.tables <- function(x) {
  
  # UP!
  quantile.t <- nrow(x) / 20
  pval.t <- qchisq(0.95, df=1)
  print(pval.t)

  g <- x$lrt.5 >= pval.t
  h <- x$lrt.1 >= pval.t
  c <- x$lrt.2 >= pval.t
#  print(x[g & h & c, 'gene_name'])

  venn.df <- data.frame(
    Human = h,
    Chimpanzee = c,
    Gorilla = g
  )
  a <- apply(venn.df,1,function(x) {paste(as.numeric(x),collapse='')})
  t <- table(a)
  print(t)
  print(t/nrow(venn.df))
#  print(names(venn.df))

  g <- x$lrt.5 <= -pval.t
  h <- x$lrt.1 <= -pval.t
  c <- x$lrt.2 <= -pval.t
  venn.df <- data.frame(
    Human = h,
    Chimpanzee = c,
    Gorilla = g
  )
  a <- apply(venn.df,1,function(x) {paste(as.numeric(x),collapse='')})
  t <- table(a)
  print(t)
#  print(t/nrow(venn.df))
#  print(names(venn.df))
#  print(x[g & h & c, 'gene_name'])

  return()


  colnames(tbl.df) <- c(
  'H', 'C', 'H+CH', 'C+CH', 'G',
  'AGA Branch', 'AGA Clade',
  'H+C', 'G+H', 'G+C',
  'bs AGA Branch', 'bs H', 'bs C', 'bs G', 'bs O',
  'bs GA Branch',
  'bs GA+marmoset'
  )


  pdf(file=pdf.f('venn.up.pdf'), width=4, height=4)
  p <- venneuler(venn.df)
  print(p)
  plot(p)
  dev.off()

  # DOWN!

  pdf(file=pdf.f('venn.down.pdf'), width=4, height=4)
  p <- venneuler(venn.df)
  plot(p)
  dev.off()
}

summary.table <- function(x) {

  x <- x[,which(grepl("(gene_name|lrt|pval)", colnames(x)))]

  tbl.df <- data.frame()

  # Get the parallel LRT_min
  x$lrt.18 <- pmin(x$lrt.1, x$lrt.2)
  x$lrt.19 <- pmin(x$lrt.1, x$lrt.5)
  x$lrt.20 <- pmin(x$lrt.2, x$lrt.5)

  x$lrt.18 <- ifelse(x$lrt.18 > 0, x$lrt.18, pmin(pmax(x$lrt.1, x$lrt.2),0))
  x$lrt.19 <- ifelse(x$lrt.19 > 0, x$lrt.19, pmin(pmax(x$lrt.1, x$lrt.5),0))
  x$lrt.20 <- ifelse(x$lrt.20 > 0, x$lrt.20, pmin(pmax(x$lrt.2, x$lrt.5),0))

  x$pval.18 <- 1 - pchisq(abs(x$lrt.18), df=1)
  x$pval.19 <- 1 - pchisq(abs(x$lrt.19), df=1)
  x$pval.20 <- 1 - pchisq(abs(x$lrt.20), df=1)

  x$pval.adj.18 <- p.adjust(x$pval.18, method='BH')
  x$pval.adj.19 <- p.adjust(x$pval.19, method='BH')
  x$pval.adj.20 <- p.adjust(x$pval.20, method='BH')

  print(x[1:5, c('gene_name', 'lrt.1', 'lrt.2', 'lrt.5', 'lrt.18', 'lrt.19', 'lrt.20', 'pval.18', 'pval.19', 'pval.20')])

  for (i in 18:20) {
    pval.key <- paste('pval.', i, sep='')
    pval.adj.key <- paste('pval.adj.', i, sep='')
    print(pval.adj.key)
    lrt.key <- paste('lrt.', i, sep='')
    x$tmp.pval <- x[, pval.key]
    x$tmp.pval.adj <- x[, pval.adj.key]
    if (i <= 10 || i >= 18) {
      x$tmp.lrt <- x[, lrt.key]
      tbl.df[i, 'Accelerated'] <- nrow(subset(x, tmp.pval < 0.05 & tmp.lrt > 0))
      tbl.df[i, 'Accelerated (FDR)'] <- nrow(subset(x, tmp.pval.adj < 0.1 & tmp.lrt > 0))
      tbl.df[i, 'Decelerated'] <- nrow(subset(x, tmp.pval < 0.05 & tmp.lrt < 0))
      tbl.df[i, 'Decelerated (FDR)'] <- nrow(subset(x, tmp.pval.adj < 0.1 & tmp.lrt < 0))
    } else {
      tbl.df[i, 'Accelerated'] <- nrow(subset(x, tmp.pval < 0.05))
      tbl.df[i, 'Accelerated (FDR)'] <- nrow(subset(x, tmp.pval.adj < 0.1))
    }
  }
  assign("summary.table", tbl.df, envir=.GlobalEnv)
  write.csv(tbl.df, file=pdf.f('summary.table.csv'))
}

get.top.genes <- function(x, stat='pval.sw.gorilla', dir=1, n=nrow(x)) {
  # FILTER: Only take genes with at least 1 substitution in each species.
#  x <- subset(x, g.ns+g.s > 0 & h.ns+h.s > 0 & c.ns+c.s > 0)

  # FILTER: at least 1 gor non-synonymous sub.
#  x <- subset(x, g.ns > 0)

  # FILTER: Only take genes with dnds > 0
#  x <- subset(x, slr_dnds > 0)

  if (stat %in% c('lrt.5','lrt.4','lrt.3')) {
    dir <- -1
  }

  x$g.f.nsyn <- x$g.ns / x$aln_length
  x$c.f.nsyn <- x$c.ns / x$aln_length
  x$h.f.nsyn <- x$h.ns / x$aln_length
  x$hc.f.nsyn <- x$hc.ns / x$aln_length

  rank.fields <- c()
  available.cols <- colnames(x)
  for (field in c(
    'lrt.5', 'pval.gr.gorilla', 'pval.sw.gorilla', 'pval.bn.gorilla', 'g.f.nsyn', 
    'lrt.4', 'pval.gr.human', 'pval.sw.human', 'pval.bn.human', 'h.f.nsyn',
    'lrt.3', 'pval.gr.chimp', 'pval.sw.chimp', 'pval.bn.chimp', 'c.f.nsyn',
    'lrt.8', 'lrt.9', 'lrt.10', 'lrt.6', 'lrt.7', 'lrt.1', 'lrt.2')) {
    rank.field <- paste('rank',field,sep='.')
    if (!(field %in% available.cols)) {
      next
    }
    if (field %in% c('g.f.nsyn','c.f.nsyn', 'h.f.nsyn', 'lrt.5', 'lrt.4', 'lrt.3', 'lrt.2', 'lrt.1', 'lrt.6', 'lrt.7', 'lrt.8', 'lrt.9', 'lrt.10')) {
      x[, rank.field] <- rank(-x[,field])
    } else {
      x[, rank.field] <- rank(x[,field])
    }
    rank.fields <- c(rank.fields, rank.field)
  }


  keep.fields <- colnames(x)
  keep.fields <- keep.fields[!grepl("lnL", keep.fields)]
  x <- x[order(x[,stat]*dir), keep.fields]
  return(head(x,n=n))
}

collapse.go.list <- function(x) {
  if (!exists("all.go.defs")) {
    all.go <- unique(unlist(x))
    defs <- getTermsDefinition(all.go, 'BP', numChar=30)
    assign("all.go.defs", defs, envir=.GlobalEnv)
  }

  x.df <- ldply(x, function(data) {
    defs <- all.go.defs[data]
    return(data.frame(
      terms=paste(data, collapse=' '),
      defs=paste(defs, collapse='; ')
    ))
  })
  return(x.df)
}

go.table <- function(x) {
  if (!exists('go.ens')) {
    go.rdata.file <- "~/scratch/gorilla/go_60.Rdata"
    load(go.rdata.file)
    assign('go.ens', go.ens, envir=.GlobalEnv)
  }
  source("~/src/greg-ensembl/scripts/go_enrichments.R")

  go.tbl <- ddply(x, .(protein_id), function(data) {
    terms <- unlist(go.ens[data$protein_id])
    defs <- getTermsDefinition(terms, 'BP', 50)
    go.ids <- paste(names(defs), collapse=' ')
    go.desc <- paste(defs, collapse='; ')
    return(data.frame(
      gene_name=data$gene_name,
      protein_id=data$protein_id,
      gene_id=data$gene_id,
      go=go.ids,
      desc=go.desc
    ))
  })
  assign("go.desc", go.tbl, envir=.GlobalEnv)
}

percent.overlap <- function(x) {

  top <- floor(nrow(y)/10)
  thresh <- qchisq(0.95, df=1)

  fu <- function(x, y) {return(length(union(x$gene_name, y$gene_name)))}
  fi <- function(x, y) {return(length(intersect(x$gene_name, y$gene_name)))}

  # X percent of genes accelerated in gorilla are accelerated in humans.
  g <- subset(x, lrt.5 >= thresh)
  gh <- subset(x, lrt.5 >= thresh & lrt.1 >= thresh)
  print(fi(g, gh) / fu(g, gh))

  # X percent of genes accelerated in gorilla are accelerated in humans.
  c <- subset(x, lrt.2 >= thresh)
  gc <- subset(x, lrt.2 >= thresh & lrt.5 >= thresh)
  print(fi(c, gc) / fu(c, gc))

  h <- subset(x, lrt.1 >= thresh)
  ch <- subset(x, lrt.1 >= thresh & lrt.2 >= thresh)
  print(fi(h, ch) / fu(h, ch))

}

parallel.numbers <- function() {
  library(RUnit)

  thresh <- qchisq(0.95, df=1)
  checkEquals(nrow(subset(y, lrt.1 > thresh)), 663, "Human branch accel. check")

  p.f <- function(lbl, number) {print(paste(lbl, number))}

  parallel.tests <- function(a, b, c, orig.thresh, lbl) {

    # Get the threshold that causes the same # of parallel genes using the overlap of the
    # independent branch models to equal the # of parallel genes using the parallel branch model.
    cur.thresh <- orig.thresh
    cur.overlap <- sum(a > cur.thresh & b > cur.thresh)
    branch.model.count <- sum(c > thresh)
    repeat {
      if (cur.overlap > branch.model.count) { break }
      cur.thresh <- cur.thresh - 0.05
      cur.overlap <- sum(pmin(a,b) > cur.thresh)
    }
    cur.overlap.genes <- y[a > cur.thresh & b > cur.thresh, ]
    branch.model.genes <- y[c > orig.thresh,]

    # Get the threshold that identifies 5% of accelerated genes using the independent   
    # branch models.
    pct.thresh <- orig.thresh
    cur.overlap <- sum(a > pct.thresh & b > pct.thresh)
    pct.target.count <- length(a) / 20
    repeat {
      if (cur.overlap > pct.target.count) { break }
      pct.thresh <- pct.thresh - 0.05
      cur.overlap <- sum(pmin(a,b) > pct.thresh)
    }

    ## Randomization test.
    random.t <- function(a, b, total, reps=100) {
      n.overlaps <- c()
      for (i in 1:reps) {
        a.sample <- base::sample(1:total, a, replace=F)
        b.sample <- base::sample(1:total, b, replace=F)
        n.overlaps[i] <- length(intersect(a.sample, b.sample))
      }
      ecdf(n.overlaps)
    }
    random.ecdf <- random.t(sum(a > orig.thresh), sum(b > orig.thresh), length(a))
    random.res <- 1 - random.ecdf(sum(pmin(a,b) > orig.thresh))

    expected.n <- (sum(a > orig.thresh) / length(a)) * (sum(b > orig.thresh) / length(b)) * length(a)
    print(expected.n)

    excess <- sum(pmin(a,b) > orig.thresh) / expected.n - 1
    print(excess)
    
    a.factor <- factor(a > orig.thresh)
    b.factor <- factor(b > orig.thresh)
    fisher.res <- fisher.test(a.factor, b.factor, alternative='greater')
    fisher.pval <- fisher.res$p.value

    overlap.n <- nrow(y[a > orig.thresh & b > orig.thresh,])
    intersect.n = length(intersect(cur.overlap.genes$data_id, branch.model.genes$data_id))

    n.a <- sum(a > orig.thresh)
    n.b <- sum(b > orig.thresh)

    return(data.frame(
      lbl = lbl,
      n.a = n.a,
      n.b = n.b,
      n.overlap = overlap.n,
      n.branch = nrow(branch.model.genes),
      pct.thresh = pct.thresh,
      ol.b.intersect = intersect.n,
      ol.b.f = intersect.n / nrow(branch.model.genes),
      random.test = random.res,
      fisher.test = fisher.pval
    ))
  }

  orig.thresh <- qchisq(0.95, df=1)

#  x.gh <- parallel.tests(y$lrt.5, y$lrt.1, y$lrt.9, orig.thresh, 'gorilla-human')
#  x.gc <- parallel.tests(y$lrt.5, y$lrt.2, y$lrt.10, orig.thresh, 'gorilla-chimp')
#  x.hc <- parallel.tests(y$lrt.1, y$lrt.2, y$lrt.8, orig.thresh, 'human-chimp')  
#  print(rbind(x.gh, x.gc, x.hc))

  tt <- qchisq(0.95, df=1)
  gene.scores <- y$lrt.7 > tt
  names(gene.scores) <- y$protein_id
  pwf <- nullp(gene.scores, bias.data=y$aln_length)

  pdf(file="pwf.pdf")
  plotPWF(pwf)
  dev.off()

  # Output the top 10 genes for each parallel combination.
  top.f <- function(scalar.score, nsyn, syn) {
    sub.y <- y
    #if (!is.na(c.field)) {
    #  sub.y$scalar.score <- pmin(sub.y[, a.field], sub.y[, b.field], sub.y[, c.field])
    #} else if (!is.na(b.field)) {
    #  sub.y$scalar.score <- pmin(sub.y[, a.field], sub.y[, b.field])
    #} else {
    #  sub.y$scalar.score <- sub.y[, a.field]
    #}

    sub.y$scalar.score <- scalar.score
    sub.y$scalar.rank <- rank(sub.y$scalar.score)
    genes <- as.integer(sub.y$scalar.rank > (nrow(sub.y) - nrow(sub.y)/5))
    names(genes) <- sub.y$protein_id
    sub.y$pwf <- pwf$pwf
    sub.y$score.corr <- sub.y$scalar.score / pwf$pwf

    sub.y$nsyn <- nsyn
    sub.y$syn <- syn

    sub.y <- subset(sub.y, nsyn > 0)

    flds.short <- c('gene_name', 'aln_length', 'slr_dnds', 'm0_dnds', 'nsyn', 'syn', 'scalar.score')
    flds <- c('gene_name', 'aln_length', 'scalar.score', 'score.corr', 'pwf', 'm0_dnds', 'slr_dnds', 'lrt.1', 'lrt.2', 'lrt.5', 'lrt.8', 'lrt.9', 'lrt.10', 'g.ns', 'g.s', 'c.ns', 'c.s', 'h.ns', 'h.s')

    x <- sub.y[order(sub.y$scalar.score, decreasing=T),]
    head(x[, flds.short], n=10)
    #x2 <- sub.y[order(sub.y$score.corr, decreasing=T),]
    #print(setdiff(x[1:20, 'gene_name'], x2[1:20, 'gene_name']))
    #rbind(head(x[, flds], n=20), head(x2[, flds], n=20))
  }

  write.csv(top.f(y$lrt.1, y$h.ns, y$h.s), file=pdf.f('top.h.csv'))
  write.csv(top.f(y$lrt.2, y$c.ns, y$c.s), file=pdf.f('top.c.csv'))
  write.csv(top.f(y$lrt.5, y$g.ns, y$g.s), file=pdf.f('top.g.csv'))
  write.csv(top.f(pmin(y$lrt.1, y$lrt.5), y$h.ns + y$g.ns, y$h.s + y$g.s), file=pdf.f('top.parallel.gh.csv'))
  write.csv(top.f(pmin(y$lrt.2, y$lrt.5), y$c.ns + y$g.ns, y$c.s + y$g.s), file=pdf.f('top.parallel.gc.csv'))
  write.csv(top.f(pmin(y$lrt.1, y$lrt.2), y$h.ns + y$c.ns, y$h.s + y$c.s), file=pdf.f('top.parallel.hc.csv'))
  write.csv(top.f(pmin(y$lrt.1, y$lrt.2, y$lrt.5), y$h.ns + y$c.ns + y$g.ns, y$h.s + y$c.s + y$g.s), file=pdf.f('top.parallel.all.csv'))
}

dump.enrichments <- function() {

  single.thresh <- 3.8415
  parallel.thresh <- 1.5

  single.stats <- c('h.val', 'g.val', 'c.val', 'aga.val', 'aga.clade.val')
  parallel.stats <- c('gh.val', 'gc.val', 'hc.val')

  source("~/src/greg-ensembl/scripts/mysql_functions.R")
  con <- connect('gj1_gorilla_go')

  sql.f <- function(stat, thresh) {
    sql.cmd <- sprintf("select n, n_sig, n_exp, go_id, descr, fis_pval, topgo_pval, goseq_wall, mean_length, genes from go where stat='%s' and t=%s and n_sig >= 5 and fis_pval < 0.05 and topgo_pval < 0.25 order by topgo_pval asc",
      stat, thresh)
    print(sql.cmd)
    db.res <- dbSendQuery(con, sql.cmd)
    res.df <- fetch(db.res, n=-1)
    write.csv(res.df, file=paste('godump', stat, 'csv', sep='.'))
  }
  
  for (stat in single.stats) {
    sql.f(stat, single.thresh)
  }

  for (stat in parallel.stats) {
    sql.f(stat, parallel.thresh)
  }
}

submit.enrichment.jobs <- function() {
  ### Single-lineage enrichments.
  vals <- c('aga.clade.val')
#  vals <- c('h.val', 'g.val', 'c.val', 'g.filt.val', 'h.filt.val', 'c.filt.val', 'aga.val')
  thresholds <- qchisq(c(0.9, 0.95, 0.97, 0.98, 0.99, 0.995), df=1)
  for (i in 1:length(vals)) {
    for(j in 1:length(thresholds)) {
      cmd <- sprintf('bsub -q normal -R "select[mem>4000] rusage[mem=4000]" -M4000000 "/software/R-2.12.0/bin/R --vanilla --args %s %.4f < enrich_script.R"',
        vals[i],
        thresholds[j]
      )
      print(cmd)
      system(cmd)
    }
  }

  return()

  ### Dual-lineage enrichments.
  vals <- c('gh.val', 'hc.val', 'gc.val', 'gh.filt.val', 'hc.filt.val', 'gc.filt.val')
  thresholds <- c(0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3)

  for (i in 1:length(vals)) {
    for(j in 1:length(thresholds)) {
      cmd <- sprintf('bsub -q normal "/software/R-2.12.0/bin/R --vanilla --args %s %.4f < enrich_script.R"',
        vals[i],
        thresholds[j]
      )
      print(cmd)
      system(cmd)
    }
  }
}

parallel.enrichments <- function(val.string, sig.t) {
  go.f <- function(scalar.scores, binary.scores, lbl) {
    ### Load the GO annotations from file.
    if (!exists('go.ens')) {
      go.rdata.file <- "~/src/greg-ensembl/projects/gorilla/go_60.Rdata"
      load(go.rdata.file)
      assign('go.ens', go.ens, envir=.GlobalEnv)
    }

    if (!exists("root.godata", envir=.GlobalEnv)) {
      scores <- y$lrt.9
      names(scores) <- y$protein_id
      geneSelFn <- function(score){return(score >= qchisq(0.95, df=1))}
      root.godata <- new("topGOdata",
        ontology = "BP",
        allGenes = scores,
        annot = annFUN.gene2GO,
        gene2GO = go.ens,
        geneSelectionFun = geneSelFn,
        nodeSize = 8,
        description = ''
      )
      assign("root.godata", root.godata, envir=.GlobalEnv)
    }
    root.godata <- get("root.godata", envir=.GlobalEnv)

    if (!exists("go.pairs", envir=.GlobalEnv)) {
      go.pairs <- load.go.annotations()
      assign("go.pairs", go.pairs, envir=.GlobalEnv)
    }
    go.pairs <- get("go.pairs", envir=.GlobalEnv)
    print(head(go.pairs))

    #go.terms <- data.frame(go_id=c('GO:0000003', 'GO:0000018', 'GO:0000070', 'GO:0006458'), stringsAsFactors=F)
    go.terms <- data.frame(go_id=usedGO(root.godata), stringsAsFactors=F)

    scores <- binary.scores
    names(scores) <- y$protein_id
    geneSelFn <- function(score){return(score == 1)}
    fields.s <- lbl
    godata.s <- paste('godata_', fields.s, sep='')
    if (!exists(godata.s, envir=.GlobalEnv)) {
      go.data <- new("topGOdata",
        ontology = "BP",
        allGenes = scores,
        annot = annFUN.gene2GO,
        gene2GO = go.ens,
        geneSelectionFun = geneSelFn,
        nodeSize = 8,
        description = ''
      )
      assign(godata.s, go.data, envir=.GlobalEnv)
    }
    go.data <- get(godata.s, envir=.GlobalEnv)
    print(paste("N sig:", sum(geneSelFn(scores))))

    y$binary.score <- binary.scores
    y$scalar.score <- scalar.scores
    sub.y <- y[y$protein_id %in% genes(go.data),]
    total.genes <- length(genes(go.data))

    pwf.scores <- as.integer(y$binary.score)
    names(pwf.scores) <- y$protein_id
    pwf <- nullp(pwf.scores, bias.data=y$aln_length)
    print(nrow(go.pairs))
    go.pairs <- subset(go.pairs, protein.id %in% genes(go.data))
    go.pairs <- subset(go.pairs, go.id %in% usedGO(go.data))
    print(nrow(go.pairs))
    print(head(pwf))
    go.corr <- goseq(pwf, gene2cat=go.pairs)
    print(head(go.corr))
    go.uncorr <- goseq(pwf, gene2cat=go.pairs, method='Hypergeometric')
    print(head(go.uncorr))

    topgo.res <- runTest(go.data, algorithm="weight01", statistic="Fisher")
    topgo.pvals <- score(topgo.res)
    topgo.ids <- names(topgo.pvals)
    topgo.df <- data.frame(id=topgo.ids, pval=as.numeric(topgo.pvals), stringsAsFactors=F)

    topgo.fis.res <- runTest(go.data, algorithm="classic", statistic="Fisher")
    topgo.fis.pvals <- score(topgo.fis.res)
    topgo.fis.ids <- names(topgo.fis.pvals)
    topgo.fis.df <- data.frame(id=topgo.fis.ids, pval=as.numeric(topgo.fis.pvals), stringsAsFactors=F)

    all.distr <- sub.y[sub.y$scalar.score >= 0, 'scalar.score']
    all.pmax <- pmax(0, sub.y$scalar.score)

    sig.indices <- geneSelFn(sub.y$binary.score)

    go.summaries <- ddply(go.terms, 'go_id', function(x) {
      cur.id <- x[1, 'go_id']
      cur.id <- as.character(cur.id)
      cur.genes <- unlist(genesInTerm(go.data, cur.id))
      n.genes <- length(cur.genes)

      row.indices <- sub.y$protein_id %in% cur.genes
      cur.rows <- sub.y[sub.y$protein_id %in% cur.genes, ]
      sig.rows <- cur.rows[geneSelFn(cur.rows$binary.score), ]
      sig.rows <- sig.rows[order(sig.rows$scalar.score, decreasing=T),]

      sig.q <- quantile(sig.rows$scalar.score, c(0.25, 0.5, 0.75))
      all.q <- quantile(cur.rows$scalar.score, c(0.25, 0.5, 0.75))

      f.in.term <- sum(row.indices) / total.genes
      f.sig <- sum(sig.indices) / total.genes
      total.n.sig <- sum(sig.indices)      
      n.expected <- f.in.term * f.sig * total.genes
      n.sig <- sum(geneSelFn(cur.rows$binary.score))

      min.n <- 1

      fisher.pval <- 1
      # Only run FET for terms with > N significant genes.
      if (n.sig >= min.n) {
        a.factor <- factor(row.indices)
        b.factor <- factor(sig.indices)
        if (length(levels(a.factor)) < 2 || length(levels(b.factor)) < 2) {
          print(cur.id)
          fisher.pval <- 1
        } else {
          fisher.res <- fisher.test(a.factor, b.factor, alternative='g')
          fisher.pval <- fisher.res$p.value
        }
      }

      hyper.pval <- 1
      if (n.sig >= min.n) {
        n.white <- sum(row.indices)
        n.black <- length(row.indices) - sum(row.indices)
        n.balls.drawn <- sum(sig.indices)
        n.white.drawn <- sum(sig.indices & row.indices)
        hyper.pval <- phyper(n.white.drawn, n.white, n.black, n.balls.drawn, lower.tail=F)
      }


      mwu.pval <- 1
      ks.pval <- 1
      contains.all <- all(row.indices, TRUE)
      if (!contains.all) {
        cur.pmax <- pmax(0, sub.y[row.indices, 'scalar.score'])
        all.pmax <- pmax(0, sub.y[!(row.indices), 'scalar.score'])

        #if (nrow(above.zero) > 0) {
        #  cur.distr <- above.zero$scalar.score
        #  ks.res <- ks.test(cur.distr, all.distr, alternative='l')
        #  ks.pval <- ks.res$p.value
        #}

        if (max(cur.pmax) > 0) {
          mwu.res <- wilcox.test(cur.pmax, all.pmax, alternative='g')
          mwu.pval <- mwu.res$p.value
        }
      }

      mean.length <- mean(cur.rows$aln_length)
      goseq.corr.pval <- go.corr[go.corr$category == cur.id, ]$over
      goseq.uncorr.pval <- go.uncorr[go.uncorr$category == cur.id, ]$over
      topgo.pval <- topgo.df[topgo.df$id == cur.id, ]$pval
      topgo.fis.pval <- topgo.fis.df[topgo.fis.df$id == cur.id, ]$pval
      descr <- getTermsDefinition(cur.id, 'BP', numChar=40)
      genes.str <- paste(sig.rows$gene_name, collapse=' ')

      cur.df <- data.frame(
        'go_id'=cur.id,
        'mean_length' = mean.length,
        'total_sig' = total.n.sig,
        'n' = sum(row.indices),
        'n_sig' = n.sig,
        'n_exp' = n.expected,
        'fis_pval' = fisher.pval,
        'hyper_pval' = hyper.pval,
        'ks_pval' = ks.pval,
        'mwu_pval' = mwu.pval,
        'goseq_wall' = goseq.corr.pval,
        'goseq_fis' = goseq.uncorr.pval,
        'topgo_pval' = topgo.pval,
        'topgo_fis' = topgo.fis.pval,
        descr=descr,
        genes=genes.str,
        stringsAsFactors=F
      )
      cur.df
    })

    adj.f <- function(x, fld) {
      adj.s <- paste(fld, '_adj', sep='')
      #b.one <- x[, fld] < 1
      large.terms <- x[, 'n'] >= 20
      x[, adj.s] <- 1
      x[large.terms, adj.s] <- p.adjust(x[large.terms, fld], method='BH')
      x
    }

    go.summaries <- adj.f(go.summaries, 'fis_pval')
    #go.summaries <- adj.f(go.summaries, 'hyper_pval')
    #go.summaries <- adj.f(go.summaries, 'ks_pval')
    go.summaries <- adj.f(go.summaries, 'mwu_pval')
    go.summaries <- adj.f(go.summaries, 'goseq_wall')
    go.summaries <- adj.f(go.summaries, 'topgo_fis')
    go.summaries
  }

  t <- qchisq(0.95, df=1)

  hc.val <- pmin(y$lrt.1, y$lrt.2)
  gh.val <- pmin(y$lrt.5, y$lrt.1)
  gc.val <- pmin(y$lrt.5, y$lrt.2)
  hc.filt.val <- pmin(y$lrt.1, y$lrt.2) - pmax(0, y$lrt.5)
  gh.filt.val <- pmin(y$lrt.5, y$lrt.1) - pmax(0, y$lrt.2)
  gc.filt.val <- pmin(y$lrt.5, y$lrt.2) - pmax(0, y$lrt.1)
  g.val <- y$lrt.5
  h.val <- y$lrt.1
  c.val <- y$lrt.2
  aga.val <- y$lrt.6
  aga.clade.val <- y$lrt.7

  g.filt.val <- y$lrt.5
  g.filt.val[y$g.ns + y$g.s < 3] <- 0
  c.filt.val <- y$lrt.2
  c.filt.val[y$c.ns + y$c.s < 3] <- 0
  h.filt.val <- y$lrt.1
  h.filt.val[y$h.ns + y$h.s < 3] <- 0

  zero <- 0

  print("#####")
  print(val.string)
  print(sig.t)
  print("#####")

  cur.val <- get(val.string)

  source("~/src/greg-ensembl/scripts/mysql_functions.R")
  con <- connect('gj1_gorilla_go')

  cur.go <- go.f(cur.val, as.integer(cur.val > sig.t), val.string)
  cur.go$t <- sig.t
  cur.go$stat <- val.string
  cur.go$label <- paste(cur.go$go_id, sig.t, val.string, sep='_')

  if (dbExistsTable(con, 'go')) {
    print(head(cur.go))
    dbUpdateVars(con, 'go', cur.go, 'label')
  } else {
    dbWriteTable(con, 'go', cur.go, row.names=F)
    dbSendQuery(con, 'ALTER TABLE go ADD UNIQUE (label(64))')
  }
}

load.go.annotations <- function() {
  if (file.exists("go.pairs.Rdata")) {
    load("go.pairs.Rdata")
    return(go.pairs)
  }

  if (!exists("root.godata", envir=.GlobalEnv)) {
    go.rdata.file <- "go_60.Rdata"
    load(go.rdata.file)
    assign('go.ens', go.ens, envir=.GlobalEnv)
    scores <- y$lrt.9
    names(scores) <- y$protein_id
    geneSelFn <- function(score){return(score >= qchisq(0.95, df=1))}
    root.godata <- new("topGOdata",
      ontology = "BP",
      allGenes = scores,
      annot = annFUN.gene2GO,
      gene2GO = go.ens,
      geneSelectionFun = geneSelFn,
      nodeSize = 5,
      description = ''
    )
    assign("root.godata", root.godata, envir=.GlobalEnv)
  }
  root.godata <- get("root.godata", envir=.GlobalEnv)

  go.terms <- usedGO(root.godata)
  go.df <- ldply(go.terms, function(x) {
    cur.term <- x
    cur.genes <- as.character(unlist(genesInTerm(root.godata, cur.term)))
    data.frame(
      protein.id = cur.genes,
      go.id = cur.term
    )
  })
  go.df
}

overlap.excess <- function(length.corr=T) {
  x <- y

  excess.f <- function(a, b, c, lbl) {
    thresholds <- c(0, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, qchisq(0.95, df=1), 4, 4.5, 5)
    #thresholds <- c(1, 1.5, 2, qchisq(0.95, df=1))

    cur.df <- data.frame(a=a, b=b, c=c)

    res.df <- data.frame()
    for (i in 1:length(thresholds)) {
      cur.t <- thresholds[i]
      orig.t <- cur.t
      cur.pval <- pchisq(cur.t, df=1, lower.tail=F)

      if (length.corr) {
        genes.a <- as.integer(a > cur.t)
        names(genes.a) <- x$protein_id
        genes.b <- as.integer(b > cur.t)
        names(genes.b) <- x$protein_id
        pwf.a <- nullp(genes.a, bias.data=x$aln_length)
        pwf.b <- nullp(genes.b, bias.data=x$aln_length)
        cur.df$pwf.a <- pwf.a$pwf
        cur.df$pwf.b <- pwf.b$pwf

        #print(head(pwf.a))
        #print(sum(a > cur.t) / length(a))

        #f.a <- sum(a > cur.t) / length(a)
        #f.b <- sum(b > cur.t) / length(b)
        #a <- a / pwf.a$pwf
        #b <- b / pwf.b$pwf
        #cur.df$a <- a
        #cur.df$b <- b
        #cur.t.a <- cur.t / f.a
        #cur.t.b <- cur.t / f.b
        cur.t.a <- cur.t
        cur.t.b <- cur.t

        #print(sum(a > cur.t) / length(a))
      }
      cur.t.a <- cur.t
      cur.t.b <- cur.t

      boot.f <- function(D, d) {
        E <- D[d,]
        obs.n <- sum(E$a > cur.t.a & E$b > cur.t.b)
        if (length.corr) {
          exp.n <- sum(E$pwf.a * E$pwf.b)
        } else {
          exp.n <- sum(E$a > cur.t.a) / nrow(E) * sum(E$b > cur.t.b) / nrow(E) * nrow(E)
        }

        if (exp.n == 0 && obs.n == 0) {
          return(0)
        } else if (exp.n == 0 && obs.n > 0) {
          return(1)
        } else {
          return(obs.n / exp.n - 1)
        }
      }
      boot.res <- boot(cur.df, boot.f, R=100)
      ci.res <- boot.ci(boot.res, type='basic', conf=0.5)
      single.res <- boot.f(cur.df, 1:nrow(cur.df))
      bt <- c(ci.res$basic[1, 4], single.res, ci.res$basic[1, 5])

      random.t <- function(x, y, total, reps=1000) {
        n.overlaps <- c()
        for (i in 1:reps) {
          if (length.corr) {
            a.sample <- base::sample(1:total, x, replace=F, prob=pwf.a$pwf)
            b.sample <- base::sample(1:total, y, replace=F, prob=pwf.b$pwf)
          } else {
            a.sample <- base::sample(1:total, x, replace=F)
            b.sample <- base::sample(1:total, y, replace=F)
          }
            n.overlaps[i] <- length(intersect(a.sample, b.sample))
        }
        ecdf(n.overlaps)
      }
      random.ecdf <- random.t(sum(a > cur.t.a), sum(b > cur.t.b), length(a))
      random.res <- 1 - random.ecdf(sum(a > cur.t.a & b > cur.t.b))
      
      obs.n <- sum(a > cur.t & b > cur.t)
      if (length.corr) {
        exp.n <- sum(pwf.a$pwf * pwf.b$pwf)
      } else {
        exp.n <- sum(a > cur.t.a) / length(a) * sum(b > cur.t.b) / length(b) * length(a)
      }
      exp.n <- as.integer(exp.n)

      res.df <- rbind(res.df, data.frame(
        lbl = lbl,
        pval = cur.pval,
        t = orig.t,
        obs.n = obs.n,
        exp.n = exp.n,
        r = bt[2],
        r.lo = bt[1],
        r.hi = bt[3],
        random.p = random.res
      ))
    }
    print(res.df)
    res.df
  }

  if (length.corr) {
    excess.df.s <- 'excess.df.corr'
  } else {
    excess.df.s <- 'excess.df.nocorr'
  }

  if (!exists(excess.df.s, envir=.GlobalEnv)) {
    dfs <- data.frame()
    dfs <- rbind(dfs, excess.f(x$lrt.5, x$lrt.2, x$lrt.1, 'Gorilla-Chimpanzee'))
    dfs <- rbind(dfs, excess.f(x$lrt.5, x$lrt.1, x$lrt.2, 'Gorilla-Human'))
    dfs <- rbind(dfs, excess.f(x$lrt.1, x$lrt.2, x$lrt.5, 'Human-Chimpanzee'))
    assign(excess.df.s, dfs, envir=.GlobalEnv)
  }
  dfs <- get(excess.df.s, envir=.GlobalEnv)
  print(dfs)

  if (length.corr) {
    pdf(file=pdf.f("overlap.excess.pdf"), width=8, height=6)
  } else {
    pdf(file=pdf.f("overlap.excess.nocorr.pdf"), width=8, height=6)
  }
  p <- ggplot(dfs, aes(x=t, y=r, colour=lbl, group=lbl, fill=lbl))
  p <- p + theme_bw()
  p <- p + geom_line()
  p <- p + geom_ribbon(aes(ymin=r.lo, ymax=r.hi), alpha=0.3, colour=NA)
  p <- p + geom_vline(x=qchisq(0.95, df=1), linetype=2, size=0.3)
  p <- p + scale_colour_discrete("Species Pair")
  p <- p + scale_fill_discrete("Species Pair")

  p <- p + scale_x_continuous("Acceleration Cutoff Threshold (LRT)")
  p <- p + scale_y_continuous("Acceleration Co-occurrence Excess")

  sub.1 <- subset(dfs, random.p <= 0.1 & random.p > 0.05)
  sub.05 <- subset(dfs, random.p <= 0.05 & random.p > 0.01)
  sub.01 <- subset(dfs, random.p <= 0.01)

  p <- p + geom_text(data=sub.1, label='*', size=4, colour='gray50', vjust=0.7)
  p <- p + geom_text(data=sub.05, label='*', size=4, colour='black', vjust=0.7)
  p <- p + geom_text(data=sub.01, label='**', size=4, colour='black', vjust=0.7)

  print(p)
  dev.off()
}

length.effect <- function() {
  x <- y
#  x <- subset(x, aln_length > 100 & aln_length < 800)  

  ind.stat <- 'aln_length'
  dep.stat <- 'lrt.5'
  x$ind.stat <- x[, ind.stat]
  x$dep.stat <- x[, dep.stat]
  x$abs.stat <- abs(x$dep.stat)

  lm.res <- lm(abs.stat ~ ind.stat, data=x)
  print(lm.res$coefficients)

  mdn <- median(x$ind.stat)
  asdf <- 1
  #x$dep.stat <- pmax(0, x$abs.stat + asdf - asdf * x$ind.stat / mdn) * sign(x$dep.stat)
  #x$dep.stat <- pmax(0, (1-asdf) * x$abs.stat + asdf * x$abs.stat * mdn / x$ind.stat) * sign(x$dep.stat)

  pdf(file="temp.pdf")
  p <- ggplot(x, aes(y=abs(dep.stat)+2.7, x=ind.stat))
  p <- p + geom_point(size=0.5, alpha=0.5)
  p <- p + scale_y_log10() + scale_x_log10()
  print(p)
  dev.off()

  thresh <- qchisq(0.95, df=1)

  n.quantiles <- 11
  length.qs <- quantile(x$ind.stat, seq(from=0, to=1, length.out=n.quantiles))
  length.qs[1] <- length.qs[1] - 0.001

  out.df <- data.frame()
  for (i in 1:(length(length.qs)-1)) {
    lo <- length.qs[i]
    hi <- length.qs[i+1]

    cur.rows <- subset(x, ind.stat > lo & ind.stat <= hi)
    n.rows <- nrow(cur.rows)
    n.sig <- nrow(subset(cur.rows, abs(dep.stat) > thresh))
    f.sig <- n.sig / n.rows

    mean.abs <- mean(abs(cur.rows$dep.stat))
    mean.lrt <- mean(cur.rows$dep.stat)

    n.below <- nrow(subset(cur.rows, dep.stat < 0))
    n.above <- nrow(subset(cur.rows, dep.stat > 0))

    out.df <- rbind(out.df, data.frame(
      lo = lo,
      hi = hi,
      mean = sprintf("%.2f",mean(cur.rows$ind.stat)),
      f.sig = f.sig,
      mean.abs = mean.abs,
      mean.lrt = mean.lrt,
      f.above = n.above / n.rows,
      f.below = n.below / n.rows
    ))
  }

  print(out.df)
}

go.sweeps <- function() {

  cur.id <- 'GO:0007265'
  len <- y$aln_length
  len <- 1

  a <- go.sweep.f(y, pmin(y$lrt.1, y$lrt.5) / len, cur.id)
  print(a)
  a <- go.sweep.f(y, y$lrt.1 / len, cur.id)
  #print(a)
  a <- go.sweep.f(y, y$lrt.2 / len, cur.id)
  #print(a)
  a <- go.sweep.f(y, y$lrt.5 / len, cur.id)
  #print(a)

}

go.sweep.f <- function(y, statistic, ids) {
  rank.cutoffs <- c(25, 50, seq(from=100, to=1500, by=100))
  out.df <- data.frame()

  if (!exists("sweep.godata", envir=.GlobalEnv)) {
    scores <- statistic
    names(scores) <- y$protein_id
    geneSelFn <- function(score){return(score >= qchisq(0.95, df=1))}
    sweep.godata <- new("topGOdata",
      ontology = "BP",
      allGenes = scores,
      annot = annFUN.gene2GO,
      gene2GO = go.ens,
      geneSelectionFun = geneSelFn,
      nodeSize = 5,
      description = ''
    )
    assign("sweep.godata", sweep.godata, envir=.GlobalEnv)
  }
  go.data <- get('sweep.godata', envir=.GlobalEnv)

  y$score <- statistic
  sub.y <- y[y$protein_id %in% genes(go.data),]
  sub.y <- sub.y[order(sub.y$score, decreasing=T), ]
  sub.y$rank <- 1:nrow(sub.y)
  total.genes <- length(genes(go.data))

  for (j in 1:length(ids)) {
    id <- ids[j]

    cur.genes <- unlist(genesInTerm(go.data, id))
    row.indices <- sub.y$protein_id %in% cur.genes
    a.factor <- factor(row.indices)  
    n.in.term <- sum(row.indices)
    all.lengths <- sub.y[row.indices, 'aln_length']

    for (i in 1:length(rank.cutoffs)) {
      cur.t <- rank.cutoffs[i]
      sig.indices <- sub.y$rank <= cur.t
      n.sig <- sum(sig.indices)
      min.n <- 1

      f.in.term <- sum(row.indices) / total.genes
      f.sig <- sum(sig.indices) / total.genes
      n.expected <- f.in.term * f.sig * total.genes
      n.expected <- sprintf("%.1f", n.expected)
      n.sig <- sum(row.indices & sig.indices)

      fisher.pval <- 1
      # Only run FET for terms with > N significant genes.
      if (n.sig >= min.n) {
        b.factor <- factor(sig.indices)
        if (length(levels(a.factor)) < 2 || length(levels(b.factor)) < 2) {
          print(cur.id)
          fisher.pval <- 1
        } else {
          fisher.res <- fisher.test(a.factor, b.factor, alternative='g')
          fisher.pval <- fisher.res$p.value
        }
      }

      sig.lengths <- sub.y[row.indices & sig.indices, 'aln_length']
      mean.len <- as.integer(mean(all.lengths))
      mean.sig.len <- as.integer(mean(sig.lengths))
      mwu.p <- 1
      if (is.finite(mean.sig.len)) {
        mwu.res <- wilcox.test(all.lengths, sig.lengths, alternative='l')
        mwu.p <- mwu.res$p.value
      }

      out.df <- rbind(out.df, data.frame(
        id = id,
        t = cur.t,
        pval = as.numeric(sprintf("%.3f",fisher.pval)),
        logp = as.numeric(sprintf("%.3f", -log(fisher.pval))),
        n.sig = n.sig,
        n.exp = n.expected,
        n.in.term = n.in.term,
        mean.len = mean.len,
        sig.len = mean.sig.len,
        len.mwu = mwu.p
      ))
    }
  }

  out.df
}

go.enrichments <- function(x) {
  if (!exists('go.ens')) {
    go.rdata.file <- "~/scratch/gorilla/go_60.Rdata"
    load(go.rdata.file)
    assign('go.ens', go.ens, envir=.GlobalEnv)
  }
  source("~/src/greg-ensembl/scripts/go_enrichments.R")

  y <- x

  top <- floor(nrow(y)/20)
  hi <- floor(nrow(y)/5)
  half <- floor(nrow(y)/2)

  bot <- nrow(y) - top

  t <- qchisq(0.95, df=1)

  do.parallel <- F
  if (do.parallel) {

    cur <- subset(y, pmin(lrt.1, lrt.5) > 1.5 & lrt.2 < pmin(lrt.1, lrt.5))
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=TRUE)
    assign('overlap.gh', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(overlap.gh, pdf.f('overlap.gh.csv'))

    cur <- subset(y, pmin(lrt.2, lrt.5) > 1.5 & lrt.1 < pmin(lrt.2, lrt.5))
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=TRUE)
    assign('overlap.gc', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(overlap.gc, pdf.f('overlap.gc.csv'))

    cur <- subset(y, pmin(lrt.1, lrt.2) > 1.5 & lrt.5 < pmin(lrt.1, lrt.2))
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=TRUE)
    assign('overlap.hc', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(overlap.hc, pdf.f('overlap.hc.csv'))
  }

  do.branch.parallel <- F
  if (do.branch.parallel) {
    cur <- subset(y, rank.lrt.8 <= top & rank.lrt.5 > top)
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=TRUE)
    assign('unique.hi.hc', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(unique.hi.hc, pdf.f('unique.hi.hc.csv'))

    cur <- subset(y, rank.lrt.9 <= top & rank.lrt.2 > top)
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=TRUE)
    assign('unique.hi.gh', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(unique.hi.gh, pdf.f('unique.hi.gh.csv'))

    cur <- subset(y, rank.lrt.10 <= top & rank.lrt.1 > top)
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=TRUE)
    assign('unique.hi.gc', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(unique.hi.gc, pdf.f('unique.hi.gc.csv'))
  }

  ### Standard enrichments.

  do.standard <- F
  if (do.standard) {
    cur <- subset(y, rank.lrt.1 <= top)
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=F)
    assign('hi.human', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(hi.human, pdf.f('hi.human.csv'))

    cur <- subset(y, rank.lrt.2 <= top)
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=F)
    assign('hi.chimpanzee', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(hi.chimpanzee, pdf.f('hi.chimpanzee.csv'))

    cur <- subset(y, rank.lrt.5 <= top)
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=F)
    assign('hi.gorilla', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(hi.gorilla, pdf.f('hi.gorilla.csv'))

    cur <- subset(y, rank.lrt.6 <= top)
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=F)
    assign('hi.aga.branch', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(hi.aga.branch, pdf.f('hi.aga.branch.csv'))

    cur <- subset(y, rank.lrt.7 <= top)
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=F)
    assign('hi.aga.clade', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(hi.aga.clade, pdf.f('hi.aga.clade.csv'))

    cur <- subset(y, rank.lrt.8 <= top)
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=F)
    assign('hi.hc', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(hi.hc, pdf.f('hi.hc.csv'))

    cur <- subset(y, rank.lrt.9 <= top)
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=F)
    assign('hi.gh', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(hi.gh, pdf.f('hi.gh.csv'))

    cur <- subset(y, rank.lrt.10 <= top)
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=F)
    assign('hi.gc', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(hi.gc, pdf.f('hi.gc.csv'))

    ### Unique enrichments.
    cur <- subset(y, rank.lrt.1 <= top & rank.lrt.2 > top & rank.lrt.5 > top)
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=TRUE)
    assign('unique.hi.human', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(unique.hi.human, pdf.f('unique.hi.human.csv'))

    cur <- subset(y, rank.lrt.2 <= top & rank.lrt.1 > top & rank.lrt.5 > top)
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=TRUE)
    assign('unique.hi.chimpanzee', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(unique.hi.chimpanzee, pdf.f('unique.hi.chimpanzee.csv'))

    cur <- subset(y, rank.lrt.5 <= top & rank.lrt.1 > top & rank.lrt.2 > top)
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=TRUE)
    assign('unique.hi.gorilla', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(unique.hi.gorilla, pdf.f('unique.hi.gorilla.csv'))

    cur <- subset(y, rank.lrt.6 <= top)
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=TRUE)
    assign('unique.hi.aga.branch', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(unique.hi.aga.branch, pdf.f('unique.hi.aga.branch.csv'))

    cur <- subset(y, rank.lrt.7 <= top & rank.lrt.1 > top & rank.lrt.2 > top & rank.lrt.5 > top)
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=TRUE)
    assign('unique.hi.aga.clade', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(unique.hi.aga.clade, pdf.f('unique.hi.aga.clade.csv'))

    cur <- subset(y, rank.lrt.8 <= top & rank.lrt.5 > top)
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=TRUE)
    assign('unique.hi.hc', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(unique.hi.hc, pdf.f('unique.hi.hc.csv'))

    cur <- subset(y, rank.lrt.9 <= top & rank.lrt.2 > top)
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=TRUE)
    assign('unique.hi.gh', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(unique.hi.gh, pdf.f('unique.hi.gh.csv'))

    cur <- subset(y, rank.lrt.10 <= top & rank.lrt.1 > top)
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=TRUE)
    assign('unique.hi.gc', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(unique.hi.gc, pdf.f('unique.hi.gc.csv'))

  }

  do.corrected = T
  if (do.corrected) {
    
  }


  if (do.extras) {
    ### Hi-masked nucs.
    cur <- subset(y, masked_nucs / aln_length > 0.2)
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=TRUE)
    assign('hi.masked.nucs', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)
    output.tbl(hi.masked.nucs, pdf.f('hi.masked.nucs.csv'))

    ### Super-unique LOW
    cur <- subset(y, rank.lrt.1 >= bot & rank.lrt.2 < bot & rank.lrt.5 < bot)
    print(nrow(cur))
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=TRUE)
    assign('lo.super.unique.human', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)   
    output.tbl(lo.super.unique.human, pdf.f('lo.super.unique.human.csv'))

    cur <- subset(y, rank.lrt.5 >= bot & rank.lrt.2 < bot & rank.lrt.1 < bot)
    print(nrow(cur))
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=TRUE)
    assign('lo.super.unique.gorilla', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)   
    output.tbl(lo.super.unique.gorilla, pdf.f('lo.super.unique.gorilla.csv'))

    cur <- subset(y, rank.lrt.2 >= bot & rank.lrt.1 < bot & rank.lrt.5 < bot)
    print(nrow(cur))
    tbl <- enrich.list(cur$protein_id, y$protein_id, include.iea=TRUE)
    assign('lo.super.unique.chimpanzee', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)   
    output.tbl(lo.super.unique.chimpanzee, pdf.f('lo.super.unique.chimpanzee.csv'))
  }

}

cum.dist.plots <- function(y, fields=c('lrt.5', 'lrt.2', 'lrt.1'), x.lim=NA) {

  cum.df <- data.frame()
  for (field in fields) {
    cur.y <- y[order(y[, field]),]
    cur.y$v <- cur.y[, field]
    cur.y <- subset(cur.y, !is.na(v))
    xs <- cur.y[, field]
    ys <- 1:nrow(cur.y) / nrow(cur.y)
    cum.df <- rbind(cum.df, data.frame(
      lbl=field,
      xs=xs,
      ys=ys
    ))
  }

  p <- ggplot(cum.df, aes(x=xs, y=ys, colour=lbl, group=lbl))
  p <- p + geom_line()
  if (!is.na(x.lim)) {
    p <- p + scale_x_continuous(limits=x.lim)
  }
  return(p)
}

top.n.plot <- function(x, dir='up') {
  x <- subset(x, !is.na(lrt.5) & !is.na(lrt.4) & !is.na(lrt.3))

  x$h.f.nsyn <- x$h.f.nsyn + x$hc.f.nsyn
  x$c.f.nsyn <- x$c.f.nsyn + x$hc.f.nsyn

  pairs.list <- list(
    'G-H LRT' = c('lrt.5', 'lrt.3'),
#    'G-H LRT2' = c('lrt.5', 'lrt.1'),
#    'G-C LRT2' = c('lrt.5', 'lrt.2')
    'G-C LRT' = c('lrt.5', 'lrt.4')
#    'H-C LRT' = c('lrt.2', 'lrt.1')
  )

  lbls <- names(pairs.list)
  df <- data.frame()
  for (i in 1:length(pairs.list)) {
    nm <- lbls[i]
    print(nm)
    cur.pair <- pairs.list[[i]]
    
    fld.a <- cur.pair[1]
    fld.b <- cur.pair[2]

    if (dir == 'down') {
      x$rnk.a <- rank(x[,fld.a])
      x$rnk.b <- rank(x[,fld.b])
    } else {
      x$rnk.a <- rank(-x[,fld.a])
      x$rnk.b <- rank(-x[,fld.b])
    }
    
    scores <- c()
    for (j in 1:nrow(x)) {
      cur.score <- sum(x$rnk.a <= j & x$rnk.b <= j)
      scores[j] <- cur.score
      if (j == 200) print(paste(nm,j,cur.score))
    }

    df <- rbind(df, data.frame(
      fld = nm,
      xx = 1:nrow(x),
      yy = scores
    ))
  }

  p <- ggplot(df, aes(x=xx, y=yy, group=fld, colour=fld))
  p <- p + theme_bw()
  p <- p + geom_line()
#  p <- p + scale_x_continuous(limits=c(0, 1000))
#  p <- p + scale_y_continuous(limits=c(0, 600))
   return(p)
}

top.n.plots <- function(x) {
  pdf(file=pdf.f("shared.n.pdf"), width=8, height=4)

  p1 <- top.n.plot(x, dir='up')
  p2 <- top.n.plot(x, dir='down')
  
  vplayout(2,2)

  print(p1, vp=subplot(1,1))
  print(p2, vp=subplot(2,1))
  
  p1 <- p1 + scale_x_continuous(limits=c(0, 500))
  p1 <- p1 + scale_y_continuous(limits=c(0, 120))
  p2 <- p2 + scale_x_continuous(limits=c(0, 500))
  p2 <- p2 + scale_y_continuous(limits=c(0, 120))

  print(p1, vp=subplot(1,2))
  print(p2, vp=subplot(2,2))
  
  dev.off()
}

shared.accel.table <- function(x, dir='up') {
  # Q: What proportion of gor accels are also hum accels for a given threshold?

  x <- subset(x, !is.na(lrt.5) & !is.na(lrt.4) & !is.na(lrt.3))

  if (dir=='down') {
    x$lrt.g <- -x$rank.lrt.5
    x$lrt.c <- -x$rank.lrt.4
    x$lrt.h <- -x$rank.lrt.3
  } else {
    x$lrt.g <- -x$rank.lrt.5
    x$lrt.c <- -x$rank.lrt.4
    x$lrt.h <- -x$rank.lrt.3
  }

  dt <- x[, c('lrt.g', 'lrt.c', 'lrt.h', 'gene_name')]
  g <- as.data.table(dt)
  print(head(g))

  all <- nrow(g)
  print("  calculating...")
  print(nrow(g))
  x.up <- data.frame()
  if (dir == 'up') {
    indices <- order(-g$lrt.g)
  } else {
    indices <- order(g$lrt.g)
  }
  for (i in 1:length(indices)) {
    cur.index <- indices[i]
    cur.threshold <- g[cur.index, lrt.g]
    print(cur.threshold)

    if (dir == 'up') {
      n.g <- nrow(g[lrt.g >= cur.threshold,])
      n.gh <- nrow(g[lrt.g >= cur.threshold & lrt.h >= cur.threshold,])
      n.gc <- nrow(g[lrt.g >= cur.threshold & lrt.c >= cur.threshold,])
    } else {
      n.g <- nrow(g[lrt.g <= cur.threshold,])
      n.gh <- nrow(g[lrt.g <= cur.threshold & lrt.h <= cur.threshold,])
      n.gc <- nrow(g[lrt.g <= cur.threshold & lrt.c <= cur.threshold,])
    }

    x.up[i, 'i'] <- i
    x.up[i, 'value'] <- cur.threshold
    x.up[i, 'n.chimp'] <- n.gc / all
    x.up[i, 'n.human'] <- n.gh / all
    x.up[i, 'f.chimp'] <- n.gc / n.g
    x.up[i, 'f.human'] <- n.gh / n.g
  }
  print("  done!")
  
  df <- data.frame(grp='GH / G_tot', yy=x.up$n.human, xx=x.up$value, xi=x.up$i)
  df <- rbind(df, data.frame(grp='GC / G_tot', yy=x.up$n.chimp, xx=x.up$value, xi=x.up$i))
  df <- rbind(df, data.frame(grp='GH / G_cur', yy=x.up$f.human, xx=x.up$value, xi=x.up$i))
  df <- rbind(df, data.frame(grp='GC / G_cur', yy=x.up$f.chimp, xx=x.up$value, xi=x.up$i))

  p <- ggplot(df, aes(x=xi, y=yy, group=grp, colour=grp))
  p <- p + geom_line(size=0.5)
  clrs <- c('#FF0000', '#0000FF', '#FFBBBB', '#BBBBFF')
  szs <- c(1, 0.5, 1, 0.5, 1, 1)
  p <- p + scale_colour_manual(values=clrs)
  p <- p + scale_size_manual(values=szs)

  if (dir == 'down') {
    pdf(file=pdf.f("shared_down.pdf"))
  } else {
    pdf(file=pdf.f("shared_up.pdf"))
  }
  print(p)
  dev.off()
}

pdf.f <- function(f) {
  return(paste("~/src/greg-ensembl/projects/gorilla/", f, sep=''))
}

ns.cum.plot <- function(x) {

  pdf(file="~/src/greg-ensembl/projects/gorilla/hgc.ns.cum.pdf", height=3, width=4)
  #x$h.f.nsyn <- x$h.f.nsyn + x$hc.f.nsyn
  #x$c.f.nsyn <- x$c.f.nsyn + x$hc.f.nsyn
  p <- cum.dist.plots(x, fields=c('g.f.nsyn', 'h.f.nsyn', 'c.f.nsyn'), x.lim=c(0, 0.04))
  p <- p + theme_bw()
  print(p)
  dev.off()

  pdf(file="~/src/greg-ensembl/projects/gorilla/hgc.lrt.cum.pdf", height=3, width=4)
  x$Gorilla <- x$lrt.5
  x$Human <- x$lrt.1
  x$Chimpanzee <- x$lrt.2
  p <- cum.dist.plots(x, fields=c('Human', 'Chimpanzee', 'Gorilla'), x.lim=c(-5, 5))
  p <- p + theme_bw()
  print(p)
  dev.off()

  pdf(file="~/src/greg-ensembl/projects/gorilla/aga.lrt.cum.pdf", height=4, width=6)
  x$AGA.Branch <- x$lrt.6
  x$AGA.Clade <- x$lrt.7
  p <- cum.dist.plots(x, fields=c('AGA.Branch', 'AGA.Clade'), x.lim=c(-10, 10))
  p <- p + theme_bw()
  print(p)
  dev.off()

  pdf(file="~/src/greg-ensembl/projects/gorilla/hc.lrt.cum.pdf", height=4, width=6)
  x$Human <- x$lrt.1
  x$Chimpanzee <- x$lrt.2
  p <- cum.dist.plots(x, fields=c('Human', 'Chimpanzee'), x.lim=c(-10, 10))
  p <- p + theme_bw()
  print(p)
  dev.off()

  pdf(file="~/src/greg-ensembl/projects/gorilla/parallel.lrt.cum.pdf", height=4, width=6)
  x$HC_Parallel <- x$lrt.8
  x$GH_Parallel <- x$lrt.9
  x$GC_Parallel <- x$lrt.10
  p <- cum.dist.plots(x, fields=c('HC_Parallel', 'GH_Parallel', 'GC_Parallel'), x.lim=c(-10, 10))
  p <- p + theme_bw()
  print(p)
  dev.off()

  df.g <- x
  df.g$value <- df.g$g.f.nsyn
  df.g$spc <- 'Gorilla'
  
  df.h <- x
  df.h$value <- df.h$h.f.nsyn
  df.h$spc <- 'Human'

  df.c <- x
  df.c$value <- df.c$c.f.nsyn
  df.c$spc <- 'Chimpanzee'
  
  df.g <- subset(df.g, h.f.nsyn <= 0.05 & g.f.nsyn <= 0.05)
  df.g <- subset(df.g, h.f.nsyn > 0 & g.f.nsyn > 0)

  l1 <- lm(g.f.nsyn ~ h.f.nsyn + 0, data=subset(df.g, m0_dnds < 0.05))
  print(l1)
  l2 <- lm(g.f.nsyn ~ h.f.nsyn + 0, data=subset(df.g, m0_dnds > 0.3))
  print(l2)

  l3 <- lm(g.f.nsyn ~ c.f.nsyn + 0, data=subset(df.g, m0_dnds < 0.05))
  print(l3)
  l4 <- lm(g.f.nsyn ~ c.f.nsyn + 0, data=subset(df.g, m0_dnds > 0.3))
  print(l4)

  data <- rbind(df.g, df.h, df.c)
  p <- ggplot(df.g, aes(x=h.f.nsyn, y=value, colour=spc, group=spc))
  p <- p + geom_point(alpha=0.2)
  pdf(file="~/src/greg-ensembl/projects/gorilla/ns.dnds.pdf", height=6, width=6)
  print(p)
  dev.off()

}


go.picked.groups <- function(x) {
  if (!exists('go.ens')) {
    go.rdata.file <- "~/scratch/gorilla/go_60.Rdata"
    load(go.rdata.file)
  }

  if (!exists("generic.godata")) {
    ids <- x$protein_id
    scores <- rep(1,times=nrow(x))
    names(scores) <- ids
    geneSelectionFun = function(score){return(score >= 1)}
    GOdata <- new("topGOdata",
      ontology = 'BP',
      allGenes = scores,
      annot = annFUN.gene2GO,
      gene2GO = go.ens,
      geneSelectionFun = geneSelectionFun,
      nodeSize = 5,
      description = 'asdf'
    )
    assign("generic.godata", GOdata, envir=.GlobalEnv)
  }
  term.groups <- get.term.groups()
  term.groups['all'] <- 'all'

  rank.fields <- c(
     'slr_dnds', 'm0_dnds', 'lrt.1', 'lrt.2', 'lrt.5', 'all'
#    'all', 'slr_dnds', 'm0_dnds', 'lrt.1', 'lrt.2', 'lrt.5', 'lrt.6', 'lrt.7', 'lrt.8', 'lrt.9', 'lrt.10'
  )

  y <- x
  for (i in 1:length(rank.fields)) {
    cur.fld <- rank.fields[i]
    if (cur.fld == 'all') {
      next
    }
    cur.rnk.fld <- paste('rank.', cur.fld, sep='')
#    y[, cur.rnk.fld] <- y[, cur.fld]
    y[, cur.rnk.fld] <- rank(y[, cur.fld], ties.method='random')
    rank.fields[i] <- cur.rnk.fld
  }  

  print(rank.fields)

  terms.df <- ldply(term.groups, function(term) {
#    print(term)
    if (term[1] != 'all') {
      genes <- getGenesForTerms(generic.godata, term)
      row.indices <- y$protein_id %in% genes
      z <- y[row.indices,]
    } else {
      z <- y
    }
    df <- data.frame(nrow=nrow(z))

    data.for.field <- function(fld, all=F) {
#      print(fld)
      if (all) {
        df <- y
      } else {
        df <- z
      }

      df$rank_tmp <- df[, fld]
      return(df)
    }

    for (field.x in rank.fields) {
      if (field.x == 'all'){
        next 
      }
#      print(field.x)
      df.x <- data.for.field(field.x)
      for (field.y in rank.fields) {
        if (field.x == field.y) {
          next
        }
        if (field.y != 'all') {
          next
        }
#        print(field.y)
        if (field.y == 'all') {
          df.y <- data.for.field(field.x, all=T)
        } else {
          df.y <- data.for.field(field.y)
        }

#        print(summary(df.x$rank_tmp))
#        print(summary(df.y$rank_tmp))

        if (nrow(z) > 0) {
          ks.up <- ks.test(df.x$rank_tmp, df.y$rank_tmp, alternative='l')
          ks.down <- ks.test(df.x$rank_tmp, df.y$rank_tmp, alternative='g')
          ws.up <- wilcox.test(df.x$rank_tmp, df.y$rank_tmp, alternative='g')
          ws.down <- wilcox.test(df.x$rank_tmp, df.y$rank_tmp, alternative='l')
          ks.up.p <- ks.up$p.value
          ks.down.p <- ks.down$p.value
          ws.up.p <- ws.up$p.value
          ws.down.p <- ws.down$p.value
        } else {
          ks.up.p <- 1
          ks.down.p <- 1
          ws.up.p <- 1
          ws.down.p <- 1
        }
        field.up <- paste(field.x, '-', field.y, '.up', sep='')
        field.down <- paste(field.x, '-', field.y, '.down', sep='')
#        df[field.up] <- ws.up.p
        df[field.down] <- ws.down.p
#        df[field.up] <- ks.up.p
#        df[field.down] <- ks.down.p
      }
    }
    return(df)
  })
  print(terms.df)
  return()
}

output.tbl <- function(tbl, file) {
  write.csv(
    orderBy(~pval.topgo, 
    subset(tbl, Annotated < 500 & Significant >= 5)), 
    file=file)
}

process.tbl <- function(genes, godata, tbl) {
  tbl <- tbl[order(tbl$pval.topgo),]
  tbl <- subset(tbl, pval.fis <= 0.05 & pval.topgo < 1)
  
  row.names <- rownames(tbl)
  for (i in 1:nrow(tbl)) {
    cur.nm <- row.names[i]
    gene.ids <- genesInTerm(godata, cur.nm)[[1]]
    scores <- scoresInTerm(godata, cur.nm)[[1]]
    gene.ids <- gene.ids[scores== 1]
    gene.names <- as.character(id.to.name(genes, gene.ids))
    if (length(gene.names) > 10) {
      gene.names <- head(gene.names, n=10)
      gene.names <- c(gene.names, '...')
    }
    tbl[i, 'genes'] <- paste(gene.names, collapse=' ')
  }
  return(tbl)
}

id.to.name <- function(x, ids) {
#  x <- x[order(x[, sort.stat]),]
  indices <- x$protein_id %in% ids
  return(x[indices, 'gene_name'])
}


dnds.figure <- function(x, dnds.field='m0_dnds', gc.field='gc_3', 
  gc.quantile.count=1, dnds.quantile.count=1, boot.reps=10) {

  x[, 'gc_tmp'] <- x[, gc.field]
  gc.q <- quantile(x$gc_tmp, c(0.33, 0.67, 1.0))
  gc.q <- c(0, gc.q)
  print(gc.q)
  gc.lo <- gc.q[-length(gc.q)]
  gc.hi <- gc.q[-1]

  if (gc.quantile.count == 1) {
    gc.lo <- min(gc.lo)
    gc.hi <- max(gc.hi)
  }

  print(gc.lo)
  print(gc.hi)

  x[, 'dnds_tmp'] <- x[, dnds.field]
  if (dnds.quantile.count == 3) {
    dnds.q <- quantile(x$dnds_tmp, c(0.33, 0.67, 1))
    dnds.labels <- c('Low', 'Medium', 'High')
  } else if (dnds.quantile.count == 2) {
    dnds.q <- c(0.7, max(x$dnds_tmp))
    dnds.labels <- c('Low', 'High')
  } else {
    dnds.q <- c(1)
    dnds.labels <- c('All genes')
  }
  dnds.q <- c(0, dnds.q)
  dnds.lo <- dnds.q[-length(dnds.q)]
  dnds.hi <- dnds.q[-1]

  x[, 'all.ns'] <- x$h.ns + x$c.ns + x$hc.ns + x$g.ns + x$hcg.ns + x$o.ns + x$hcgo.ns + x$m.ns + x$hcgom.ns + x$r.ns
  x[, 'all.s'] <- x$h.s + x$c.s + x$hc.s + x$g.s + x$hcg.s + x$o.s + x$hcgo.s + x$m.s + x$hcgom.s + x$r.s

  fld.pairs <- list(
    'Human' = c('h.ns', 'h.s'),
    'Chimpanzee' = c('c.ns', 'c.s'),
    'H/C ancestor' = c('hc.ns', 'hc.s'),
    'Gorilla' = c('g.ns', 'g.s'),
    'H/C/G ancestor' = c('hcg.ns', 'hcg.s'),
    'Orangutan' = c('o.ns', 'o.s'),
    'H/C/G/O ancestor' = c('hcgo.ns', 'hcgo.s'),
    'Macaque' = c('m.ns', 'm.s'),
    'H/C/G/O/M ancestor' = c('hcgom.ns', 'hcgom.s'),
    'Marmoset' = c('r.ns', 'r.s'),
    'All' = c('all.ns', 'all.s')
  )
  fld.names <- names(fld.pairs)

  comb.df <- data.frame()
  for (i in 1:length(dnds.lo)) {
    dlo <- dnds.lo[i]
    dhi <- dnds.hi[i]
    x1 <- subset(x, dnds_tmp > dlo & dnds_tmp <= dhi)
    print(nrow(x1))

    dnds.labels[i] <- paste(dnds.labels[i], ' (', base::format(dlo, digits=2), ' - ', base::format(dhi, digits=2), ')', '\n n=', nrow(x1), sep='')

    for (j in 1:length(gc.lo)) {
      gclo <- gc.lo[j]
      gchi <- gc.hi[j]
     
      x2 <- subset(x1, gc_tmp > gclo & gc_tmp <= gchi)
      gc.lbl <- paste('GC content', ': ', base::format(gclo, digits=2), ' - ', base::format(gchi, digits=2), sep='')

      for (k in 1:length(fld.pairs)) {
        cur.pair <- fld.pairs[[k]]
        cur.nm <- fld.names[k]
        ns.fld <- cur.pair[1]
        s.fld <- cur.pair[2]

        all.ns <- sum(x2[, ns.fld])
        all.s <- sum(x2[, s.fld])

        #x2 <- x2[x2[, s.fld] > 0 & x2[, ns.fld] > 0,]

        sum.ratio <- function(D, d) {
          E <- D[d,]
          return(sum(E[, ns.fld]) / sum(E[, s.fld]) * .29 / .71)
        }

        boot.ratio <- boot(x2, sum.ratio, R=boot.reps)
        ci.ratio <- boot.ci(boot.ratio, type='basic')
        d.mean <- sum.ratio(x2, 1:nrow(x2))
        bt <- c(ci.ratio$basic[1, 4], d.mean, ci.ratio$basic[1, 5])
        print(paste(cur.nm, format(bt, digits=3), collapse=' '))

        real.mean.dnds <- mean(x2$dnds_tmp)

        cur.df <- data.frame(
          dnds = i,
          codon.m0.dnds = real.mean.dnds,
          n.genes = nrow(x2),
          tot.subs = sum(x[, ns.fld]) + sum(x[, s.fld]),
          n.ns = sum(x2[, ns.fld]),
          n.s = sum(x2[, s.fld]),
          gc = gc.lbl,
          species = cur.nm,
          mean = bt[2],
          lo = max(0, bt[1]),
          hi = bt[3]
        )
        if (cur.nm == 'All') {
          cur.df$tot.subs <- 0
        }        
        comb.df <- rbind(comb.df, cur.df)
      }
    }
  }
  
  print(dnds.labels)
  comb.df$dnds <- as.factor(comb.df$dnds)
  comb.df$max.subs <- max(comb.df$tot.subs)
  print(comb.df)

  p <- ggplot(comb.df, aes(
    x=species,
    xmin=as.numeric(species)-0.3,
    xmax=as.numeric(species)+0.3,    
    ymin=lo, ymax=hi, y=mean,
    fill=dnds
  ))
  p <- p + geom_hline(aes(yintercept=codon.m0.dnds), colour='gray', alpha=0.7, linetype='dashed')
  p <- p + geom_rect(width=0.75, alpha=0.3, position='identity')
  p <- p + geom_rect(fill='black', aes(ymin=mean-0.001, ymax=mean+0.001), width=0.75, alpha=0.6, position='identity')

  p <- p + geom_rect(data=comb.df,
    aes(xmin=as.numeric(species)-0.3, xmax=as.numeric(species)+0.3,
      ymin=-0.2, ymax=-0.2+(tot.subs/max.subs)*0.1),
    fill='gray', width=0.75, alpha=1.0, position='identity')

    breaks <- c(-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
    labels <- c(0, comb.df[1, 'max.subs'], 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)

  p <- p + scale_y_continuous("dN/dS Approximation", breaks=breaks, labels=labels)
#  p <- p + scale_y_continuous("dN/dS Approximation", limits=c(-0.2, max(comb.df$hi)))
  p <- p + scale_x_discrete("Species / ancestral lineage")

  if (dnds.quantile.count > 1) {
    p <- p + scale_fill_manual("dN/dS bin", breaks=1:length(dnds.labels),
    labels = dnds.labels,
  #    labels=c("Low", "Medium", "High"),
      values=c('blue', 'black', 'red')
    )
    p <- p + scale_colour_manual("dN/dS bin", breaks=1:length(dnds.labels),
    labels = dnds.labels,
  #    labels=c("Low", "Medium", "High"),
      values=c('blue', 'black', 'red')
    )
  }
  
  p <- p + facet_grid(. ~ gc)

  p <- p + theme_bw()
  p <- p + opts(
    axis.text.x = theme_text(angle=90, hjust=1)
  )
  if (gc.quantile.count == 1) {
    pdf(file=pdf.f("dnds.figure.pdf"), width=5, height=4)
  } else {
    pdf(file=pdf.f("dnds.figure.pdf"), width=10, height=4)
  }
  print(p)
  dev.off()

  assign('dnds.figure.df', comb.df, envir=.GlobalEnv)
}

alt.dnds.calcs <- function() {

  load(pdf.f("nofilter_merged.Rdata"))
  alt.dnds.calc(merged, prefix='nofilter_')
  rm(merged, envir=.GlobalEnv)

  load(pdf.f("nofilter_merged.Rdata"))
  alt.dnds.calc(merged, prefix='nofilter_all_')
  rm(merged, envir=.GlobalEnv)
  
  load(pdf.f("merged.Rdata"))
  alt.dnds.calc(merged, prefix='')
  rm(merged, envir=.GlobalEnv)

  load(pdf.f("merged.Rdata"))
  alt.dnds.calc(merged, prefix='all_')  
  rm(merged, envir=.GlobalEnv)


  # Doing X chromosome on its own is pretty much stupid.
  #load(pdf.f("merged.Rdata"))
  #merged <- merge(merged, ggenes[, c('data_id', 'chr_name')])
  #print(table(merged$chr_name))
  #alt.dnds.calc(subset(merged, chr_name=='X'), prefix='x_')
}

alt.dnds.calc <- function(x, boot.reps=10, prefix='') {

    print("  coding subs")
    #x <- subset(x, nongap_mammal > 20)
#    x <- subset(x, confidence > 0.9 & nongap_count==6 & nongap_mammal > 25)
    nsyn <- !is.na(x$mut_nsyn) & x$mut_nsyn==1
    syn <- !is.na(x$mut_nsyn) & x$mut_nsyn==0
    x[, 'h.ns'] <- as.numeric(x$taxon_id == 9606) & nsyn
    x[, 'h.s'] <- as.numeric(x$taxon_id == 9606) & syn
    x[, 'c.ns'] <- as.numeric(x$taxon_id == 9598) & nsyn
    x[, 'c.s'] <- as.numeric(x$taxon_id == 9598) & syn
    x[, 'g.ns'] <- as.numeric(x$taxon_id == 9593) & nsyn
    x[, 'g.s'] <- as.numeric(x$taxon_id == 9593) & syn
    print("  hcg")
    x[, 'hc.ns'] <- as.numeric(x$taxon_id == 1234) & nsyn
    x[, 'hc.s'] <- as.numeric(x$taxon_id == 1234) & syn
    x[, 'hcg.ns'] <- as.numeric(x$taxon_id == 207598) & nsyn
    x[, 'hcg.s'] <- as.numeric(x$taxon_id == 207598) & syn
    x[, 'o.ns'] <- as.numeric(x$taxon_id == 9600) & nsyn
    x[, 'o.s'] <- as.numeric(x$taxon_id == 9600) & syn
    print("  hcgo")
    x[, 'hcgo.ns'] <- as.numeric(x$taxon_id == 9604) & nsyn
    x[, 'hcgo.s'] <- as.numeric(x$taxon_id == 9604) & syn
    x[, 'm.ns'] <- as.numeric(x$taxon_id == 9544) & nsyn
    x[, 'm.s'] <- as.numeric(x$taxon_id == 9544) & syn
    x[, 'hcgom.ns'] <- as.numeric(x$taxon_id == 376913) & nsyn
    x[, 'hcgom.s'] <- as.numeric(x$taxon_id == 376913) & syn
    x[, 'r.ns'] <- as.numeric(x$taxon_id == 9483) & nsyn
    x[, 'r.s'] <- as.numeric(x$taxon_id == 9483) & syn
    print("  all")
    x[, 'all.ns'] <- nsyn
    x[, 'all.s'] <- syn
    drop.cols <- c('aln_length', 'mut_nsyn', 'taxon_id')
    for (i in 1:length(drop.cols)) {
      x[, drop.cols[i]] <- NULL
    }
    x$pattern <- as.factor(x$pattern)
    coded.subs <- x
    save(coded.subs, file=pdf.f(paste(prefix, 'coded.subs.Rdata', sep='')))
    assign('coded.subs', coded.subs, envir=.GlobalEnv)
    print(nrow(coded.subs))
    print("  done!")

  x <- coded.subs

  print("  filtering...")
  print(nrow(x))
  x <- subset(x, is.na(nongap_count) | nongap_count==6)
  print("  done!")
  print(nrow(x))

  below.q <- quantile(coded.subs$lrt_stat, c(0, 0.05, 0.33, 0.67, 0.98, 1))

  lrt.lo <- c(-209.746, -42.855, -22.819, -10.092, 1.729)  
  lrt.hi <- c(-42.855, -22.819, -10.092, 1.729, 118.981)

  if (prefix == 'all_' || prefix == 'nofilter_all_') {
    lrt.lo <- c(-500)
    lrt.hi <- c(500)
  }

  names(lrt.lo) <- NULL
  names(lrt.hi) <- NULL

  print(lrt.lo)
  print(lrt.hi)

  fld.pairs <- list(
    'Human' = c('h.ns', 'h.s'),
    'Chimpanzee' = c('c.ns', 'c.s'),
    'H/C ancestor' = c('hc.ns', 'hc.s'),
    'Gorilla' = c('g.ns', 'g.s'),
    'H/C/G ancestor' = c('hcg.ns', 'hcg.s'),
    'Orangutan' = c('o.ns', 'o.s'),
    'H/C/G/O ancestor' = c('hcgo.ns', 'hcgo.s'),
    'Macaque' = c('m.ns', 'm.s'),
    'H/C/G/O/M ancestor' = c('hcgom.ns', 'hcgom.s'),
    'Marmoset' = c('r.ns', 'r.s'),
    'All' = c('all.ns', 'all.s')
  )
  fld.names <- names(fld.pairs)

  lrt.labels <- c()
  comb.df <- data.frame()
  for (i in 1:length(lrt.lo)) {
    lo <- lrt.lo[i]
    hi <- lrt.hi[i]

    lbl <- paste(lo, hi, sep=" - ")
    lrt.labels[i] <- lbl
    print(lbl)

    x2 <- x[x$lrt_stat > lo & x$lrt_stat < hi, ]

    for (k in 1:length(fld.pairs)) {
      cur.pair <- fld.pairs[[k]]
      cur.nm <- fld.names[k]
      ns.fld <- cur.pair[1]
      s.fld <- cur.pair[2]

      sum.ratio <- function(D, d=-1) {
        if (d != -1) {
          E <- D[d,]
          return(sum(E[, ns.fld]) / sum(E[, s.fld]) * .29 / .71)
        } else {
          return(sum(D[, ns.fld]) / sum(D[, s.fld]) * .29 / .71)
        }
      }

      do.boot = T
      if (do.boot) {
        boot.ratio <- boot(x2, sum.ratio, R=boot.reps)
        ci.ratio <- boot.ci(boot.ratio, type='basic')
        d.mean <- sum.ratio(x2, 1:nrow(x2))
        bt <- c(ci.ratio$basic[1, 4], d.mean, ci.ratio$basic[1, 5])
        print(bt)
      } else {
        d.mean <- sum.ratio(x2)
        bt <- c(d.mean, d.mean, d.mean)
      }

      cur.df <- data.frame(
        lrt = i,
        lrt.lo = lo,
        lrt.hi = hi,
        n.sites = nrow(x2),
        tot.subs = sum(x2[, ns.fld]) + sum(x2[, s.fld]),
        n.ns = sum(x2[, ns.fld]),
        n.s = sum(x2[, s.fld]),
        species = cur.nm,
        mean = bt[2],
        lo = max(0, bt[1]),
        hi = bt[3]
      )
      comb.df <- rbind(comb.df, cur.df)
    }
  }

#  assign("dnds.df", comb.df, envir=.GlobalEnv)
  print(comb.df)
  write.csv(comb.df, file=pdf.f(paste(prefix, "dnds.csv", sep='')), row.names=F)
}

alt.dnds.figs <- function() {
  alt.dnds.fig(prefix='')
  alt.dnds.fig(prefix='all_')
  alt.dnds.fig(prefix='nofilter_')
  alt.dnds.fig(prefix='nofilter_all_')
}

alt.dnds.fig <- function(prefix='') {
  dnds.df <- read.csv(pdf.f(paste(prefix, 'dnds.csv', sep='')), stringsAsFactors=F)
  print(head(dnds.df))

  x <- dnds.df
  x$max.subs <- max(x$tot.subs)

  n <- length(unique(x$lrt))
  x <- subset(x, !(species %in% c('H/C/G/O/M ancestor', 'Marmoset')))

  x$species <- as.character(x$species)
  x$species <- factor(x$species, labels=ancestor.order(), levels=ancestor.order())
  print(unique(x$species))

  p.df <- function(df) {
    lo <- df[1, 'lrt.lo']
    hi <- df[1, 'lrt.hi']
    lbl <- paste(lo, ' - ', hi, sep='')

    p <- ggplot(df, aes(
      x=species,
      xmin=as.numeric(species)-0.3,
      xmax=as.numeric(species)+0.3,    
      ymin=lo, ymax=hi, y=mean
    ))
    p <- p + geom_rect(width=0.75, alpha=0.3, position='identity')
     p <- p + geom_segment(colour='black', aes(x=as.numeric(species)-0.3, xend=as.numeric(species)+0.3, y=mean, yend=mean), size=0.5)
#    p <- p + scale_y_continuous("dN/dS Estimate", limits=c(0, max(df$hi)))
    p <- p + scale_x_discrete("Species / ancestral lineage")
    p <- p + theme_bw()
    p <- p + opts(
      axis.text.x = theme_text(angle=90, hjust=1),
      title = paste('LRT', lbl)
    )
    return(p)
  }

  print(n)
  pdf(file=pdf.f(paste(prefix, "alt.dnds.pdf", sep='')), width=5*n, height=5)
  vplayout(n, 1)
  for (i in 1:n) {
    cur.df <- subset(x, lrt == i)
    p <- p.df(cur.df)
    print(p, vp=subplot(i, 1))
  }
  dev.off()
}

combined.dnds.figs <- function() {

  combined.dnds.fig(prefix='')
  combined.dnds.fig(prefix='all_')
  combined.dnds.fig(prefix='nofilter_')
  combined.dnds.fig(prefix='nofilter_all_')

}

combined.dnds.fig <- function(prefix='') {

  subs.df <- read.csv(file=paste(prefix, 'dnds.csv', sep=''))
  paml.df <- read.csv(file=paste(prefix, 'paml_dnds.csv', sep=''))

  subs.df$group <- 'Substitutions'
  paml.df$group <- 'PAML'
  subs.df$x_off <- -1
  paml.df$x_off <- 1

  keep.flds <- c('species', 'lrt', 'mean', 'lo', 'hi', 'group', 'x_off')
  x <- rbind(subs.df[, keep.flds], paml.df[, keep.flds])

  x <- subset(x, !(species %in% c('H/C/G/O/M ancestor', 'Marmoset')))
  x$species <- as.character(x$species)
  x$species <- factor(x$species, labels=ancestor.order(), levels=ancestor.order())
  print(unique(x$species))

  x$dx <- 0.2
  x$wd <- 0.25/2

  lrt.lbls <- lrt.quantiles()
  y.lim <- NA

  if (length(grep('all', prefix)) > 0) {
    lrt.lbls <- c('[0, 1]')
    y.lim <- c(0.15, 0.3)
  }

  p.df <- function(df) {
    lrt.bin <- df[1, 'lrt']
    lrt.lbl <- lrt.lbls[lrt.bin]

    p <- ggplot(df, aes(
      x=species,
      xmin=as.numeric(species)-wd + x_off*dx,
      xmax=as.numeric(species)+wd + x_off*dx,
      ymin=lo, ymax=hi, y=mean,
      fill=group
    ))
    p <- p + geom_rect(width=0.33, alpha=0.3, colour='black')
    p <- p + geom_segment(colour='black', aes(
      x=as.numeric(species)-wd + x_off*dx,
      xend=as.numeric(species)+wd + x_off*dx,
      y=mean, yend=mean), size=0.5)
#    p <- p + scale_y_continuous('', limits=c(0, max(df$hi)))
    if (!is.na(y.lim)) {
      p <- p + scale_y_continuous('', limits=y.lim)
    } else {
      p <- p + scale_y_continuous('')
    }
    p <- p + scale_x_discrete('')
    p <- p + scale_fill_manual('', values=c('blue', 'red'))
    p <- p + theme_bw()
    p <- p + opts(
      axis.text.x = theme_text(angle=90, hjust=1),
      title = lrt.lbl,
      legend.position = 'none'
    )
    return(p)
  }

  n <- length(unique(x$lrt))
  pdf(file=pdf.f(paste(prefix, "comb.dnds.pdf", sep='')), width=4*n, height=4)
  vplayout(n, 1)
  for (i in 1:n) {
    cur.df <- subset(x, lrt == i)
    p <- p.df(cur.df)
    print(p, vp=subplot(i, 1))
  }
  dev.off()  
}

paml.xy.fig <- function() {
  source("~/lib/greg-ensembl/scripts/mysql_functions.R")
  dbname <- 'gj1_dnds'
  con <- connect(dbname)

  res.df <- dbReadTable(con, 'results')

  res.df <- subset(res.df, data_id != -500)

  spc.to.lbl <- c(
    'h' = 'Human',
    'c' = 'Chimpanzee',
    'ch' = 'H/C ancestor',
    'g' = 'Gorilla',
    'cgh' = 'H/C/G ancestor',
    'o' = 'Orangutan',
    'cgho' = 'H/C/G/O ancestor',
    'm' = 'Macaque'
  )
  nms <- names(spc.to.lbl)

  col_prefix <- ''

  data.ids <- sort(res.df[, 'data_id'])
  plot.df <- data.frame()
  for (i in 1:length(data.ids)) {
    data.id <- data.ids[i]
    for (j in 1:length(nms)) {
      cur.row <- res.df[i, ]
      dnds.lbl <- paste(nms[j], '_dnds', sep='')
      se.lbl <- paste(nms[j], '_dnds_se', sep='')

      dnds.lbl <- paste(col_prefix, dnds.lbl, sep='')
      se.lbl <- paste(col_prefix, se.lbl, sep='')
      dnds.val <- cur.row[, dnds.lbl]
      dnds.se <- cur.row[, se.lbl]

      dnds.h <- cur.row[, 'h_dnds']
      dnds.h.se <- cur.row[, 'h_dnds_se']

      if(is.na(dnds.se)) {
        dnds.se <- 0
      }

      cur.df <- data.frame(
        species = spc.to.lbl[j],
        data.id = data.id,
        x = dnds.h,
        y = dnds.val / dnds.h,
        y_lo = (dnds.val-dnds.se) / (dnds.h),
        y_hi = (dnds.val+dnds.se) / (dnds.h)
      )
      plot.df <- rbind(plot.df, cur.df)
    }
  }
  print(plot.df)

  plot.df$species <- with(plot.df, reorder(species, y, function(x){-x[1]}))
  
  p <- ggplot(plot.df, aes(x=x, y=y, colour=species))
  p <- p + theme_bw()
  p <- p + geom_line()
  p <- p + geom_point()
  p <- p + geom_errorbar(aes(ymin=y_lo, ymax=y_hi), width=0.1, alpha=0.5)
  p <- p + scale_colour_discrete("Species")
   p <- p + scale_x_log10("dN/dS in human")
   p <- p + scale_y_continuous("dN/dS relative to human")

  pdf(file=pdf.f("paml.xy.pdf"), width=6, height=5)
  print(p)
  dev.off()
  

}

paml.dnds.figs <- function() {
  .paml.dnds.fig(prefix='all_')
  .paml.dnds.fig(prefix='nofilter_all_')
  .paml.dnds.fig(prefix='nofilter_')
  .paml.dnds.fig(prefix='')
}

.paml.dnds.fig <- function(prefix='') {
  source("~/src/greg-ensembl/scripts/mysql_functions.R")
  dbname <- 'gj1_dnds'
  con <- connect(dbname)

  res.df <- dbReadTable(con, 'results')

  if (length(grep('all', prefix)) > 0) {
    res.df <- subset(res.df, data_id == -500)
  } else {
    res.df <- subset(res.df, data_id != -500)
  }

  col_prefix <- ''
  if (length(grep('nofilter', prefix)) > 0) {
    col_prefix <- 'nofilt_'
  }

  spc.to.lbl <- c(
    'h' = 'Human',
    'c' = 'Chimpanzee',
    'ch' = 'H/C ancestor',
    'g' = 'Gorilla',
    'cgh' = 'H/C/G ancestor',
    'o' = 'Orangutan',
    'cgho' = 'H/C/G/O ancestor',
    'm' = 'Macaque',
#    'cghmo' = 'H/C/G/O/M ancestor',
#    'r' = 'Marmoset',
    'm0' = 'All'
  )
  nms <- names(spc.to.lbl)

  data.ids <- sort(res.df[, 'data_id'])
  plot.df <- data.frame()
  for (i in 1:length(data.ids)) {
    data.id <- data.ids[i]
    for (j in 1:length(nms)) {
      cur.row <- res.df[i, ]
      dnds.lbl <- paste(nms[j], '_dnds', sep='')
      se.lbl <- paste(nms[j], '_dnds_se', sep='')

      dnds.lbl <- paste(col_prefix, dnds.lbl, sep='')
      se.lbl <- paste(col_prefix, se.lbl, sep='')
      dnds.val <- cur.row[, dnds.lbl]
      dnds.se <- cur.row[, se.lbl]

      if(is.na(dnds.se)) {
        dnds.se <- 0
      }

      cur.df <- data.frame(
        species = spc.to.lbl[j],
        data.id = data.id,
        lrt = i,
        lrt.lbl = lrt.quantiles()[i],
        mean = dnds.val,
        lo = dnds.val - dnds.se,
        hi = dnds.val + dnds.se
      )
      plot.df <- rbind(plot.df, cur.df)
    }
  }
  print(plot.df)
  write.csv(plot.df, file=pdf.f(paste(prefix, 'paml_dnds.csv', sep='')), row.names=F)

  p.df <- function(df) {
    lbl <- df[1, 'data.id']

    p <- ggplot(df, aes(
      x=species,
      xmin=as.numeric(species)-0.3,
      xmax=as.numeric(species)+0.3,    
      ymin=lo, ymax=hi, y=mean
    ))
    p <- p + theme_bw()
    p <- p + geom_rect(width=0.75, alpha=0.3, position='identity')
#    p <- p + geom_rect(fill='black', aes(ymin=mean-0.0002, ymax=mean+0.0002), width=0.75, alpha=0.6, position='identity')
     p <- p + geom_segment(colour='black', aes(x=as.numeric(species)-0.3, xend=as.numeric(species)+0.3, y=mean, yend=mean), size=0.5)
     p <- p + scale_y_continuous('')
#    p <- p + scale_y_continuous('', limits=c(0, max(df$hi)))
#    p <- p + scale_y_continuous('', limits=c(0, 1.8))
    p <- p + scale_x_discrete('')
    p <- p + opts(
      axis.text.x = theme_text(angle=90, hjust=1),
      title = df[1, 'lrt.lbl']
    )
    return(p)
  }

  n <- length(unique(plot.df$data.id))
  data.ids <- unique(plot.df$data.id)
  pdf(file=pdf.f(paste(prefix, "paml.dnds.pdf", sep='')), width=5*n, height=5)
  vplayout(n, 1)
  for (i in 1:length(data.ids)) {
    cur.id <- data.ids[i]
    cur.df <- subset(plot.df, data.id == cur.id)
    p <- p.df(cur.df)
    print(p, vp=subplot(i, 1))
  }
  dev.off()
}

lrt.quantiles <- function() {
  return(c(
    '[0, 0.05]',
    '[0.05, 0.33]',
    '[0.33, 0.67]',
    '[0.67, 0.98]',
    '[0.98, 1]'))
}

ancestor.order <- function() {
  fld.pairs <- list(
    'Human' = c('h.ns', 'h.s'),
    'Chimpanzee' = c('c.ns', 'c.s'),
    'H/C ancestor' = c('hc.ns', 'hc.s'),
    'Gorilla' = c('g.ns', 'g.s'),
    'H/C/G ancestor' = c('hcg.ns', 'hcg.s'),
    'Orangutan' = c('o.ns', 'o.s'),
    'H/C/G/O ancestor' = c('hcgo.ns', 'hcgo.s'),
    'Macaque' = c('m.ns', 'm.s'),
#    'H/C/G/O/M ancestor' = c('hcgom.ns', 'hcgom.s'),
#    'Marmoset' = c('r.ns', 'r.s'),
    'All' = c('all.ns', 'all.s')
  )
  return(names(fld.pairs))
}

pop.size.calc <- function(dnds_a, dnds_b) {
  w <- 1e50
  n <- 1

  root.a <- uniroot(function(x) 2*x / (1 - exp(-2*x)) - dnds_a, lower=-w, upper=w)
  root.b <- uniroot(function(x) 2*x / (1 - exp(-2*x)) - dnds_b, lower=-w, upper=w)
  
  #print(root.a$root)
  #print(root.b$root)
  return(root.b$root / root.a$root)

  dnds.f <- function(x) 2*n*x / (1 - exp(-2*n*x))
  xs <- seq(from=-2, to=2, length.out=100)
  print(xs)
  ys <- dnds.f(xs)
  print(ys)
  pdf(pdf.f('temp.pdf'))
  plot(xs, ys, ylim=c(0, 3))
  dev.off()

}

lrt.cumsum <- function(dir='down') {
  x <- coded.subs

  if (dir == 'up') {
    x <- x[order(x$lrt_stat),]
  } else {
    x <- x[order(-x$lrt_stat),]
  }  

  x$cum.h.ns <- cumsum(x$h.ns)
  x$cum.h.s <- cumsum(x$h.s) + 1
  x$cum.c.ns <- cumsum(x$c.ns)
  x$cum.c.s <- cumsum(x$c.s) + 1
  x$cum.hc.ns <- cumsum(x$hc.ns)
  x$cum.hc.s <- cumsum(x$hc.s) + 1
  x$cum.g.ns <- cumsum(x$g.ns)
  x$cum.g.s <- cumsum(x$g.s) + 1
  x$cum.hcg.ns <- cumsum(x$hcg.ns)
  x$cum.hcg.s <- cumsum(x$hcg.s) + 1
  x$cum.o.ns <- cumsum(x$o.ns)
  x$cum.o.s <- cumsum(x$o.s) + 1
  x$cum.hcgo.ns <- cumsum(x$hcgo.ns)
  x$cum.hcgo.s <- cumsum(x$hcgo.s) + 1
  x$cum.m.ns <- cumsum(x$m.ns)
  x$cum.m.s <- cumsum(x$m.s) + 1
  x$cum.hcgom.ns <- cumsum(x$hcgom.ns)
  x$cum.hcgom.s <- cumsum(x$hcgom.s) + 1
  x$cum.r.ns <- cumsum(x$r.ns)
  x$cum.r.s <- cumsum(x$r.s) + 1
  indices <- seq(from=1, to=nrow(x), by=50)
  x <- x[indices,]
  x$index <- 1:nrow(x)
  x$f <- .29 / .71

  if (dir == 'up') {
    negative <- sum(x$lrt_stat < -3.99)
    neutral <- sum(x$lrt_stat < 0)
    positive <- sum(x$lrt_stat < 3.99)
  } else {
    negative <- sum(x$lrt_stat > -3.99)
    neutral <- sum(x$lrt_stat > 0)
    positive <- sum(x$lrt_stat > 3.99)
  }

  plot.df <- data.frame()
  species <- c('h', 'c', 'g', 'o', 'm', 'r')
  for (i in 1:length(species)) {
    cur.s <- species[i]
    ns.fld <- paste('cum.', cur.s, '.ns', sep='')
    s.fld <- paste('cum.', cur.s, '.s', sep='')
    plot.df <- rbind(plot.df, data.frame(
      lbl = cur.s,
      ns = x[, ns.fld],
      s = x[, s.fld],
      max.ns = max(x[, ns.fld]),
      max.s = max(x[, s.fld]),
      f = x[, 'f'],
      index = x[, 'index'],
      lrt_stat = x[, 'lrt_stat']
    ))
  }

  pdf(pdf.f('cum.dnds.pdf'))
  p <- ggplot(plot.df, aes(x=index, y=ns / s * f, colour=lbl))

  p <- p + geom_vline(xintercept=negative, colour='blue', alpha=0.5, linetype='dashed')
  p <- p + geom_vline(xintercept=neutral, colour='gray', alpha=0.5, linetype='dashed')
  p <- p + geom_vline(xintercept=positive, colour='red', alpha=0.5, linetype='dashed')

  p <- p + geom_line()
#  p <- p + scale_x_continuous(limits=c(1, nrow(plot.df)/40))
#  p <- p + xlim(1, floor(nrow(plot.df)/5))
#  p <- p + ylim(0, 0.2)
#   p <- p + coord_cartesian(xlim=c(1, nrow(plot.df)/2))

  print(p)
  dev.off()
}

ils.genes <- function() {
  x <- get('y', envir=.GlobalEnv)

  dnds.q <- quantile(x$m0_dnds, c(0, 0.2, 0.4, 0.6, 0.8, 1.0))

  x$dnds <- x[, 'm0_dnds']
  for (i in 2:length(dnds.q)) {
    dnds.lo <- dnds.q[i-1]
    dnds.hi <- dnds.q[i]
    dnds.sub <- subset(x, dnds > dnds.lo & dnds <= dnds.hi)
    x[x$dnds >= dnds.lo & x$dnds <= dnds.hi, 'dnds.label'] <- 
      paste(base::format(dnds.lo, digits=2),
            base::format(dnds.hi, digits=2),
            sep=' -\n')
  }

  print(unique(x$dnds.label))
  dnds <- ddply(x, .(dnds.label), function(data) {
    hg.bt <- smean.cl.boot(chunk.b.gene(data))
    cg.bt <- smean.cl.boot(chunk.c.gene(data))

    x.v <- chunk.x.gene(data)
    x.v <- as.numeric(x.v > 0)
    print(summary(x.v))
    x.v <- subset(x.v, !is.na(x.v) & !is.infinite(x.v))
    x.bt <- smean.cl.boot(x.v)

    return(data.frame(
      mean.dnds = mean(data$dnds),
      hg = hg.bt[1],
      hg.lo = hg.bt[2],
      hg.hi = hg.bt[3],
      cg = cg.bt[1],
      cg.lo = cg.bt[2],
      cg.hi = cg.bt[3],
      x = x.bt[1],
      x.lo = x.bt[2],
      x.hi = x.bt[3]
    ))
  })

  dnds <- dnds[order(dnds$mean.dnds),]
  dnds$dnds.label <- factor(dnds$dnds.label, levels = dnds$dnds.label, ordered=T)
  print(as.numeric(dnds$dnds.label))

  print(dnds$dnds.label)

  dnds$value.label <- 'Human-Gorilla'
  dnds$lo <- dnds$hg.lo
  dnds$hi <- dnds$hg.hi
  dnds$mid <- dnds$hg
  dnds2 <- dnds
  dnds2$value.label <- 'Chimpanzee-Gorilla'
  dnds2$lo <- dnds$cg.lo
  dnds2$hi <- dnds$cg.hi
  dnds2$mid <- dnds$cg
  dnds3 <- dnds
  dnds3$value.label <- 'ILS Density'
  dnds3$lo <- dnds$x.lo
  dnds3$hi <- dnds$x.hi
  dnds3$mid <- dnds$x

  print(dnds3$mean.dnds)

  pdf(pdf.f("ils.genes.pdf"), width=5, height=4)
  p <- ggplot(dnds3, aes(x=as.numeric(dnds.label), fill=mean.dnds))

#  p <- p + geom_rect(data=dnds, aes(
#    ymin=lo, ymax=hi,
#    xmin=as.numeric(dnds.label) - .33,
#    xmax=as.numeric(dnds.label) + .33)
#  )
#  p <- p + geom_line(data=dnds, colour='blue', aes(x=as.numeric(dnds.label),y=mid))

#  p <- p + geom_rect(data=dnds2, aes(
#    ymin=lo, ymax=hi,
#    xmin=as.numeric(dnds.label) + .33,
#    xmax=as.numeric(dnds.label) + .05)
#  )
#  p <- p + geom_line(data=dnds2, colour='gray', aes(x=as.numeric(dnds.label),y=mid))

  p <- p + geom_rect(data=dnds3, alpha=0.3, aes(
    ymin=lo, ymax=hi,
    xmin=as.numeric(dnds.label) + .33,
    xmax=as.numeric(dnds.label) - .33)
  )
  p <- p + geom_line(data=dnds3, colour='black', aes(x=as.numeric(dnds.label), y=mid))

  p <- p + scale_fill_gradient("dN/dS", low="blue", high="red")
#  p <- p + scale_fill_manual("ILS type", values=c('#AAAAAA', '#AAAAAA'))
  p <- p + scale_x_continuous("dN/dS Bin",
    breaks=c(1:length(labels(dnds$dnds.label))),
    labels=levels(dnds$dnds.label)
  )
  p <- p + scale_y_continuous("ILS density", limits=c(0, 0.35))
#  p <- p + opts(title="ILS densities vs dN/dS")
   p <- p + theme_bw()
  print(p)
  dev.off()

  print(dnds)
}

ils.densities <- function() {
  pdf(pdf.f("ils.regions.quantile.pdf"), width=5, height=4)
  p <- ils.density.plot('exon_dist', 'm0_dnds', bin_method='quantile', relative=F)
#  p <- p + opts(title="ILS density vs Exon Distance, Mammal dN/dS")
  print(p)
  dev.off()  

  return()

  pdf(pdf.f("ils.regions.window.pdf"))
  p <- ils.density.plot('exon_dist', 'slr_dnds', bin_method='window', relative=F)
  p <- p + opts(title="ILS density vs Exon Distance, Mammal dN/dS")
  print(p)
  dev.off()  

  pdf("ils_exon_mammal.pdf")
  p <- ils.density.plot('exon_dist', 'gene_slr', bin_method='window', relative=F)
  p <- p + opts(title="ILS density vs Exon Distance, Mammal dN/dS")
  print(p)

  p <- ils.density.plot('exon_dist', 'gene_slr', bin_method='quantile', relative=F)
  p <- p + opts(title="ILS density vs Exon Distance, Mammal dN/dS")
  print(p)

  p <- ils.density.plot('exon_dist', 'gene_slr', bin_method='window', relative=T)
  p <- p + opts(title="Relative ILS density vs Exon Distance, Mammal dN/dS")
  print(p)

  p <- ils.density.plot('exon_dist', 'gene_slr', bin_method='quantile', relative=T)
  p <- p + opts(title="Relative ILS density vs Exon Distance, Mammal dN/dS")
  print(p)

  dev.off()
  return()

}

ils.density.plot <- function(distance_field, dnds_field, bin_method='window', relative=T) {
  if (!exists("orig.chunks")) {
    print("  fetching chunks from DB")
    assign("dbname", 'gj1_ils', envir=.GlobalEnv)
    source("~/src/greg-ensembl/scripts/collect_sitewise.R")
    orig.chunks <- get.vector(con, 'select * from chunks_10')
    assign("orig.chunks", orig.chunks, envir=.GlobalEnv)
    print("  done!")
  }
  chunks <- orig.chunks

  y <- get("y", envir=.GlobalEnv)
  chunks <- merge(chunks, y[, c('gene_name', dnds_field)])

  chunks$dist <- abs(chunks[, distance_field])
  chunks$dnds <- chunks[, dnds_field]

  chunks <- subset(chunks, !(chr_name %in% c('MT', 'X', 'Y')))
#  chunks <- subset(chunks, dnds > -1)
  chunks <- subset(chunks, !is.na(dist))

  assign('cur.chunks', chunks, envir=.GlobalEnv)

  x.labels <- c()
  if (bin_method == 'quantile') {
    q.mids <- c()
    qs <- quantile(chunks$dist,c(0, 0.2, 0.4, 0.6, 0.8, 1.0),na.rm=T)
    print(qs)
    for (i in 1:length(qs)) {
      bin.hi <- qs[i]
      if (i == 1) {
        bin.lo <- min(chunks$dist)
      } else {
        bin.lo <- qs[i-1]
      }
      print(paste(i,bin.lo,bin.hi))
      if (bin.lo == 0 && bin.hi == 0) {
        chunks[chunks$dist == 0 , 'x_bin'] <- i
      } else {
        chunks[chunks$dist > bin.lo & chunks$dist <= bin.hi , 'x_bin'] <- i
      }
      q.mids[i] <- paste(as.integer(bin.lo), as.integer(bin.hi), sep="\n")
    }
    chunks$x_bin <- as.factor(chunks$x_bin)
    x.labels <- q.mids
    print(x.labels)
  } else {  
     win.width <- 20000
     chunks$x_bin <- 'all'
     chunks$x_bin <- floor(chunks$dist / win.width) * win.width
     chunks <- subset(chunks, dist < 200000)
     chunks$x_bin <- as.factor(chunks$x_bin)

     x.labels <- paste(
       sprintf("%.2g", as.numeric(levels(chunks$x_bin))),
       sprintf("%.2g", as.numeric(levels(chunks$x_bin))+win.width),
       sep="\n"
     )
     print(x.labels)
  }

  dnds.subset <- subset(chunks, dnds > -1)
  dnds.qs <- quantile(dnds.subset$dnds, c(0.25, 0.75))
  print(dnds.qs)

  lo.dnds <- chunks$dnds >= 0 & chunks$dnds <= dnds.qs[1]
  hi.dnds <- chunks$dnds >= dnds.qs[2]
  chunks$dnds_bin <- NA
  chunks$x_off <- 0
  chunks[lo.dnds, 'dnds_bin'] <- 'lo'
  chunks[lo.dnds, 'x_off'] <- -.15
  chunks[hi.dnds, 'dnds_bin'] <- 'hi'   
  chunks[hi.dnds, 'x_off'] <- .15
  c.copy <- chunks
  c.copy$dnds_bin <- 'all'
  c.copy$x_off <- 0
  chunks <- rbind(chunks, c.copy)
  chunks <- subset(chunks, !is.na(dnds_bin))
  chunks$dnds_bin <- factor(chunks$dnds_bin, levels=c('all', 'lo', 'hi'), labels=c('All', 'Low', 'High'))

  print(table(chunks$x_bin))
  print(table(chunks$dnds_bin))

  chunks$value <- chunk.x(chunks)
  chunks$value <- as.numeric(chunks$value > 0)

  f.df <- function(x) {
    x <- subset(x, !is.na(value) & !is.infinite(value))
    bt <- smean.cl.boot(x$value)

    #stderr <- function(x) sqrt(var(x)/length(x))
    #err.v <- stderr(x$value)

    return(data.frame(
      n = nrow(x),
      v = bt[1],
      v.lo = bt[2],
      v.hi = bt[3],
      x_off = x[1, 'x_off']
    ))
  }
  print("  ddplying")
  res.df <- ddply(chunks, c('x_bin','dnds_bin'), f.df)

  if (relative) {
    # Normalize by the mean value in 'all' windows for each x bin
    print("  normalizing")
    x_bins <- unique(res.df$x_bin)
    new.res.df <- data.frame()
    for (i in 1:length(x_bins)) {
      cur.bin <- x_bins[i]
      cur.df <- subset(res.df, x_bin==cur.bin)
     cur.df.all <- subset(cur.df, dnds_bin == 'All')
      # The mean value across all windows.
       all.y.value <- cur.df.all$v
  
      # Normalize all dnds bins to that value.
      cur.df$v <- cur.df$v / all.y.value
      cur.df$v.hi <- cur.df$v.hi / all.y.value
      cur.df$v.lo <- cur.df$v.lo / all.y.value
      new.res.df <- rbind(new.res.df, cur.df)
    }
    res.df <- new.res.df
    all.df <- subset(res.df, dnds_bin == 'All')
    res.df <- subset(res.df, dnds_bin != 'All')
  } else {
    all.df <- subset(res.df, dnds_bin == 'All')
    res.df <- subset(res.df, dnds_bin != 'All')
  }

  print(res.df)

  print("  plotting")
  p <- ggplot(res.df, aes(colour=dnds_bin, fill=dnds_bin))
    p <- p + geom_rect(data=all.df, fill='darkgray', colour='darkgray', alpha=1, aes(
      xmin=as.numeric(x_bin) - .3,
      xmax=as.numeric(x_bin) + .3,
      ymin=v.lo,
      ymax=v.hi
    ))
    p <- p + geom_rect(data=res.df, alpha=0.3, size=0, colour=NA,
      aes(
        xmin=as.numeric(x_bin)+ x_off - .15,
        xmax=as.numeric(x_bin)+ x_off + .15,
        ymin=v.lo,
        ymax=v.hi
      )
    )
  p <- p + geom_line(data=res.df, aes(x=as.numeric(x_bin),y=v))
  p <- p + scale_colour_manual("dN/dS Bin", values=c('black', 'blue', 'red'))
  p <- p + scale_fill_manual("dN/dS Bin", values=c('black', 'blue', 'red'))
  p <- p + scale_y_continuous("ILS density", limits=c(0, 0.35))
  p <- p + scale_x_continuous(name="Distance Bin", breaks=c(1:length(x.labels)), labels=x.labels)
  p <- p + theme_bw()
  return(p)
}

chunk.a <- function(x) {
  ((x$n_100 + x$n_010)/2 + x$n_110) / x$n_001
}

chunk.b <- function(x) {
  (x$n_101+x$n_011) / ((x$n_100 + x$n_010)/2 + x$n_110 + x$n_001)
}

chunk.x <- function(x) {
  #x$n_110 <- pmax(x$n_110, 1)
  #(x$n_101 + x$n_011) / (x$n_110)
  (x$n_101 + x$n_011) - x$n_110
}

chunk.x.gene <- function(x) {
  #x$ptrn_110 <- pmax(x$ptrn_110, 1)
  #(x$ptrn_101 + x$ptrn_011) / (x$ptrn_110)
  (x$ptrn_101 + x$ptrn_011) - x$ptrn_110
}

chunk.a.gene <- function(x) {
  ((x$ptrn_100)) / (x$ptrn_001+1)
}

chunk.b.gene <- function(x) {
  x$ptrn_110 <- pmax(x$ptrn_110, 1)
  (x$ptrn_101) / (x$ptrn_110)
}

chunk.c.gene <- function(x) {
  x$ptrn_110 <- pmax(x$ptrn_110, 1)
  (x$ptrn_011) / (x$ptrn_110)
}

chunk.c <- function(x) {
  (x$n_101+x$n_011)/2 / (x$n_110+1)
}

subplot <- function(x, y) viewport(layout.pos.col=x, layout.pos.row=y)
vplayout <- function(x, y) {
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(y,x)))
}

format.numeric.df <- function(x,digits=3) {  
  nums <- unlist(lapply(x,is.numeric))
  nums <- nums[nums == TRUE]
  for (col in names(nums)) {
    x[,col] <- base::format(x[,col], digits=digits)
  }
  return(x)
}