library(plyr)
library(ggplot2)
library(boot)
library(Hmisc)
source("~/src/greg-ensembl/scripts/go_enrichments.R")

dbname <- 'gj1_gorilla'

pst <- function(...) { return(paste(..., sep=''))}

gene.sub.counts <- function(x) {
  x <- subset(x, !is.na(mut_nsyn))
  return(data.frame(
    g.ns = nrow(subset(x, mut_nsyn==1 & taxon_id %in% c(9593))),
    g.s = nrow(subset(x, mut_nsyn==0 & taxon_id %in% c(9593))),
    h.ns = nrow(subset(x, mut_nsyn==1 & taxon_id %in% c(9606))),
    h.s = nrow(subset(x, mut_nsyn==0 & taxon_id %in% c(9606))),
    c.ns = nrow(subset(x, mut_nsyn==1 & taxon_id %in% c(9598))),
    c.s = nrow(subset(x, mut_nsyn==0 & taxon_id %in% c(9598))),
    hc.ns = nrow(subset(x, mut_nsyn==1 & taxon_id %in% c(1234))),
    hc.s = nrow(subset(x, mut_nsyn==0 & taxon_id %in% c(1234))),
    hcg.ns = nrow(subset(x, mut_nsyn==1 & taxon_id %in% c(5678))),
    hcg.s = nrow(subset(x, mut_nsyn==0 & taxon_id %in% c(5678)))
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

get.pvals <- function() {
  # Collect genes and PAML inferred substitutions from the SQL tables.
  if (!exists("genes")) {
    source("~/src/greg-ensembl/scripts/collect_sitewise.R")
    print("  loading genes")
    ggenes <- get.vector(con, 'select * from genes')
    ggenes$slr_dnds <- as.numeric(ggenes$slr_dnds)
    ggenes$m0_dnds <- as.numeric(ggenes$m0_dnds)
    assign("ggenes", ggenes, envir=.GlobalEnv)
  }
  if (!exists("subs")) {
    source("~/src/greg-ensembl/scripts/collect_sitewise.R")
    print("  loading subs")
    subs <- get.vector(con, 'select * from subs')

    hc.sub.indices <- grep("ENSP0.*ENSPTRP0.*",subs$leaves_beneath)
    subs[hc.sub.indices, 'taxon_id'] <- 1234
    hcg.sub.indices <- grep("ENSGGOP0.*ENSP0.*ENSPTRP0.*",subs$leaves_beneath)
    subs[hcg.sub.indices, 'taxon_id'] <- 5678
    assign("subs", subs, envir=.GlobalEnv)
  }
  if (!exists("sites")) {
    source("~/src/greg-ensembl/scripts/collect_sitewise.R")
    print("  loading sites")
    sites <- get.vector(con, 'select data_id, aln_position as aln_pos, omega_upper, omega from sites')
    assign("sites", sites, envir=.GlobalEnv)
  }
  if (!exists("subs.sites")) {
    print("  merging subs & sites")
    s <- subset(subs, taxon_id %in% c(9593, 9606, 9598, 1234))
    subs.sites <- merge(sites, s, by=c('data_id', 'aln_pos'), all.x=TRUE)
    print(nrow(subs.sites))
    assign("subs.sites", subs.sites, envir=.GlobalEnv)
  }
  if (!exists("merged")) {
    print("  merging genes & subs")
    merged <- merge(ggenes[, c('data_id', 'aln_length', 'm0_dnds', 'slr_dnds')], 
      subs.sites, by=c('data_id'))
    assign("merged", merged, envir=.GlobalEnv)
    save(merged, file=pdf.f('subs.sites.Rdata'))
  }
  g <- ggenes
  s <- subs.sites

  muts.per.codon <- nrow(subset(s, taxon_id==9593 & (mut_nsyn==1 | mut_syn==1))) / nrow(s)

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

  counts <- ddply(merged, .(data_id), gene.sub.counts)
  res <- merge(res, counts)

  res.merged <- merge(ggenes, res)
  assign("res.merged", res.merged, envir=.GlobalEnv)

  x <- merge(res.merged, stats.all)
  y <- get.top.genes(x, stat='lrt.5', n=nrow(x))
  assign('y', y, envir=.GlobalEnv)
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

  for (cur.type in c('all', 'noC', 'noH')) {
    cur.stats <- subset(g, select=c('data_id', 'gene_name', 'protein_id', 'm0_dnds', 'slr_dnds', 'data_prefix', 'aln_length'))
    print(str(cur.stats))
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
    for (i in 11:17) {    
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

    cur.stats <- cur.stats[order(-cur.stats$lrt.5),]
    assign(paste("stats.", cur.type, sep=''), cur.stats, envir=.GlobalEnv)
  }

  cols <- colnames(stats.all)
  pval.fields <- cols[grepl('pval',cols)]
  lrt.fields <- cols[grepl('lrt',cols)]

  for (prefix in c('noC', 'noH')) {
    cur.stats <- get(paste('stats.', prefix,sep=''))
    for (field in c(pval.fields, lrt.fields)) {
      prefix.lbl <- paste(prefix, '.', field, sep='')
      stats.all[, prefix.lbl] <- cur.stats[, field]
    }
  }
  assign('stats.all', stats.all, envir=.GlobalEnv)
  save(stats.all, file=pdf.f("stats.all.Rdata"))
}

all.results.table <- function() {
  get.lrts()
#  get.pvals()
  x <- merge(res.merged, stats.all)
  y <- get.top.genes(x, stat='lrt.5', n=nrow(x))

  write.csv(y, file="~/src/greg-ensembl/projects/gorilla/table.csv",row.names=F)
  assign('y', y, envir=.GlobalEnv)
  save(y, file=pdf.f("y.Rdata"))
}

diff.datasets <- function() {
  compare.top <- function(x, y, fld) {
    mg <- merge(x, y, by=c('data_id','gene_name'))
    mg <- subset(mg, !is.na(lrt.5.x) & !is.na(lrt.5.y))
    f.x <- paste(fld,'.x',sep='')
    f.y <- paste(fld,'.y',sep='')
    mg.x <- mg[order(-mg[,f.x]),]
    mg.y <- mg[order(-mg[,f.y]),]
    
    print(head(subset(mg, lrt.5.x > 10 & lrt.5.y < 10)))
  }

  compare.top(stats.noH, stats.noC, 'lrt.5')
}

compare.pval.methods <- function(x) {
#  x <- subset(x, lrt.5 < 0)

  print("LRT vs NS fraction")
  print(cor(x$rank.lrt.5, x$rank.g.f.nsyn, method='spearman'))

  print("LRT vs Binomial")
  print(cor(x$rank.lrt.5, x$rank.pval.bn.gorilla, method='spearman'))

  print("LRT vs Sitewise")
  print(cor(x$rank.lrt.5, x$rank.pval.sw.gorilla, method='spearman'))

  print("NS fraction vs Binomial")
  print(cor(x$rank.g.f.nsyn, x$rank.pval.bn.gorilla, method='spearman'))

  print("NS fraction vs Sitewise")
  print(cor(x$rank.g.f.nsyn, x$rank.pval.sw.gorilla, method='spearman'))

  print("Binomial vs Sitewise")
  print(cor(x$rank.pval.bn.gorilla, x$rank.pval.sw.gorilla, method='spearman'))

  print("Grantham vs Sitewise")
  print(cor(x$rank.pval.gr.gorilla, x$rank.pval.sw.gorilla, method='spearman'))

  print("Grantham vs Binomial")
  print(cor(x$rank.pval.gr.gorilla, x$rank.pval.bn.gorilla, method='spearman'))

  print("Grantham vs LRT")
  print(cor(x$rank.pval.gr.gorilla, x$rank.lrt.5, method='spearman'))
}

submit.jobs <- function(rows) {
  l_ply(rows$job_id, 
    function(x){
      system(paste('bsub -q normal runWorker.pl -url mysql://ensadmin:ensembl@ens-research:3306/gj1_gorilla -debug 1 -job_id ',x,sep=''))
    }
  )
}

top.genes.set <- function(x, n=5, dir='up') {
  comb.df <- data.frame()

  if (dir == 'up') {
    x <- subset(x, g.ns > 0)
    vars <- c('rank.lrt.5', 'rank.g.f.nsyn', 'rank.pval.bn.gorilla', 'rank.pval.sw.gorilla', 'rank.pval.gr.gorilla')
  } else {
    x <- subset(x, g.s > 0 & c.s > 0 & h.s > 0)
    vars <- c('rank.lrt.5', 'rank.lrt.7', 'rank.lrt.6')
  }

  for (var in vars) {
    fields <- c('gene_name', 'job_id', 'data_prefix', 'slr_dnds', 'm0_dnds', 'aln_length', 'g.f.nsyn','c.f.nsyn', 'h.f.nsyn',
      'rank.lrt.5', 'rank.g.f.nsyn', 'rank.pval.bn.gorilla', 'rank.pval.sw.gorilla', 'rank.pval.gr.gorilla',
      'lrt.5', 'lrt.7', 'lrt.6',
      'lrt.14', 'pval.adj.14')
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
    write.csv(comb.df, file=pdf.f("top.genes.csv"), row.names=F)
  } else {
    write.csv(comb.df, file=pdf.f("bot.genes.csv"), row.names=F)
  }
  return(comb.df)
}

summary.table <- function(x) {
  x <- subset(x, !is.na(lrt.5) & !is.na(lrt.1) & !is.na(lrt.2))

  tbl.df <- data.frame()

  for (i in 1:17) {
    pval.key <- paste('pval.', i, sep='')
    pval.adj.key <- paste('pval.adj.', i, sep='')
    print(pval.adj.key)
    lrt.key <- paste('lrt.', i, sep='')
    x$tmp.pval <- x[, pval.key]
    x$tmp.pval.adj <- x[, pval.adj.key]
    if (i <= 10) {
      x$tmp.lrt <- x[, lrt.key]
      tbl.df[1, i] <- nrow(subset(x, tmp.pval < 0.05 & tmp.lrt > 0))
      tbl.df[2, i] <- nrow(subset(x, tmp.pval.adj < 0.1 & tmp.lrt > 0))
      tbl.df[3, i] <- nrow(subset(x, tmp.pval < 0.05 & tmp.lrt < 0))
      tbl.df[4, i] <- nrow(subset(x, tmp.pval.adj < 0.1 & tmp.lrt < 0))
    } else {
      tbl.df[1, i] <- nrow(subset(x, tmp.pval < 0.05))
      tbl.df[2, i] <- nrow(subset(x, tmp.pval.adj < 0.1))
    }
  }

  colnames(tbl.df) <- c(
  'H', 'C', 'H+CH', 'C+CH', 'G',
  'AGA Branch', 'AGA Clade',
  'H+C', 'G+H', 'G+C',
  'bs AGA Branch', 'bs H', 'bs C', 'bs G', 'bs O',
  'bs GA Branch',
  'bs GA+marmoset'
  )

  library(venneuler)
  
  # UP!
  pval.t <- 0.05
  g <- x$pval.5 < pval.t & x$lrt.5 > 0
  h <- x$pval.1 < pval.t & x$lrt.1 > 0
  c <- x$pval.2 < pval.t & x$lrt.2 > 0

  print(x[g & h & c,])

  venn.df <- data.frame(
    Human = h,
    Chimpanzee = c,
    Gorilla = g
  )
  a <- apply(venn.df,1,function(x) {paste(as.numeric(x),collapse='')})
  t <- table(a)
  print(t)
  print(t/nrow(venn.df))
  print(names(venn.df))

  pdf(file=pdf.f('venn.up.pdf'), width=4, height=4)
  p <- venneuler(venn.df)
  print(p)
  plot(p)
  dev.off()

  # DOWN!
  pval.t <- 0.05
  g <- x$pval.5 < pval.t & x$lrt.5 < 0
  h <- x$pval.1 < pval.t & x$lrt.1 < 0
  c <- x$pval.2 < pval.t & x$lrt.2 < 0
  venn.df <- data.frame(
    Human = h,
    Chimpanzee = c,
    Gorilla = g
  )
  a <- apply(venn.df,1,function(x) {paste(as.numeric(x),collapse='')})
  t <- table(a)
  print(t)
  print(t/nrow(venn.df))
  print(names(venn.df))

  pdf(file=pdf.f('venn.down.pdf'), width=4, height=4)
  p <- venneuler(venn.df)
  plot(p)
  dev.off()
}

get.top.genes <- function(x, stat='pval.sw.gorilla', dir=1, n=20) {
  # FILTER: Only take genes with at least 1 substitution in each species.
#  x <- subset(x, g.ns+g.s > 0 & h.ns+h.s > 0 & c.ns+c.s > 0)

  # FILTER: at least 1 gor non-synonymous sub.
#  x <- subset(x, g.ns > 0)

  # FILTER: Only take genes with dnds > 0
  x <- subset(x, slr_dnds > 0)

  if (stat %in% c('lrt.5','lrt.4','lrt.3')) {
    dir <- -1
  }

  x$g.f.nsyn <- x$g.ns / x$aln_length
  x$c.f.nsyn <- x$c.ns / x$aln_length
  x$h.f.nsyn <- x$h.ns / x$aln_length
  x$hc.f.nsyn <- x$hc.ns / x$aln_length

  rank.fields <- c()
  for (field in c(
    'lrt.5', 'pval.gr.gorilla', 'pval.sw.gorilla', 'pval.bn.gorilla', 'g.f.nsyn', 
    'lrt.4', 'pval.gr.human', 'pval.sw.human', 'pval.bn.human', 'h.f.nsyn',
    'lrt.3', 'pval.gr.chimp', 'pval.sw.chimp', 'pval.bn.chimp', 'c.f.nsyn',
    'lrt.6', 'lrt.7', 'lrt.1', 'lrt.2')) {
    rank.field <- paste('rank',field,sep='.')
    if (field %in% c('g.f.nsyn','c.f.nsyn', 'h.f.nsyn', 'lrt.5', 'lrt.4', 'lrt.3', 'lrt.6', 'lrt.7')) {
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

go.enrichments <- function(x) {
  if (!exists('go.ens')) {
    go.rdata.file <- "~/scratch/gorilla/go_60.Rdata"
    load(go.rdata.file)
  }
  source("~/src/greg-ensembl/scripts/go_enrichments.R")

  y <- x
  y$rank.g.sw <- y[, 'rank.pval.sw.gorilla']
  y$rank.g.gr <- y[, 'rank.pval.gr.gorilla']
  y$rank.g.bn <- y[, 'rank.pval.bn.gorilla']
  y$rank.g.lrt <- y[, 'rank.lrt.5']
  y$rank.h <- y[, 'rank.pval.sw.human']
  y$rank.c <- y[, 'rank.pval.sw.chimp']
  y$rank.aga.branch <- y[, 'rank.lrt.6']
  y$rank.aga.clade <- y[, 'rank.lrt.7']

  all.gor <- subset(y)

  hi.gor <- subset(y, rank.g.sw <= 250)
  tbl <- enrich.list(hi.gor$protein_id, all.gor$protein_id, include.iea=TRUE)
  assign('hi.gor.sw', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)

  hi.gor <- subset(y, rank.g.gr <= 250)
  tbl <- enrich.list(hi.gor$protein_id, all.gor$protein_id, include.iea=TRUE)
  assign('hi.gor.gr', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)

  hi.gor <- subset(y, rank.g.bn <= 250)
  tbl <- enrich.list(hi.gor$protein_id, all.gor$protein_id, include.iea=TRUE)
  assign('hi.gor.bn', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)

  hi.gor <- subset(y, rank.g.lrt <= 250)
  tbl <- enrich.list(hi.gor$protein_id, all.gor$protein_id, include.iea=TRUE)
  assign('hi.gor.lrt', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)

  hi <- subset(y, rank.aga.branch <= 250)
  tbl <- enrich.list(hi$protein_id, all.gor$protein_id, include.iea=TRUE)
  assign('hi.aga.branch', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)

  lo <- subset(y, rank.aga.branch >= nrow(y) - 250)
  tbl <- enrich.list(lo$protein_id, all.gor$protein_id, include.iea=TRUE)
  assign('lo.aga.branch', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)

  hi <- subset(y, rank.aga.clade <= 250)
  tbl <- enrich.list(hi$protein_id, all.gor$protein_id, include.iea=TRUE)
  assign('hi.aga.clade', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)

  lo <- subset(y, rank.aga.clade >= nrow(y) - 250)
  tbl <- enrich.list(lo$protein_id, all.gor$protein_id, include.iea=TRUE)
  assign('lo.aga.clade', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)

  ### Parallel enrichments.
  rank.threshold <- 800

  parallel.hg <- subset(y,
    rank.g.lrt <= rank.threshold & rank.h <= rank.threshold & rank.c > rank.threshold
  )
  print(nrow(parallel.hg))
  tbl <- enrich.list(parallel.hg$protein_id, all.gor$protein_id, include.iea=TRUE)
  assign('parallel.hg', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)

  parallel.cg <- subset(y,
    rank.g.lrt <= rank.threshold & rank.c <= rank.threshold & rank.h > rank.threshold
  )
  print(nrow(parallel.cg))
  tbl <- enrich.list(parallel.cg$protein_id, all.gor$protein_id, include.iea=TRUE)
  assign('parallel.cg', process.tbl(x, last.godata, tbl), envir=.GlobalEnv)

}

cum.dist.plots <- function(y, fields=c('lrt.5', 'lrt.4', 'lrt.3'), x.lim=NA) {

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
    x$lrt.g <- -x$lrt.5
    x$lrt.c <- -x$lrt.2
    x$lrt.h <- -x$lrt.1
  } else {
    x$lrt.g <- -x$pval.sw.gorilla
    x$lrt.c <- -x$pval.sw.chimp
    x$lrt.h <- -x$pval.sw.human
  }

  g <- x$lrt.g
  c <- x$lrt.c
  h <- x$lrt.h

  #x.up <- x[order(x$lrt.g),]
  #x.up <- subset(x, lrt.g > 0)

  all <- nrow(x.up)
  print("  calculating...")
  print(nrow(x.up))
  for (i in 1:nrow(x.up)) {
    cur.threshold <- x.up[i, 'lrt.g']
    #cur.sub <- subset(x.up, lrt.g >= cur.threshold)
    cur.sub <- x.up$lrt.g >= cur.threshold
    gor.total <- sum(cur.sub)
    x.up[i, 'i'] <- all - i
    n.chimp <- sum(cur.sub & x.up$lrt.c >= cur.threshold)
    n.human <- sum(cur.sub & x.up$lrt.h >= cur.threshold)
    x.up[i, 'n.chimp'] <- n.chimp / all
    x.up[i, 'n.human'] <- n.human / all
    x.up[i, 'f.chimp'] <- n.chimp / gor.total
    x.up[i, 'f.human'] <- n.human / gor.total
  }
  print("  done!")
  
  df <- data.frame(grp='GH / G_tot', yy=x.up$n.human, xx=x.up$lrt.g, xi=x.up$i)
  df <- rbind(df, data.frame(grp='GC / G_tot', yy=x.up$n.chimp, xx=x.up$lrt.g, xi=x.up$i))
  df <- rbind(df, data.frame(grp='GH / G_cur', yy=x.up$f.human, xx=x.up$lrt.g, xi=x.up$i))
  df <- rbind(df, data.frame(grp='GC / G_cur', yy=x.up$f.chimp, xx=x.up$lrt.g, xi=x.up$i))

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
  pdf(file="~/src/greg-ensembl/projects/gorilla/ns.cum.pdf", height=3, width=4)
  x$h.f.nsyn <- x$h.f.nsyn + x$hc.f.nsyn
  x$c.f.nsyn <- x$c.f.nsyn + x$hc.f.nsyn
  p <- cum.dist.plots(x, fields=c('g.f.nsyn', 'h.f.nsyn', 'c.f.nsyn'), x.lim=c(0, 0.04))
  p <- p + theme_bw()
  print(p)
  dev.off()

  pdf(file="~/src/greg-ensembl/projects/gorilla/lrt.cum.pdf", height=4, width=6)
  x$Gorilla <- x$lrt.5
  x$Human <- x$lrt.3
  x$Chimpanzee <- x$lrt.4
  p <- cum.dist.plots(x, fields=c('Human', 'Chimpanzee', 'Gorilla'), x.lim=c(-10, 10))
  p <- p + theme_bw()
  print(p)
  dev.off()

  pdf(file="~/src/greg-ensembl/projects/gorilla/lrt2.cum.pdf", height=4, width=6)
  x$AGA.Branch <- x$lrt.6
  x$AGA.Clade <- x$lrt.7
  p <- cum.dist.plots(x, fields=c('AGA.Branch', 'AGA.Clade'), x.lim=c(-10, 10))
  p <- p + theme_bw()
  print(p)
  dev.off()

  pdf(file="~/src/greg-ensembl/projects/gorilla/lrt3.cum.pdf", height=4, width=6)
  x$Human <- x$lrt.1
  x$Chimpanzee <- x$lrt.2
  p <- cum.dist.plots(x, fields=c('Human', 'Chimpanzee'), x.lim=c(-10, 10))
  p <- p + theme_bw()
  print(p)
  dev.off()

  pdf(file="~/src/greg-ensembl/projects/gorilla/lrt4.cum.pdf", height=4, width=6)
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
      nodeSize = 8,
      description = 'asdf'
    )
    assign("generic.godata", GOdata, envir=.GlobalEnv)
  }

  term.groups <- get.term.groups()

  term.groups['all'] <- 'all'

  rank.fields <- c(
    'all',
#    'rank.pval.sw.gorilla', 'rank.pval.sw.human', 'rank.pval.sw.chimp',
#    'rank.lrt.5', 'rank.lrt.4', 'rank.lrt.3'
     'slr_dnds', 'm0_dnds'
  )

  terms.df <- ldply(term.groups, function(term) {
    if (term != 'all') {
      genes <- getGenesForTerms(generic.godata, term)
      row.indices <- y$protein_id %in% genes
      z <- y[row.indices,]
    } else {
      z <- y    
    }
    df <- data.frame(nrow=nrow(z))

    data.for.field <- function(fld, all=F) {
      if (all) {
        df <- y
      } else {
        df <- z
      }
      if (grepl('rank', fld)) {
        df$rank_tmp <- -df[, fld]
      } else {
        df$rank_tmp <- df[, fld]
      }
      return(df)
    }

    for (field.x in rank.fields) {
      if (field.x == 'all'){
        next 
      }
      df.x <- data.for.field(field.x)
      for (field.y in rank.fields) {
        if (field.x == field.y) {
          next
        }
        if (field.y == 'all') {
          df.y <- data.for.field(field.x, all=T)
        } else {
          df.y <- data.for.field(field.y)
        }

        if (nrow(z) > 0) {
          ks.up <- ks.test(df.x$rank_tmp, df.y$rank_tmp, alternative='l')
          ks.down <- ks.test(df.x$rank_tmp, df.y$rank_tmp, alternative='g')
          ws <- wilcox.test(df.x$rank_tmp, df.y$rank_tmp, alternative='g')
          ks.up.p <- ks.up$p.value
          ks.down.p <- ks.down$p.value
          ws.p <- ws$p.value
        } else {
          ks.up.p <- 1
          ks.down.p <- 1
          ws.p <- 1
        }
        field.up <- paste(field.x, '-', field.y, '.up', sep='')
        field.down <- paste(field.x, '-', field.y, '.down', sep='')
        df[field.up] <- ks.up.p
        df[field.down] <- ks.down.p
      }
    }
    return(df)
  })
  print(terms.df)
  return()

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
  indices <- x$protein_id %in% ids
  return(x[indices, 'gene_name'])
}


ils.densities <- function() {
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
    assign("dbname", 'gj1_ils', envir=.GlobalEnv)
    source("~/src/greg-ensembl/scripts/collect_sitewise.R")
    orig.chunks <- get.vector(con, 'select * from chunks_10')
    assign("orig.chunks", orig.chunks, envir=.GlobalEnv)
  }
  chunks <- orig.chunks

  chunks$dist <- abs(chunks[, distance_field])
  chunks$dnds <- chunks[, dnds_field]

  chunks <- subset(chunks, !(chr_name %in% c('MT', 'X', 'Y')))
#  chunks <- subset(chunks, dnds > -1)
  chunks <- subset(chunks, !is.na(dist))

  x.labels <- c()
  if (bin_method == 'quantile') {
    q.mids <- c()
    qs <- quantile(chunks$dist,c(0.2, 0.4, 0.6, 0.8, 1.0),na.rm=T)
    print(qs)
    for (i in 1:length(qs)) {
      bin.hi <- qs[i]
      if (i == 1) {
        bin.lo <- min(chunks$dist)
      } else {
        bin.lo <- qs[i-1]
      }
      print(paste(i,bin.lo,bin.hi))
      chunks[chunks$dist >= bin.lo & chunks$dist <= bin.hi , 'x_bin'] <- i
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

  #chunks$value <- chunk.a(chunks)
  chunks$value <- chunk.b(chunks)
  #chunks$value <- chunk.c(chunks)

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
  p <- p + scale_y_continuous("Scaled ILS density")
  p <- p + scale_x_continuous(name="Distance Bin", breaks=c(1:length(x.labels)), labels=x.labels)
  return(p)
}

chunk.a <- function(x) {
  ((x$n_100 + x$n_010)/2 + x$n_110) / x$n_001
}

chunk.b <- function(x) {
  (x$n_101+x$n_011) / ((x$n_100 + x$n_010)/2 + x$n_110 + x$n_001)
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