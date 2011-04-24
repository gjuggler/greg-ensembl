library(plyr)
library(ggplot2)
library(boot)
library(Hmisc)

dbname <- 'gj1_gorilla'

pst <- function(...) { return(paste(..., sep=''))}

get.pvals <- function() {
  source("~/src/greg-ensembl/scripts/collect_sitewise.R")
  if (!exists("genes")) {
    genes <- get.vector(con, 'select * from genes')
    assign("genes", genes, envir=.GlobalEnv)
  }
  if (!exists("subs")) {
    subs <- get.vector(con, 'select * from subs')
    assign("subs", subs, envir=.GlobalEnv)
  }

  g <- genes
  s <- subs

  print("  merging genes & subs")
  merged <- merge(g[, c('data_id', 'aln_length', 'm0_dnds', 'slr_dnds')], s, by=c('data_id'))
  print("  binomial testing...")
  binom.test <- ddply(merged, .(data_id), function(x) {
    aln.length <- x[1, 'aln_length']
    cur.dnds <- x[1, 'slr_dnds']
    data.id <- x[1, 'data_id']
    y <- subset(x, mut_nsyn==1)

    if (nrow(subset(x, taxon_id==9598)) == 0) {
      has.chimp <- FALSE
    } else {
      has.chimp <- TRUE
    }

    hc.subs <- subs[grep("ENSP0.*,ENSPTRP0.*$", y$leaves_beneath),]
    h.subs <- subset(y, taxon_id==9606)
    c.subs <- subset(y, taxon_id==9598)
    g.subs <- subset(y, taxon_id==9593)
    n.hc <- nrow(hc.subs)
    n.h <- nrow(h.subs)
    n.c <- nrow(c.subs)
    n.g <- nrow(g.subs)

    if (has.chimp) {
      p.nsyn <- (n.hc + (n.h + n.c)/2) / aln.length
    } else {
      p.nsyn <- (n.hc + n.h) / aln.length
    }
    
    mut.rate <- 0.5e-9
    mut.per.site <- mut.rate * 10e6
    p.per.site <- mut.per.site * cur.dnds
    ns.sites <- aln.length * 3 * 2/3

    res <- binom.test(n.g, ns.sites, p.per.site, "greater")
    pval <- res$p.value
    data.frame(
      mut.per.site=mut.per.site,
      p.per.site=p.per.site,
      binom.pval=pval,
      gor.ns = nrow(subset(x, mut_nsyn==1 & taxon_id==9593)),
      gor.s = nrow(subset(x, mut_syn==1 & taxon_id==9593))
    )
  })
  print(head(binom.test))
  g <- merge(g, binom.test, by=c('data_id'))
  print(head(g))
  print("  done!")  

  for (cur.type in c('all', 'noC', 'noH')) {
    cur.stats <- subset(g, select=c('data_id', 'gene_name', 'protein_id', 'm0_dnds', 'slr_dnds', 'data_prefix', 'gor.ns', 'gor.s', 'binom.pval', 'p.per.site', 'aln_length'))
    print(str(cur.stats))
    # Branch models.
    for (i in 1:10) {
      lrt.key <- paste('lrt','.',i,sep='')
      pval.key <- paste('pval','.',i,sep='')
      print(paste(cur.type, lrt.key))

      null.lnl <- g[, pst(cur.type, '_', '0_lnL')]
      null.omega <- g[, pst(cur.type, '_', '0_w_0')]
      
      exp.lnl <- g[, pst(cur.type, '_', i, '_lnL')]
      w.b.key <- pst(cur.type, '_', i, '_w_1')
      exp.omega.a <- g[, pst(cur.type, '_', i, '_w_0')]
      exp.omega.b <- exp.omega.a
      if (w.b.key %in% colnames(genes)) {
        exp.omega.b <- g[, w.b.key]
      }
      signed.lrt <- 2 * pmax(0, exp.lnl - null.lnl) * sign(exp.omega.b - exp.omega.a)
      pval <- 1 - pchisq( 2 * (exp.lnl - null.lnl), df=1)

      cur.stats[, lrt.key] <- signed.lrt
      cur.stats[, pval.key] <- pval
    }

    # Branch-sites models.
    for (i in 11:17) {    
      lrt.key <- paste('lrt','.',i,sep='')
      pval.key <- paste('pval','.',i,sep='')
      print(paste(cur.type, lrt.key))

      alt.lnl <- g[, pst(cur.type, '_', i, '_', 'alt_lnL')]
      null.lnl <- g[, pst(cur.type, '_', i, '_', 'null_lnL')]

      alt.lnl <- as.numeric(alt.lnl)
      null.lnl <- as.numeric(null.lnl)

      lrt <- 2 * pmax(0, alt.lnl - null.lnl)
      pval <- 1 - pchisq( 2 * (alt.lnl - null.lnl), df=1)

      cur.stats[, lrt.key] <- lrt
      cur.stats[, pval.key] <- pval
    }

    cur.stats <- cur.stats[order(-cur.stats$lrt.5),]
    assign(paste("stats.", cur.type, sep=''), cur.stats, envir=.GlobalEnv)
  }
}

subs.test <- function() {
  if (!exists('subs')) {
    source("~/src/greg-ensembl/scripts/collect_sitewise.R")
    subs <- get.vector(con, 'select * from subs')
    assign("subs", subs, envir=.GlobalEnv)
  }

  print(nrow(subs))
  x <- subs
#  x <- subset(subs, mut_syn == 1)
  print(nrow(x))

  x <- merge(genes[, c('data_id', 'slr_dnds', 'm0_dnds')], x, by=c('data_id'))

  x <- subset(x, m0_dnds >= 0.2)

  h.subs <- subset(x, taxon_id==9606)
  c.subs <- subset(x, taxon_id==9598)
  g.subs <- subset(x, taxon_id==9593)

  print(nrow(h.subs))
  print(nrow(c.subs))
  print(nrow(g.subs))

  hg <- merge(h.subs, g.subs, by=c('data_id', 'aln_pos'), suffixes=c('.h','.g'))
  cg <- merge(c.subs, g.subs, by=c('data_id', 'aln_pos'), suffixes=c('.c','.g'))
  hg.per.gene <- ddply(hg, .(data_id), nrow)
  print(hg.per.gene$V1)
  
  print(mean(x$slr_dnds))
  print(median(hg$slr_dnds.g))
  print(median(cg$slr_dnds.g))
}

ils.densities <- function() {
  pdf("ils_exon_mammal.pdf")
  p <- ils.density.plot('exon_dist', 'gene_slr', bin_method='window', relative=F)
  p <- p + opts(title="ILS density vs Exon Distance, Mammal dN/dS")
  print(p)

  p <- ils.density.plot('exon_dist', 'gene_slr', bin_method='window', relative=T)
  p <- p + opts(title="Relative ILS density vs Exon Distance, Mammal dN/dS")
  print(p)

  p <- ils.density.plot('exon_dist', 'gene_slr', bin_method='quantile', relative=F)
  p <- p + opts(title="ILS density vs Exon Distance, Mammal dN/dS")
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

  chunks$a <- chunk.a(chunks)
  chunks$b <- chunk.b(chunks)

  f.df <- function(x) {
    x <- subset(x, !is.na(a) & !is.na(b) & !is.infinite(a) & !is.infinite(b))

    bt.a <- smean.cl.boot(x$a)
    bt.b <- smean.cl.boot(x$b)

    #stderr <- function(x) sqrt(var(x)/length(x))
    #err.a <- stderr(x$a)
    #err.b <- stderr(x$b)

    return(data.frame(
      n = nrow(x),
      a = bt.a[1],
      a.lo = bt.a[2],
      a.hi = bt.a[3],
      b = bt.b[1],
      b.lo = bt.b[2],
      b.hi = bt.b[3],
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
       all.y.value <- cur.df.all$b
  
      # Normalize all dnds bins to that value.
      cur.df$b <- cur.df$b / all.y.value
      cur.df$b.hi <- cur.df$b.hi / all.y.value
      cur.df$b.lo <- cur.df$b.lo / all.y.value
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
      ymin=b.lo,
      ymax=b.hi
    ))
    p <- p + geom_rect(data=res.df, alpha=0.3, size=0, colour=NA,
      aes(
        xmin=as.numeric(x_bin)+ x_off - .15,
        xmax=as.numeric(x_bin)+ x_off + .15,
        ymin=b.lo,
        ymax=b.hi
      )
    )
  p <- p + geom_line(data=res.df, aes(x=as.numeric(x_bin),y=b))
  #p <- p + geom_errorbar(data=res.df, size=0.5, alpha=0.8, aes(x=x_bin, ymin=b.lo, ymax=b.hi))
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