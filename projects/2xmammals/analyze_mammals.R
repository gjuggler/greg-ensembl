library(plyr)
library(ggplot2)
source("~/src/greg-ensembl/scripts/mysql_functions.R")

do.all <- function() {
  
}

get.genes <- function() {
  if (!exists('genes.df', envir=.GlobalEnv)) {
    con <- connect('gj1_2x_63'); 
    df <- dbGetQuery(con, 'select * from genes')
    assign('genes.df', df, envir=.GlobalEnv)
  }
  df <- get('genes.df', envir=.GlobalEnv)
  df
}

get.genes.split <- function() {
  if (!exists('genes.split.df', envir=.GlobalEnv)) {
    genes <- get.genes()

    out.df <- data.frame()
    clades <- c('p', 'g', 'l', 'm', 'e')
    
    for (i in 1:4) {
      clade <- clades[i]
      bl.s <- paste(clade, '_slr_total_length', sep='')
      sites.s <- paste(clade, '_sites', sep='')
      dnds.s <- paste(clade, '_slr_dnds', sep='')
      ecdf.f <- ecdf(genes[, dnds.s])
      cum <- ecdf.f(genes[, dnds.s])
      cur.genes <- genes
      cur.genes$parameter_set_id <- i
      cur.genes$clade <- clade
      cur.genes$dnds <- genes[, dnds.s]
      cur.genes$sites <- genes[, sites.s]
      cur.genes$bl <- genes[, bl.s]
      cur.genes$cum_dnds <- cum
      out.df <- rbind(out.df, cur.genes)
    }
    assign('genes.split.df', out.df, envir=.GlobalEnv)
  }
  df <- get('genes.split.df', envir=.GlobalEnv)
  df
}

get.sites <- function() {
  if (!exists('sites.df', envir=.GlobalEnv)) {
    con <- connect('gj1_2x_63'); 
    df <- dbGetQuery(con, 'select * from sites limit 4000000')

    # Remove all gaps and 'single char' sites
    df <- subset(df, !(note %in% c('all_gaps', 'single_char')))

    # Add p-values
    df[, 'pval'] <- 1 - pchisq(abs(df[, 'lrt_stat']), 1)
    
    # Add p-set factor.
    pset.levels <- c('Primates', 'Glires', 'Laurasiatheria', 'Eutheria', 'Mammalia')
    df$clade <- factor(df$parameter_set_id, levels=1:5, labels=pset.levels)

    # Add bl quantiles grouped by pset.
    df <- ddply(df, .(parameter_set_id), function(x) {
      fn <- ecdf(x$nongap_bl)
      x$bl_z <- fn(x$nongap_bl)
      x$bl_factor <- round_any(x$bl_z, 0.2)
      x
    })
    df$bl_factor <- as.factor(df$bl_factor)

    # Add lrt z-scores
    df <- ddply(df, .(parameter_set_id), function(x) {
      fn <- ecdf(x$lrt_stat)
      x$lrt_z <- fn(x$lrt_stat)
      x
    })

    # Try a two-sided LRT z-score.
    df <- ddply(df, .(parameter_set_id), function(x) {
      abv.i <- x$lrt_stat >= 0
      abv <- subset(x, lrt_stat >= 0)
      blw <- subset(x, lrt_stat < 0)
      fn.abv <- ecdf(abv$lrt_stat)
      fn.blw <- ecdf(blw$lrt_stat)
      abv.z <- fn.abv(x$lrt_stat)
      blw.z <- fn.blw(x$lrt_stat)
      
      lrt_zz <- c()
      lrt_zz[abv.i] <- .5 + abv.z[abv.i] / 2
      lrt_zz[!abv.i] <- blw.z[!abv.i] / 2
      x$lrt_zz <- lrt_zz
      print(summary(lrt_zz))
      x
    })

    assign('sites.df', df, envir=.GlobalEnv)
  }
  df <- get('sites.df', envir=.GlobalEnv)
  df
}

summarize.sites <- function() {
  data <- get.sites()

  summ <- ddply(data, .(parameter_set_id), function(sites) {
    print(quantile(sites$nongap_bl))
    print(quantile(sites$ncod))

    sub.above.one <- subset(sites, omega > 1)
    pos.sites <- subset(sub.above.one, pval < 0.05)
    pos.domains <- subset(pos.sites, !is.na(pfam_domain))
  
    n.total.domains <- nrow(unique(sites[,c('data_id','pfam_domain')]))
    n.total.domain.types <- length(unique(sites$pfam_domain))
    n.total.genes <- length(unique(sites$data_id))

    n.pos.domains <- nrow(unique(pos.domains[,c('data_id','pfam_domain')]))
    n.pos.domain.types <- length(unique(pos.domains$pfam_domain))
    n.pos.genes <- length(unique(pos.sites$data_id))

    n.pos.sites <- nrow(pos.sites)
    n.sites <- nrow(sites)  
    f.pos.sites <- n.pos.sites / n.sites

    pos.sites <- subset(sites, omega > 1 & pval < 0.05)
    neg.sites <- subset(sites, omega < 1 & pval < 0.05)
    neutral.sites <- subset(sites, pval > 0.05)

    `f.less.0.5` <- nrow(subset(sites, omega < 0.5)) / nrow(sites)
    `f.less.1` <- nrow(subset(sites, omega < 1)) / nrow(sites)
    `f.gt.1` <- nrow(subset(sites, omega > 1)) / nrow(sites)
    `f.gt.1.5` <- nrow(subset(sites, omega > 1.5)) / nrow(sites)

    df.out <- data.frame(
      n.sites = n.sites,
      'f.less.0.5' = `f.less.0.5`,
      'f.less.1' = `f.less.1`,
      'f.gt.1' = `f.gt.1`,
      'f.gt.1.5' = `f.gt.1.5`,
  
      n.total.domains = n.total.domains,
      n.total.domain.types = n.total.domain.types,
      n.total.genes = n.total.genes,
    
      n.pos.domains = n.pos.domains,
      n.pos.domain.types = n.pos.domain.types,
      n.pos.genes = n.pos.genes,

      pos.n = n.pos.sites,
      pos.f = nrow(pos.sites)/n.sites,
      neg.f = nrow(neg.sites)/n.sites,
      neutral.f = nrow(neutral.sites)/n.sites
    )
    df.out
  })
  print(summ)
}

plot.sites.qc <- function() {
  sites <- get.sites()
  sites <- subset(sites, !(note %in% c('all_gaps', 'single_char')))

  sites.bl <- sites
  sites.bl$nongap_bl <- pmin(sites.bl$nongap_bl, 20)
  p.bl <- ggplot(sites.bl, aes(x=nongap_bl))
  p.bl <- p.bl + geom_histogram(binwidth=0.2) + facet_grid(clade ~ .)

  sites.lrt <- sites
  sites.lrt$lrt_stat <- pmin(sites.lrt$lrt_stat, 10)
  sites.lrt$lrt_stat <- pmax(sites.lrt$lrt_stat, -50)
  p.lrt <- ggplot(sites.lrt, aes(x=lrt_stat))
  p.lrt <- p.lrt + geom_histogram(binwidth=0.2) + facet_grid(clade ~ .)

  pdf(file="mammals_qc_histograms.pdf", width=12, height=6)
  vplayout(2, 1)
  print(p.bl, vp=subplot(1,1))
  print(p.lrt, vp=subplot(2,1))
  dev.off()

}

plot.subclade.distributions <- function() {
  sites <- get.sites()
  sites <- filter.bl.sites(sites)

  custom.merge <- function(s.a, s.b, xtras) {
    sel.cols <- c(xtras, 'data_id', 'aln_position', 'clade')
    s.a <- s.a[, c('lrt_factor', sel.cols)]
    s.b <- s.b[, sel.cols]
    merged.s <- merge(s.a, s.b, by=c('data_id', 'aln_position'))
    merged.s
  }

  if (!exists('out.df', envir=.GlobalEnv)) {
    # Separate sites by Eutherian LRT constraint, and plot distributions
    # of sub-clade dN/dS values.
    euth.s <- subset(sites, clade == "Eutheria")  
    euth.factor <- as.factor(cut_interval(euth.s$lrt_z, 6))
    euth.s$lrt_factor <- euth.factor
    euth.dnds <- dlply(euth.s, .(lrt_factor), function(x) {
      mean(x$omega)
    })
    lbls <- paste(levels(euth.s$lrt_factor), "\n", sprintf("%.2g",euth.dnds))
    euth.s$lrt_factor <- factor(euth.s$lrt_factor, levels=levels(euth.s$lrt_factor), labels=lbls)
    print(summary(euth.s$lrt_factor))

    tbl.df <- data.frame()
    out.df <- data.frame()
    for (i in 1:3) {
      mgd <- custom.merge(euth.s, subset(sites, parameter_set_id==i), 'omega')
      mgd.cum <- ddply(mgd, .(lrt_factor), function(x) {
        ecdf.f <- ecdf(x$omega.y)
        x$cum <- ecdf.f(x$omega.y)
        x
      })
      print(head(mgd.cum, n=2))
      out.df <- rbind(out.df, mgd.cum)
    }
    assign('out.df', out.df, envir=.GlobalEnv)
  }

  out.df$clade.y <- factor(out.df$clade.y)
  out.df$omega.y <- pmin(out.df$omega.y, 3)

  tbl.df <- ddply(out.df, .(clade.y, lrt_factor), function(x) {
    blw.one <- 
    strong.neg <- nrow(subset(x, omega.y < 0.1)) / nrow(x)
    abv.one <- nrow(subset(x, omega.y > 1)) / nrow(x)
    data.frame(
      `blw.1` = nrow(subset(x, omega.y < 1)) / nrow(x),
      `blw.5` = nrow(subset(x, omega.y < 0.5)) / nrow(x),
      `blw.25` = nrow(subset(x, omega.y < 0.25)) / nrow(x),
      `abv.1` = nrow(subset(x, omega.y > 1)) / nrow(x),
      `abv.2` = nrow(subset(x, omega.y > 2)) / nrow(x)
    )
  })
  tbl.df <- tbl.df[order(tbl.df[, 'lrt_factor'], tbl.df[, 'clade.y']), ]
  print(tbl.df)

  p <- ggplot(out.df, aes(x=omega.y, fill=clade.y))
  p <- p + geom_histogram(binwidth=0.1, position="dodge")
  p <- p + scale_fill_discrete("Clade")
  p <- p + scale_y_log10("Count")
  p <- p + scale_x_continuous("dN/dS")
  p <- p + facet_grid(lrt_factor ~ ., scales="free_y")
  
  p <- p + opts(
    strip.text.y = theme_text(angle=0, hjust=0, size=12)
  )
  
  pdf(file="mammals_subclade_histograms.pdf", width=8, height=6)
  print(p)
  dev.off()

  p <- ggplot(out.df, aes(x=omega.y, y=cum, colour=clade.y))
  p <- p + geom_line()
  p <- p + scale_colour_discrete("Clade")
  p <- p + scale_y_continuous("Cumulative Sites", limits=c(0, 1))
  p <- p + scale_x_continuous("dN/dS")
  p <- p + facet_grid(. ~ lrt_factor)
  
  png(file="mammals_subclade_cum_histograms.png", width=1200, height=600)
  print(p)
  dev.off()
}

plot.sites.clade.correlations <- function() {
  sites <- get.sites()

  plot.lrt <- function(pset.a, pset.b) {
    mrgd <- merge.clades(sites, pset.a, pset.b, 'lrt_stat')
    clade.x <- mrgd[1, 'clade.x']
    clade.y <- mrgd[1, 'clade.y']
    p <- ggplot(mrgd, aes(x=lrt_stat.x, y=lrt_stat.y))
    p <- p + geom_point(size=0.25, alpha=0.5)
    p <- p + geom_abline(slope=1, intercept=0, linetype='dashed')
 
    scale.lim <- c(-20, 10)
    if (pset.a == 4 || pset.a == 5) {
      scale.lim <- c(-70, 10)
    } 
    p <- p + scale_x_continuous(clade.x, limits=scale.lim)
    p <- p + scale_y_continuous(clade.y, limits=scale.lim)
    p
  }

  plot.omega <- function(pset.a, pset.b) {
    mrgd <- merge.clades(sites, pset.a, pset.b, 'omega')
    clade.x <- mrgd[1, 'clade.x']
    clade.y <- mrgd[1, 'clade.y']
    p <- ggplot(mrgd, aes(x=omega.x+0.01, y=omega.y+0.01))
    p <- p + geom_point(size=0.25, alpha=0.5)
    p <- p + geom_abline(slope=1, intercept=0, linetype='dashed')
 
    scale.lim <- log10(c(0.01, 3))
    p <- p + scale_x_log10(clade.x)
    p <- p + scale_y_log10(clade.y)
    p
  }

  plot.z <- function(pset.a, pset.b) {
    mrgd <- merge.clades(sites, pset.a, pset.b, 'lrt_zz')
    clade.x <- mrgd[1, 'clade.x']
    clade.y <- mrgd[1, 'clade.y']
    p <- ggplot(mrgd, aes(x=lrt_zz.x, y=lrt_zz.y))
    p <- p + geom_point(size=0.25, alpha=0.5)
    p <- p + geom_abline(slope=1, intercept=0, linetype='dashed')
 
    p <- p + scale_x_continuous(clade.x)
    p <- p + scale_y_continuous(clade.y)
    p
  }

  plot.old.z <- function(pset.a, pset.b) {
    mrgd <- merge.clades(sites, pset.a, pset.b, 'lrt_z')
    clade.x <- mrgd[1, 'clade.x']
    clade.y <- mrgd[1, 'clade.y']
    p <- ggplot(mrgd, aes(x=lrt_z.x, y=lrt_z.y))
    p <- p + geom_point(size=0.25, alpha=0.5)
    p <- p + geom_abline(slope=1, intercept=0, linetype='dashed')
 
    p <- p + scale_x_continuous(clade.x)
    p <- p + scale_y_continuous(clade.y)
    p
  }

  plot.bl <- function(pset.a, pset.b) {
    mrgd <- merge.clades(sites, pset.a, pset.b, 'nongap_bl')
    clade.x <- mrgd[1, 'clade.x']
    clade.y <- mrgd[1, 'clade.y']
    p <- ggplot(mrgd, aes(x=nongap_bl.x, y=nongap_bl.y))
    p <- p + geom_point(size=0.25, alpha=0.5)
    p <- p + geom_abline(slope=1, intercept=0, linetype='dashed')
 
    scale.lim <- c(0, 5)
    if (pset.a == 4 || pset.a == 5) {
      scale.lim <- c(0, 15)
    } 
    p <- p + scale_x_continuous(clade.x, limits=scale.lim)
    p <- p + scale_y_continuous(clade.y, limits=scale.lim)
    p
  }
  
  plot.corrs.f <- function(inner.function, filename) {
    const <- 200
    png(file=filename, width=4*3*const, height=1*3*const)
    vplayout(4, 1)
    print(inner.function(1, 2), vp=subplot(1, 1))
    print(inner.function(1, 3), vp=subplot(2, 1))
    print(inner.function(2, 3), vp=subplot(3, 1))
    print(inner.function(4, 5), vp=subplot(4, 1))
    dev.off()
  }

  plot.mamm.corrs.f <- function(inner.function, filename) {
    const <- 200
    png(file=filename, width=3*3*const, height=1*3*const)
    vplayout(4, 1)
    print(inner.function(1, 4), vp=subplot(1, 1))
    print(inner.function(2, 4), vp=subplot(2, 1))
    print(inner.function(3, 4), vp=subplot(3, 1))
    print(inner.function(5, 4), vp=subplot(4, 1))
    dev.off()
  }


#   plot.corrs.f(plot.z, "mammals_corrs_z.png")
#   plot.mamm.corrs.f(plot.z, "mammals_corrs_z_euth.png")

#   plot.corrs.f(plot.old.z, "mammals_corrs_oldz.png")
#   plot.mamm.corrs.f(plot.old.z, "mammals_corrs_oldz_euth.png")

   plot.mamm.corrs.f(plot.lrt, "mammals_corrs_lrt_euth.png")
#  plot.corrs.f(plot.lrt, "mammals_corrs_lrt.png")
#  plot.corrs.f(plot.omega, "mammals_corrs_omega.png")
#  plot.corrs.f(plot.bl, "mammals_corrs_bl.png")

}

plot.sites.lrt.omega <- function() {
  sites <- get.sites()

  sites$lrt_stat <- pmin(sites$lrt_stat, 20)
  sites$lrt_stat <- pmax(sites$lrt_stat, -60)
  sites$omega <- pmin(sites$omega, 20)

  library(RColorBrewer)
  p <- ggplot(sites, aes(y=lrt_stat, x=omega+0.01, colour=as.factor(bl_factor)))
  p <- p + geom_point(size=1, alpha=0.5)
  p <- p + scale_y_continuous("LRT statistic")
  p <- p + scale_colour_discrete("Branch length quantile")
  p <- p + scale_x_log10("dN/dS")
  p <- p + facet_wrap(~ clade)

  png(file="mammals_lrt_omega_scatter.png", width=1200, height=1200)
  print(p)
  dev.off()
}

plot.sites.boxplots <- function() {
  sites <- get.sites()
  sites <- subset(sites, !(note %in% c('all_gaps', 'single_char')))

  sites$nongap_bl <- pmin(sites$nongap_bl, 7)
  sites$lrt_stat <- pmin(sites$lrt_stat, 10)
  sites$lrt_stat <- pmax(sites$lrt_stat, -60)

  sites$nongap_y <- round_any(sites$nongap_bl, 2)

  p <- ggplot(sites, aes(x=factor(nongap_y), y=lrt_stat))
  p <- p + geom_boxplot(aes(fill=factor(clade)), outlier.shape=NA)

  pdf(file="mammals_bl_boxplots.pdf")
  print(p)
  dev.off()
}

plot.sites.histograms <- function() {
  sites <- get.sites()

  # Subset the data
  sites <- subset(sites, !(note %in% c('all_gaps', 'single_char')))
  sites <- subset(sites, nongap_bl > 0.5)
  sites <- subset(sites, is.na(random))
  zero.sites <- subset(sites, omega == 0)
  sites <- subset(sites, omega > 0)

  sites$omega <- pmin(sites$omega, 3)
  print(head(sites, n=3))
  
  p <- ggplot(sites, aes(x=omega, colour=lrt_stat))
  p <- p + geom_histogram(binwidth=0.02)
  p <- p + scale_colour_gradient2("LRT statistic", low='blue', mid='white', high='red', limits=c(-4,4))
  p <- p + facet_grid(clade ~ .)

  q <- ggplot(zero.sites, aes(x=note))
  q <- q + geom_bar()
  q <- q + facet_grid(clade ~ .)

  pdf(file="mammals_histogram.pdf", width=8, height=12)
  vplayout(4, 1)
  print(q, vp=subplot(1,1))
  print(p, vp=subplot(2:4,1))
  dev.off()
}

plot.genes.dnds <- function() {
  genes <- get.genes()

  out.df <- data.frame()
  for (clade in c('p', 'g', 'l', 'm')) {
    fld.s <- paste(clade, '_slr_dnds', sep='')
    ecdf.f <- ecdf(genes[, fld.s])
    cum <- ecdf.f(genes[, fld.s])
    out.df <- rbind(out.df, data.frame(
      clade = clade,
      dnds = genes[, fld.s],
      cum = cum
    ))
  }

  p <- ggplot(out.df, aes(x=dnds, fill=clade))
  p <- p + geom_histogram(binwidth=0.1, position="dodge")
  pdf(file="mammals_genes_dnds.pdf")
  print(p)
  dev.off()

  p <- ggplot(out.df, aes(x=dnds, y=cum, colour=clade))
  p <- p + geom_line()
  pdf(file="mammals_genes_cum_dnds.pdf")
  print(p)
  dev.off()
}

plot.genes.maxlrt <- function() {
  genes <- get.genes.split()
  sites <- get.sites()
  sites <- filter.bl.sites(sites)

  max.lrt <- ddply(sites, .(data_id, parameter_set_id), function(x) {
    x <- subset(x, omega < 99)
    data.frame(
      max_lrt = max(x$lrt_stat)
    )
  })
  print(head(max.lrt))
  genes.max <- merge(genes, max.lrt)

  genes.clade <- ddply(genes.max, .(clade), summarize,
    max_lrt = max_lrt,
    ecdf = ecdf(max_lrt)(max_lrt)
  )

  p <- ggplot(genes.clade, aes(x=max_lrt, y=ecdf, colour=clade))
  p <- p + geom_line()
  pdf(file="mammals_genes_maxlrt.pdf")
  print(p)
  dev.off()
}

merge.clades <- function(sites, pset.a, pset.b, xtra) {
  sel.cols <- c(xtra, 'data_id', 'aln_position', 'clade')
  s.a <- subset(sites, parameter_set_id == pset.a, select=sel.cols)
  s.b <- subset(sites, parameter_set_id == pset.b, select=sel.cols)

  merged.s <- merge(s.a, s.b, by=c('data_id', 'aln_position'))
  merged.s
}

filter.bl.sites <- function(sites) {
  param.sets <- sort(unique(sites$parameter_set_id)  )
  for (pset in param.sets) {
    pset.bls <- sites[sites$parameter_set_id == pset, 'nongap_bl']
    qnts <- quantile(pset.bls, probs=c(0, 0.1, 0.5, 0.9, 1))
    sites <- subset(sites, !(parameter_set_id == pset & nongap_bl < qnts[2]))
    sites <- subset(sites, !(parameter_set_id == pset & nongap_bl > qnts[4]))
  }
  sites
}

filter.random.sites <- function(sites) {
  subset(sites, is.na(random))
}

filter.const.sites <- function(sites) {
  subset(sites, !(note %in% c('constant')))
}

filter.pfam <- function(sites) {
  subset(sites, !is.na(pfam_domain))
}

subplot <- function(x, y) viewport(layout.pos.col=x, layout.pos.row=y)
vplayout <- function(x, y) {
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(y,x)))
}
