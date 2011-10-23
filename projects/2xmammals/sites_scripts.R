uname  <- Sys.getenv("USER")
if (uname == 'gj1') {
  source("~/src/greg-ensembl/projects/2xmammals/analyze_mammals.R")
  source("~/src/greg-ensembl/projects/2xmammals/analyze_genes.R")
  source("~/src/greg-ensembl/projects/2xmammals/calc_recomb.R")
  source("~/src/greg-ensembl/scripts/xtable_utils.R")
} else {
  source("~/lib/greg-ensembl/projects/2xmammals/analyze_mammals.R")
  source("~/lib/greg-ensembl/projects/2xmammals/analyze_genes.R")
  source("~/lib/greg-ensembl/projects/2xmammals/calc_recomb.R")
  source("~/lib/greg-ensembl/scripts/xtable_utils.R")
}

args <- commandArgs(trailingOnly=T)

main <- function() {
  if (!is.na(args[1])) {
    fn.name <- args[1]
    print(paste("Calling", fn.name))

    # Get the LSF_JOBINDEX environment variable
    jobindex <- Sys.getenv('LSB_JOBINDEX')
    if (jobindex != '' && as.integer(jobindex) > 0) {
      args <- c(args, as.integer(jobindex))
    }
    print(args)

    error.f <- function(e) {
      print("###ERROR###")
      out.df <- data.frame(
        tstmp = paste(timestamp(), fn.name, collapse=' '),
        fn_name = as.character(fn.name),
        fn_args = paste(args, collapse=' '),
        error = as.character(e)
      )
      con <- connect(db())
      write.or.update(out.df, 'errors', con, 'tstmp')
      disconnect(con)
    }

    tryCatch(
      do.call(fn.name, as.list(args[2:length(args)])),
      error = error.f
    )
  }
}

collect_genes <- function(pset, filter='default', subset.index=NULL, test=F) {
  pset <- as.integer(pset)
  if (!is.null(subset.index)) {
    subset.index <- as.integer(subset.index)
  }
  test <- as.logical(test)

  orig.genes.f <- paste(scratch(), "genes_orig.Rdata", sep='')
  if (file.exists(orig.genes.f)) {
    print("Loading genes...")
    load(orig.genes.f)
  } else {
    con <- connect(db()); 
    cmd <- sprintf("select * from genes")
    genes <- dbGetQuery(con, cmd)
    disconnect(con)
    save(genes, file=orig.genes.f)
  }

  print("Processing...")
  process.genes(genes, pset, filter=filter, subset.index=subset.index, test=test)
  print("  done!")
}

collect_clusters <- function(taxon_id, test=F) {
  taxon_id <- as.integer(taxon_id)
  test <- as.logical(test)

  process.cluster(taxon_id, win.size=15, test=test)
}

collect_sites <- function(pset, filter='default', test=F) {
  pset <- as.integer(pset)
  test <- as.logical(test)

  limit.s <- ''
  if (test) {
    limit.s <- 'limit 500000'
  }

  sites.f <- get.sites.filename(pset=pset, filter=filter, test=test)
  if (file.exists(sites.f)) {
    print("processed sites file already exists -- not processing again!")
    return()
  }

  recomb.filters <- c('recomb_10k', 'recomb_100k', 'recomb_1mb')

  orig.sites.f <- get.sites.filename(pset=pset, filter='orig', test=test)
  if (file.exists(orig.sites.f)) {
    if (!(filter %in% recomb.filters) ) {
      load(orig.sites.f)
    }
  } else {
    if (filter != 'orig') {
      stop("Original sites haven't been gathered yet!")
    }
    if (test) {
      print("  limiting sites count for testing")
    }
    con <- connect(db()); 
    not.null.s <- 'and data_id is not null and aln_position is not null and parameter_set_id is not null'
    cmd <- sprintf("select * from sites where parameter_set_id=%d %s %s", pset, not.null.s, limit.s)
    sites <- dbGetQuery(con, cmd)
    disconnect(con)
    save(sites, file=orig.sites.f)
  }

  if (filter == 'orig') {
    return()
  }

  if (!file.exists(sites.f)) {
    if (filter %in% recomb.filters) {
      default.f <- get.sites.filename(pset=pset, filter='default', test=test)
      load(default.f)
    } else {
      sites <- process.sites(sites)
    }

    sites <- switch(filter,
      none = sites,
      default = filter.default(sites),
      stringent = filter.stringent(sites),
      pfam = filter.pfam(sites),
      pfam_stringent = filter.pfam.stringent(sites),
      clusters = filter.clusters(sites),
      clusters_inverse = filter.clusters(sites, return.inverse=T),
      recomb_10k = add.recomb.gc(sites, width=1e4, test=test),
      recomb_100k = add.recomb.gc(sites, width=1e5, test=test),
      recomb_1mb = add.recomb.gc(sites, width=1e6, test=test)
    )
    print("Saving file...")
    save(sites, file=sites.f)
  } else {
    print("  sites file already exists -- not processing again!")
  }
  print("  done!")
}

summary_table <- function(pset, test=F, filter='default') {
  test <- as.logical(test)
  pset <- as.numeric(pset)
  sites <- get.pset.sites(pset, filter=filter, test=test)
  df <- summarize.sites(sites, filter=filter)
  con <- connect(db())
  write.or.update(df, 'summaries', con, 'label')
  disconnect(con)
}

collect_gc <- function(window.width=10000) {
  window.width <- as.integer(window.width)
  get.gc.content(width=window.width)
}

poisson_sim <- function(dnds, neutral_thresh=1) {
  dnds <- as.numeric(dnds)
  neutral_thresh <- as.numeric(neutral_thresh)
  do.poisson.sim(dnds, neutral_thresh)
}

recomb_calc <- function(pset, filter='recomb', test=F) {
  pset <- as.integer(pset)
  test <- as.logical(test)
  calc.recomb(pset, filter=filter, test=test)
}

recomb_rate_calc <- function(sex, width, chr) {
  width <- as.integer(width)
  calc.recomb.rate(sex, width, chr)
}

genes_enrichment <- function(pset, filter, method, excl.iea) {
  pset <- as.integer(pset)
  excl.iea <- as.logical(excl.iea)

  calc.genes.enrichment(pset, filter, method, excl.iea)
}

cumulative_calc <- function(pset, filter, sort.f, direction, test=F) {
  pset <- as.integer(pset)
  test <- as.logical(test)
  cumulative.calc(pset, filter, sort.f, direction, test)
}

bl_pos_sel_breakdown <- function(pset, test=F) {
  pos_sel_summary <- function(x) {
    n_abv <- nrow(subset(x, omega > 1))
    n_blw <- nrow(subset(x, omega < 1))
    n_pos <- nrow(subset(x, lrt_stat > qchisq(0.95, df=1)))
    n <- nrow(x)
    med_ncod <- median(x$ncod)
    ncod_qs <- as.integer(quantile(x$ncod, c(0.25, 0.5, 0.75)))
    bl_qs <- quantile(x$nongap_bl, c(0.25, 0.5, 0.75))
    data.frame(
      bl_lo = bl_qs[1],
      bl_med = bl_qs[2],
      bl_hi = bl_qs[3],
      ncod_lo = ncod_qs[1],
      ncod_med = ncod_qs[2],
      ncod_hi = ncod_qs[3],
      below_one = n_blw / n * 100,
      above_one = n_abv / n * 100,
      pos_f = n_pos / n * 100
    )
  }

  sites <- get.pset.sites(6, test=test)
  sites$bl_quantile <- round_any(sites$bl_z, 0.1, f=ceiling)
  df <- ddply(sites, .(bl_quantile), pos_sel_summary)
  df <- cbind(data.frame(species='', stringsAsFactors=F), df)
  df[1, 'species'] <- pset.to.alias(6)

  sites <- get.pset.sites(1, test=test)
  sites$bl_quantile <- round_any(sites$bl_z, 0.1, f=ceiling)
  df2 <- ddply(sites, .(bl_quantile), pos_sel_summary)
  df2 <- cbind(data.frame(species='', stringsAsFactors=F), df2)
  df2[1, 'species'] <- pset.to.alias(1)

  xt <- xtable(df)
  xt <- color.columns(xt, c('above_one', 'below_one', 'pos_f'), log=T)

  xt2 <- xtable(df2)
  xt2 <- color.columns(xt2, c('above_one', 'below_one', 'pos_f'), log=T)

  xt.blank <- xt[1, ]
  xt.blank[, names(xt.blank)] <- ''

  xt.f <- scratch.f("bl_pos_sel_breakdown.txt")
  print.latex(xt, xt.f)
  print.latex(xt.blank, xt.f, append=TRUE)
  print.latex(xt2, xt.f, append=TRUE)

  #print(df)
}

fit_distr <- function(pset,
  filter=c('default', 'stringent', 'pfam'),
  distr=c('lnorm', 'gamma', 'beta', 'exp'), 
  use=c('ci', 'imputed', 'omega'),
  test=F
) {
  pset <- as.numeric(pset)
  distr = distr[1]
  use = use[1]
  test <- as.logical(test)
  sites <- get.pset.sites(pset, filter=filter, test=test)
  sites <- subset(sites, select=c('data_id', 'omega_lower', 'omega_upper', 'omega', 'parameter_set_id'))

  print("Fitting...")
  for (i in 1:50) {
#    con <- connect(db())
#    df <- dbGetQuery(con, sprintf("select * from fitdistr where 
#      pset=%s and filter='%s' and dist='%s' and use_type='%s' and i=%d", pset, filter, distr, use, i)
#    )
#    dbDisconnect(con)
#    if (nrow(df) == 0) {
      fit.sites(sites, distr=distr, filter=filter, use=use, i=i, write.to.table=T)
#    }
  }
}

pset_correlation <- function(filter, pset.a, pset.b, test=F) {
  test <- as.logical(test)
  pset.a <- as.integer(pset.a)
  pset.b <- as.integer(pset.b)

  x <- get.merged.sites(pset.a, pset.b,
    xtra=c('omega', 'lrt_stat', 'lrt_z', 'omega_lower', 'omega_upper'),
    test=test,
    filter=filter
  )
  
  # Only use sites with \omg less than 1 and greater than 0 for the PCA calculations.
  lm.max <- 1
  lm.min <- 0
  lm.sub <- subset(x, omega.x > lm.min & omega.y > lm.min & omega.x < lm.max & omega.y < lm.max)

  # The general approach here is to use PCA to identify the rotation of coordinates that
  # produces the greatest explanation of variance in the dataset. So we do that (using weights
  # correlating to the inverse of the confidence intervals) and then use the slope of the PCA
  # as the slope of correlation between the two datasets.
  require(aroma.light)
  require(boot)
  # see http://rss.acs.unt.edu/Rdoc/library/aroma.light/html/wpca.matrix.html
  # and http://books.google.co.uk/books?id=KIHuSXyhawEC&pg=PA340&lpg=PA340&dq=orthogonal+regression+pca&source=bl&ots=d4OnbDGR1m&sig=Vp_72VdOSYz9a5demOusToz4fTs&hl=en&ei=-s98TsmtHKmU0QWZ9OXwDw&sa=X&oi=book_result&ct=result&resnum=10&ved=0CHMQ6AEwCQ#v=onepage&q&f=false
  # for examples...
  iwpca.slope <- function(data, indices, center) {
    xs <- data[indices, 'omega.x']
    ys <- data[indices, 'omega.y']
    ws <- data[indices, 'wts']
    r <- wpca(cbind(xs,ys), center=center, w=ws)
    r.slope <- r$vt[1,2] / r$vt[1,1]
    return(r.slope)
  }
  iwpca.slope.int <- function(xs, ys, ws=NULL, center=T) {
    r <- wpca(cbind(xs,ys), center=center, w=ws)
    r.slope <- r$vt[1,2] / r$vt[1,1]
    r.int <- r$xMean[2] - r.slope * r$xMean[1]
    return(c(r.slope, r.int))
  }

  # Let's use the mean confidence interval width at each site to weight the PCA calculation.
  w.x <- lm.sub$omega_upper.x - lm.sub$omega_lower.x
  w.y <- lm.sub$omega_upper.y - lm.sub$omega_lower.y
  wt <- (w.x + w.y) / 2
  wt <- pmax(0.025, wt)
  wt <- pmin(20, wt)
  wts <- 1 / wt
  lm.sub$wts <- wts

  # Do a PCA without re-centering the data... this forces the 'root' of the best fit line to
  # be at the origin.
  zero.pca <- iwpca.slope.int(lm.sub$omega.x, lm.sub$omega.y, ws=wts, center=F)
  # Do the PCA with re-centering, allowing the line to sit anywhere.
  nonzero.pca <- iwpca.slope.int(lm.sub$omega.x, lm.sub$omega.y, ws=wts, center=T)

  # Run a number of bootstraps of the PCA calculation... we'll use this to get a confidence
  # interval on the slope of the PCA.
  print("  bootstrapping...")
  n.t <- 100
  if (test) {
    n.t <- 10
  }
  boot.data <- subset(lm.sub, select=c('omega.x', 'omega.y', 'wts'))
  zero.reps <- boot(lm.sub, iwpca.slope, R=n.t, center=F)
  nonzero.reps <- boot(lm.sub, iwpca.slope, R=n.t, center=T)
  zero.ci <- boot.ci(zero.reps, type='basic', conf=0.95)
  nonzero.ci <- boot.ci(nonzero.reps, type='basic', conf=0.95)

  print("  plotting...")
  clr.rmp <- colorRampPalette(c('white', 'blue', 'red'), bias=0.5)
  n.lvls <<- 10
  clr.val <- clr.rmp(n.lvls+1)
  brks <- 1 : (n.lvls)

  lm.f <- paste(scratch(), "omega_cor_", pset.a, '_', pset.b, ".pdf", sep='')
  p <- ggplot(lm.sub, aes(x=omega.x, y=omega.y))
  p <- p + theme_bw()
  p <- p + stat_density2d(geom='tile', aes(
    fill=factor(round((..density.. - min(..density..)) / diff(range(..density..)) * 10))
    ),
    contour=FALSE)
  p <- p + scale_fill_manual("Density", values=clr.val, breaks=brks)
  p <- p + geom_abline(slope=1, intercept=0, colour='black', linetype='dashed')
  # Only plot the upper and lower slopes from the bootstrap of the 'zero' PCA 
  p <- p + geom_abline(slope=zero.ci$basic[1,4], intercept=0, colour='black')
  p <- p + geom_abline(slope=zero.ci$basic[1,5], intercept=0, colour='black')
  p <- p + scale_x_continuous(pset.to.alias(pset.a), limits=c(lm.min, lm.max), expand=c(0,0))
  p <- p + scale_y_continuous(pset.to.alias(pset.b), limits=c(lm.min, lm.max), expand=c(0,0))
  p <- p + coord_equal()
  pdf(file=lm.f, width=5, height=4)
  print.ggplot(p)
  dev.off()    
  print("Done!")

  # Calculate some rather dry correlations between the datasets. Always use spearman.
  omega.cor <- cor(x$omega.y, x$omega.x, method='spearman', use='complete.obs')
  lrt.cor <- cor(x$lrt_stat.y, x$lrt_stat.x, method='spearman', use='complete.obs')
  sub.omega.cor <- cor(lm.sub$omega.y, lm.sub$omega.x, method='spearman', use='complete.obs')
  # This looks at the correlation for sites within the plotted range (0 < \omgml < 1)
  sub.lrt.cor <- cor(lm.sub$lrt_stat.y, lm.sub$lrt_stat.x, method='spearman', use='complete.obs')

  # Look at how many sites had non-overlapping 95% CIs.
  # X upper below Y lower
  x.blw.y <- x$omega_upper.x < x$omega_lower.y
  y.blw.x <- x$omega_upper.y < x$omega_lower.x  
  blw.abv.ratio <- sum(x.blw.y) / (sum(y.blw.x) + 1)

  # Calculate, for each site, the fraction of A's CI that are above B.
  # In other words, calculate the length of x's by the length of y's
  #        yyyyyyyyyy
  #               xxx
  #        |--------|
  #   |----------|
  #   aaaaa
  #   bbbbbbbbbbbb
  # If none of A's CI is above B, then a zero is contributed to the score.

  w.x <- x$omega_upper.x - x$omega_lower.x
  any.abv <- x$omega_upper.x > x$omega_upper.y
  w.abv <- (x$omega_upper.x - x$omega_upper.y) / w.x
  x$f.above <- 0
  x[any.abv, 'f.above'] <- w.abv[any.abv]
  median.abv <- median(x$f.above)
  mean.abv <- mean(x$f.above)  

  any.below <- x$omega_lower.x < x$omega_lower.y
  w.below <- (x$omega_lower.y - x$omega_lower.x) / w.x
  x$f.below <- 0
  x[any.below, 'f.below'] <- w.below[any.below]
  median.below <- median(x$f.below)
  mean.below <- mean(x$f.below)

  # OK, let's do a t-test on the sites between 0 and 1.
  t.res <- t.test(lm.sub$omega.x, lm.sub$omega.y, paired=T)
  t.pval <- t.res$p.value
  t.est <- t.res$estimate

  tz.res <- t.test(x$lrt_z.x, x$lrt_z.y, paired=T)
  tz.pval <- tz.res$p.value
  tz.est <- tz.res$estimate

  con <- connect(db())
  out.df <- data.frame(
    label = paste(pset.a, pset.b, filter, test, sep=' '),
    pset_a = pset.a,
    pset_b = pset.b,
    filter = filter,

    n_merged = nrow(x),
    n_pca = nrow(lm.sub),

    omega_cor = omega.cor,
    lrt_cor = lrt.cor,
    sub_omega_cor = sub.omega.cor,
    sub_lrt_cor = sub.lrt.cor,

    pca_slope = nonzero.pca[1],
    pca_lo = nonzero.ci$basic[1, 4],
    pca_hi = nonzero.ci$basic[1, 5],
 
    pcaz_slope = zero.pca[1],
    pcaz_lo = zero.ci$basic[1, 4],
    pcaz_hi = zero.ci$basic[1, 5],

    a_below_b = sum(x.blw.y),
    b_below_a = sum(y.blw.x),
    below_above_ratio = blw.abv.ratio,

    med_f_above = median.abv,
    mean_f_above = mean.abv,
    med_f_below = median.below,
    mean_f_below = mean.below,

    t_test_pval = t.pval,
    t_test_est = t.est, 
    z_t_test_pval = tz.pval,
    z_t_test_est = tz.est
  )

  ## Store this in the 'correlations' table.
  if (!test) {
    write.or.update(out.df, 'correlations', con, 'label')
  } else {
    print(out.df)
  }

}

parallel_sites <- function(filter, pset.a, pset.b, test=F) {
  pset.a <- as.integer(pset.a)
  pset.b <- as.integer(pset.b)
  test <- as.logical(test)
  x <- get.merged.sites(pset.a, pset.b,
    xtra=c('omega', 'lrt_stat', 'lrt_z', 'lrt_q', 'lrt_qq'), 
    filter=filter, 
    test=test
  ) 

  con <- connect(db())

  # Here, we're looking at the excess overlap between positively and
  #  negatively selected sites at a variety of significance
  #  thresholds, and also comparing the use of LRT p-values and
  #  percentile-based q-values for detecting such overlap. The idea is
  #  that q-values probably work better (esp. for identifying overlap
  #  in purifying sites) because of the different pop. sizes and
  #  efficacies of purifying selection. For positive selection... who knows?
  threshs <- c(0.01, 0.25, 0.5, 0.75, 0.8, 0.9, 0.95, 0.99)
  for (i in 1:length(threshs)) {
    thresh <- threshs[i]
    for (type in c('p', 'z', 'q', 'qq')) {
      out.df <- data.frame(
        label = paste(pset.a, pset.b, type, thresh, test, sep=' '),
        pset_a = pset.a,
        pset_b = pset.b,
        filter = filter,
        test = test,
        thresh = thresh,
        thresh_type = type,
        n_sites = nrow(x),
        pos_pval = NA,
        pos_est = NA,
        pos_lo = NA,
        pos_hi = NA,
        pos_n_ab = NA,
        pos_n_a = NA,
        pos_n_b = NA,
        neg_pval = NA,
        neg_est = NA,
        neg_lo = NA,
        neg_hi = NA,
        neg_n_ab = NA,
        neg_n_a = NA,
        neg_n_b = NA
      )

      if (type == 'p') {
        p_thresh <- threshs[i]
        neg.a <- factor(x$lrt_stat.x < -qchisq(p_thresh, df=1))
        neg.b <- factor(x$lrt_stat.y < -qchisq(p_thresh, df=1))
        pos.a <- factor(x$lrt_stat.x > qchisq(p_thresh, df=1))
        pos.b <- factor(x$lrt_stat.y > qchisq(p_thresh, df=1))
      } else if (type == 'z') {
        z_thresh <- threshs[i]
        neg.a <- factor(x$lrt_z.x < -qnorm(z_thresh))
        neg.b <- factor(x$lrt_z.y < -qnorm(z_thresh))
        pos.a <- factor(x$lrt_z.x > qnorm(z_thresh))
        pos.b <- factor(x$lrt_z.y > qnorm(z_thresh))
      } else if (type == 'qq') {
        q_thresh <- threshs[i] / 2
        neg.a <- factor(x$lrt_qq.x < 0.5 - q_thresh)
        neg.b <- factor(x$lrt_qq.y < 0.5 - q_thresh)
        pos.a <- factor(x$lrt_qq.x >= 0.5 + q_thresh)
        pos.b <- factor(x$lrt_qq.y >= 0.5 + q_thresh)
      } else if (type == 'q') {
        q_thresh <- threshs[i]
        neg.a <- factor(x$lrt_q.x < q_thresh)
        neg.b <- factor(x$lrt_q.y < q_thresh)
        pos.a <- factor(x$lrt_q.x > q_thresh)
        pos.b <- factor(x$lrt_q.y > q_thresh)
      }

      neg.n.a <- sum(neg.a == TRUE)
      neg.n.b <- sum(neg.b == TRUE)
      if (neg.n.a > 0 && neg.n.b > 0) {
        neg.res <- fisher.test(neg.a, neg.b, alternative='two.sided')
        neg.pval <- neg.res$p.value
        neg.est <- neg.res$estimate
        neg.ci <- neg.res$conf.int
        neg.ci[2] <- pmin(neg.ci[2], 200)
        neg.n.ab <- sum(neg.a == TRUE & neg.b == TRUE)
        out.df$neg_pval <- neg.pval
        out.df$neg_est <- neg.est
        out.df$neg_lo <- neg.ci[1]
        out.df$neg_hi <- neg.ci[2]
        out.df$neg_n_ab <- neg.n.ab
        out.df$neg_n_a <- neg.n.a
        out.df$neg_n_b <- neg.n.b
      }

      
      pos.n.a <- sum(pos.a == TRUE)
      pos.n.b <- sum(pos.b == TRUE)
      if (pos.n.a > 0 && pos.n.b > 0) {
        pos.res <- fisher.test(pos.a, pos.b, alternative='two.sided')
        pos.pval <- pos.res$p.value
        pos.est <- pos.res$estimate
        pos.ci <- pos.res$conf.int
        pos.ci[2] <- pmin(pos.ci[2], 200)
        pos.n.ab <- sum(pos.a == TRUE & pos.b == TRUE)
        out.df$pos_pval <- pos.pval
        out.df$pos_est <- pos.est
        out.df$pos_lo <- pos.ci[1]
        out.df$pos_hi <- pos.ci[2]
        out.df$pos_n_ab <- pos.n.ab
        out.df$pos_n_a <- pos.n.a
        out.df$pos_n_b <- pos.n.b
      }

      if (!test) {
        write.or.update(out.df, 'parallel_sites', con, 'label')
      } else {
        write.or.update(out.df, 'parallel_sites', con, 'label')
        print(out.df)
      }
    }
  }
  disconnect(con)

}

add_fill_density <- function(x, p, legend=F) {
  p <- p + theme_bw()
  p <- p + geom_point(stat='sum', aes(colour=..prop..), size=1.3)
  p <- p + scale_colour_gradientn("Density", colour=c('white', rgb(0.5, 0.5, 1), 'red'), trans='log')
  p <- p + opts(
    legend.position='none'
  )
  return(p)
}

plot_global_distribution <- function(pset, filter='default', test=F) {
  pset <- as.integer(pset)
  test <- as.logical(test)
  sites <- get.pset.sites(pset, filter=filter, test=test)
  p.list <- plot.global.distribution(sites, test=test)

  test.f <- ifelse(test, '_test', '')  
  pdf.f <- scratch.f(sprintf("global_dist_%s_%s%s.pdf", pset, filter, test.f))
  
  pdf(file=pdf.f, width=7, height=2)
  vplayout(6, 1)
  print(p.list$zero, vp=subplot(1:2,1))
  print(p.list$hist, vp=subplot(3:6,1))
  dev.off()
}

plot_sites_scatters <- function(pset, test=F) {
  sites <- get.pset.sites(6, filter='orig', test=T)
  sites <- process.sites(sites)
  sites <- filter.default(sites)
  sites$ci <- sites$omega_upper - sites$omega_lower
  x <- subset(sites, select=c('lrt_stat', 'omega', 'ci'))

  r.f <- function(x, a, b, c) {
    x$lrt_stat <- round_any(x$lrt_stat, a*2.2)
    x$omega <- round_any(x$omega, b*2.2)
    x$ci <- round_any(x$ci, c*2.2)
    x
  }

  print("  p")
  x <- subset(sites, lrt_stat > -50 & lrt_stat < 5 & omega < 2)
  x <- r.f(x, 0.4, 0.03, 0.06)
  p <- ggplot(x, aes(x=lrt_stat, y=omega))
  p <- add_fill_density(x, p)
  p <- p + xlab("Signed LRT") + ylab("Omega")

  print("  q")
  x <- subset(sites, lrt_stat > -50 & lrt_stat < 5 & ci >= 0 & ci < 5)
  x <- r.f(x, 0.4, 0.03, 0.06)
  q <- ggplot(x, aes(x=lrt_stat, y=ci))
  q <- add_fill_density(x, q)
  q <- q + xlab("Signed LRT") + ylab("95% CI Width")

  print("  z")
  x <- subset(sites, omega < 1.5 & ci < 5)
  x <- r.f(x, 0.15, 0.015, 0.03)
  z <- ggplot(x, aes(x=omega, y=ci))
  z <- add_fill_density(x, z)
  z <- z + xlab("Omega") + ylab("95% CI Width")

  sites <- get.pset.sites(1, filter='orig', test=T)
  sites <- process.sites(sites)
  sites <- filter.default(sites)
  sites$ci <- sites$omega_upper - sites$omega_lower
  sites <- subset(sites, select=c('lrt_stat', 'omega', 'ci'))

  print("  p")
  x <- subset(sites, lrt_stat > -20 & lrt_stat < 5 & omega < 2)
  x <- r.f(x, 0.15, 0.02, 0.06)
  p2 <- ggplot(x, aes(x=lrt_stat, y=omega))
  p2 <- add_fill_density(x, p2)
  p2 <- p2 + xlab("Signed LRT (w > 0)") + ylab("Omega")

  print("  q")
  x <- subset(sites, lrt_stat > -20 & lrt_stat < 5 & ci >= 0 & ci < 5)
  x <- r.f(x, 0.15, 0.02, 0.06)
  q2 <- ggplot(x, aes(x=lrt_stat, y=ci))
  q2 <- add_fill_density(x, q2)
  q2 <- q2 + xlab("Signed LRT") + ylab("95% CI Width")

  print("  z")
  x <- subset(sites, omega < 1.5 & ci < 5)
  x <- r.f(x, 0.15, 0.015, 0.03)
  z2 <- ggplot(x, aes(x=omega, y=ci))
  z2 <- add_fill_density(x, z2)
  z2 <- z2 + xlab("Omega (0 > w > 1)") + ylab("95% CI Width")

  n <- 3
  ff <- scratch.f(paste("sites_scatters.pdf", sep=''))
  svg(file=scratch.f("sites_scatters_orig.svg"), width=n*3, height=3*2)
  #pdf(file=ff, width=n*3, height=3 * 2)
  vplayout(n, 2)
  print.ggplot(p, vp=subplot(1, 1))
  print.ggplot(q, vp=subplot(2, 1))
  print.ggplot(z, vp=subplot(3, 1))
  print.ggplot(p2, vp=subplot(1, 2))
  print.ggplot(q2, vp=subplot(2, 2))
  print.ggplot(z2, vp=subplot(3, 2))
  dev.off()
}

plot_qc_histograms <- function(pset, filter='default', test=F) {
  pset <- as.integer(pset)
  test <- as.logical(test)

  sites <- get.pset.sites(pset, filter=filter, test=test)
  
  qc.fields <- list(
    'omega' = c(-0.2, 3),
    'lrt_stat' = c(-50, 10),
    'ncod' = c(0, 40),
    'nongap_bl' = c(0, 15),
    'note' = c('Constant', 'Synonymous', 'Nonsynonymous'),
    'random' = c('Nonrandom', 'Random')
  )

  p.list <- list()
  for (i in 1:length(qc.fields)) {
    fld <- names(qc.fields)[i]
    rng <- qc.fields[[i]]
    print(fld)
    print(rng)
    sites$cur.val <- sites[, fld]
    cur.df <- subset(sites, select='cur.val')

    if (fld == 'random') {
      cur.df$cur.val <- ifelse(is.na(cur.df$cur.val), 'Nonrandom', 'Random')
      cur.df$cur.val <- factor(cur.df$cur.val, levels=rng, labels=rng)
    }
    if (fld == 'note') {
      cur.df$cur.val <- ifelse(is.na(cur.df$cur.val), 'nonsynonymous', as.character(cur.df$cur.val))
      cur.df$cur.val <- factor(cur.df$cur.val, levels=c('constant', 'synonymous', 'nonsynonymous'), labels=rng)
    }

    p <- ggplot(cur.df, aes(x=cur.val))
    p <- p + theme_bw()
    if (any(is.numeric(rng))) {
      if (pset %in% c(1, 2, 3, 4)) {
        if (fld == 'nongap_bl') {
          rng <- c(0, 5)
        } else if (fld == 'lrt_stat') {
          rng <- c(-20, 5)
        }
      }
      if (test) {
        binw <- diff(rng) / 30
      } else {
        binw <- diff(rng) / 100
      }
      if (fld == 'ncod') {
        binw <- 1
        if (pset %in% c(1, 2, 3, 4, 7, 8)) {
          rng <- c(0, 10)
        }        
      }

      cur.df$cur.val <- pmin(cur.df$cur.val, rng[2])
      cur.df$cur.val <- pmax(cur.df$cur.val, rng[1])
      rng[2] <- rng[2] * 1.05
      rng[1] <- rng[1] * 0.95
      p <- p + geom_histogram(binwidth=binw)
      p <- p + scale_x_continuous(fld, limits=rng)
    } else {
      p <- p + geom_histogram()
      p <- p + scale_x_discrete(fld)
      p <- p + opts(
        axis.text.x = theme_blank()
      )
    }

    p <- p + opts(
      axis.title.y = theme_blank(),
      axis.title.x = theme_blank()
    )
    p.list[[i]] <- p
  }

  edge.len <- 3
  n <- length(p.list)
  test.s <- ifelse(test, '_test', '')
  file.s <- scratch.f(sprintf("qc_hist_%s_%s%s.pdf", pset, filter, test.s))
  pdf(file=file.s,
    height=edge.len,
    width=edge.len * n
  )
  vplayout(n, 1)
  for (i in 1:n) {
    print(paste("  printing ", i, sep=''))
    cur.p <- p.list[[i]]
    print.ggplot(cur.p, vp=subplot(i, 1))
  }
  dev.off()
}

main()