source("~/src/greg-ensembl/projects/2xmammals/analyze_mammals.R")
args <- commandArgs(trailingOnly=T)

main <- function() {
  if (!is.na(args[1])) {
    fn.name <- args[1]
    print(paste("Calling", fn.name))
    do.call(fn.name, as.list(args[2:length(args)]))
  }
}

collect_sites <- function(pset, test=F) {
  pset <- as.integer(pset)
  test <- as.logical(test)

  limit.s <- ''
  if (test) {
    print("  limiting sites count for testing")
    limit.s <- 'limit 500000'
  }

  orig.sites.f <- paste(scratch(), "sites_orig_", pset, ".Rdata", sep='')
  if (!test && file.exists(orig.sites.f)) {
    load(orig.sites.f)
  } else {
    con <- connect(db()); 
    cmd <- sprintf("select * from sites where parameter_set_id=%d %s", pset, limit.s)
    sites <- dbGetQuery(con, cmd)
  }

  if (!test && !file.exists(orig.sites.f)) {
    save(sites, file=orig.sites.f)
  }

  sites.f <- paste(scratch(), "sites_", pset, ".Rdata", sep='')
  if (test) {
    sites.f <- paste(scratch(), "sites_test_", pset, ".Rdata", sep='')
  }
  if (!file.exists(sites.f)) {
    sites <- process.sites(sites)
    save(sites, file=sites.f)
  } else {
    print("  sites file already exists -- not processing again!")
  }
  print("  done!")
}

summary_table <- function(pset, test=F, filter='') {
  sites <- get.pset.sites(pset, test=test)

  df <- summarize.sites(sites, filter=filter)

  con <- connect(db())
  write.or.update(df, 'summaries', con, 'name')
  dbDisconnect(con)
}

pset_correlation <- function(pset.a, pset.b, test=F) {
  mgd <- get.merged.sites(pset.a, pset.b, xtra=c('omega', 'lrt_stat', 'lrt_zz', 'omega_lower', 'omega_upper'), test=test)
  
  lm.max <- 1
  lm.min <- 0
  lm.sub <- subset(mgd, omega.x > lm.min & omega.y > lm.min & omega.x < lm.max & omega.y < lm.max)

  iwpca.slope <- function(data, indices) {
    xs <- data[indices, 'omega.x']
    ys <- data[indices, 'omega.y']
    r <- wpca(cbind(xs,ys), center=center, w=ws)
    r.slope <- r$vt[1,2] / r$vt[1,1]
  }
  iwpca.slope.int <- function(xs, ys, ws=NULL, center=T) {
    r <- wpca(cbind(xs,ys), center=center, w=ws)
    r.slope <- r$vt[1,2] / r$vt[1,1]
    r.int <- r$xMean[2] - r.slope * r$xMean[1]
    return(c(r.slope, r.int))
  }



  library(nlme)

  w.x <- lm.sub$omega_upper.x - lm.sub$omega_lower.x
  w.y <- lm.sub$omega_upper.y - lm.sub$omega_lower.y
  lm.sub$wt <- (w.x + w.y) / 2
  lm.sub$wt <- pmax(0.05, lm.sub$wt)
  lm.sub$wt <- pmin(10, lm.sub$wt)

  lm.res <- lm(omega.y ~ omega.x, data=lm.sub, weights=I(1/wt))
  coefs <- coef(lm.res)
  lm.coef <- coefs[2]
  lm.intercept <- coefs[1]
  lm.sum <- summary(lm.res)
  lm.rsquared <- lm.sum$r.squared
  lm.conf <- confint(lm.res)

  lm0.res <- lm(omega.y ~ omega.x - 1, data=lm.sub, weights=I(wt))
  coefs <- coef(lm0.res)
  lm0.coef <- coefs[1]
  lm0.sum <- summary(lm0.res)
  lm0.rsquared <- lm0.sum$r.squared
  lm0.conf <- confint(lm0.res)

  omega.cor.p <- cor(mgd$omega.y, mgd$omega.x, method='pearson', use='complete.obs')
  omega.cor.s <- cor(mgd$omega.y, mgd$omega.x, method='spearman', use='complete.obs')
  lrt.cor.p <- cor(mgd$lrt_stat.y, mgd$lrt_stat.x, method='pearson', use='complete.obs')
  lrt.cor.s <- cor(mgd$lrt_stat.y, mgd$lrt_stat.x, method='spearman', use='complete.obs')

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
  p <- p + geom_abline(slope=1, intercept=0, colour='black', linetype='solid')

  xs <- seq(from=lm.min, to=lm.max, length.out=100)
  lines.df <- data.frame(omega.x = xs)
  lines.df$omega.y <- stats::predict(lm.res, newdata=lines.df)
  p <- p + geom_line(data=lines.df, colour='orange', linetype='dashed')
  lines.df <- data.frame(omega.x = xs)
  lines.df$omega.y <- stats::predict(lm0.res, newdata=lines.df)
  p <- p + geom_line(data=lines.df, colour='black', linetype='dashed')

  p <- p + scale_x_continuous(pset.to.alias(pset.a), limits=c(lm.min, lm.max), expand=c(0,0))
  p <- p + scale_y_continuous(pset.to.alias(pset.b), limits=c(lm.min, lm.max), expand=c(0,0))
  p <- p + coord_equal()
  pdf(file=lm.f, width=5, height=5)
  print.ggplot(p)
  dev.off()  
  
  print("Done!")

  con <- connect(db())

  z_threshs <- c(0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.49, 0.495)
  for (i in 1:length(z_threshs)) {
    z_thresh <- z_threshs[i]
    neg.a <- factor(mgd$lrt_zz.x < 0.5 - z_thresh)
    neg.b <- factor(mgd$lrt_zz.y < 0.5 - z_thresh)
    pos.a <- factor(mgd$lrt_zz.x >= 0.5 + z_thresh)
    pos.b <- factor(mgd$lrt_zz.y >= 0.5 + z_thresh)

    neg.fis.res <- fisher.test(neg.a, neg.b, alternative='two.sided')
    neg.fis.pval <- neg.fis.res$p.value
    neg.fis.est <- neg.fis.res$estimate
    neg.ci <- neg.fis.res$conf.int
    neg.ci[2] <- pmin(neg.ci[2], 99)
    neg.fis.n <- sum(neg.a == TRUE & neg.b == TRUE)

    pos.fis.res <- fisher.test(pos.a, pos.b, alternative='two.sided')
    pos.fis.pval <- pos.fis.res$p.value
    pos.fis.est <- pos.fis.res$estimate
    pos.ci <- pos.fis.res$conf.int
    pos.ci[2] <- pmin(pos.ci[2], 99)
    pos.fis.n <- sum(pos.a == TRUE & pos.b == TRUE)

    out.df <- data.frame(
      psets = paste(pset.a, pset.b, z_thresh, sep=' '),
      pset_a = pset.a,
      pset_b = pset.b,

      omega_cor_p = omega.cor.p,
      omega_cor_s = omega.cor.s,
      lrt_cor_p = lrt.cor.p,
      lrt_cor_s = lrt.cor.s,

      lm_coef = lm.coef,
      lm_rsquared = lm.rsquared,
      lm_lo = lm.conf[1],
      lm_hi = lm.conf[2],

      lm0_coef = lm0.coef,
      lm0_rsquared = lm0.rsquared,
      lm0_lo = lm0.conf[1],
      lm0_hi = lm0.conf[2],

      z_thresh = z_thresh,    
      pos_fis_pval = pos.fis.pval,
      pos_fis_est = pos.fis.est,
      pos_fis_lo = pos.ci[1],
      pos_fis_hi = pos.ci[2],
      pos_fis_n = pos.fis.n,

      neg_fis_pval = neg.fis.pval,
      neg_fis_est = neg.fis.est,
      neg_fis_lo = neg.ci[1],
      neg_fis_hi = neg.ci[2],
      neg_fis_n = neg.fis.n
    )
    write.or.update(out.df, 'correlations', con, 'psets')
  }
  dbDisconnect(con)

}

plot_global_distribution <- function(pset, test=F) {
  sites <- get.pset.sites(pset, test=test)
  
  sites <- filter.default(sites)
  p.list <- plot.sites.histograms(sites, test=test)
  
  pdf.f <- scratch.f(paste('global_dist_', pset, '.pdf', sep=''))
  
  pdf(file=pdf.f, width=6, height=3)
  vplayout(5, 1)
  print(p.list$zero, vp=subplot(1,1))
  print(p.list$hist, vp=subplot(2:5,1))
  dev.off()
}

plot_qc_histograms <- function(pset, test=F, post.filter=F) {
  sites <- get.pset.sites(pset, test=test)
  
  post.filter <- as.logical(post.filter)
  if (post.filter) {
    sites <- filter.default(sites)
  }

  qc.fields <- list(
  'omega' = c(-0.2, 3),
  'lrt_stat' = c(-50, 10),
  'ncod' = c(0, 40),
  'nongap_bl' = c(0, 15),
  'note' = NA,
  'random' = NA
  )

  p.list <- list()
  for (i in 1:length(qc.fields)) {
    fld <- names(qc.fields)[i]
    rng <- qc.fields[[i]]
    print(fld)
    print(rng)
    sites$cur.val <- sites[, fld]
    cur.df <- subset(sites, select='cur.val')
    p <- ggplot(cur.df, aes(x=cur.val))
    if (any(!is.na(rng))) {
      if (test) {
        binw <- diff(rng) / 30
      } else {
        binw <- diff(rng) / 100
      }
      if (fld == 'ncod') {
        binw <- 1
      }
      cur.df$cur.val <- pmin(cur.df$cur.val, rng[2])
      cur.df$cur.val <- pmax(cur.df$cur.val, rng[1])
      p <- p + geom_histogram(binwidth=binw)
      p <- p + scale_x_continuous(fld, limits=rng)
    } else {
      p <- p + geom_histogram()
      p <- p + scale_x_discrete(fld)
      p <- p + opts(
        axis.text.x = theme_text(angle=90, hjust=1)
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
  file.s <- scratch.f(paste('qc_hist_', pset, '_', as.character(post.filter), '.pdf', sep=''))
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