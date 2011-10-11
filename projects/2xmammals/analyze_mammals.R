library(R.oo)
library(ape)
library(plyr)
source("~/src/greg-ensembl/scripts/mysql_functions.R")
source("~/src/greg-ensembl/projects/2xmammals/analyze_filters.R")
source("~/src/greg-ensembl/projects/orthologs/analyze_orthologs.R")
source("~/src/greg-ensembl/scripts/xtable_utils.R")
library(ggplot2)

db <- function() {
   'gj1_2x_63_alt'
}

scratch.f <- function(f) {
  paste(scratch(), f, sep='')
}

project.f <- function(f) {
  paste("~/src/greg-ensembl/projects/2xmammals/", f, sep='')
}

scratch <- function() {
  "~/scratch/gj1_2x_63_alt/current/data/"
}

get.genes <- function() {
  if (!exists('genes.df', envir=.GlobalEnv)) {
    con <- connect(db()); 
    df <- dbGetQuery(con, 'select * from genes')
    dbDisconnect(con)
    assign('genes.df', df, envir=.GlobalEnv)
  }
  df <- get('genes.df', envir=.GlobalEnv)
  df
}

get.genes.split <- function() {
  if (!exists('genes.split.df', envir=.GlobalEnv)) {
    genes <- get.genes()

    pset.df <- pset.df()
    clades <- pset.df$char

    out.df <- data.frame()
    
    for (i in 1:length(clades)) {
      clade <- clades[i]
      bl.s <- paste(clade, '_slr_total_length', sep='')
      sites.s <- paste(clade, '_slr_sites', sep='')
      dnds.s <- paste(clade, '_slr_dnds', sep='')
      mpl.s <- paste(clade, '_slr_mean_path', sep='')
      leaves.s <- paste(clade, '_leaf_count', sep='')
      ecdf.f <- ecdf(genes[, dnds.s])
      cum <- ecdf.f(genes[, dnds.s])

      cur.genes <- genes
      cur.genes$pset <- i
      cur.genes$clade <- clade
      cur.genes$n_leaves <- genes[, leaves.s]
      cur.genes$dnds <- genes[, dnds.s]
      cur.genes$n_sites <- genes[, sites.s]
      cur.genes$bl <- genes[, bl.s]
      cur.genes$mpl <- genes[, mpl.s]
      cur.genes$cum_dnds <- cum
      out.df <- rbind(out.df, cur.genes)
    }
    assign('genes.split.df', out.df, envir=.GlobalEnv)
  }
  df <- get('genes.split.df', envir=.GlobalEnv)
  df
}

get.pset.sites <- function(pset, filter='default', test=F) {
  sites.f <- get.sites.filename(pset, filter, test)
  if (file.exists(sites.f)) {
    load(sites.f)
  } else {
    stop("Sites file does not exist!")
  }
  sites
}

get.sites.filename <- function(pset, filter='orig', test=F) {
  test.str <- ''
  if (test) {
    test.str <- '_test'
  }
  return(sprintf("%ssites_%d_%s%s.Rdata", scratch(), pset, filter, test.str))
}

bsub.pset.function <- function(fn_name, psets=NA, queue='normal', mem=4, extra.args='', ...) {
  if (any(is.na(psets))) {
    df <- pset.df()
    psets <- df$pset_id
  }

  for (i in 1:length(psets)) {
    pset <- psets[i]
    cur.xtra <- paste(pset, extra.args, sep=' ')
    bsub.function(fn_name, queue, mem, extra.args=cur.xtra, ...)
  }
}

bsub.function <- function(fn_name, queue='normal', mem=4, extra.args='', jobarray=NULL, jobarray_id=fn_name) {
  array_s = ''
  if (!is.null(jobarray)) {
    array_s <- paste('-J ', jobarray_id, '[1-', jobarray, ']', sep='')
  }
  args_s <- paste(fn_name, ' ', extra.args, sep='')
  cmd <- sprintf('bsub -q %s -R "select[mem>%d000] rusage[mem=%d000]" -M%d000000 %s "/software/R-2.13.0/bin/R --vanilla --args %s < sites_scripts.R"',
    queue, mem, mem, mem,
    array_s, args_s
  )
  print(cmd)
  system(cmd)
}

bsub.collect.sites <- function(...) {
  bsub.pset.function('collect_sites', queue='normal', mem=18, ...)
}

bsub.collect.sites.filtered <- function(...) {
  for (filter in c('none', 'default', 'stringent', 'pfam', 'pfam_stringent', 'clusters', 'clusters_inverse')) {
    bsub.pset.function('collect_sites', extra.args=filter, queue='normal', mem=18, ...)
  }
}

bsub.collect.genes <- function(...) {
  bsub.pset.function('collect_genes', queue='normal', mem=18, ...)
}

bsub.summaries <- function(...) {
  filters <- c('none', 'default', 'stringent', 'pfam', 'clusters_inverse')
  for (i in 1:length(filters)) {
    bsub.pset.function('summary_table', queue='normal', mem=12, extra.args=paste('F ', filters[i], sep=''), ...)
  }
}

bsub.plots <- function() {
  plots <- c('plot_global_distribution')
  for (i in 1:length(plots)) {
    bsub.pset.function(plots[i], queue='normal', mem=18)
  }

  plots <- c('plot_qc_histograms')
  for (i in 1:length(plots)) {
    bsub.function(plots[i], queue='normal', mem=18, extra.args='1 default F')
    bsub.function(plots[i], queue='normal', mem=18, extra.args='6 default F')
    bsub.function(plots[i], queue='normal', mem=18, extra.args='1 none F')
    bsub.function(plots[i], queue='normal', mem=18, extra.args='6 none F')
    bsub.function(plots[i], queue='normal', mem=18, extra.args='1 stringent F')
    bsub.function(plots[i], queue='normal', mem=18, extra.args='6 stringent F')
  }

}

bsub.fits <- function() {
  df <- pset.df()
  psets <- df$pset_id

  for (i in 1:length(psets)) {
    for (j in c('lnorm', 'gamma', 'beta', 'exp', 'weibull')) {
      for (k in c('ci', 'omega')) {
        jobarray_id <- paste(i, j, k, sep='_')
        bsub.function('fit_distr', mem=8, extra.args=paste(i, j, k, 'F', sep=' '), jobarray=100, jobarray_id=jobarray_id)
      }
    }
  }
}

write.fits.table <- function() {
  # One row per pset, showing mean AIC / mean / F > 1 for each model fit
    con <- connect(db())
    df <- dbGetQuery(con, 'select * from fitdistr')
    dbDisconnect(con)

    df.sub <- subset(df)

    xt.df <- ddply(df.sub, .(pset, use_type), function(x) {
      cur.df <- data.frame(
        pset = x[1, 'pset'],
        use_type = x[1, 'use_type']
      )
      for (distr_type in c('lnorm', 'gamma', 'exp', 'beta', 'weibull')) {
        distr.row <- subset(x, dist == distr_type)
        distr.df <- data.frame(
          mean = median(distr.row$mean),
          f_above_1 = median(distr.row$f_above_1) * 100
        )
        colnames(distr.df) <- paste(distr_type, colnames(distr.df), sep='_')
        cur.df <- cbind(cur.df, distr.df)
      }
      cur.df
    })

    xt.df$pset <- pset.to.alias(xt.df$pset)
    xt.df <- xt.df[order(xt.df$use_type),]
    xt.df$use_type <- factor(xt.df$use_type, levels=c('omega', 'ci', ''), labels=c('\\omgml', '\\ci', ''))

    df1 <- subset(xt.df, use_type == '\\omgml')
    df2 <- subset(xt.df, use_type == '\\ci')
    xt1 <- xtable(df1)
    xt2 <- xtable(df2)
    clr.cols <- colnames(xt1)
    clr.cols <- setdiff(clr.cols, c('pset', 'use_type'))
    xt1 <- color.columns(xt1, clr.cols, log=T)
    xt2 <- color.columns(xt2, clr.cols, log=T)

    xt1[-1, 'use_type'] <- ''
    xt2[-1, 'use_type'] <- ''
    ff <- scratch.f("distribution_fits.txt")
    print.latex(xt1, ff)
    cat("\\midrule", file=ff, append=T)
    print.latex(xt2, ff, append=T)


    # Now sort by the best-fit model type for each dataset (only using CIs) and write a table
    # with the AIC and best fit parameters
    df.sub <- subset(df, use_type == 'ci')
    
    df.sub$dist <- dist.to.factor(df.sub$dist)

    mean.df <- ddply(df.sub, .(pset, dist), function(x) {
      data.frame(
        aic = median(x$aic),
        param_1 = median(x$est_1),
        param_2 = median(x$est_2)
      )
    })

    mean.df <- ddply(mean.df, .(pset), function(x) {
      min.aic <- min(x$aic)
      max.aic <- max(x$aic)
      x <- x[order(x$aic), ]
      x$aic_diff <- c(0,diff(x$aic))
      x
    })

    mean.df <- mean.df[with(mean.df, order(pset, aic)),]
    mean.df[ (1:nrow(mean.df)) %% 5 != 1, 'pset'] <- ''
    mean.df$pset <- pset.to.alias(mean.df$pset)
    mean.df$pset <- ifelse(is.na(mean.df$pset), '', mean.df$pset)
    mean.df[mean.df$dist == 'Exponential', 'param_2'] <- NA
    mean.df <- subset(mean.df, select=c('pset', 'dist', 'aic', 'aic_diff', 'param_1', 'param_2'))

    print(mean.df)
    xt <- xtable(mean.df)
    xt <- color.columns(xt, c('aic_diff'))

    ff <- scratch.f("distribution_params.txt")
    print.latex(xt, file=ff)


    # Finally, plot a figure for some species groups showing the bootstrap distributions of mean omega and % > 1.
    con <- connect(db())
    df <- dbGetQuery(con, 'select * from fitdistr')
    dbDisconnect(con)

    df$dist <- dist.to.factor(df$dist)

    df.sub <- subset(df, use_type == 'ci')
    df.sub$pset <- pset.to.alias(df.sub$pset, factors=T)
    df.sub$pset <- factor(df.sub$pset, levels=df.sub$pset, labels=paste(" ", df.sub$pset, sep=''))

    p <- ggplot(df.sub, aes(y=-as.numeric(dist), x=mean, colour=dist, fill=dist))
    p <- p + theme_bw()
    p <- p + geom_point(aes(y=jitter(-as.numeric(dist), factor=2)), size=0.3)
    p <- p + scale_colour_brewer("Distribution")
    p <- p + scale_fill_brewer("Distribution")
    p <- p + scale_x_continuous("Mean Omega", limits=c(0.15, 0.45))
    p <- p + scale_y_continuous("Distribution", breaks=c(0))
    p <- p + facet_grid(pset ~ .)
    p <- p + opts(
      axis.ticks.y = theme_blank(),
      axis.text.y = theme_blank(),
      axis.text.x = theme_text(angle=90, hjust=1, size=8),
      strip.text.y = theme_text(angle=0, hjust=0, size=10),
      strip.background = theme_blank()
    )
    
    ff <- scratch.f("distribution_mean.svg")
    svg(file=ff, width=8, height=4)
    print.ggplot(p)
    dev.off()

    df.sub <- subset(df, use_type == 'ci' & dist != 'beta')
    df.sub$pset <- pset.to.alias(df.sub$pset, factors=T)
    df.sub$pset <- factor(df.sub$pset, levels=df.sub$pset, labels=paste(" ", df.sub$pset, sep=''))

    p <- ggplot(df.sub, aes(y=-as.numeric(dist), x=f_above_1 * 100, colour=dist, fill=dist))
    p <- p + theme_bw()
    p <- p + geom_point(aes(y=jitter(-as.numeric(dist), factor=2)), size=0.3)
    p <- p + scale_colour_brewer("Distribution")
    p <- p + scale_fill_brewer("Distribution")
    p <- p + scale_x_continuous("% Omega > 1", limits=c(0 * 100, 0.07 * 100))
    p <- p + scale_y_continuous("Distribution", breaks=c(0))
    p <- p + facet_grid(pset ~ .)
    p <- p + opts(
      axis.ticks.y = theme_blank(),
      axis.text.y = theme_blank(),
      axis.text.x = theme_text(angle=90, hjust=1, size=8),
      strip.text.y = theme_text(angle=0, hjust=0, size=10),
      strip.background = theme_blank()
    )

    ff <- scratch.f("distribution_above_one.svg")
    svg(file=ff, width=8, height=4)
    print.ggplot(p)
    dev.off()
}

bsub.corrs <- function() {
  df <- pset.df()
  psets <- df$pset_id

  for (i in 1:length(psets)) {
    for (j in 1:length(psets)) {
      if (i == j) {
        next()  
      }
      bsub.function('pset_correlation', extra.args=paste(i, j, 'F', sep=' '), mem=11)
    }
  }
}

bsub.parallel <- function() {
  df <- pset.df()
  psets <- df$pset_id

  for (i in 1:length(psets)) {
    for (j in 1:length(psets)) {
      if (i == j) {
        next()  
      }
      bsub.function('parallel_sites', extra.args=paste(i, j, 'F', sep=' '), mem=11)
    }
  }
}

process.genes <- function(genes, pset, filter='default', subset.index=NULL, test=F) {
  gc(T)

  # Get the sites corresponding to this pset.
  sites <- get.pset.sites(pset, test=test)

  abv.one <- sites$lrt_stat > 0
  sites[, 'pos.pval'] <- 1
  if (any(abv.one)) {
    sites[abv.one, 'pos.pval'] <- sites[abv.one, 'pval'] / 2
  }
  if (any(!abv.one)) {
    sites[!abv.one, 'pos.pval'] <- 0.5 + (1 - sites[!abv.one, 'pval']) / 2
  }
  sites[, 'pos.bh'] <- p.adjust(sites[, 'pos.pval'], method='BH')
  print(summary(sites$pos.pval))

  gc(T)

  # Filter sites.
  if (filter == 'stringent') {
    sites <- filter.stringent(sites)
  } else {
    sites <- filter.default(sites)
  }

  cum.z <- function(x) {
    ecdf.f <- ecdf(x)
    ecdf.f(x)
  }

  # Extract the parameter-set specific columns from the genes df
  pset.char <- pset.to.short(pset)
  bl.s <- paste(pset.char, '_slr_total_length', sep='')
  sites.s <- paste(pset.char, '_slr_sites', sep='')
  dnds.s <- paste(pset.char, '_slr_dnds', sep='')
  mpl.s <- paste(pset.char, '_slr_mean_path', sep='')
  leaves.s <- paste(pset.char, '_leaf_count', sep='')
  # Calculate z-scores for dnds and branch length values.
  dnds.z <- cum.z(genes[, dnds.s])
  bl.z <- cum.z(genes[, bl.s])
  gene.cols <- subset(genes, select=c(
    'data_id', 'aln_length', 'chr_end', 'chr_start', 'chr_name', 'chr_strand',
    'data_prefix', 'dup_species_count', 'dup_species_list', 'gc_3', 'gc_cds', 'gc_genomic',
    'gene_name', 'job_id', 'leaf_count', 'masked_ids', 'masked_nucs', 'paml_dnds', 'paml_mean_path',
    'paml_total_length', 'ref_gene_description', 'ref_gene_id', 'ref_protein_id', 'ref_taxon_id',
    'skipped_clade_count', 'skipped_clades', 'tree_mean_path', 'tree_total_length')
  )
  gene.cols$slr_bl <- genes[, bl.s]
  gene.cols$slr_mpl <- genes[, mpl.s]
  gene.cols$slr_sites <- genes[, sites.s]
  gene.cols$slr_dnds <- genes[, dnds.s]
  gene.cols$slr_dnds_z <- dnds.z
  gene.cols$slr_bl_z <- bl.z
  genes <- gene.cols

  genes$parameter_set_id <- pset
  genes$pset_char <- pset.char

  gc(T)

  # Take a subset of data_ids to process gene-wise
  if (!is.null(subset.index)) {
    n.indices <- 100

    all.ids <- sort(unique(sites$data_id))
    genes.per.index <- floor(length(all.ids) / n.indices)

    cur.lo <- (subset.index-1) * genes.per.index + 1
    cur.hi <- cur.lo + (genes.per.index - 1)
    if (subset.index == n.indices) {
      cur.hi <- length(all.ids)
    }
    print(paste(cur.lo, cur.hi))
    cur.ids <- all.ids[cur.lo:cur.hi]
    sites <- subset(sites, data_id %in% cur.ids)
    rm(all.ids)
    rm(cur.ids)
  }

  print(sprintf("%d genes to process", length(unique(sites$data_id))))

  # ddply each gene, extract the pset SLR values and analyze the sites.
  con <- connect(db())
  d_ply(sites, .(data_id), function(x) {
    cur.gene <- subset(genes, data_id == x[1, 'data_id'])
    if (nrow(cur.gene) == 0) {
      print(paste("Gene not found for data_id", x[1, 'data_id']))
      return()
    }

    sites.cols <- process.gene.sites(cur.gene, x)
    cur.df <- cbind(cur.gene, sites.cols)

    cur.df$filter <- filter
    cur.df$label <- paste(cur.df$filter, cur.df$parameter_set_id, cur.df$data_id, collapse=' ')
  
    write.or.update(cur.df, 'genes_sites', con, 'label')
  })

  dbDisconnect(con)
}

# Given a row corresponding to a gene, summarize the sitewise data w.r.t. pos and neg selection.
process.gene.sites <- function(gene, sites) {
  print(sprintf("%d sites for gene %s", nrow(sites), gene$gene_name))

  # omg_ml > 1 and < 1
  n.abv.1 <- nrow(subset(sites, omega > 1))
  n.blw.1 <- nrow(subset(sites, omega < 1))

  n.pos.20 <- nrow(subset(sites, omega > 1 & pval < 0.2))
  n.neg.20 <- nrow(subset(sites, omega < 1 & pval < 0.2))
  n.ntr.20 <- nrow(subset(sites, pval > 0.2))

  n.pos.10 <- nrow(subset(sites, omega > 1 & pval < 0.1))
  n.neg.10 <- nrow(subset(sites, omega < 1 & pval < 0.1))
  n.ntr.10 <- nrow(subset(sites, pval > 0.1))

  n.pos.05 <- nrow(subset(sites, omega > 1 & pval < 0.05))
  n.neg.05 <- nrow(subset(sites, omega < 1 & pval < 0.05))
  n.ntr.05 <- nrow(subset(sites, pval > 0.05))

  n.pos.01 <- nrow(subset(sites, omega > 1 & pval < 0.01))
  n.pos.bh <- nrow(subset(sites, omega > 1 & pos.bh < 0.05))

  sites$hoch.pval <- p.adjust(sites$pval, method='hochberg')
  n.pos.fwer <- nrow(subset(sites, hoch.pval < 0.05 & omega > 1))

  # Use the one-tailed p-values to do BH adjustment.
  sites$bh.pval <- p.adjust(sites$pos.pval, method='BH')
  n.pos.gene.bh.10 <- nrow(subset(sites, bh.pval < 0.10))
  n.pos.gene.bh.05 <- nrow(subset(sites, bh.pval < 0.05))

  print(sprintf("SLR: 05 %d hoch %d  ME: 05 %d hoch %d",
    nrow(subset(sites, type %in% c('positive1', 'positive2', 'positive3', 'positive4'))),
    nrow(subset(sites, type %in% c('positive3', 'positive4'))),
    n.pos.05,
    n.pos.fwer
  ))

  omgs <- pmin(sites$omega, 2)
  mean.omg <- mean(omgs)
  rm(omgs)

  sorted.omgs <- sort(sites$omega)
  n.quart <- nrow(sites) / 10
  mid.omgs <- sorted.omgs[floor(.5*n.quart):floor(9.5*n.quart)]
  mean.mid.omega <- mean(mid.omgs)
  rm(sorted.omgs)
  rm(mid.omgs)

  # Get some estimate of the lrt_stat autocorrelation using a Mantel test
  library(ade4)
  aln.dists <- dist(sites$aln_position)
  lrt.dists <- dist(sites$lrt_stat)
  #mantel.res <- mantel.test(aln.dists, lrt.dists, nperm=100)
  #mantel.res <- mantel.rtest(aln.dists, lrt.dists, nrepet = 999)
  #alt.mantel <- mantel(lrt.dists ~ aln.dists, nperm=1000)
  #print(str(alt.mantel))
  # Use the fast mantel test implemented in C.
  mantel.res <- mantel.randtest(aln.dists, lrt.dists, nrepet = 1000)
  mantel.cor <- mantel.res$obs
  mantel.p <- mantel.res$pvalue

  # Also calculate Moran's I
  aln.dists <- as.matrix(dist(sites$aln_position))
  aln.dists.inv <- 1/aln.dists
  diag(aln.dists.inv) <- 0
  moran.res <- Moran.I(sites$lrt_stat, aln.dists.inv)
  moran.dir <- moran.res$observed - moran.res$expected
  moran.p <- moran.res$p.value

  # Combine p-values using Fisher's method.
  fis.stat <- 2 * sum(-log(sites$pos.pval))
  fis.n <- 2 * nrow(sites)
  fis.p <- pchisq(fis.stat, df=fis.n, lower=F)

  # Combine the top N, where N=ceiling(length/10), p-values using Fisher's method.
  #blw.50 <- subset(sites, pos.pval < 0.5)
  #blw.10 <- subset(sites, pos.pval < 0.1)
  #top.fis.stat <- 2 * sum(-log(top.sites$pos.pval))
  #top.fis.n <- 2 * nrow(top.sites)
  #top.fis.p <- pchisq(top.fis.stat, df=top.fis.n, lower=F)

  # Use the weighted TPM method from
  # Zaykin et al. 2002
  source("~/src/greg-ensembl/projects/2xmammals/Tpmw.r")
  weights <- sites$ncod
  pvals <- sites$pos.pval
  loops <- 1000
  tpm.05 <- Tpmw(pvals, weights, 0.05, loops=loops)
  tpm.10 <- Tpmw(pvals, weights, 0.1, loops=loops)
  tpm.20 <- Tpmw(pvals, weights, 0.2, loops=loops)
  tpm.50 <- Tpmw(pvals, weights, 0.5, loops=loops)

  n <- nrow(sites)
  out.df <- data.frame(
    n_sites = nrow(sites),
    f_sites = nrow(sites) / gene$aln_length,
    mean_ncod = mean(sites$ncod),

    mean_omega = mean.omg,
    mean_mid_omega = mean.mid.omega,

    f_abv_1 = n.abv.1 / n,
    f_blw_1 = n.blw.1 / n,
    
    f_pos_20 = n.pos.20 / n,
    f_neg_20 = n.neg.20 / n,
    f_ntr_20 = n.ntr.20 / n,
    f_pos_10 = n.pos.10 / n,
    f_neg_10 = n.neg.10 / n,
    f_ntr_10 = n.ntr.10 / n,
    f_pos_05 = n.pos.05 / n,
    f_neg_05 = n.neg.05 / n,
    f_ntr_05 = n.ntr.05 / n,

    n_pos_05 = n.pos.05,
    n_pos_01 = n.pos.01,
    n_pos_bh = n.pos.bh,
    n_pos_fwer = n.pos.fwer,
    n_pos_gene_bh_10 = n.pos.gene.bh.10,
    n_pos_gene_bh_05 = n.pos.gene.bh.05,

    f_pos_01 = n.pos.01 / n,
    f_pos_bh = n.pos.bh / n,
    f_pos_fwer = n.pos.fwer / n,
    f_pos_gene_bh_10 = n.pos.gene.bh.10 / n,
    f_pos_gene_bh_05 = n.pos.gene.bh.05 / n,

    fis_p = fis.p,
    tpm_05 = tpm.05,
    tpm_10 = tpm.10,
    tpm_20 = tpm.20,
    tpm_50 = tpm.50,

    mantel_cor = mantel.cor,
    mantel_p = mantel.p,
    moran_dir = moran.dir,
    moran_p = moran.p
  )
  out.df
}

process.sites <- function(df, do.quantiles=T) {
  gc(F)

  # Two-tailed p-value.
  df[, 'pval'] <- 1 - pchisq(abs(df[, 'lrt_stat']), 1)
 
  # One-tailed p-value for pos-sel.
  abv.one <- df$lrt_stat > 0
  # Positive sites: sites with \omgml > 1 get p-value from 0 to 0.5,
  # sites with \omgml < 1 get p-value from 0.5 to 1.
  df[, 'pos.pval'] <- 1
  df[abv.one, 'pos.pval'] <- df[abv.one, 'pval'] / 2
  df[!abv.one, 'pos.pval'] <- 0.5 + (1 - df[!abv.one, 'pval']) / 2
  rm(abv.one)
  gc(F)

  # Do BH adjustment on separate one-tailed p-values.
  df$pos.bh <- p.adjust(df$pos.pval, method='BH')
  df$neg.bh <- p.adjust(1 - df$pos.pval, method='BH')
  gc(F)
    
  # Add p-set factor.
  df$clade <- factor(pset.to.alias(df[1, 'parameter_set_id']))
  gc(F)

  df$note <- factor(df$note)
  df$chr_name <- factor(df$chr_name)
  df$exon_position <- factor(df$exon_position)
  df$pfam_domain <- factor(df$pfam_domain)
  df$type <- factor(df$type)
  df$random <- factor(df$random)
  gc(F)              

  if (do.quantiles) {
    print("  doing quantiles...")
    # Add bl quantiles.
    gc(F)
    fn <- ecdf(df$nongap_bl)
    ng.bl <- pmin(df$nongap_bl, max(df$nongap_bl - 1e-10))
    df$bl_z <- fn(ng.bl)
    df$bl_factor <- round_any(df$bl_z, 0.1, ceiling)
    df$bl_factor <- as.factor(df$bl_factor)
    rm(fn)
    rm(ng.bl)
    gc(F)

    # Add lrt_stat z-scores
    fn <- ecdf(df$lrt_stat)
    df$lrt_q <- fn(df$lrt_stat)
    rm(fn)
    gc(F)

    # lrt_z is a true z-score.
    df$lrt_z <- (df$lrt_stat - mean(df$lrt_stat)) / sd(df$lrt_stat)

    # Try a two-sided LRT z-score.
    abv.i <- df$lrt_stat >= 0
    abv <- subset(df, lrt_stat >= 0, select='lrt_stat')
    blw <- subset(df, lrt_stat < 0, select='lrt_stat')
    fn.abv <- ecdf(abv$lrt_stat)
    fn.blw <- ecdf(blw$lrt_stat)
    rm(abv)
    rm(blw)
    abv.z <- fn.abv(df$lrt_stat)
    blw.z <- fn.blw(df$lrt_stat)
    gc(F)                    

    lrt_zz <- c()
    lrt_zz[abv.i] <- .5 + abv.z[abv.i] / 2
    lrt_zz[!abv.i] <- blw.z[!abv.i] / 2
    df$lrt_qq <- lrt_zz
    rm(lrt_zz)
    rm(abv.i)
    rm(abv.z)
    rm(blw.z)
    rm(fn.abv)
    rm(fn.blw)
    gc(F)              
  }
  df
}

test.summary <- function() {
  sites <- get.pset.sites(1, filter='orig', test=T)
  sites <- process.sites(sites)
  summarize.sites(sites, filter='test')
}

summarize.sites <- function(sites, filter='') {
    pset <- sites[1, 'parameter_set_id']
    pset.alias <- pset.to.alias(pset)

    lbl <- paste(filter, ' ', pset.alias, sep='')

    n.sites <- nrow(sites)  

    sub.above.one <- subset(sites, omega > 1)
    pos.sites <- subset(sub.above.one, pval < 0.05)
    pos.domains <- subset(pos.sites, !is.na(pfam_domain))  

    n.total.domains <- nrow(unique(sites[,c('data_id','pfam_domain')]))
    n.total.domain.types <- length(unique(sites$pfam_domain))
    n.total.genes <- length(unique(sites$data_id))

    n.pos.domains <- nrow(unique(pos.domains[,c('data_id','pfam_domain')]))
    n.pos.domain.types <- length(unique(pos.domains$pfam_domain))
    n.pos.genes <- length(unique(pos.sites$data_id))

    `f.less.0.5` <- nrow(subset(sites, omega < 0.5)) / nrow(sites)
    `f.less.1` <- nrow(subset(sites, omega < 1)) / nrow(sites)
    `f.gt.1` <- nrow(subset(sites, omega > 1)) / nrow(sites)
    `f.gt.1.5` <- nrow(subset(sites, omega > 1.5)) / nrow(sites)

    n.pos.sites <- nrow(pos.sites)
    pos.1 <- nrow(subset(sub.above.one, pval < 0.01))
    pos.5 <- nrow(subset(sub.above.one, pval < 0.05))
    pos.10 <- nrow(subset(sub.above.one, pval < 0.1))

    pos.bh.05 <- nrow(subset(sites, pos.bh < 0.05))
    pos.bh.10 <- nrow(subset(sites, pos.bh < 0.10))

    n.const <- nrow(subset(sites, note == 'constant'))
    n.syn <- nrow(subset(sites, note == 'synonymous'))

    pos.f5 <- nrow(subset(sites, omega > 1 & pval < 0.05))
    neg.f5 <- nrow(subset(sites, omega < 1 & pval < 0.05))
    neutral.f5 <- nrow(subset(sites, pval > 0.05))

    pos.f10 <- nrow(subset(sites, omega > 1 & pval < 0.1))
    neg.f10 <- nrow(subset(sites, omega < 1 & pval < 0.1))
    neutral.f10 <- nrow(subset(sites, pval > 0.1))

    med.bl <- median(sites$nongap_bl)
    mean.bl <- mean(sites$nongap_bl)
    sd.bl <- sd(sites$nongap_bl)
    med.ncod <- median(sites$ncod)

    omega.cap <- 3
    capped.omegas <- pmin(omega.cap, sites$omega)
    med.omega <- median(capped.omegas)
    mean.omega <- mean(capped.omegas)
    sd.omega <- sd(capped.omegas)
    max.omega <- max(capped.omegas)
    mean.abv <- mean(capped.omegas[capped.omegas > 0])
    sd.abv <- sd(capped.omegas[capped.omegas > 0])
    mean.blw <- mean(capped.omegas[capped.omegas < 1])
    sd.blw <- sd(capped.omegas[capped.omegas < 1])

    sorted.omgs <- sort(sites$omega)
    n.quart <- nrow(sites) / 10
    mid.omgs <- sorted.omgs[c(floor(.5*n.quart):floor(9.5*n.quart))]
    mean.mid.omega <- mean(mid.omgs)

    df.out <- data.frame(
      pset = pset,
      label = lbl,
      filter = filter,

      n.sites = n.sites,
      f.const = n.const / n.sites,
      f.syn = n.syn / n.sites,

      med.bl = med.bl,
      med.ncod = med.ncod,
      mean.omega = mean.omega,
      omega.cap = omega.cap,
      mean.mid.omega = mean.mid.omega,

      mean.bl = mean.bl,
      sd.bl = sd.bl,
      mean.abv = mean.abv,
      sd.abv = sd.abv,
      sd.omega = sd.omega,
      mean.blw = mean.blw,
      sd.blw = sd.blw,

      'f.less.0.5' = `f.less.0.5`,
      'f.less.1' = `f.less.1`,
      'f.gt.1' = `f.gt.1`,
      'f.gt.1.5' = `f.gt.1.5`,

      pos.10 = pos.10,
      pos.05 = pos.5,
      pos.01 = pos.1,
      pos.bh.05 = pos.bh.05,
      pos.bh.10 = pos.bh.10,

      pos.f05 = pos.f5 / n.sites,
      neg.f05 = neg.f5 / n.sites,
      neutral.f05 = neutral.f5 /n.sites,

      pos.f10 = pos.f10 / n.sites,
      neg.f10 = neg.f10 / n.sites,
      neutral.f10 = neutral.f10 / n.sites
    )
    df.out
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

plot.global.distribution <- function(sites, test=F) {
  test <- as.logical(test)

  sites$note <- as.character(sites$note)

  # Subset the data
  all.sites <- sites
  zero.sites <- subset(all.sites)
  sites <- subset(all.sites, omega > 0)

  print(table(zero.sites$note, useNA='ifany'))
  print(table(sites$note, useNA='ifany'))

  max.omega <- 3
  dx <- 0.02

  sites$omega <- pmin(sites$omega, max.omega)
  #print(head(sites, n=3))

  max.y <- 10 * 1000 * 1000 * 0.04
  if (test) {
    max.y <- nrow(all.sites) * 0.05
  }

  # Also plot the cdf of sites in the main histogram.
  cum.dx <- 0.005
  omega.f <- ecdf(all.sites$omega)
  lo.f <- ecdf(all.sites$omega_lower)
  hi.f <- ecdf(all.sites$omega_upper)
  omega.vals <- unique(round_any(all.sites$omega, cum.dx))
  lo.vals <- unique(round_any(all.sites$omega_lower, cum.dx))
  hi.vals <- unique(round_any(all.sites$omega_upper, cum.dx))
  omega.df <- data.frame(omega=omega.vals, cum=omega.f(omega.vals) * max.y)
  lo.df <- data.frame(omega=lo.vals, cum=lo.f(lo.vals) * max.y)
  hi.df <- data.frame(omega=hi.vals, cum=hi.f(hi.vals) * max.y)

  p <- ggplot(sites, aes(x=omega))
  p <- p + theme_bw()
  p <- p + geom_histogram(binwidth=dx, colour=NA, fill='black')

  clr <- rgb(0.2, 0.2, 0.8)
  sz <- 1.2
  ap <- 0.8
  p <- p + geom_line(data=omega.df, aes(y=cum), colour=clr, size=sz, alpha=ap)
  p <- p + geom_line(data=lo.df, aes(y=cum), colour=clr, linetype='dashed', size=sz, alpha=ap)
  p <- p + geom_line(data=hi.df, aes(y=cum), colour=clr, linetype='dashed', size=sz, alpha=ap)

  frmt <- function(x, ...) {
    sprintf("%.1e", x)
  }

  p <- p + scale_x_continuous(limits=c(0, max.omega))
  p <- p + scale_y_continuous(limits=c(0, max.y), formatter=frmt)
  p <- p + opts(
    axis.title.y = theme_blank(),
    axis.title.x = theme_blank(),
    axis.text.x = theme_blank()
  )

  max.zero.y <- 6 * 1000 * 1000
  if (test) {
    max.zero.y <- max(table(zero.sites$note, useNA='ifany'))
  }
  print(table(zero.sites$note, useNA='ifany'))
  zero.sites$note <- ifelse(is.na(zero.sites$note) | zero.sites$note == '', 'Non-synonymous', zero.sites$note)
  zero.sites$note <- factor(as.character(zero.sites$note), levels=c('constant', 'synonymous', 'Non-synonymous'),
    labels=c('Constant', 'Synonymous', 'Non-synonymous')
  )
  zero.sites$omega <- 0

  print(table(zero.sites$note, useNA='ifany'))

  q <- ggplot(zero.sites, aes(x=note, fill=note))
  q <- q + theme_bw()
  q <- q + geom_histogram(binwidth=0.04, width=0.5, colour='black')
  q <- q + scale_fill_manual(values=c('white', 'gray', 'black'))
  q <- q + scale_y_continuous(limits=c(0, max.zero.y * 1.05), formatter=frmt, expand=c(0,0))
  q <- q + opts(
    legend.position = 'none',
    axis.title.y = theme_blank(),
    axis.title.x = theme_blank(),
    axis.text.x = theme_blank(),
    axis.text.y = theme_blank(),
    strip.text.x = theme_blank(),
    panel.grid.minor = theme_blank(),
    panel.grid.major = theme_blank()
  )

  list(
    zero=q,
    hist=p
  )
}

plot.subclade.corrs <- function() {
  # First run the bsub.corrs() method to populate the correlations table!
  con <- connect(db()); 
  all.df <- dbGetQuery(con, 'select * from parallel_sites')
  dbDisconnect(con)

  # temporary fix: parse pset a and b from the label.
  lbl <- all.df$label
  pset.a <- substr(lbl, 1, 1)
  pset.b <- substr(lbl, 3, 3)
  all.df$pset_a <- as.integer(pset.a)
  all.df$pset_b <- as.integer(pset.b)

  min.y <- min(c(all.df$pos_est, all.df$pos_lo, all.df$neg_lo))
  max.y <- max(c(all.df$pos_est, all.df$pos_hi, all.df$neg_hi))
  min.y <- 0
  max.y <- 1600

  print(min.y)
  print(max.y)

  df.z <- data.frame()
  df.p <- data.frame()

  subclade.plots <- list(
    c(1, 2, 3, 4),
    c(2, 1, 3, 4),
    c(3, 1, 2, 4),
    c(4, 1, 2, 3)
  )

  mamm.plots <- list(
    c(1, 5, 6, 8),
    c(2, 5, 6, 8),
    c(3, 5, 6, 8)
  )

  df.plot <- function(df, type) {
    p <- ggplot(df, aes(x=as.numeric(thresh), y=pos_est, fill=pset_b, colour=pset_b, linetype=pset_b))
    p <- p + theme_bw()

    df$pset_b <- pset.to.alias(df$pset_b, factors=T)

    n.values <- length(levels(df$thresh))
    p <- p + geom_vline(xintercept=2, colour=gray(0.3))
    p <- p + geom_vline(xintercept=n.values-1, colour=gray(0.3))
    p <- p + geom_vline(xintercept=1, colour=gray(0.6))
    p <- p + geom_vline(xintercept=n.values, colour=gray(0.6))

    p <- p + geom_line(size=1)
    #p <- p + geom_ribbon(aes(ymin=pos_lo, ymax=pos_hi), alpha=0.3, colour=NA)
    p <- p + scale_fill_discrete("Other Group")
    p <- p + scale_colour_discrete("Other Group")
    p <- p + scale_linetype_discrete("Other Group")

    if (type == 'z') {
       x.lbl <- "Sitewise Z-score Threshold"
    } else if (type == 'p') {
       x.lbl <- "Signed LRT Threshold"
    } else if (type == 'z2') {
       x.lbl <- "Sitewise Z2 Threshold"
    }
  
    p <- p + scale_x_continuous(x.lbl, expand=c(0.05,0),
      breaks=1:n.values,
      labels=sprintf("%.2f", as.numeric(levels(df$thresh)))
    )
    #p <- p + scale_y_continuous("Odds Ratio", limits=c(10^0.5, 200))
    p <- p + scale_y_log10("Odds Ratio")#, limits=c(10^0.5, 500))
    #p <- p + coord_trans(y="log10")
    p <- p + opts(
      title=df[1, 'pset_a'],
      axis.text.x = theme_text(angle=90, size=8, hjust=1)
    )
    p
  } 

  do.plots.todo <- function(plots.list, type) {
    df.collect <- data.frame()
    for (i in 1:length(plots.list)) {
      cur.psets <- plots.list[[i]]
      print(cur.psets)
      cur.df <- subset(all.df, pset_a == cur.psets[1] & pset_b %in% cur.psets[-1])

      print(type)
      df <- subset(cur.df, thresh_type == type)
  
      df$pset_a <- pset.to.alias(df$pset_a, factors=T)
      df$pset_b <- pset.to.alias(df$pset_b, factors=T)

      if (type != 'z2') {
        df2 <- df
        df2$thresh <- - df2$thresh
        df2$pos_est <- df2$neg_est
        df2$pos_lo <- df2$neg_lo
        df2$pos_hi <- df2$neg_hi
        df2$pos_n_a <- df2$neg_n_a
        df2$pos_n_b <- df2$neg_n_b
        df2$pos_n_ab <- df2$neg_n_ab
        df <- rbind(df, df2)
      }

      df$i <- paste(i, df[1, 'pset_a'])

      if (type == 'p') {
        df$thresh <- sign(df$thresh) * qchisq(abs(df$thresh), df=1)
        df$thresh <- factor(df$thresh)
      } else {
        df$thresh <- factor(df$thresh)
      }
      df.collect <- rbind(df.collect, df)
    }
    df.collect
  }

  df <- do.plots.todo(subclade.plots, 'p')
  p.all <- df.plot(df, 'p')
  p.all <- p.all + facet_grid(. ~ i)
  pdf(file=scratch.f("subclade.corr.all.p.pdf"), width=9, height=2.5)
  print.ggplot(p.all)
  dev.off()

  df <- do.plots.todo(subclade.plots, 'z')
  p.all <- df.plot(df, 'z')
  p.all <- p.all + facet_grid(. ~ i)
  pdf(file=scratch.f("subclade.corr.all.z.pdf"), width=9, height=2.5)
  print.ggplot(p.all)
  dev.off()

  df <- do.plots.todo(subclade.plots, 'z2')
  p.all <- df.plot(df, 'z2')
  p.all <- p.all + facet_grid(. ~ i)
  pdf(file=scratch.f("subclade.corr.all.z2.pdf"), width=9, height=2.5)
  print.ggplot(p.all)
  dev.off()

  df <- do.plots.todo(mamm.plots, 'p')
  p.all <- df.plot(df, 'p')
  p.all <- p.all + facet_grid(. ~ i)
  pdf(file=scratch.f("subclade.corr.mamm.p.pdf"), width=9, height=2.5)
  print.ggplot(p.all)
  dev.off()

  df <- do.plots.todo(mamm.plots, 'z')
  p.all <- df.plot(df, 'z')
  p.all <- p.all + facet_grid(. ~ i)
  pdf(file=scratch.f("subclade.corr.mamm.z.pdf"), width=9, height=2.5)
  print.ggplot(p.all)
  dev.off()

  df <- do.plots.todo(mamm.plots, 'z2')
  p.all <- df.plot(df, 'z2')
  p.all <- p.all + facet_grid(. ~ i)
  pdf(file=scratch.f("subclade.corr.mamm.z2.pdf"), width=9, height=2.5)
  print.ggplot(p.all)
  dev.off()

}

get.merged.sites <- function(pset.a, pset.b, xtra, test=F, filter.fn=filter.default) {
  pset.a <- as.integer(pset.a)
  pset.b <- as.integer(pset.b)
  test <- as.logical(test)

  print("  getting a...")
  sites.a <- filter.fn(get.pset.sites(pset.a, test=test))
  print("  getting b...")
  sites.b <- filter.fn(get.pset.sites(pset.b, test=test))

  sel.cols <- c(xtra, 'data_id', 'aln_position', 'clade')
  sites.a <- subset(sites.a, select=sel.cols)
  sites.b <- subset(sites.b, select=sel.cols)

  print(paste("  merging ", pset.a, ' - ', pset.b, "...", sep=''))
  merged.s <- merge(sites.a, sites.b, by=c('data_id', 'aln_position'))
  merged.s
}

filter.default <- function(sites) {
  sites <- filter.branchlength(sites)
  sites <- filter.random(sites)
  sites <- filter.ncod(sites)
  sites <- filter.omega(sites)
  sites
}

test.filters <- function() {
  sites <- get.pset.sites(1, filter='orig', test=T)
  sites <- process.sites(sites)
  print(str(sites))
  print(paste("orig: ",nrow(sites)))
  x <- filter.default(sites)
  print(paste("default: ",nrow(x)))
  x <- filter.default(filter.default(sites))
  print(paste("double default: ",nrow(x)))
  x <- filter.stringent(sites)
  print(paste("stringent: ",nrow(x)))
  x <- filter.pfam(sites)
  print(paste("pfam: ",nrow(x)))
  x <- filter.clusters(sites)
  print(paste("clusters: ",nrow(x)))
}

filter.stringent <- function(sites) {
  sites <- filter.ncod.stringent(sites)
  sites <- filter.gc(sites)
  sites <- filter.dups(sites)
  sites <- filter.clusters(sites)
  sites <- filter.default(sites)

  sites
}

filter.clusters <- function(sites, return.inverse=F) {
  sites <- filter.default(sites)

  # Collect data_id / aln_positions of alignment segments w/ bad windows.
  con <- connect(db())
  cmd <- sprintf("select * from win_baddies")
  baddies <- dbGetQuery(con, cmd)
  disconnect(con)

  # Create a 'bad.df' data frame indicating data_id / aln_pos
  #  sites contained in clusters above the 99.9th percentile for any
  #  given species or ancestral node.
  win.size <- baddies[1, 'win_size']
  bad.data_ids <- c()
  bad.aln_pos <- c()
  for (i in 0:(win.size-1)) {
    d_id <- baddies$data_id
    aln_pos <- baddies$aln_position + i
    bad.data_ids <- c(bad.data_ids, d_id)
    bad.aln_pos <- c(bad.aln_pos, aln_pos)
  }
  bad.df <- data.frame(
    data_id = bad.data_ids,
    aln_position = bad.aln_pos
  )

  sum.df <- ddply(sites, .(data_id), function(x) {
    data.frame(
      n_sites = nrow(x)
    )
  })

  # Merge the bad data_id/aln_positions with the genes data frame
  # to identify genes with a high percentage of bad sites.
  #pp("bad.df: %d", nrow(bad.df))
  #genes <- get.genes()
  #genes <- subset(genes, select=c('data_id', 'aln_length'))
  #pp("genes: %d", nrow(genes))
  #mgd <- merge(genes, bad.df)
  #pp("merged: %d", nrow(mgd))
  
  mgd <- merge(sum.df, bad.df)

  bad.summaries <- ddply(mgd, .(data_id), function(x) {
    sum.bad <- nrow(x)
    n.bad.positions <- length(unique(x$aln_position))
    n.total <- x[1, 'n_sites']
    data.frame(
      n.total = n.total,
      bad.total = sum.bad,
      n.bad.positions = n.bad.positions,
      f.bad.positions = n.bad.positions / n.total
    )
  })
  top.q <- quantile(bad.summaries$f.bad.positions, 0.9)
  bad.genes <- subset(bad.summaries, f.bad.positions > top.q)
  bad.gene.ids <- bad.genes$data_id

  rm(bad.summaries)
  rm(bad.genes)
  rm(mgd)

  bad.df <- bad.df[!duplicated(bad.df),]

  # Merge the cluster_sites data frame with the sites data frame...
  # There should be no loss here, so double-check that.  Remove any
  # sites that should be filtered based on high NS subs within
  # 15-codon windows.
  bad.str <- paste(bad.df$data_id, bad.df$aln_position, sep=' ')
  #print(sprintf("%d total bad sites!", length(bad.str)))
  rm(bad.df)
  sites.str <- paste(sites$data_id, sites$aln_position)

  bad.sites <- sites.str %in% bad.str
  bad.genes <- sites$data_id %in% bad.gene.ids

  #print(sprintf("Removed %d bad sites!", sum(bad.sites)))
  if (return.inverse) {
    sites[bad.sites | bad.genes,]
  } else {
    sites[!(bad.sites | bad.genes),]
  }
}

filter.test <- function() {
  sites <- get.pset.sites(1, test=T)
  sites <- filter.stringent(sites)
}

get.cluster <- function(taxon_id, win.size=15, test=F) {
  cluster.f <- scratch.f(paste('win_', taxon_id, win.size, '_', as.character(test), '.Rdata', sep=''))
  if (!file.exists(cluster.f)) {
    con <- connect(db())
    null.str <- 'and taxon_id is not null and aln_position is not null and data_id is not null'
    if (test) {
      cmd <- sprintf("select * from windows_%d where taxon_id=%d %s limit 100000;", win.size, taxon_id, null.str)
    } else {
      cmd <- sprintf("select * from windows_%d where taxon_id=%d %s;", win.size, taxon_id, null.str)
    }
    clusters <- dbGetQuery(con, cmd)
    disconnect(con)
    save(clusters, file=cluster.f)
  }
  load(cluster.f)
  clusters
}

process.cluster <- function(taxon_id, win.size=15, test=F) {
  win <- get.cluster(taxon_id, win.size=win.size, test=test)

  is.leaf <- win[1, 'is_leaf']

  s.e2 <- quantile(win$n_s_subs, 1 - 1e-2)
  s.e3 <- quantile(win$n_s_subs, 1 - 1e-3)
  s.e4 <- quantile(win$n_s_subs, 1 - 1e-4)
  ns.e2 <- quantile(win$n_ns_subs, 1 - 1e-2)
  ns.e3 <- quantile(win$n_ns_subs, 1 - 1e-3)
  ns.e4 <- quantile(win$n_ns_subs, 1 - 1e-4)

  # Find the number of nsyn subs with > than the synonymous 99th percentile
  n.ns.se3 <- nrow(subset(win, n_ns_subs > s.e3))
  n.ns.se4 <- nrow(subset(win, n_ns_subs > s.e4))
  n.ns.nse3 <- nrow(subset(win, n_ns_subs > ns.e3))
  n.ns.nse4 <- nrow(subset(win, n_ns_subs > ns.e4))

  s.tbl <- table(win$n_s_subs)
  ns.tbl <- table(win$n_ns_subs)
  s.counts.str <- paste(s.tbl, collapse=' ')
  ns.counts.str <- paste(ns.tbl, collapse=' ')

  ns.windows <- subset(win, n_ns_subs > 0)
  conf <- weighted.mean(ns.windows$mean_ns_confidence, ns.windows$n_ns_subs)

  clust.df <- data.frame(
    label = paste(taxon_id, win.size, test, sep=' '),
    taxon_id = taxon_id,
    win_size = win.size,
    s_counts = s.counts.str,
    ns_counts = ns.counts.str,
    is_leaf = is.leaf,
    n_clusters = nrow(win),
    mean_conf = conf,
    syn_e2 = s.e2,
    syn_e3 = s.e3,
    syn_e4 = s.e4,
    nsyn_e2 = ns.e2,
    nsyn_e3 = ns.e3,
    nsyn_e4 = ns.e4,
    n_ns_se3 = n.ns.se3,
    n_ns_se4 = n.ns.se4,
    n_ns_nse3 = n.ns.nse3,
    n_ns_nse4 = n.ns.nse4
  )

  con <- connect(db())
  write.or.update(clust.df, 'win_thresholds', con, 'label')
  disconnect(con)

  # Output all the bad windows...
  bad.wins <- subset(win, n_ns_subs > ns.e3)
  bad.df <- subset(bad.wins, select=c('data_id', 'aln_position', 'taxon_id', 'n_ns_subs', 'is_leaf'))
  bad.df$label <- paste(bad.df$data_id, bad.df$taxon_id, bad.df$aln_position)
  bad.df$win_size <- win.size
  
  con <- connect(db())
  write.or.update(bad.df, 'win_baddies', con, 'label')
  disconnect(con)
}

bsub.collect.clusters <- function() {
  con <- connect(db())
  cmd <- sprintf("select distinct(taxon_id) from windows_15 limit 75;")
  df <- dbGetQuery(con, cmd)
  disconnect(con)
  tx.ids <- df$taxon_id

  for (id in tx.ids) {
    bsub.function('collect_clusters', extra.args=paste(id, 'F', sep=' '), mem=6)
  }
}

filter.dups <- function(sites) {
  genes <- get.genes()
  bad.genes <- subset(genes, dup_species_count > 10)
  bad.ids <- bad.genes$data_id
  subset(sites, !(data_id %in% bad.ids))
}


filter.gc <- function(sites) {
  # Remove sites in genes with the top and bottom 10% GC content.
  # Get the GC content from the genes vector.
  genes <- get.genes()

  gc.fn <- ecdf(genes$gc_3)
  genes$gc_cum <- gc.fn(genes$gc_3)
  bad.genes <- subset(genes, gc_cum < .10 | gc_cum > .90)
  bad.ids <- bad.genes$data_id

  subset(sites, !(data_id %in% bad.ids))
}

filter.omega <- function(sites) {
  subset(sites, omega < omega_upper & omega < 50)
}

filter.branchlength <- function(sites) {
  subset(sites, bl_z > 0.1)
}

filter.random <- function(sites) {
  subset(sites, is.na(random))
}

filter.ncod <- function(sites) {
  subset(sites, ncod > 3)
}

filter.ncod.stringent <- function(sites) {
  # Only allow sites with at least 75% of the max ncod.
  max.ncod <- max(sites$ncod)
  ncod.thresh <- floor(max.ncod * 0.75)
  subset(sites, ncod >= ncod.thresh)
}

filter.pfam <- function(sites) {
  sites <- filter.default(sites)
  subset(sites, !is.na(pfam_domain))
}

filter.pfam.stringent <- function(sites) {
  sites <- filter.stringent(sites)
  subset(sites, !is.na(pfam_domain))
}

pset.df <- function(factors=F) {
  map.list <- list(
    '1' = 'Primates',
    '2' = 'Glires',
    '3' = 'Laurasiatheria',
    '4' = 'Atlantogenata',
    '5' = 'Eutheria',
    '6' = 'Mammalia',
    '7' = 'Sparse Glires',
    '8' = 'Sparse Mammalia',
    '9' = 'HQ Mammalia',
    '10' = 'HMRD'
  )
  char.list <- list(
    '1' = 'p',
    '2' = 'g',
    '3' = 'l',
    '4' = 'a',
    '5' = 'e',
    '6' = 'm',
    '7' = 'sg',
    '8' = 'sm',
    '9' = 'h',
    '10' = 'f'
  )
  df <- data.frame(
    pset_id = as.numeric(names(map.list)),
    pset = names(map.list),
    label = unlist(map.list),
    char = unlist(char.list),
    stringsAsFactors=factors
  )
  if (factors) {
    df$label <- factor(df$pset, levels=names(map.list), labels=unlist(map.list))
  }
  df
}

dist.to.factor <- function(dists) {
  factor(dists, levels=c('lnorm', 'gamma', 'weibull', 'beta', 'exp'), labels=c('Lognormal', 'Gamma', 'Weibull', 'Beta', 'Exponential'))
}

filter.to.alias <- function(filters, factors=F) {
  fct <- factor(filters,
    levels=c('none', 'default', 'stringent', 'pfam', 'pfam_stringent', 'clusters', 'clusters_inverse'),
    labels = c(
      'None',
      'Default',
      'Stringent',
      'Pfam',
      'Pfam Stringent',
      'Clusters',
      '(WCS)'
    )
  )
  if (factors) {
    fct
  } else {
    as.character(fct)    
  }
}

pset.to.alias <- function(psets, factors=F) {
  map <- pset.df(factors=factors)
  match.inds <- match(psets, map$pset)
  map$label[match.inds]
}

pset.to.short <- function(psets, factors=F) {
  map <- pset.df(factors=factors)
  match.inds <- match(psets, map$pset)
  map$char[match.inds]
}


write.subsets.table <- function() {
  source("~/src/greg-ensembl/projects/phylosim/PhyloSimPlots.R")

  taxids.list <- list(
    '1' = c(9478, 9483, 9544, 9593, 9598, 9601, 9606, 30608,
  30611, 61853),
    '2' = c(9978, 9986, 10020, 10090, 10116, 10141, 43179),
    '3' = c(9365, 9615, 9646, 9685, 9739, 9796, 9823, 9913, 30538,
  42254, 59463, 132908),
    '4' = c(9358, 9361, 9371, 9785, 9813),
    '5' = c(9358, 9361, 9365, 9371, 9478, 9483, 9544, 9593, 9598,
  9601, 9606, 9615, 9646, 9685, 9739, 9785, 9796, 9813, 9823, 9913,
  9978, 9986, 10020, 10090, 10116, 10141, 30538, 30608, 30611, 37347,
  42254, 43179, 59463, 61853, 132908),
    '6' = c(9258, 9315, 9358, 9361, 9365, 9371, 9478, 9483, 9544,
  9593, 9598, 9601, 9606, 9615, 9646, 9685, 9739, 9785, 9796, 9813,
  9823, 9913, 9978, 9986, 10020, 10090, 10116, 10141, 13616, 30538,
  30608, 30611, 37347, 42254, 43179, 59463, 61853, 132908),
    '7' = c(10090, 10116, 10020, 43179, 10141),
    '8' = c(9258, 9315, 9361, 9785, 9606, 10090, 9615)
  )

  aliases <- llply(taxids.list, function(x) {
    paste(sort(taxid.to.alias(x)), collapse=', ')
  })
  counts <- llply(taxids.list, function(x) {
    length(x)
  })

  tree <- read.tree(file="~/src/greg-ensembl/ensembl-compara-63/scripts/pipeline/species_tree_blength.nh")
  all.tips <- tree$tip.label

  get.subset.length <- function(x, total.length=F) {
    x <- taxid.to.alias(x, binomials=T)
    x <- gsub(' ', '_', x)
    remove.tips <- setdiff(all.tips, x)
    remaining.tree <- drop.tip(tree, remove.tips)

    if (total.length) {
      ret.val <- sum(remaining.tree$edge.length)
    } else {
      ret.val <- mean.path.length(remaining.tree)
    }
    ret.val
  }
  mpls <- llply(taxids.list, get.subset.length)
  sizes <- llply(taxids.list, get.subset.length, total.length=T)

  genes.df <- get.genes.split()
  pset.mpls <- dlply(genes.df, .(pset), function(x) {
    median(x$mpl, na.rm=T)
  })
  pset.bls <- dlply(genes.df, .(pset), function(x) {
    median(x$bl, na.rm=T)
  })

  pop.size <- as.integer(round(c(
    20000, # Primates
    230000, # Glires
    34410, # Laurasiatheria
    30000, # Atlantogenata
    110000, # Eutheria
    120000, # Mammalia
    230000, # Sparse G
    120000 # Sparse M
  )))

  out.df <- data.frame(
#    'Pset' = names(taxids.list),
    'Name' = pset.to.alias(names(taxids.list)),
    'Species' = paste("\\tiny{", "(", unlist(counts), ") ", unlist(aliases), "}", sep=""),
    'Ne' = as.integer(pop.size),
    'MPL' = unlist(mpls),
    'Total' = unlist(sizes),
    'slrMPL' = 0,
    'slrTotal' = 0
  )
  out.df[names(pset.mpls), 'slrMPL'] <- unlist(pset.mpls)
  out.df[names(pset.bls), 'slrTotal'] <- unlist(pset.bls)

  xt <- xtable(out.df)
  xt <- color.column(xt, 'Ne')
  xt <- color.column(xt, 'MPL')
  xt <- color.column(xt, 'Total')
  xt <- color.column(xt, 'slrMPL')
  xt <- color.column(xt, 'slrTotal')
  print.latex(xt, "species_set_summary.txt")

  #print(out.df)
}

write.summary.tables <- function() {
  con <- connect(db()); 
  summaries <- dbGetQuery(con, 'select * from summaries')
  disconnect(con)

  summaries <- cbind(data.frame(
    lbl=pset.to.alias(summaries$pset),
    pset.fact=pset.to.alias(summaries$pset, factors=T),
    flt=filter.to.alias(summaries$filter),
    flt.fact=filter.to.alias(summaries$filter, factors=T)
  ), summaries)

  df <- subset(summaries, pset %in% c(1, 6))
  df$crap <- df$filter == 'clusters_inverse'
  df <- df[order(df$crap, df$pset.fact, df$flt.fact),]
  df[which(df$pset == 1)[-1], 'lbl'] <- NA
  df[which(df$pset == 6)[-1], 'lbl'] <- NA
  df[df$crap, 'lbl'] <- pset.to.alias(df[df$crap, 'pset'])

  first.split <- sum(!df$crap) / 2
  sec.split <- sum(!df$crap)
  .write.summary.table(df, prefix='filter_summaries', split.at.indices=c(first.split, sec.split))

  df <- subset(summaries, filter == 'default')
  df <- df[order(df$med_bl),]
  df[which(df$filter == 'default')[-1], 'flt'] <- NA
  .write.summary.table(df, prefix='pset_summaries')
}

.write.summary.table <- function(df, prefix, split.at.indices=c()) {
  df$lbl <- gsub('_', '\\\\_', df$lbl)

  # Do some formatting...
  df$f_nsyn <- 1 - (df$f_syn + df$f_const)
  df$n_sites_s <- sprintf("%.2e", df$n_sites)
  to.pct <- function(df, pct.columns) {
    for (cl in pct.columns) {
      df[, cl] <- df[, cl] * 100
    }
    df
  }
  pct.cols <- c('f_const', 'f_syn', 'f_nsyn', 'pos_f05', 'neg_f05', 'neutral_f05', 'pos_f10', 'neg_f10', 'neutral_f10',
    'f_less_0_5', 'f_less_1', 'f_gt_1', 'f_gt_1_5')
  df <- to.pct(df, pct.cols)

  df$med_ncod <- as.integer(df$med_ncod)
  df$pos_10 <- as.integer(df$pos_10)
  df$pos_05 <- as.integer(df$pos_05)
  df$pos_01 <- as.integer(df$pos_01)
  df$pos_bh_05 <- as.integer(df$pos_bh_05)

  df$f_pos_10 <- df$pos_10 / df$n_sites
  df$f_pos_05 <- df$pos_05 / df$n_sites
  df$f_pos_01 <- df$pos_01 / df$n_sites
  df$f_pos_bh_05 <- df$pos_bh_05 / df$n_sites
  df <- to.pct(df, c('f_pos_10', 'f_pos_05', 'f_pos_01', 'f_pos_bh_05'))

  first.tbl.cols <- c('lbl', 'flt', 'n_sites_s', 'f_const', 'f_syn',
  'f_nsyn', 'med_ncod', 'med_bl', 'mean_bl', 'sd_bl', 'mean_omega',
  'sd_omega', 'f_less_0_5', 'f_less_1', 'f_gt_1', 'f_gt_1_5')

  second.tbl.cols <- c('lbl', 'flt', 'pos_10', 'f_pos_10', 'pos_05',
  'f_pos_05', 'pos_01', 'f_pos_01', 'pos_bh_05', 'f_pos_bh_05', 'neg_f10',
  'neutral_f10', 'pos_f10', 'neg_f05', 'neutral_f05', 'pos_f05')
    
  first.t <- subset(df, select=first.tbl.cols)
  second.t <- subset(df, select=second.tbl.cols)

  ccs <- c('med_bl', 'mean_bl', 'sd_bl', 'med_ncod', 'mean_omega', 'sd_omega', 'f_const',
  'f_syn', 'f_nsyn', 'f_less_0_5', 'f_less_1', 'f_gt_1', 'f_gt_1_5',
  'f_pos_10', 'f_pos_05', 'f_pos_01', 'f_pos_bh_05', 
  'neg_f05', 'neutral_f05', 'pos_f05', 'neg_f10', 'neutral_f10', 'pos_f10')

  if (length(split.at.indices) > 0) {
    xt1 <- xtable(first.t)
    xt2 <- xtable(second.t)

    ind <- which(colnames(xt2) == 'f_pos_bh_05') + 1
    attr(xt2, 'digits')[ind] <- 3
    
    xt1.list <- list()
    xt2.list <- list()
    cur.lo <- 1
    split.at.indices <- c(split.at.indices, nrow(xt1))
    for (i in 1:length(split.at.indices)) {
      cur.hi <- split.at.indices[i]
      xt1.top <- xt1[cur.lo:cur.hi,]
      xt2.top <- xt2[cur.lo:cur.hi,]
      skip.clr <- F
      if (i == 3) {
        skip.clr <- T
      }
      xt1.top <- color.columns(xt1.top, intersect(ccs,first.tbl.cols), log=F, skip.coloring=skip.clr)
      xt2.top <- color.columns(xt2.top, intersect(ccs,second.tbl.cols), log=F, skip.coloring=skip.clr)

      xt1.list[[i]] <- xt1.top
      xt2.list[[i]] <- xt2.top

      cur.lo <- cur.hi+1
    }
  } else {
    xt1 <- xtable(first.t)
    xt1 <- color.columns(xt1, intersect(ccs,first.tbl.cols), log=T)
    xt1.list <- list(xt1)
    xt2 <- xtable(second.t)

    ind <- which(colnames(xt2) == 'f_pos_bh_05') + 1
    attr(xt2, 'digits')[ind] <- 3

    xt2 <- color.columns(xt2, intersect(ccs,second.tbl.cols), log=T)
    xt2.list <- list(xt2)
  }
  
  add.pct <- function(xt, colb) {
    paste('(', trim(colb), ')', sep='')
  }

  for (i in 1:length(xt2.list)) {
    x <- xt2.list[[i]]
    x[, 'f_pos_10'] <- add.pct(x, x$f_pos_10)
    x[, 'f_pos_05'] <- add.pct(x, x$f_pos_05)
    x[, 'f_pos_01'] <- add.pct(x, x$f_pos_01)
    x[, 'f_pos_bh_05'] <- add.pct(x, x$f_pos_bh_05)
    xt2.list[[i]] <- x
  }

  ff1 <- scratch.f(paste(prefix, "_1.txt", sep=''))
  ff2 <- scratch.f(paste(prefix, "_2.txt", sep=''))

  for (i in 1:length(xt1.list)) {
    appnd <- (i > 1)
    print.latex(xt1.list[[i]], ff1, append=appnd)
    if (i != length(xt1.list)) {
      cat("\\midrule", file=ff1, append=T)
    }
  }
  for (i in 1:length(xt2.list)) {
    appnd <- (i > 1)
    print.latex(xt2.list[[i]], ff2, append=appnd)
    if (i != length(xt2.list)) {
      cat("\\midrule", file=ff2, append=T)
    }
  }
}

do.poisson.sim <- function(dnds, neutral_thresh=1) {
  n <- 1000 * 1000
  print(sprintf("  %.2f", dnds))
  for (bl in seq(from=0.05, to=10, length.out=20)) {
    #print(sprintf("    %.2f", bl))
    for (distr_shape in c('const_bl_no_inv', 'constant', 'gamma_2', 'gamma_0.5')) {
      #print(sprintf("      %s", distr_shape))

      if (distr_shape == 'const_bl_no_inv') {
        bl.dist <- bl
      } else {
        # Model a normally-distributed branchlength.
        bl.dist <- rnorm(n=n, mean=bl, sd=bl * .67)
        # Cut off the BL distribution at a 0.1 quantile.
        bl.dist <- pmax(quantile(bl.dist, probs=c(0.1)), bl.dist)
      }

      # Use constant DNDS or sample from a distribution
      if (distr_shape %in% c('constant', 'const_no_inv', 'const_bl_no_inv')) {
        dnds.dist <- dnds
      } else if (distr_shape %in% c('weibull_0.5', 'weibull_0.8', 'weibull_1')) {
        weibull.shape <- 1
        if (distr_shape == 'weibull_0.5') {
          weibull.shape <- 0.5
        } else if (distr_shape == 'weibull_0.8') {
          weibull.shape <- 0.8
        }
        weibull.scale <- dnds / gamma(1 + 1 / weibull.shape)
        dnds.dist <- rweibull(n, shape=weibull.shape, scale=weibull.scale)
      } else if (distr_shape %in% c('gamma_1', 'gamma_0.5', 'gamma_2')) {
        gamma.shape <- 1
        if (distr_shape == 'gamma_0.5') {
          gamma.shape <- 0.5
        } else if (distr_shape == 'gamma_2') {
          gamma.shape <- 2
        }
        gamma.scale <- dnds / gamma.shape
        dnds.dist <- rgamma(n, shape=gamma.shape, scale=gamma.scale)
      }
      #print(mean(dnds.dist))

      # Consider a 'neutral_thresh' which is the quantile above which to push sites     
      # to neutral.
      if (neutral_thresh < 1) {
        # Figure out the dnds level at the given quantile.
        dnds.q <- quantile(dnds.dist, probs=neutral_thresh)
        dnds.dist[dnds.dist > dnds.q] <- 1
      }

      # Sample synonymous and non-synonymous subs as a poisson process
      syn.subs <- rpois(n, lambda=bl.dist)
      nsyn.subs <- rpois(n, lambda=bl.dist * dnds.dist)

      # Add a %-age of invariant sites.
      f.invariant <- 0.065
      if (distr_shape == 'const_no_inv') {
        f.invariant <- 0
      }
      invariant.sites <- runif(n=n, min=0, max=1) <= f.invariant
      syn.subs[invariant.sites] <- 0
      nsyn.subs[invariant.sites] <- 0

      # Count constant, nsyn, and syn sites.
      const.sites <- sum(syn.subs == 0 & nsyn.subs == 0)
      syn.sites <- sum(syn.subs > 0 & nsyn.subs == 0)
      nsyn.sites <- sum(nsyn.subs > 0)

      # Calculate sitewise dN/dS estimates
      nsyn.syn <- syn.subs > 0 & nsyn.subs > 0
      nsyn.no.syn <- syn.subs == 0 & nsyn.subs > 0
      syn.no.nsyn <- nsyn.subs == 0 & syn.subs > 0
      const <- syn.subs == 0 & nsyn.subs == 0

      # Synonymous, no nsyn: zero
      syn.dnds <- nsyn.subs[syn.no.nsyn] / syn.subs[syn.no.nsyn]
      # Constant: zero
      const.dnds <- nsyn.subs[const]
      # Nsyn, no syn: put at dnds=1
      nsyn.dnds <- rep(1, times=sum(nsyn.no.syn))
      # Syn & nsyn: take ratio
      both.dnds <- nsyn.subs[nsyn.syn] / syn.subs[nsyn.syn]
          
      all.dnds <- c(syn.dnds, const.dnds, nsyn.dnds, both.dnds)
          
      # Clip to max value.
      omega.top.q <- 2
      all.dnds <- pmin(all.dnds, omega.top.q)

      mean.dnds <- mean(all.dnds)
      sd.dnds <- sd(all.dnds)
      mean.nz <- mean(all.dnds[all.dnds > 0])
      sd.nz <- sd(all.dnds[all.dnds > 0])

      mean.blw <- mean(all.dnds[all.dnds < 1])
      sd.blw <- sd(all.dnds[all.dnds < 1])

      f_lt_0_5 <- sum(all.dnds < 0.5) / length(all.dnds)
      f_lt_1 <- sum(all.dnds < 1) / length(all.dnds)
      f_gt_1 <- sum(all.dnds > 1) / length(all.dnds)
      f_gt_1_5 <- sum(all.dnds > 1.5) / length(all.dnds)

      lbl <- paste(dnds, bl, distr_shape, neutral_thresh, sep=" ")

      out.df <- data.frame(
        lbl = lbl,
        neutral_thresh = neutral_thresh,
        dnds = dnds,
        bl = bl,
        distr_shape = distr_shape,
        n_const = const.sites,
        n_syn = syn.sites,
        n_nsyn = nsyn.sites,
        f_const = const.sites / n,
        f_syn = syn.sites / n,
        f_nsyn = nsyn.sites / n,
        omega_top_q = omega.top.q,
        mean_dnds = mean.dnds,
        mean_nz = mean.nz,
        mean_blw = mean.blw,
        sd_dnds = sd.dnds,
        sd_nz = sd.nz,
        sd_blw = sd.blw,
        f_lt_0_5 = f_lt_0_5,
        f_lt_1 = f_lt_1,
        f_gt_1 = f_gt_1,
        f_gt_1_5 = f_gt_1_5
      )

      con <- connect(db())
      write.or.update(out.df, 'poisson', con, 'lbl')
      dbDisconnect(con)
    }
  }
}

bsub.poisson.sims <- function() {
  # Trying to see if the proportion of synonymous sites is invariant to the branch
  # length of the species set being investigated...
  # Simple model: we have a syn rate, a (lower) nsyn rate, and observations.
  # We observe after time T, and classify into three site types: constant, syn, nsyn

  for (dnds in seq(from=0.1, to=0.7, by=0.1)) {
    neutral_thresh <- 1
#    for (neutral_thresh in c(0.8, 0.9, 1)) {
      bsub.function('poisson_sim', mem=8, extra.args=paste(dnds, neutral_thresh))
#    }
  }
}

plot.poisson.sims <- function() {
  con <- connect(db()); 
  out.df <- dbGetQuery(con, 'select * from poisson')
  dbDisconnect(con)

  plot.vals <- c('f_const', 'f_syn', 'f_nsyn', 'mean_dnds', 'f_lt_1', 'f_gt_1')
  plot.lbls <- c('Const.', 'Syn.', 'Nsyn', 'Mean dN/dS', '< 1', '> 1')

  comb.df <- data.frame()
  for (val in plot.vals) {
    cur.df <- out.df
    cur.df$val.lbl <- val
    cur.df$val <- cur.df[, val]
    if (val == 'f_lt_1') {
      cur.df$val <- cur.df$val - 0.6
    } else if (val == 'f_gt_1') {
      cur.df$val <- cur.df$val*3 + 0.3
    }
    comb.df <- rbind(comb.df, cur.df)
  }
  comb.df$pset <- NA

  comb.emp <- data.frame()

  empirical.vals <- 
  con <- connect(db()); 

  empirical.df <- dbGetQuery(con, 'select pset, mean_bl as bl,
    f_const, f_syn, (1 - f_const - f_syn) as f_nsyn, mean_omega as
    mean_dnds, sd_omega as sd_dnds, mean_abv, sd_abv, mean_blw,
    sd_blw, f_less_0_5 as f_lt_0_5, f_less_1 as f_lt_1, f_gt_1_5,
    f_gt_1 from summaries where filter="default" order by pset' )

  dbDisconnect(con)
  empirical.df$pset <- pset.to.short(empirical.df$pset)
  empirical.df$dnds <- 1
  print(head(empirical.df))

  for (val in plot.vals) {
    cur.df <- empirical.df
    cur.df$val.lbl <- val
    cur.df$val <- cur.df[, val]
    if (val == 'f_lt_1') {
      cur.df$val <- cur.df$val - 0.6
    } else if (val == 'f_gt_1') {
      cur.df$val <- cur.df$val*3 + 0.3
    }
    comb.emp <- rbind(comb.emp, cur.df)
  }

  comb.emp$val.lbl <- factor(comb.emp$val.lbl, levels=plot.vals, labels=plot.lbls)
  comb.df$val.lbl <- factor(comb.df$val.lbl, levels=plot.vals, labels=plot.lbls)
  comb.df$sim.lbl <- paste(comb.df$distr_shape, comb.df$neutral_thresh)

  # Turn the continuous dnds into a factor for better control of the color gradient.
  comb.df$dnds <- factor(comb.df$dnds)
  clr.f <- colorRampPalette(c('blue', 'red'), bias=1)
  clrs <- clr.f(length(unique(comb.df$dnds)))

  sub.df <- subset(comb.df, neutral_thresh == 1 & distr_shape %in% c('constant', 'gamma_2'))

  p <- ggplot(sub.df, aes(x=bl, y=val, group=dnds))
  p <- p + theme_bw()
  p <- p + geom_line(alpha=0.25, aes(colour=dnds))
  p <- p + scale_colour_manual("Mean dN/dS", values=clrs)
  p <- p + geom_text(data=comb.emp, aes(x=bl*1.05, label=pset), size=2.5, hjust=0, vjust=0.5)
  p <- p + geom_point(data=comb.emp, size=0.5)
  p <- p + scale_y_log10("Measured Value")
  p <- p + scale_x_log10("Branch Length")
  p <- p + coord_cartesian(xlim=c(0.5, 10), ylim=c(0.05, 2))
  p <- p + facet_grid(sim.lbl ~ val.lbl)

  pdf(file=scratch.f("poisson_sims.pdf"), width=10, height=5)
  print.ggplot(p)
  dev.off()

  sub.df <- subset(comb.df, neutral_thresh == 1)

  p <- ggplot(sub.df, aes(x=bl, y=val, group=dnds))
  p <- p + theme_bw()
  p <- p + geom_line(alpha=0.25, aes(colour=dnds))
  p <- p + scale_colour_manual("Mean dN/dS", values=clrs)
  p <- p + geom_point(data=comb.emp, size=1)
  p <- p + scale_y_continuous("Measured Value")
  p <- p + scale_x_continuous("Branch Length")
  p <- p + facet_grid(distr_shape ~ val.lbl)
  pdf(file=scratch.f("poisson_nolog.pdf"), width=9, height=5)
  print.ggplot(p)
  dev.off()

}

plot.gamma.dnds <- function() {
  n <- 20000
  out.df <- data.frame()
  for (shape in c('gamma 2', 'gamma 0.5')) {
    for (dnds in seq(from=0.1, to=0.7, by=0.1)) {
      if (shape == 'gamma 2') {
        gamma.shape <- 2
      } else {
        gamma.shape <- 0.5
      }
      gamma.scale <- dnds / gamma.shape
      dnds.dist <- rgamma(n, shape=gamma.shape, scale=gamma.scale)
      out.df <- rbind(out.df, data.frame(
        shape = shape,
        dnds = dnds,
        val = dnds.dist       
      ))
    }
  }
  # Turn the continuous dnds into a factor for better control of the color gradient.
  out.df$dnds <- factor(out.df$dnds)
  clr.f <- colorRampPalette(c('blue', 'red'), bias=1)
  clrs <- clr.f(length(unique(out.df$dnds)))

  p <- ggplot(out.df, aes(x=val, colour=dnds, group=dnds))
  p <- p + theme_bw()
  p <- p + geom_line(stat='density', alpha=0.7)
  p <- p + scale_colour_manual("Mean dN/dS", values=clrs)
  p <- p + scale_x_continuous("dN/dS", limits=c(-0.1, 3))
  p <- p + scale_y_continuous("Density")
  p <- p + coord_cartesian(xlim=c(0, 3))
  p <- p + facet_grid(shape ~ .)
  pdf(file=scratch.f("poisson_dnds.pdf"), width=5, height=4)
  print.ggplot(p)
  dev.off()
}


bl.breakdown.plot <- function() {
  con <- connect(db()); 
  df <- dbGetQuery(con, 'select * from summaries where filter="default" order by pset')
  dbDisconnect(con)

  df$pop_size <- as.integer(round(c(
    20000, # Primates
    230000, # Glires
    34410, # Laurasiatheria
    30000, # Atlantogenata
    110000, # Eutheria
    120000, # Mammalia
    230000, # Sparse G
    120000 # Sparse M
  )))

  plt.vals <- c('neg_f05', 'neg_f10', 'f_less_1', 'pos_f05', 'pos_f10', 'f_gt_1')
  out.df <- reshape(df,
    varying=plt.vals, 
    times=plt.vals,
    v.names='val', 
    timevar='lbl',
    idvar='pset',
    direction='l'
  )

  out.df$lbl <- factor(out.df$lbl, levels=plt.vals,
    labels = c('Purifying (5% FPR)', 'Purifying (10% FPR)', 'Omega < 1',
      'Positive (5% FPR)', 'Positive (10% FPR)', 'Omega > 1')
  )

  out.df$type <- 1
  out.df[as.integer(out.df$lbl) >= 4, 'type'] <- 2

  clrs <- c(
    rgb(0.3, 0.3, 1),
    rgb(0.5, 0.5, 1),
    rgb(0.7, 0.7, 1),
    rgb(1, 0.3, 0.3),
    rgb(1, 0.5, 0.5),
    rgb(1, 0.7, 0.7)
  )

  p <- ggplot(out.df, aes(x=med_bl, y=val, colour=lbl))
  p <- p + theme_bw()
  p <- p + geom_point()
  p <- p + geom_line()

  p <- p + scale_colour_manual(values=clrs)
  p <- p + scale_x_continuous("Median BL")
  p <- p + scale_y_continuous("Fraction of sites")
  p <- p + facet_grid(type ~ ., scale="free_y")

  pdf(file='test.pdf', width=10, height=4)
  vplayout(2, 1)
  print.ggplot(p, vp=subplot(1,1))

  p <- ggplot(out.df, aes(x=pop_size, y=val, colour=lbl))
  p <- p + theme_bw()
  p <- p + geom_point()
  p <- p + geom_line()

  p <- p + scale_colour_manual(values=clrs)
  p <- p + scale_x_log10("Ne")
  p <- p + scale_y_continuous("Fraction of sites")
  p <- p + facet_grid(type ~ ., scale="free_y")

  print.ggplot(p, vp=subplot(2,1))
  dev.off()


}

sites.bl.effect.correlations <- function() {
  con <- connect(db()); 
  df <- dbGetQuery(con, 'select * from summaries where filter="default" order by pset')
  dbDisconnect(con)

  # Data from http://www.plosone.org/article/info:doi/10.1371/journal.pone.0004396
  # Gray whale estimate from http://www.nature.com/nrg/journal/v10/n3/fig_tab/nrg2526_T1.html
  # Low elephant pop size from http://www.nature.com/hdy/journal/v84/n3/full/6886740a.html
  # Evidence for larger ancestral mammalian pop. size: http://www.indiana.edu/~lynchlab/PDF/Lynch173.pdf
  df$pop_size <- as.integer(round(c(
    20000, # Primates
    230000, # Glires
    34410, # Laurasiatheria
    30000, # Atlantogenata
    110000, # Eutheria
    120000, # Mammalia
    230000, # Sparse G
    120000 # Sparse M
  )))
   
  cor.f <- function(df, fld_a) {
    df$cur_val <- df[, fld_a]
    m <- 's'
    bl.cor <- cor.test(df$cur_val, df$med_bl, method=m)
    ps.cor <- cor.test(df$cur_val, df$pop_size, method=m)
    omg.cor <- cor.test(df$cur_val, df$mean_omega, method=m)
    syn.cor <- cor.test(df$cur_val, df$f_syn, method=m)
    nsyn.cor <- cor.test(df$cur_val, 1 - (df$f_syn+df$f_const), method=m)
    const.cor <- cor.test(df$cur_val, df$f_const, method=m)
    blw.cor <- cor.test(df$cur_val, df$f_less_1, method=m)
    blw5.cor <- cor.test(df$cur_val, df$f_less_0_5, method=m)
    neg.cor <- cor.test(df$cur_val, df$neg_f10, method=m)
    pos.cor <- cor.test(df$cur_val, df$pos_f10, method=m)
    neut.cor <- cor.test(df$cur_val, df$neutral_f10, method=m)
    cor.l <- list(
      'Branch Length' = bl.cor,
      '\\Ne' = ps.cor,
      'Mean \\omgml' = omg.cor,
      'Constant' = const.cor,
      'Synonymous' = syn.cor,
      'Nsynonoymous' = nsyn.cor,
      '\\omgml $<$1' = blw.cor,
      '\\omgml $<$0.5' = blw5.cor,
      '\\psfive' = pos.cor,
      '\\nsfive' = neg.cor
    )    

    out.df <- data.frame()
    for (i in c(1:length(cor.l))) {
      cur.cor <- cor.l[[i]]
      cur.lbl <- names(cor.l)[i]
      cur.df <- data.frame(
        b = cur.lbl,
        'Correlation' = cur.cor$estimate,
        'Pval' = cur.cor$p.value
      )
      out.df <- rbind(out.df, cur.df)
    }
    out.df
  }

  df.bl <- cor.f(df, 'med_bl')
  df.pop <- cor.f(df, 'pop_size')
  df.pop$b <- NULL
  colnames(df.pop) <- c('Correlation 2', 'Pval 2')

  out.df <- cbind(df.bl, df.pop)
  print(out.df)

  xt <- xtable(out.df)
  xt <- color.columns(xt, c('Pval', 'Pval 2'), low=rgb(0.3, 0.3, 1), high='white')
  print.latex(xt, file=scratch.f("summary.correlations.txt"))  
}

subplot <- function(x, y) viewport(layout.pos.col=x, layout.pos.row=y)
vplayout <- function(x, y) {
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(y,x)))
}


pp <- function(...) {
  print(sprintf(...))
}