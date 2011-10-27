library(R.oo)
library(ape)
library(plyr)
uname  <- Sys.getenv("USER")
if (uname == 'gj1') {
  source("~/src/greg-ensembl/scripts/mysql_functions.R")
  source("~/src/greg-ensembl/projects/2xmammals/analyze_filters.R")
  source("~/src/greg-ensembl/projects/orthologs/analyze_orthologs.R")
  source("~/src/greg-ensembl/scripts/xtable_utils.R")
} else {
  source("~/lib/greg-ensembl/scripts/mysql_functions.R")
  source("~/lib/greg-ensembl/projects/2xmammals/analyze_filters.R")
  source("~/lib/greg-ensembl/projects/orthologs/analyze_orthologs.R")
  source("~/lib/greg-ensembl/scripts/xtable_utils.R")
}
library(ggplot2)
library(boot)
library(ppcor)
library(GenomicRanges)

db <- function() {
   'gj1_2x_63_alt'
}

ifebi <- function(a, b) {
  ifelse(Sys.getenv("USER") == 'greg', return(a), return(b))
}

is.ebi <- function() {
  Sys.getenv("USER") == 'greg'
}

scratch.f <- function(f) {
  paste(scratch(), f, sep='')
}

project.f <- function(f) {
  ifebi(
    paste("~/src/greg-ensembl/projects/2xmammals/", f, sep=''),
    paste("~/lib/greg-ensembl/projects/2xmammals/", f, sep='')
  )
}

scratch <- function() {
  "~/scratch/gj1_2x_63_alt/current/data/"
}

get.genes <- function() {
  if (!exists('genes.df', envir=.GlobalEnv)) {
    genes.f <- scratch.f("genes_orig.Rdata")
    if (file.exists(genes.f)) {
      load(genes.f)
      df <- genes
    } else {
      con <- connect(db()); 
      df <- dbGetQuery(con, 'select * from genes')
      dbDisconnect(con)
    }
    assign('genes.df', df, envir=.GlobalEnv)
  }
  df <- get('genes.df', envir=.GlobalEnv)
  df
}

get.genes.split <- function() {
  split.genes.f <- scratch.f("genes.split.Rdata")
  if (!file.exists(split.genes.f)) {
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
    #assign('genes.split.df', out.df, envir=.GlobalEnv)
    genes.split.df <- out.df
    save(genes.split.df, file=split.genes.f)
  }
  load(split.genes.f)
  genes.split.df
}

get.pset.sites <- function(pset, filter='default', test=F, with.primate.subs=F) {
  if (with.primate.subs) {
    orig.filter <- filter
    filter <- paste(filter, 'primate_subs', sep='_')
  }
  sites.f <- get.sites.filename(pset, filter, test)
  if (file.exists(sites.f)) {
    load(sites.f)
  } else {
    if (with.primate.subs) {
      print("Primate Subs file doesn't exist... collecting subs and saving!")
      orig.sites.f <- get.sites.filename(pset, orig.filter, test=test)
      if (!file.exists(orig.sites.f)) {
        print(pset)
        print(filter)
        print(sites.f)
        print(orig.sites.f)
        stop("Original sites file for primates subs does not exist!")
      }
      load(orig.sites.f)
      sites <- add.primate.subs(sites, test=test)
      save("sites", file=sites.f)
    } else {
      print(pset)
      print(filter)
      print(sites.f)
      stop("Sites file does not exist!")
    }
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
  cmd <- ifebi(
    sprintf('bsub -q research -M%d000 %s -o "/homes/greg/scratch/lsf_logs/%s_%%J_%%I.txt" "R-2.12.0 --vanilla --args %s < sites_scripts.R"',
      mem, array_s, fn_name, args_s
    ),
    sprintf('bsub -q %s -R "select[mem>%d000] rusage[mem=%d000]" -M%d000000 %s "/software/R-2.13.0/bin/R --vanilla --args %s < sites_scripts.R"',
      queue, mem, mem, mem,
      array_s, args_s
    )
  )
  print(cmd)
  system(cmd)
}

bsub.plot.pvals <- function(...) {
   bsub.function('plot_pvals', mem=6, extra.args='1')
   bsub.function('plot_pvals', mem=6, extra.args='2')
   bsub.function('plot_pvals', mem=6, extra.args='6')
}

bsub.collect.sites <- function(...) {
#  bsub.pset.function('collect_sites', queue='normal', mem=18, ...)
   bsub.function('collect_sites', queue='normal', mem=12, extra.args='2 default')
   bsub.function('collect_sites', queue='normal', mem=12, extra.args='2 default')
}

bsub.collect.sites.filtered <- function(...) {
  bsub.function('collect_sites', mem=18, extra.args='1 clusters_inverse')
#  bsub.function('collect_sites', mem=18, extra.args='1 dups_inverse')
#  bsub.function('collect_sites', mem=18, extra.args='1 stringent')

#  bsub.function('collect_sites', mem=18, extra.args='6 clusters_inverse')
#  bsub.function('collect_sites', mem=18, extra.args='6 dups_inverse')
#  bsub.function('collect_sites', mem=18, extra.args='6 stringent')

  for (filter in c('default', 'stringent', 'pfam', 'pfam_stringent', 'clusters_inverse', 'dups_inverse')) {
#    bsub.pset.function('collect_sites', extra.args=filter, mem=18, ...)
     
  }
}

bsub.collect.sites.recomb <- function(...) {
  bsub.pset.function('collect_sites', extra.args='recomb_10k', queue='hugemem', mem=24, ...)
  bsub.pset.function('collect_sites', extra.args='recomb_100k', queue='hugemem', mem=24, ...)
  #bsub.pset.function('collect_sites', extra.args='recomb_1mb', queue='hugemem', mem=24, ...)
}

bsub.collect.gc <- function(...) {
  fn <- 'collect_gc'
  bsub.function(fn, mem=18, extra.args='10000', queue='normal')
  bsub.function(fn, mem=18, extra.args='100000', queue='normal')
  bsub.function(fn, mem=18, extra.args='1000000', queue='normal')
}

bsub.recomb.calc <- function(...) {
  fn <- 'recomb_calc'
  mem <- 18
  queue <- 'hugemem'
  for (filter in c('recomb_10k', 'recomb_100k', 'recomb_1mb')) {
    bsub.function(fn, mem=mem, extra.args=sprintf("1 %s F",filter), queue=queue)
    bsub.function(fn, mem=mem, extra.args=sprintf("2 %s F",filter), queue=queue)
    bsub.function(fn, mem=mem, extra.args=sprintf("3 %s F",filter), queue=queue)
    bsub.function(fn, mem=mem, extra.args=sprintf("6 %s F",filter), queue=queue)
  }
}

bsub.collect.genes <- function(...) {
  n.parallel.jobs <- 100

  df <- pset.df()
  psets <- df$pset_id
  for (pset in c(7, 8, 9, 10)) {
    for (filter in c('default', 'stringent', 'pfam')) {
#  for (pset in 1:length(psets)) {
#    for (filter in c('default', 'stringent', 'pfam')) {
      jobarray_id <- paste('genes', pset, filter, sep='_')
      bsub.function('collect_genes', mem=5,
        extra.args=paste(pset, filter, n.parallel.jobs, 'FALSE', sep=' '),
        jobarray=n.parallel.jobs,
        jobarray_id=jobarray_id
      )
    }
  }
}

bsub.summaries <- function(...) {
  filters <- c('none', 'default', 'stringent', 'pfam', 'pfam_stringent', 'clusters_inverse', 'dups_inverse')
  for (i in 1:length(filters)) {
    bsub.pset.function('summary_table', mem=12, extra.args=filters[i], ...)
  }
}

bsub.plots <- function() {
#  plots <- c('plot_global_distribution')
#  for (i in 1:length(plots)) {
#    bsub.pset.function(plots[i], queue='normal', mem=18)
#  }

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
#  for (i in c(1)) {
#    for (j in c('default')) {
#      for (k in c('lnorm')) {
#        for (l in c('ci')) {

  con <- connect(db())
  for (pset in 1:length(psets)) {
    for (filter in c('default', 'stringent', 'pfam')) {
      for (dist in c('lnorm', 'gamma', 'beta', 'exp', 'weibull')) {
        for (use_type in c('ci', 'omega', 'ci_noconstant', 'omega_noconstant')) {
          df <- dbGetQuery(con, sprintf("select * from fitdistr where 
            pset=%s and filter='%s' and dist='%s' and use_type='%s' and i=50", pset, filter, dist, use_type
          ))
          if (nrow(df) == 0) {
            bsub.function('fit_distr', mem=6, extra.args=paste(pset, filter, dist, use_type, 'FALSE', sep=' '))
          }
        }
      }
    }
  }
  dbDisconnect(con)
}

bsub.corrs <- function() {
  df <- pset.df()
  psets <- df$pset_id

  con <- connect(db())
  for (filter in c('default', 'stringent')) {
    for (i in 1:length(psets)) {
      for (j in 1:length(psets)) {
        if (i == j) {
          next()  
        }

#        df <- dbGetQuery(con, sprintf("select * from correlations where filter='%s' and pset_a=%s and pset_b=%s", filter, i, j))
#        if (nrow(df) == 0) {
          print(paste(filter, i, j, nrow(df)))
          bsub.function('pset_correlation', extra.args=paste(filter, i, j, 'F', sep=' '), mem=18)
#        }
      }
    }
  }
  dbDisconnect(con)
}

bsub.parallel <- function() {
  df <- pset.df()
  psets <- df$pset_id

  con <- connect(db())
  for (filter in c('default', 'stringent')) {
    for (i in 1:length(psets)) {
      for (j in 1:length(psets)) {
        if (i == j) {
          next()  
        }
        df <- dbGetQuery(con, sprintf("select * from parallel_sites where filter='%s' and pset_a=%s and pset_b=%s", filter, i, j))
        if (nrow(df) == 0) {
          bsub.function('parallel_sites', extra.args=paste(filter, i, j, 'F', sep=' '), mem=12)
        }                                         
      }
    }
  }
  dbDisconnect(con)
}

plot.pval.example <- function(pset=1, test=T) {
  sites <- get.pset.sites(pset, filter='default', test=test)
  sites <- subset(sites, select=c('omega', 'pval', 'pos.pval'))

  n.cols <- 4
  bw <- 0.02

  out.f <- scratch.f(sprintf("pval_example_%d.pdf", pset))
  pdf(file=out.f, width=3 * n.cols, height=2.5)
  vplayout(n.cols, 1)

  p <- ggplot(sites, aes(x=pval, y=..count..+10))
  p <- p + geom_histogram(binwidth=bw)
#  p <- p + scale_y_log10()
  p <- p + scale_x_continuous("Two-tailed")
  p <- generic.opts(p)
  p <- p + opts(
    axis.title.y = theme_blank()  
  )
  print.ggplot(p, vp=subplot(1, 1))

  sub.neg <- subset(sites, omega < 1)
  p <- ggplot(sub.neg, aes(x=pval, y=..count..+1))
  p <- p + geom_histogram(binwidth=bw)
#  p <- p + scale_y_log10()
  p <- p + scale_x_continuous("Two-tailed (w < 1)")
  p <- generic.opts(p)
  p <- p + opts(
    axis.title.y = theme_blank()
  )
  print.ggplot(p, vp=subplot(2, 1))

  sub.pos <- subset(sites, omega > 1)
  p <- ggplot(sub.pos, aes(x=pval, y=..count..+1))
  p <- p + geom_histogram(binwidth=bw)
#  p <- p + scale_y_log10()
  p <- p + scale_x_continuous("Two-tailed (w > 1)")
  p <- generic.opts(p)
  p <- p + opts(
    axis.title.y = theme_blank()
  )
  print.ggplot(p, vp=subplot(3, 1))

  n.at.first <- nrow(subset(sites, pos.pval <= bw))
  print(n.at.first)

  p <- ggplot(sites, aes_string(x="pos.pval", y=paste("pmin(", n.at.first, ", ..count..+1)", sep='')))
  p <- p + geom_histogram(binwidth=bw)
#  p <- p + scale_y_log10()
  p <- p + scale_x_continuous("One-tailed")
  p <- generic.opts(p)
  p <- p + opts(
    axis.title.y = theme_blank()  
  )
  print.ggplot(p, vp=subplot(4, 1))

  dev.off()
}

write.fits.tables <- function() {
  write.fits.table('ci')
  write.fits.table('omega')
}

write.fits.table <- function(use_t='ci') {
  # One row per pset, showing mean AIC / mean / F > 1 for each model fit
    con <- connect(db())
    df <- dbGetQuery(con, 'select * from fitdistr')
    dbDisconnect(con)
    #load(scratch.f("fitdistr.Rdata"))

    df.sub <- subset(df, filter == 'stringent')

    xt.df <- ddply(df.sub, .(pset, use_type), function(x) {
      cur.df <- data.frame(
        pset = x[1, 'pset'],
        filter = x[1, 'filter'],
        use_type = x[1, 'use_type']
      )
      for (distr_type in c('lnorm', 'gamma', 'exp', 'beta', 'weibull')) {
        distr.row <- subset(x, dist == distr_type)
        if (nrow(distr.row) > 0) {
          distr.df <- data.frame(
            mean = median(distr.row$mean),
            f_above_1 = median(distr.row$f_above_1) * 100
          )
          colnames(distr.df) <- paste(distr_type, colnames(distr.df), sep='_')
          cur.df <- cbind(cur.df, distr.df)
        }
      }
      cur.df
    })

    xt.df <- cbind(data.frame(
      flt.fact=filter.to.alias(xt.df$filter, factors=T),
      use.type.fact=use.type.to.alias(xt.df$use_type, factors=T),
      pset.fact=pset.to.alias(xt.df$pset, factors=T)
      ), xt.df)

    xt.df <- xt.df[order(xt.df$use.type.fact, xt.df$pset.fact),]
   
    xt.df[duplicated(xt.df$use.type.fact), 'use.type.fact'] <- ''

    rcs <- c('pset', 'filter', 'use_type')

    df1 <- subset(xt.df, use_type == 'omega')
    df2 <- subset(xt.df, use_type == 'ci')
    df1[duplicated(df1$flt.fact), 'flt.fact'] <- ''
    df2[duplicated(df2$flt.fact), 'flt.fact'] <- ''

    xt1 <- xtable(remove.cols(df1, rcs))
    xt2 <- xtable(remove.cols(df2, rcs))

    clr.cols <- colnames(xt1)
    clr.cols <- setdiff(clr.cols, c('pset.fact', 'flt.fact', 'use.type.fact'))

    xt1 <- color.columns(xt1, clr.cols, log=T)
    xt2 <- color.columns(xt2, clr.cols, log=T)

    ff <- scratch.f("fitdistr_fits_all.txt")
    print.latex(xt1, ff)
    cat("\\midrule", file=ff, append=T)
    print.latex(xt2, ff, append=T)

    # Now sort by the best-fit model type for each dataset (only using CIs) and write a table
    # with the AIC and best fit parameters
    df.sub <- subset(df, use_type == use_t)
    
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

    #print(mean.df)
    xt <- xtable(mean.df)
    xt <- color.columns(xt, c('aic_diff'))

    ff <- scratch.f(sprintf("fitdistr_params_%s.txt", use_t))
    print.latex(xt, file=ff)

    # Finally, plot a figure for some species groups showing the bootstrap distributions of mean omega and % > 1.
    con <- connect(db())
    df <- dbGetQuery(con, 'select * from fitdistr')
    dbDisconnect(con)

    print("Plotting...")

    df$dist.fact <- dist.to.factor(df$dist)
    df$filter.fact <- filter.to.alias(df$filter, factors=T)

    df.sub <- subset(df, use_type == use_t & filter == 'stringent')

    df.sub$pset.fact <- pset.to.alias(df.sub$pset, factors=T)
    df.sub$pset.filter <- paste(" ", df.sub$pset.fact, sep=' ')

    df.sub <- df.sub[order(df.sub$pset.fact, df.sub$filter.fact), ]

    df.sub$pset.filter <- factor(df.sub$pset.filter, levels=unique(df.sub$pset.filter), labels=unique(df.sub$pset.filter))

    # For each parameter set, find the mean-across-genes dN/dS and plot that too.
    # [ TODO ]

    p <- ggplot(df.sub, aes(x=mean, fill=dist.fact))
    p <- p + theme_bw()
    p <- p + geom_histogram(binwidth=0.002)
    p <- p + scale_fill_discrete("Distribution")
    p <- p + scale_x_continuous("Mean Omega", limits=c(0, 0.5))
    p <- p + scale_y_continuous("Filter", breaks=c(0))
    p <- p + facet_grid(pset.filter ~ ., scales="free")
    p <- p + opts(
      axis.ticks.y = theme_blank(),
      axis.text.y = theme_blank(),
      axis.text.x = theme_text(angle=90, hjust=1),
      strip.text.y = theme_text(angle=0, hjust=0),
      strip.background = theme_blank()
    )
    
    ff <- scratch.f(sprintf("fitdistr_means_%s.pdf", use_t))
    pdf(file=ff, width=8, height=8)
    print.ggplot(p)
    dev.off()

    df.sub <- subset(df.sub, dist != 'beta')

    p <- ggplot(df.sub, aes(x=f_above_1 * 100, fill=dist))
    p <- p + theme_bw()
    p <- p + geom_histogram(binwidth=0.05)
    p <- p + scale_fill_discrete("Distribution")
    p <- p + scale_x_continuous("% Omega > 1", limits=c(0 * 100, 0.07 * 100))
    p <- p + scale_y_continuous("Filter", breaks=c(0))
    p <- p + facet_grid(pset.filter ~ ., scales="free")
    p <- p + opts(
      axis.ticks.y = theme_blank(),
      axis.text.y = theme_blank(),
      axis.text.x = theme_text(angle=90, hjust=1),
      strip.text.y = theme_text(angle=0, hjust=0),
      strip.background = theme_blank()
    )

    ff <- scratch.f(sprintf("fitdistr_above_one_%s.pdf", use_t))
    pdf(file=ff, width=8, height=8)
    print.ggplot(p)
    dev.off()
}

process.genes <- function(genes, pset, filter='default', subset.index=NULL, n.indices=0, test=F) {
  # Get the sites corresponding to this pset.
  print("Getting sites...")
  sites <- get.pset.sites(pset, filter=filter, test=test)

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

  # Calculate quantiles and z-scores for dnds and branch length values.
  dnds.q <- cum.z(genes[, dnds.s])
  bl.q <- cum.z(genes[, bl.s])
  # Try a log-transformed dN/dS z-score
  log.dnds <- log(pmax(genes[, dnds.s], 0.001))
  log.mean <- mean(log.dnds, na.rm=T)
  log.sd <- sd(log.dnds, na.rm=T)
  dnds.z <- (log.dnds - log.mean) / log.sd

  log.bl <- log(pmax(genes[, bl.s], 0.001))
  log.mean <- mean(log.bl, na.rm=T)
  log.sd <- sd(log.bl, na.rm=T)
  bl.z <- (log.bl - log.mean) / log.sd

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
  gene.cols$slr_dnds_q <- dnds.q
  gene.cols$slr_bl_q <- bl.q
  gene.cols$slr_dnds_z <- dnds.z
  gene.cols$slr_bl_z <- bl.z
  genes <- gene.cols

  genes$parameter_set_id <- pset
  genes$pset_char <- pset.char

  # Take a subset of data_ids to process gene-wise
  process.all <- FALSE
  if (!is.null(subset.index) && !process.all) {
    all.ids <- sort(unique(sites$data_id))
    print(sprintf("%d total genes", length(all.ids)))
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
  } else {
    # Shuffle the order of the data_ids
    all.ids <- sort(unique(sites$data_id))
    all.ids <- sample(all.ids,length(all.ids))
    keep.ids <- all.ids[1:floor(length(all.ids)/5)]
    sites <- subset(sites, data_id %in% keep.ids)
  }

  print(sprintf("%d genes to process", length(unique(sites$data_id))))

  con <- connect(db())
  df <- dbGetQuery(con, 
    sprintf("select * from genes_sites where filter='%s' and pset=%s",
      filter, pset)
  )
  disconnect(con)
  existing.data.ids <- df$data_id

  # ddply each gene, extract the pset SLR values and analyze the sites.
  d_ply(sites, .(data_id), function(x) {
    cur.gene <- subset(genes, data_id == x[1, 'data_id'])
    if (nrow(cur.gene) == 0) {
      print(paste("Gene not found for data_id", x[1, 'data_id']))
      return()
    }

    already.exists <- cur.gene$data_id %in% existing.data.ids
    if (!already.exists) {
      sites.cols <- process.gene.sites(cur.gene, x)
      cur.df <- cbind(cur.gene, sites.cols)

      cur.df$filter <- filter
      cur.df$label <- paste(cur.df$filter, cur.df$parameter_set_id, cur.df$data_id, collapse=' ')
      cur.df$pset <- cur.df$parameter_set_id  
      cur.df$subset.index <- subset.index

      con <- connect(db())
      if (!test) {
        write.or.update(cur.df, 'genes_sites', con, 'label')

        #existing.df <- dbGetQuery(con, 
        #  sprintf("select * from genes_sites where filter='%s' and pset=%s",
        #  filter, pset)
        #)
        #existing.data.ids <- existing.df$data_id
      }
      disconnect(con)
    } else {
      print(sprintf("Gene %s already exists in table!", cur.gene$gene_name))
    }
  })

}

# Given a row corresponding to a gene, summarize the sitewise data w.r.t. pos and neg selection.
process.gene.sites <- function(gene, sites) {
  print(sprintf("%d sites for gene %s", nrow(sites), gene$gene_name))
  
  if (is.ebi()) {
    source("~/lib/greg-ensembl/projects/2xmammals/Tpmw.r")
  } else {
    source("~/src/greg-ensembl/projects/2xmammals/Tpmw.r")
  }

  # Get GC and recombination rate at various resolutions.
  chr_name <- gene$chr_name
  chr_start <- gene$chr_start 
  chr_end <- gene$chr_end
  recombM_10k <- NA_real_
  recombM_100k <- NA_real_
  recombM_1mb <- NA_real_
  recombF_10k <- NA_real_
  recombF_100k <- NA_real_
  recombF_1mb <- NA_real_
  gc_10k <- NA_real_
  gc_100k <- NA_real_
  gc_1mb <- NA_real_
  if (!is.na(chr_start) && length(grep("ENSP0.*", gene$ref_protein_id)) > 0) {
    cur.range <- GRanges(
      seqnames = chr_name,
      ranges = IRanges(chr_start, chr_end),
      strand = '*'
    )
    for (width in c(1e4, 1e5, 1e6)) {
      recomb.ranges <- get.recomb.ranges(width=width)
      recomb.sub <- subsetByOverlaps(recomb.ranges, cur.range)
      gc.ranges <- get.gc.ranges(width=width)
      gc.sub <- subsetByOverlaps(gc.ranges, cur.range)
      mean.m <- NA_real_
      mean.f <- NA_real_
      mean.gc <- NA_real_
      elmeta <- elementMetadata(recomb.sub)
      if (length(recomb.sub) > 0) {
        print(elementMetadata(recomb.sub)$recombM)
        mean.m <- mean(elementMetadata(recomb.sub)$recombM)
        mean.f <- mean(elementMetadata(recomb.sub)$recombF)
      }
      if (length(gc.sub) > 0) {
        mean.gc <- mean(elementMetadata(gc.sub)$gc)
      }
      if (width == 1e4) {
        recombM_10k <- mean.m
        recombF_10k <- mean.f
        gc_10k <- mean.gc
      } else if (width == 1e5) {
        recombM_100k <- mean.m
        recombF_100k <- mean.f
        gc_100k <- mean.gc
      } else if (width == 1e6) {
        recombM_1mb <- mean.m
        recombF_1mb <- mean.f
        gc_1mb <- mean.gc
      }
    }
  }
  print(gene)

  # omg_ml > 1 and < 1
  n.abv.1 <- nrow(subset(sites, omega > 1))
  n.blw.1 <- nrow(subset(sites, omega < 1))

  n.pos.50 <- nrow(subset(sites, omega > 1 & pval < 0.5))
  n.neg.50 <- nrow(subset(sites, omega < 1 & pval < 0.5))
  n.ntr.50 <- nrow(subset(sites, pval > 0.5))

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

  n.pos.fwer <- 0
  n.pos.gene.bh.10 <- 0
  n.pos.gene.bh.05 <- 0
  
  if (nrow(sites > 1)) {
    sites$hoch.pval <- p.adjust(sites$pval, method='hochberg')
    n.pos.fwer <- nrow(subset(sites, hoch.pval < 0.05 & omega > 1))

    # Use the one-tailed p-values to do BH adjustment.
    sites$bh.pval <- p.adjust(sites$pos.pval, method='BH')
    n.pos.gene.bh.10 <- nrow(subset(sites, bh.pval < 0.10))
    n.pos.gene.bh.05 <- nrow(subset(sites, bh.pval < 0.05))
  }

  omgs <- pmin(sites$omega, 5)
  mean.omg <- mean(omgs)
  rm(omgs)

  geom.mean.omg <- exp(mean(log(pmax(sites$omega, 0.01))))

  sorted.omgs <- sort(sites$omega)
  n.quart <- nrow(sites) / 10
  mid.omgs <- sorted.omgs[floor(.5*n.quart):floor(9.5*n.quart)]
  mean.mid.omega <- mean(mid.omgs)
  rm(sorted.omgs)
  rm(mid.omgs)

  # Get some estimate of the lrt_stat autocorrelation using a Mantel test
  #print("  mantel test")
  library(ade4)
  aln.dists <- dist(sites$aln_position)
  lrt.dists <- dist(sites$lrt_stat)
  #mantel.res <- mantel.test(aln.dists, lrt.dists, nperm=100)
  #mantel.res <- mantel.rtest(aln.dists, lrt.dists, nrepet = 999)
  #alt.mantel <- mantel(lrt.dists ~ aln.dists, nperm=1000)
  #print(str(alt.mantel))
  # Use the fast mantel test implemented in C.
  mantel.cor <- 0
  mantel.p <- 1
  if (nrow(sites) > 10) {
    mantel.res <- mantel.randtest(aln.dists, lrt.dists, nrepet = 500)
    mantel.cor <- mantel.res$obs
    mantel.p <- mantel.res$pvalue
  }

  # Also calculate Moran's I
  #print("  moran's I")
  moran.dir <- 0
  moran.p <- 1
  if (nrow(sites) > 10) {
    aln.dists <- as.matrix(dist(sites$aln_position))
    aln.dists.inv <- 1/aln.dists
    diag(aln.dists.inv) <- 0
    moran.res <- Moran.I(sites$lrt_stat, aln.dists.inv)
    moran.dir <- moran.res$observed - moran.res$expected
    moran.p <- moran.res$p.value
  }

  # Combine p-values using Fisher's method.
  fis.stat <- 2 * sum(-log(sites$pos.pval))
  fis.n <- 2 * nrow(sites)
  fis.p <- pchisq(fis.stat, df=fis.n, lower=F)

  # Combine the top N, where N=ceiling(length/10), p-values using Fisher's method.
  #blw.50 <- subset(sites, pos.pval < 0.5)
  #blw.10 <- subset(sites, pos.pval < 0.1)

  # Use the weighted TPM method from
  # Zaykin et al. 2002
  #print("  tpm p-value")
  weights <- sites$ncod
  pvals <- sites$pos.pval
  loops <- 2500
  tpm.05 <- Tpmw(pvals, weights, 0.05, loops=loops)
  tpm.10 <- Tpmw(pvals, weights, 0.1, loops=loops)
  tpm.20 <- Tpmw(pvals, weights, 0.2, loops=loops)
  tpm.50 <- Tpmw(pvals, weights, 0.5, loops=loops)

  # Fit the distribution of sites.
  #print("  fit distr")
  lnorm.df <- fit.sites(sites, distr='lnorm', use='ci', i=0, write.to.table=F, boot.sample=F)
  gamma.df <- fit.sites(sites, distr='gamma', use='ci', i=0, write.to.table=F, boot.sample=F)
  exp.df <- fit.sites(sites, distr='exp', use='ci', i=0, write.to.table=F, boot.sample=F)

  l.aic <- lnorm.df$aic
  l.mean <- lnorm.df$mean
  l.abv <- lnorm.df$f_above_1
  g.aic <- gamma.df$aic
  g.mean <- gamma.df$mean
  g.abv <- gamma.df$f_above_1
  e.aic <- exp.df$aic
  e.mean <- exp.df$mean
  e.abv <- exp.df$f_above_1

  print("  output")
  n <- nrow(sites)
  out.df <- data.frame(
    n_sites = nrow(sites),

    f_sites = nrow(sites) / gene$aln_length,
    mean_ncod = mean(sites$ncod),

    fit_l_aic = l.aic,
    fit_l_mean = l.mean,
    fit_l_abv = l.abv,
    fit_g_aic = g.aic,
    fit_g_mean = g.mean,
    fit_g_abv = g.abv,
    fit_e_aic = e.aic,
    fit_e_mean = e.mean,
    fit_e_abv = e.abv,

    recombM_10k = recombM_10k,
    recombM_100k = recombM_100k,
    recombM_1mb = recombM_1mb,
    recombF_10k = recombF_10k,
    recombF_100k = recombF_100k,
    recombF_1mb = recombF_1mb,
    gc_10k = gc_10k,
    gc_100k = gc_100k,
    gc_1mb = gc_1mb,

    mean_omega = mean.omg,
    mean_mid_omega = mean.mid.omega,
    mean_omega_geom = geom.mean.omg,

    f_abv_1 = n.abv.1 / n,
    f_blw_1 = n.blw.1 / n,
    
    f_pos_50 = n.pos.50 / n,
    f_neg_50 = n.neg.50 / n,
    f_ntr_50 = n.ntr.50 / n,
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

w.to.s <- function(ws) {
  # Go from dnds to s
  dnds.fn <- function(x) {
    rt <- uniroot(function(y) y / (1 - exp(-y)) - x, lower=-1e50, upper=1e50)
    rt$root
  }

  dnds.x <- seq(from=0, to=5, length.out=101)
  dnds.x <- 1.1^(-100:10)
  s.vals <- c()
  for (i in 2:length(dnds.x)) {
    mid.val <- (dnds.x[i-1] + dnds.x[i]) / 2
    s.vals <- c(s.vals, dnds.fn(mid.val))
  }

  w.f <- cut(ws, dnds.x)
  s.bins <- as.integer(table(w.f))
  ss <- rep(s.vals, times=s.bins)
  ss <- ss + rnorm(n=length(ss), mean=0, sd=mean(abs(ss))/20)
  ss
}

plot.dfe.cum <- function(style='a') {
  small.vals <- c(
    -1e-1, 
    -5e-2, -2e-2, -1e-2,
    -5e-3, -2e-3, -1e-3,
    -5e-4, -2e-4, -1e-4,
    -5e-5, -2e-5, -1e-5,
    -5e-6, -2e-6, -1e-6,
    0,
    1e-6, 2e-6, 5e-6,
    1e-5, 2e-5, 5e-5,
    1e-4, 2e-4, 5e-4
  )

  df.fn <- function(lbl) {
    print(lbl)
    dfe.res <- get.dfe(lbl, plot=F)
    # S = 2*Ne*s
    # s = S / (2*Ne)
    pop.size <- dfe.res$pop.size
    xs <- dfe.res$s / (pop.size*2)
    ecdf.f <- ecdf(xs)

    if (lbl %in% c('boyko08_normal', 'boyko08_lnorm', 'boyko08_gamma', 'eyrewalker06')) {
      ltype <- 'Segregating'
    } else {
      ltype <- 'Fixed'
    }

    if (style == 'a') {
      xpoints <- seq(from=-1e-5, to=1e-5, length.out=200)
      cum.vals <- ecdf.f(xpoints)
      cum.df <- data.frame(xx=xpoints, yy=cum.vals, lbl=lbl, ltype=ltype)
    } else if (style == 'b') {
      xpoints <- seq(from=-0.1, to=1e-4, length.out=200)
      cum.vals <- ecdf.f(xpoints)
      cum.df <- data.frame(xx=xpoints, yy=cum.vals, lbl=lbl, ltype=ltype)
    } else {
      cum.vals <- ecdf.f(small.vals)
      cum.df <- data.frame(xx=as.integer(factor(small.vals)), yy=cum.vals, lbl=lbl, ltype=ltype)
    }
    cum.df
  }

  plot.df <- data.frame()
  plot.lbls <- c('boyko08_lnorm', 'boyko08_gamma', 'eyrewalker06',
    'boyko08_normal', 'nielsen03', 'me_mamms', 'me_primates',
    'me_rodents', 'me_laur')

  for (lbl in plot.lbls) {
    plot.df <- rbind(plot.df, df.fn(lbl))
  }

  plot.df$lbl <- factor(plot.df$lbl, levels=plot.lbls, labels=plot.lbls)

  library(RColorBrewer)
  p <- ggplot(plot.df, aes(x=xx, y=yy, colour=lbl))
  p <- p + theme_bw()
  p <- p + geom_line(size=1)
#  p <- p + scale_colour_discrete("Distribution")
  plot.clrs <- c(brewer.pal(5, "Set1")[1:4], brewer.pal(5, "Set1"))
  p <- p + scale_colour_manual("Distribution", values=plot.clrs)
  p <- p + facet_grid(. ~ ltype)
  lbl <- sprintf("%.0e", small.vals)
  if (style == 'c') {
    p <- p + scale_x_continuous("Selection Coefficient", breaks=1:length(small.vals), labels=sprintf("%.0e", small.vals))
  }
  if (style == 'a') {
    p <- p + scale_y_continuous("Cumulative Probability", limits=c(0.7, 1))
  } else {
    p <- p + scale_y_continuous("Cumulative Probability", limits=c(0, 1))
  }
  p <- p + opts(
    axis.text.x = theme_text(angle=90, hjust=1, size=8),
    legend.position = c(0.8, 0.5)
  )

  pdf(file=scratch.f(sprintf("dfe_cum_%s.pdf", style)), width=11, height=6)
  print.ggplot(p)
  dev.off()
}

cum.dfe.plots <- function() {
  plot.dfe.cum('a')
  plot.dfe.cum('b')
  plot.dfe.cum('c')
}

get.dfes <- function() {
  get.dfe('nielsen03')
  get.dfe('boyko08_normal')
  get.dfe('boyko08_lnorm')
  get.dfe('boyko08_gamma')
  get.dfe('eyrewalker06')
  get.dfe('me_mamms')
  get.dfe('me_mamms_gamma')
  get.dfe('me_primates')
  get.dfe('me_rodents')
  get.dfe('me_laur')
}

get.ne.diff <- function(a, b) {
  dfe.a <- get.dfe(a, plot=F)
  dfe.b <- get.dfe(b, plot=F)

  mean.a <- mean(dfe.a$s)
  mean.b <- mean(dfe.b$s)

  return(mean.a/mean.b)
}

get.dfe <- function(ref='nielsen03', plot=T) {
  n <- 100000

  ss <- switch(ref, 
    nielsen03 = rnorm(n, -1.72*2, 0.72*2),
    boyko08_normal = rnorm(n, -38.5, 28.6),
    boyko08_lnorm = -rlnorm(n, 5.02, 5.94),
    boyko08_gamma = -rgamma(n, shape=0.184, scale=8200),
    boyko08_exp = -rexp(n, 0.0365),
    eyrewalker06 = -rgamma(n, shape=0.23, scale=425/0.23),
    me_mamms_gamma = w.to.s(rgamma(n, 0.75, 4.30)),
    me_mamms = w.to.s(rlnorm(n, -2.51, 1.26)),
    me_primates = w.to.s(rlnorm(n, -1.15, 0.66)),
    me_rodents = w.to.s(rlnorm(n, -1.78, 0.60)),
    me_laur = w.to.s(rlnorm(n, -1.80, 0.77))
  )

  # Default primate pop. size
  pop.size <- 30000
  if (ref == 'me_primates') {
    pop.size <- 70000
  } else if (ref %in% c('me_mamms', 'me_mamms_gamma')) {
    pop.size <- 150000
  } else if (ref == 'me_rodents') {
    pop.size <- 230000
  } else if (ref == 'me_laur') {
    pop.size <- 130000
  } else if (ref == 'nielsen03') {
    pop.size <- 50000
  }
  if (ref %in% c('boyko08_exp', 'boyko08_gamma', 'boyko08_normal', 'boyko08_lnorm')) {
    pop.size <- 25636
    ss <- pmax(-25636, ss)
    ss <- pmin(25636, ss)
  }

  # Go from S to dnds
  w.fn <- function(x) {
    x / (-expm1(-x))
  }
  ws <- w.fn(ss)

  cum.p <- function(xs, is.s=F, is.omega=F) {
    if (is.omega) {
      xs <- pmin(xs, 1.5)
    }
    bw <- diff(range(xs))/200

    xs.tmp <- xs
    if (ref %in% c('boyko08_lnorm', 'boyko08_gamma', 'boyko08_normal', 'eyrewalker06')) {
      xs.tmp <- pmax(xs.tmp, -5000)
    }
    ct <- cut(xs.tmp, seq(from=range(xs)[1], to=range(xs)[2], by=bw), include.lowest=T)
    maxval <- max(table(ct))

    ecdf.f <- ecdf(xs)

    if (is.s) {
      small.vals <- c(
        -1e-1, 
        -5e-2, -2e-2, -1e-2,
        -5e-3, -2e-3, -1e-3,
        -5e-4, -2e-4, -1e-4,
        -5e-5, -2e-5, -1e-5,
        -5e-6, -2e-6, -1e-6,
        0,
        1e-6, 2e-6, 5e-6,
        1e-5, 2e-5, 5e-5,
        1e-4, 2e-4, 5e-4
      )
      cum.vals <- ecdf.f(small.vals * pop.size)
      cum.df <- data.frame(xx=as.integer(factor(small.vals)), yy=cum.vals)
      #print(cum.df)
      p <- ggplot(cum.df, aes(x=xx, y=yy))
      p <- p + theme_bw()
      p <- p + geom_line()
      lbl <- sprintf("%.0e", small.vals)
      p <- p + scale_x_continuous("Selection Coefficient", breaks=1:length(small.vals), labels=sprintf("%.0e", small.vals))
      p <- p + scale_y_continuous("", limits=c(0, 1))
      p <- p + opts(
        axis.text.x = theme_text(angle=90, hjust=1, size=7)
      )
      return(p)
    } else {

      unq.xs <- unique(round_any(xs, bw/5, f=ceiling))
      cumx <- ecdf.f(unq.xs) * maxval
      unq.df <- data.frame(xx=unq.xs, cum=cumx)
      df.x <- data.frame(xx=xs, cum=0)

      if (!is.omega && ref %in% 
        c('boyko08_lnorm', 'boyko08_gamma', 'boyko08_normal', 'eyrewalker06')
      ) {
        #print("Capping xx")
        df.x$xx <- pmax(df.x$xx, -5000)
      }

      p <- ggplot(df.x, aes(x=xx))
      p <- p + theme_bw()
      if (is.omega) {
        p <- p + geom_vline(xintercept=1, colour=gray(0.5), linetype='dashed')
      } else {
        p <- p + geom_vline(xintercept=0, colour=gray(0.5), linetype='dashed')
      }

      p <- p + geom_histogram(binwidth=bw, colour=NA, fill='black')
      p <- p + geom_line(data=unq.df, aes(x=xx, y=cum), colour=rgb(0.2, 0.2, 0.8), size=1)
      p <- p + scale_y_continuous("")
      if (ref %in% c('me_primates', 'me_mamms', 'me_mamms_gamma', 'me_rodents', 'me_laur', 'nielsen03')) {
        if (is.omega) {
          p <- p + scale_x_continuous("dN/dS", limits=c(0, 1.5))
        } else {
          p <- p + scale_x_continuous("Scaled Selection Coefficient", limits=c(-10, 2))
        }
      } else {
        if (is.omega) {
          p <- p + scale_x_continuous("dN/dS", limits=c(-0.1, 1.1))
        } else {
          #print(bw)
          p <- p + scale_x_continuous("Scaled Selection Coefficient", c(-5200, 100))
        }
      }
      p <- p + opts(
        axis.text.x = theme_blank()
      )
      p
    }
  }

  if (plot) {
    out.f <- scratch.f(sprintf("dfe_%s.pdf", ref))
    pdf(file=out.f, width=8, height=2.5)
    vplayout(2, 1)
    print.ggplot(cum.p(ws, is.omega=T), vp=subplot(1, 1))
    print.ggplot(cum.p(ss), vp=subplot(2, 1))
    dev.off()
  } else {
    return(list(
      w=ws,
      s=ss,
      pop.size=pop.size
    ))
  }
}


fit.sites <- function(sites, distr, use, i=NA, boot.sample=T, write.to.table=T, filter='default') {
  library(dfoptim)

  # Use a slightly customized version of fitdistrplus.
  if (is.ebi()) {
    lapply(dir("~/lib/greg-ensembl/projects/2xmammals/fitdistr", full.name=T), source)
  } else {
    lapply(dir("~/src/greg-ensembl/projects/2xmammals/fitdistr", full.name=T), source)
  }

  pset <- sites[1, 'parameter_set_id']

  small.value <- 0.001
  big.value <- 5
  fix.arg <- NULL

  sites <- sites[, c('data_id', 'omega_lower', 'omega_upper', 'omega', 'note')]
  colnames(sites) <- c('data_id', 'left', 'right', 'omega', 'note')

  # Make sure the upper values are within a reasonable range.
  sites$left <- pmin(sites$left, big.value)
  sites$right <- pmin(sites$right, big.value * 2)

  sites$left <- pmax(sites$left, small.value)
  sites$right <- pmax(sites$right, small.value * 2)

  if (boot.sample) {
    print(nrow(sites))
    data.ids <- unique(sites$data_id)
    use.n <- floor(length(data.ids) * 0.5)
    ids <- data.ids[sample(1:length(data.ids), size=use.n, replace=T)]
    sites <- sites[sites$data_id %in% ids,]
    print(nrow(sites))
  }
  sites <- subset(sites, !is.na(left) & !is.na(right))

  # Create a list of starting values depending on the distribution.
  param.min = sqrt(.Machine$double.eps)
  start <- NULL
  lower <- -Inf
  upper <- Inf
  if (distr == 'beta') {
    start <- list(1.7, 2.5)
  }
  if (distr == 'weibull') {
    start <- list(0.5, 0.5)
  }
  if (distr == 'gamma') {
    start <- list(1.5, 5)
  }
  if (distr == 'lnorm') {
    start <- list(-5,4)
  }
  if (distr == 'exp') {
    start <- list(3)
  }

  do.fit <- function() {
    fit <- NA
    if (use=='imputed' || use == 'omega' || use == 'omega_noconstant') {
      # For imputed values or \omgml estimates, use 'fitdist' to fit the distribution.
      if (use == 'imputed') {
        vals <- (sites$omega + sites$right + sites$left) / 3
      } else if (use == 'omega') {
        vals <- sites$omega
      } else if (use == 'omega_noconstant') {
        sites <- filter.constant(sites)
        vals <- sites$omega
      }

      # Cap all omegas to the 'large value'
      vals <- pmax(vals, small.value)
      vals <- pmin(vals, big.value)

      if (distr == 'beta') {
        # Beta must be between 0 and 1, so cap values
        vals <- pmin(vals, 1 - small.value)
      }

      fit <- fitdist(vals, distr, lower=lower, upper=upper, start=start)
      if (!is.na(fit) && is.null(fit$aic)) {
  	fit$aic = AIC(fit)
      }
    } else {
      if (distr == 'beta') {
        sites <- subset(sites, left < 1 - small.value)
        sites$right <- pmin(sites$right, 1 - small.value)
        sites <- subset(sites, right > 0 + small.value)
        sites$left <- pmax(sites$left, 0 + small.value)
      } else if (distr == 'exp') {
      } else if (distr == 'gamma') {
        #sites$left <- sites$left * 10
        #sites$right <- sites$right * 10
      }

      if (use == 'ci_noconstant') {
        sites <- filter.constant(sites)
      }
      
      sites <- subset(sites, select=c('left', 'right'))

      if (length(start) > 1) {
        fit <- fitdistcens(sites, distr, lower=lower, upper=upper, start=start, custom.optim=nmk, fix.arg=fix.arg)
      } else {
        fit <- fitdistcens(sites, distr, lower=lower, upper=upper, start=start, fix.arg=fix.arg)
      }
    }
    fit
  }

  out.df <- data.frame(
    label = paste(pset, distr, filter, use, i),
    pset = pset,
    filter = filter,
    dist = distr,
    use_type = use,
    i = i,
    n_sites = nrow(sites),
    
    error = 1,
    error_str = '',

    aic = -1.00,
    fit_str = '',
    'est_1' = -1.00,
    'est_2' = -1.00,
    mean = -1,
    sd = -1,
    'f_above_1' = 0,
    'f_above_1_5' = 0,
    'f_below_0_5' = 0
  )

  if (nrow(sites) < 10) {
    return(out.df)
  }

  error.f = function(e) {
    print("### Error ###")
    print(e)
    out.df$error <- 1
    out.df$error_str <- as.character(e)
    if (write.to.table) {
      con <- connect(db())
      write.or.update(out.df, 'fitdistr', con, 'label')
      dbDisconnect(con)  
    }
    return(NA)
  }
  #print("  doing fit...")
  fit <- tryCatch(do.fit(), error=error.f)
  
  #print("  done!")
  
  if (any(!is.na(fit))) {
    # Calculate the mean of the parameterized distribution.
    fn.str <- paste('r', distr, sep='')

    if (length(fit$estimate) > 1) {
      sampled.values <- do.call(fn.str, list(10 * 1000 * 1000, fit$estimate[1], fit$estimate[2]))
    } else {
      sampled.values <- do.call(fn.str, list(10 * 1000 * 1000, fit$estimate[1]))    
    }

    # Limit the sampled values to between small.value and big.value
    sampled.values <- pmax(sampled.values, small.value)
    sampled.values <- pmin(sampled.values, big.value)
    
    out.df$mean <- mean(sampled.values)
    out.df$sd <- sd(sampled.values)
    out.df$f_above_1 <- sum(sampled.values > 1) / length(sampled.values)
    out.df$f_above_1_5 <- sum(sampled.values > 1.5) / length(sampled.values)
    out.df$f_below_0_5 <- sum(sampled.values < 0.5) / length(sampled.values)

    out.df$error <- 0
    out.df$aic <- fit$aic
#    out.df$lnl <- fit$loglik
    out.df$fit_str <- as.character(fit)[1]
    out.df$est_1 <- fit$estimate[1]
    if (length(fit$estimate) > 1) {
      out.df$est_2 <- fit$estimate[2]
    }
    if (write.to.table) {
      con <- connect(db())
      write.or.update(out.df, 'fitdistr', con, 'label')
      dbDisconnect(con)  
    } else {
      print(out.df)
    }
  }
  out.df
}


process.sites <- function(df, do.quantiles=T) {
  gc(F)
  print("Processing sites...")

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

    # Add lrt_stat quantiles and z-scores
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

summarize.sites <- function(sites, filter='', 
  short.output=F, 
  recomb.gc=F, 
  primate.subs=F,
  do.bootstrap=F
) {
  if (primate.subs) {
    sub.sites <- subset(sites, !is.na(mut_ws))
    #print(nrow(sites))
    sites <- sites[!duplicated(sites[, c('data_id', 'aln_position')]),]
    #print(nrow(sites))
  }
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
    rm(pos.domains)
    n.pos.genes <- length(unique(pos.sites$data_id))

    `f.less.0.5` <- nrow(subset(sites, omega < 0.5)) / nrow(sites)
    `f.less.1` <- nrow(subset(sites, omega < 1)) / nrow(sites)
    `f.gt.1` <- nrow(subset(sites, omega > 1)) / nrow(sites)
    `f.gt.1.5` <- nrow(subset(sites, omega > 1.5)) / nrow(sites)

    n.pos.sites <- nrow(pos.sites)
    rm(pos.sites)
    pos.1 <- nrow(subset(sub.above.one, pval < 0.01))
    pos.5 <- nrow(subset(sub.above.one, pval < 0.05))
    pos.10 <- nrow(subset(sub.above.one, pval < 0.1))
    rm(sub.above.one)

    gc(F)

    pos.bh.05 <- nrow(subset(sites, pos.bh < 0.05))
    pos.bh.10 <- nrow(subset(sites, pos.bh < 0.10))

    n.const <- nrow(subset(sites, note == 'constant'))
    n.syn <- nrow(subset(sites, note == 'synonymous'))

    pos.f5 <- nrow(subset(sites, omega > 1 & pval < 0.05))
    neg.f5 <- nrow(subset(sites, omega < 1 & pval < 0.05))
    neutral.f5 <- nrow(subset(sites, pval > 0.05))

    #print(head(sites$pval))
    pos.f10 <- nrow(subset(sites, omega > 1 & pval < 0.1))
    neg.f10 <- nrow(subset(sites, omega < 1 & pval < 0.1))
    neutral.f10 <- nrow(subset(sites, pval > 0.1))

    med.bl <- median(sites$nongap_bl)
    mean.bl <- mean(sites$nongap_bl)
    sd.bl <- sd(sites$nongap_bl)
    med.ncod <- median(sites$ncod)

    omega.cap <- 5
    capped.omegas <- pmin(omega.cap, sites$omega)

    med.omega <- median(capped.omegas)
    mean.omega <- mean(capped.omegas)
    sd.omega <- sd(capped.omegas)
    max.omega <- max(capped.omegas)
    mean.abv <- mean(capped.omegas[capped.omegas > 0])
    sd.abv <- sd(capped.omegas[capped.omegas > 0])
    mean.blw <- mean(capped.omegas[capped.omegas < 1])
    sd.blw <- sd(capped.omegas[capped.omegas < 1])
    rm(capped.omegas)

    sorted.omgs <- sort(sites$omega)
    n.quart <- nrow(sites) / 10
    mid.omgs <- sorted.omgs[c(floor(.5*n.quart):floor(9.5*n.quart))]
    mean.mid.omega <- mean(mid.omgs)
    rm(sorted.omgs)

    boot.mean.fn <- function(D, d) {
      E <- D[d]
      mean(E)
    }
    boot.fraction.fn <- function(D, d) {
      E <- D[d]
      sum(E) / length(E)
    }
    pos.ci.lo <- 0
    pos.ci.hi <- 1
    neg.ci.lo <- 0
    neg.ci.hi <- 1
    if (do.bootstrap) {
      sites$tmp <- sites$omega > 1 & sites$pval < 0.1
      pos.res <- boot(sites$tmp, boot.fraction.fn, R=100)
      sites$tmp <- sites$omega < 1 & sites$pval < 0.1
      neg.res <- boot(sites$tmp, boot.fraction.fn, R=100)
      ci.res <- boot.ci(pos.res, type='basic', conf=0.95) 
      if (!is.null(ci.res$basic)) {
        pos.ci.lo <- ci.res$basic[1, 4]
        pos.ci.hi <- ci.res$basic[1, 5]
      }
      ci.res <- boot.ci(neg.res, type='basic', conf=0.95) 
      if (!is.null(ci.res$basic)) {
        neg.ci.lo <- ci.res$basic[1, 4]
        neg.ci.hi <- ci.res$basic[1, 5]
      }
    }

    mean.lrt <- mean(sites$lrt_stat)
    lrt.qs <- quantile(sites$lrt_stat, c(0.25, 0.5, 0.75))

    df.out <- data.frame(
      pset = pset,
      label = lbl,
      filter = filter,

      n.sites = n.sites,
      n.genes = n.total.genes,
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

      mean.lrt = mean.lrt,
      lrt.25 = lrt.qs[1],
      lrt.50 = lrt.qs[2],
      lrt.75 = lrt.qs[3],

      'f.less.0.5' = `f.less.0.5`,
      'f.less.1' = `f.less.1`,
      'f.gt.1' = `f.gt.1`,
      'f.gt.1.5' = `f.gt.1.5`,

      pos.10 = pos.10,
      pos.05 = pos.5,
      pos.01 = pos.1,
      pos.bh.10 = pos.bh.10,
      pos.bh.05 = pos.bh.05,

      pos.f05 = pos.f5 / n.sites,
      neg.f05 = neg.f5 / n.sites,
      neutral.f05 = neutral.f5 /n.sites,

      pos.f10 = pos.f10 / n.sites,
      neg.f10 = neg.f10 / n.sites,
      neutral.f10 = neutral.f10 / n.sites,
      
      pos.f10.lo = pos.ci.lo,
      pos.f10.hi = pos.ci.hi,
      neg.f10.lo = neg.ci.lo,
      neg.f10.hi = neg.ci.hi
    )

    if (short.output) { 
      out.flds <- c('label', 'n.sites', 'med.bl', 'f.const', 
        'neg.f10', 'f.less.0.5', 'pos.f10', 'f.gt.1.5')
      df.out <- subset(df.out,
        select=out.flds
      )
    }

    if (recomb.gc) {
      avg <- (sites$recombM + sites$recombF) / 2
      mean.recombAvg <- mean(avg)
      mean.recombM <- mean(sites$recombM, na.rm=T)
      mean.recombF <- mean(sites$recombF, na.rm=T)
      mean.gc <- mean(sites$gc, na.rm=T)
      df.out <- cbind(df.out, data.frame(
        mean.recombAvg = mean.recombAvg,
        mean.recombM = mean.recombM,
        mean.recombF = mean.recombF,
        mean.gc = mean.gc
      ))
    }
    
    if (primate.subs) {
      calc.subs <- function(x, prefix) {
        mean.gc <- mean(x$gc, na.rm=T)
        f.subs <- nrow(x) / nrow(sites)
        f.ws <- sum(x$mut_ws) / nrow(sites) / (1 - mean.gc)
        f.sw <- sum(x$mut_sw) / nrow(sites) / mean.gc
        f.sw.nocpg <- sum(x$mut_sw_nocpg) / nrow(sites) / mean.gc
        gc.star.nocpg <- f.ws / (f.sw.nocpg + f.ws)
        gc.star <- f.ws / (f.sw + f.ws)
        n.nsyn <- sum(x$mut_nsyn)
        n.syn <- sum(1 - x$mut_nsyn)
        f.nsyn <- n.nsyn / (n.nsyn + n.syn + 1)

        out.x <- data.frame(
          f.subs = f.subs,
          f.ws = f.ws,
          f.sw = f.sw.nocpg,
          f.sw.nocpg = f.sw.nocpg,
          f.nsyn = f.nsyn,
          gc.star = gc.star
        )
        colnames(out.x) <- paste(prefix, colnames(out.x), sep='.')
        out.x
      }

      df.out <- cbind(df.out, calc.subs(sub.sites, 'p'))
      df.out <- cbind(df.out, calc.subs(subset(sub.sites, taxon_id==9606), 'h'))
      df.out <- cbind(df.out, calc.subs(subset(sub.sites, mut_nsyn==0), 'syn'))
    }

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

get.merged.sites <- function(pset.a, pset.b, xtra, test=F, filter='default') {
  print("  getting a...")
  sites.a <- get.pset.sites(pset.a, filter=filter, test=test)
  print("  getting b...")
  sites.b <- get.pset.sites(pset.b, filter=filter, test=test)

  sel.cols <- c(xtra, 'data_id', 'aln_position')
  sites.a <- subset(sites.a, select=sel.cols)
  sites.b <- subset(sites.b, select=sel.cols)

  print(paste("  merging ", pset.a, ' - ', pset.b, "...", sep=''))
  merged.s <- merge(sites.a, sites.b, by=c('data_id', 'aln_position'))
  merged.s
}

get.recomb.rates <- function() {
  final.f <- scratch.f("recombRate.Rdata")
  if (!file.exists(final.f)) {
    url <- 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/recombRate.txt.gz'
    out.f <- scratch.f("recombRate.txt.gz")
    out.txt <- scratch.f("recombRate.txt")
    if (!file.exists(out.f) && !file.exists(out.txt)) {
      system(sprintf("wget %s -O %s", url, out.f))
      system(sprintf("gunzip %s", out.f))
    }

    recomb.df <- read.table(out.txt, header=F)
    colnames(recomb.df) <- c(
      'chr', 'start', 'end',
      'name',
     'decodeAvg', 'decodeF', 'decodeM',
       'marshfieldAvg', 'marshfieldF', 'marshfieldM',
      'genethonAvg', 'genethonF', 'genethonM'
    )
    save(recomb.df, file=final.f)
  }
  load(final.f)
  recomb.df
}

summarize.quantiles <- function(sites, field, quantiles=c(0, seq(from=0.05, to=0.95, length.out=5), 1), ...) {
  # Remove NA values in the given field.
  sites <- sites[!is.na(sites[, field]),]
  
  q.values <- quantile(sites[, field], probs=quantiles, names=F)
  ecdf.f <- ecdf(sites[, field])
  sites$cum <- ecdf.f(sites[, field])
  sites$grps <- cut(sites$cum, quantiles)

  out.df <- ddply(sites, .(grps), summarize.sites, short.output=T, recomb.gc=T, ...)
  out.df <- cbind(data.frame(lb=out.df[1, 'label']), out.df)
  out.df <- cbind(data.frame(fld=field), out.df)
  out.df[duplicated(out.df$lb), 'lb'] <- ''
  out.df[duplicated(out.df$fld), 'fld'] <- ''
  out.df
}

summarize.matrix <- function(sites, a, b, quantiles=seq(from=0, to=1, by=0.2), filter='recomb', plot=T, store=T, test=F) {
  # Remove NA values in the given fields.
  sites <- sites[!is.na(sites[, a]),]
  sites <- sites[!is.na(sites[, b]),]

  q.values <- unique(quantile(sites[, a], probs=quantiles, names=F))
  sites$grps_a <- cut(sites[, a], q.values, include.lowest=T)
  q.values <- unique(quantile(sites[, b], probs=quantiles, names=F))
  sites$grps_b <- cut(sites[, b], q.values, include.lowest=T)

  out.df <- ddply(sites, .(grps_a, grps_b), summarize.sites, short.output=F, recomb.gc=T, primate.subs=T)
  pset <- sites[1, 'parameter_set_id']
  rm(sites)
  gc(F)
  out.df$field_a <- a
  out.df$field_b <- b
  out.df$lvl_a <- as.integer(out.df$grps_a)
  out.df$lvl_b <- as.integer(out.df$grps_b)
  out.df$label <- paste(a, b, out.df$grps_a, out.df$grps_b, pset, filter, test, sep='_')
  out.df$test <- test

  if (store) {
    con <- connect(db())
    write.or.update(out.df, 'matrix_summaries', con, 'label')
    disconnect(con)
  }

}

plot.matrices <- function() {
  for (pset in c(1, 2, 3, 6)) {
    for (filter in c('recomb_10k', 'recomb_100k', 'recomb_1mb')) {
      plot.matrix(pset, filter, 'recombM', 'gc')
      plot.matrix(pset, filter, 'recombF', 'gc')
    }
  }
}

plot.matrix <- function(pset, filter, field_a, field_b) {
  con <- connect(db())
  cmd <- sprintf("select * from matrix_summaries where label like '%%%s%%' and pset=%d and field_a='%s' and field_b='%s'", filter, pset, field_a, field_b)
  out.df <- dbGetQuery(con, cmd)
  disconnect(con)

  out.df <- out.df[order(out.df$lvl_a, out.df$lvl_b),]
  print(out.df$n_sites)
  #print(out.df$lvl_a)
  #print(out.df$lvl_b)

  out.df$grps_a <- factor(out.df$lvl_a, levels=unique(out.df$lvl_a), labels=unique(out.df$grps_a))
  out.df$grps_b <- factor(out.df$lvl_b, levels=unique(out.df$lvl_b), labels=unique(out.df$grps_b))
  #print(out.df$grps_a)
  #print(out.df$grps_b)

  plot.flds <- c(
    'p.f.subs', 'med.bl', 'f.const', 'pos.f10', 'f.gt.1.5', 'neg.f10', 'f.less.0.5', 'p.f.ws', 'p.f.sw.nocpg',
    'p.f.nsyn', 'p.gc.star', 'h.f.ws', 'h.f.nsyn', 'syn.f.ws', 'syn.f.sw.nocpg'
  )
  plot.flds <- gsub('\\.', '_', plot.flds)
  plot.lbls <- c()
  plot.df <- data.frame()
  for (i in 1:length(plot.flds)) {
    fld <- plot.flds[i]
    print(fld)
    df1 <- out.df
    df1$val <- rescale(df1[, fld])
    plot.lbls[i] <- sprintf("%s (%.3f - %.3f)", fld, min(df1[, fld]), max(df1[, fld]))
    df1$fld <- fld
    plot.df <- rbind(plot.df, df1)
  }

  plot.df$fld <- factor(plot.df$fld, levels=plot.flds, labels=plot.lbls)

  p <- ggplot(plot.df, aes(x=grps_a, y=grps_b, fill=val))
  p <- p + theme_bw()
  p <- p + geom_tile()
  p <- p + scale_fill_gradientn("Value", colour=c('white', rgb(0.5, 0.5, 1), 'red'))
  p <- p + scale_colour_gradientn("Value", colour=c('white', rgb(0.5, 0.5, 1), 'red'))
#  p <- p + geom_point(aes(size=n_sites), colour='black', alpha=0.3)
  p <- p + scale_size(to=c(0.5, 3))
  p <- p + scale_x_discrete(field_a)
  p <- p + scale_y_discrete(field_b)
  p <- p + facet_wrap(~ fld, ncol=6)
  p <- p + coord_equal()
  p <- p + opts(
    axis.text.x = theme_text(angle=90, hjust=1)
  )

  plot.f <- scratch.f(sprintf("matrix_%s_%s_%s_%s.pdf", pset, filter, field_a, field_b))
  pdf(file=plot.f, height=10, width=18)
  print.ggplot(p)
  dev.off()
}

bsub.cumulative.calc <- function() {
  con <- connect(db())
  for (pset in c(3)) {
    for (filter in c('recomb_100k')) {
      for (direction in c('window')) {
        for (sort.f in c('lrt_stat', 'recombM', 'recombAvg', 'gc', 'nongap_bl')) {
          df <- dbGetQuery(con, sprintf("select * from top_sites_cum where 
            pset=%s and filter='%s' and direction='%s' and sort_f='%s' and i=20", pset, filter, direction, sort.f
          ))
          if (nrow(df) == 0) {
            xtra <- sprintf("%d %s %s %s %s", pset, filter, sort.f, direction, FALSE)
            bsub.function('cumulative_calc', extra.args=xtra, mem=14)
          }
        }
      }
    }
  }
  dbDisconnect(con)
}

# One final recombination test -- look at most strongly positively-selected sites,
# and do a cumulative plot of their M/F recombination, GC content, and WS rates.
cumulative.calc <- function(
  pset,
  filter='recomb_100k',
  sort.f = 'lrt_stat',
  direction='window',
  test=T
) {
  sites <- get.pset.sites(pset, filter=filter, test=test, with.primate.subs=T)
  sites$recombAvg <- (sites$recombF + sites$recombM) / 2

  decrease <- ifelse(direction == 'pos', T, F)
  sites <- sites[order(sites[, sort.f], decreasing=decrease, na.last=NA),]
  sites <- sites[!duplicated(sites[, c('data_id', 'aln_position')]),]
  sites <- sites[!is.na(sites[, sort.f]),]

  if (sort.f %in% c('recombM', 'recombF', 'recombAvg')) {
    sites <- sites[sites[, sort.f] > 0, ]
  }

  # Now, step through each few sites and re-calculate the cumulative summary.
  summary.f <- function(lo, hi, i, x.value) {
    print(paste(lo, hi))
    x.sub <- sites[lo:hi,]
    df.out <- summarize.sites(x.sub, filter=filter, short.output=F, recomb.gc=T, primate.subs=T)

    do.boot <- F
    if (do.boot) {
      # Generate bootstraps and add 95% CIs.
      boot.fn <- function(D, d) {
        boot.sub <- D[d, ]
        res <- summarize.sites(boot.sub, filter=filter, short.output=T, recomb.gc=T, primate.subs=T)
        vals <- c(res$mean.recombM, res$mean.recombF, res$mean.gc, res$p.f.ws, res$p.f.sw)
        vals
      }
      print("Bootstrapping...")
      boot.res <- boot(x.sub, boot.fn, R=50)
      rM.ci <- boot.ci(boot.res, type='basic', conf=0.95, index=1)$basic
      print(rM.ci)
      rF.ci <- boot.ci(boot.res, type='basic', conf=0.95, index=2)$basic
      gc.ci <- boot.ci(boot.res, type='basic', conf=0.95, index=3)$basic
      ws.ci <- boot.ci(boot.res, type='basic', conf=0.95, index=4)$basic
      sw.ci <- boot.ci(boot.res, type='basic', conf=0.95, index=5)$basic
      df.out$mean_recombM_lo <- rM.ci[4]
      df.out$mean_recombM_hi <- rM.ci[5]
      df.out$mean_recombF_lo <- rF.ci[4]
      df.out$mean_recombF_hi <- rF.ci[5]
      df.out$mean_gc_lo <- gc.ci[4]
      df.out$mean_gc_hi <- gc.ci[5]
      df.out$p_f_ws_lo <- ws.ci[4]
      df.out$p_f_ws_hi <- ws.ci[5]
      df.out$p_f_sw_lo <- sw.ci[4]
      df.out$p_f_sw_hi <- sw.ci[5]
    }

    df.out$min_lrt <- min(x.sub$lrt_stat)
    df.out$max_lrt <- max(x.sub$lrt_stat)

    df.out$min_sort_f <- min(x.sub[, sort.f])
    df.out$max_sort_f <- max(x.sub[, sort.f])
    df.out$mean_sort_f <- mean(x.sub[, sort.f])

    df.out$label <- paste(pset, df.out$label, direction, sort.f, lo, hi, sep=' ')

    df.out$filter <- filter
    df.out$pset <- pset
    df.out$lo <- lo
    df.out$hi <- hi
    df.out$i <- i
    df.out$x.value <- x.value
    df.out$direction <- direction
    df.out$sort.f <- sort.f

    con <- connect(db())
    write.or.update(df.out, 'top_sites_cum', con, 'label')
    disconnect(con)
  }

  if (direction == 'window') {
    n <- nrow(sites)
    seq.by <- ifelse(test, floor(n/20), floor(n/20))
    indices <- seq(from=seq.by, to=n, by=seq.by)
  } else {
    n <- nrow(sites) / 10
    seq.len <- ifelse(test, 20, 20)
    indices <- floor(exp(seq(from=floor(log(n/100)), to=log(n), length.out=seq.len)))
    print(indices)
  }

  for (i in 1:length(indices)) {
    if (direction == 'window') {
      summary.f(indices[i] - seq.by, indices[i], i, indices[i] - seq.by/2)
    } else {
      summary.f(1, indices[i], i, indices[i])
    }
  }
}

download.table <- function(tbl) {
  filename <- scratch.f(sprintf("%s.Rdata", tbl))
  con <- connect(db())
  cmd <- sprintf("select * from %s", tbl)
  df <- dbGetQuery(con, cmd)
  disconnect(con)
  save(df, file=filename)
}

cumulative.plots <- function() {
  fn <- function(pset, sort, fields) {
    cumulative.plot(pset=pset, sort.f=sort, fields=fields)
  }

  # LRT stat
  for (i in c(1, 2, 3, 6)) {
    fn(i, 'lrt_stat', 'subs')
  }
  fn(1, 'lrt_stat', 'recomb_gc')

  # Branch Length
  for (i in c(1, 6)) {
    fn(i, 'nongap_bl', 'slr')
    fn(i, 'nongap_bl', 'recomb_gc')
  }

  # GC
  for (i in c(1, 2, 3, 6)) {
    fn(i, 'gc', 'slr')
  }
  fn(1, 'gc', 'subs')
  fn(1, 'gc', 'recomb_gc')

  # Recombination
  for (i in c(1, 2, 3, 6)) {
    fn(i, 'recombAvg', 'slr')
  }
  fn(1, 'recombAvg', 'subs') 
  fn(1, 'recombAvg', 'recomb_gc') 

}

cumulative.plot <- function(pset=pset, 
  filter='recomb_100k', 
  sort.f='lrt_stat', 
  direction='window',
  fields = ''
) {
  con <- connect(db())
  cmd <- sprintf("select * from top_sites_cum where pset=%d and
  filter='%s' and direction='%s' and sort_f='%s'", pset, filter, direction, sort.f)
  df <- dbGetQuery(con, cmd)
  disconnect(con)

  df$logmean_lrt <- log(pmax(0.1, abs(df$mean_lrt)))
  df$nsyn_f_subs <- df$p_f_subs - df$syn_f_subs
  df$f_nsyn <- 1 - (df$f_const + df$f_syn)
  df$f_sub <- 1 - (df$f_const)
  df$f_pos_01 <- df$pos_01 / df$n_sites
  df$log_recombM <- log(df$mean_recombM)
  df$log_gc <- log(df$mean_gc)

  possible.fields <- list(
    mean_recombM = 'Male Recombination',
    mean_recombF = 'Female Recombination',
    log_recombM = 'Male Recombination (log)',
    mean_gc = 'GC Content',
    log_gc = 'GC Content (log)',
    mean_bl = 'Mean Nongap Branch Length',
    p_f_subs = 'Primate Substitution Rate',
    syn_f_subs = 'Primate Synonymous Substitutions Per Codon',
    nsyn_f_subs = 'Primate Nonsynonymous Substitutions Per Codon',
    f_const = 'Fraction Constant Sites',
    f_syn = 'Fraction  Synonymous Sites',
    f_nsyn = 'Fraction Nonsynonymous Sites',
    f_sub = 'Fraction Non-Constant Sites',
    p_f_ws = 'Primate W-S Substitution Rate',
    p_f_sw = 'Primate S-W Substitution Rate',
    syn_f_ws = 'Synonymous W-S Substitution Rate',
    syn_f_sw = 'Synonymous S-W Substitution Rate',
    neg_f10 = 'Fraction NSCs (10% FPR)',
    pos_f10 = 'Fraction PSCs (10% FPR)',
    f_pos_01 = 'Fraction PSCs (1% FPR)',
    f_gt_1 = 'Fracton Sites w>1',
    f_gt_1_5 = 'Fraction Sites w>1.5',
    f_less_0_5 = 'Fraction Sites w<0.5',
    mean_lrt = 'Mean LRT Statistic',
    logmean_lrt = 'Log-Mean LRT Statistic',
    p_gc_star = 'Primate GC*',
    syn_gc_star = 'Synonymous GC*',
    h_gc_star = 'Human GC*',
    mean_mid_omega = 'Mean Omega',
    p_f_nsyn = 'Primate Nsyn Fraction'
  )

  if (fields == 'subs') {
    if ( pset == 1) {
      plot.flds <- c('f_const', 'f_syn', 'f_nsyn', 'syn_f_subs', 'nsyn_f_subs')
    } else {
      plot.flds <- c('f_const', 'f_syn', 'f_nsyn')
    }
  } else if (fields == 'recomb_gc') {
    plot.flds <- c('mean_recombM', 'mean_recombF', 'mean_gc')
  } else if (fields == 'slr') {
    plot.flds <- c('neg_f10', 'pos_f10', 'f_gt_1_5', 'f_less_0_5')
  } else {
    if (sort.f == 'lrt_stat') {
      plot.flds <- unique(c('mean_recombM', 'mean_recombF', 'mean_gc', 'syn_f_subs', 'f_const', 'syn_gc_star'))
    } else if (sort.f == 'gc') {
      plot.flds <- c('f_const', 'p_f_subs', 'syn_f_subs', 'nsyn_f_subs')
    } else if (sort.f == 'recombM') {
      plot.flds <- unique(c('mean_gc', 'syn_f_subs', 'f_const', 'pos_f10', 'neg_f10' ))
    } else if (sort.f == 'nongap_bl') {
      plot.flds <- unique(c('mean_recombM', 'mean_gc', 'syn_f_subs', 'f_const', 'pos_f10', 'neg_f10' ))
    }
  }

  if (fields == 'subs' && sort.f != 'gc') {
    y.str <- 'Value'
  } else {
    y.str <- 'Value'
  }

  plot.lbls <- possible.fields[plot.flds]
  
  add.fld <- function(src, dest, fld) {
    df <- src
    df$fld <- fld
    df$val <- df[, fld]
    df$orig.val <- df$val
    if (direction %in% c('neg', 'pos')) {
      df$val <- df$val / df[nrow(df), 'val']
    } else {
      if (fields == 'subs' && sort.f != 'gc') {

      } else {
#        df$val <- df$val / mean(df[, 'val'])
#         df$val <- (df$val - min(df$val)) / diff(range(df$val))
      }
    }
    dest <- rbind.fill(dest, df)
  }
  plot.df <- data.frame()
  for (fld in plot.flds) {
    plot.df <- add.fld(df, plot.df, fld)
  }


#  if (sort.f == 'lrt_stat' && fields == 'subs') {
#    plot.df$val <- pmin(2, plot.df$val)
#    plot.df$val <- pmax(0, plot.df$val)
#  } else if (sort.f == 'recombM' && fields == 'slr') {
#    plot.df$val <- pmin(1.5, plot.df$val)
#    plot.df$val <- pmax(0.5, plot.df$val)    
#  } else {
#    plot.df$val <- pmin(2, plot.df$val)
#    plot.df$val <- pmax(0.5, plot.df$val)
#  }

  plot.df$fld <- factor(plot.df$fld, levels=plot.flds, labels=plot.lbls)
  plot.df$fld <- factor(plot.df$fld)

  x.df <- plot.df[!duplicated(plot.df$i),]
  x.val <- switch(sort.f,
    gc = 'mean_gc',
    recombM = 'mean_recombM',
    recombAvg = 'mean_recombAvg',
    recombF = 'mean_recombF',
    lrt_stat = 'lrt_50',
    nongap_bl = 'mean_bl'
  )
  x.df <- x.df[order(x.df[, x.val]),]
  cuts <- floor(seq(from=1, to=nrow(x.df), length.out=5))
  #print(cuts)
  cuts <- 1:nrow(x.df) %in% cuts
  x.vals <- x.df[cuts, 'hi']
  #print(x.vals)
  cut.vals <- x.df[cuts, x.val]
  #print(cut.vals)
  cut.vals <- sprintf("%.3f", cut.vals)
  #print(cut.vals)

  sort.lbls <- list(
    gc = 'GC Content',
    recombM = 'Male Recombination Rate',
    recombAvg = 'Sex-Averaged Recombination Rate',
    recombF = 'Female Recombination Rate',
    lrt_stat = 'LRT Statistic',
    nongap_bl = 'Non-gap Branch Length'
  )
  sort.f.lbl <- sort.lbls[sort.f]

  p <- ggplot(plot.df, aes(x=hi, y=val, colour=fld, linetype=fld))

  p <- p + stat_smooth(size=0.5, colour='gray', alpha=0.2)
  p <- p + geom_line(alpha=1, size=1.2)

  if (sort.f == 'lrt_stat') {
    neutral.index <- max(which(x.df$lrt_50 < 0))
    neutral.xval <- x.df[neutral.index, 'hi']
    p <- p + geom_vline(xintercept=neutral.xval, colour=gray(0.5), linetype='dashed')
  }

  x.str <- sprintf("%s", sort.f.lbl)
  p <- p + scale_x_continuous(x.str, breaks=x.vals, labels=cut.vals)
  p <- p + scale_y_continuous(y.str)
  p <- p + scale_colour_discrete("Measurement")
  p <- p + scale_linetype_discrete("Measurement")
  p <- p + facet_grid(fld ~ ., scales="free_y")
  p <- generic.opts(p)
  p <- p + opts(
    axis.text.x = theme_text(angle=90, hjust=1),
    legend.position = 'none',
    strip.text.y = theme_blank(),
    strip.background = theme_blank()
  )
  out.f <- sprintf("top_cum_%s_%s_%s_%s_%s.pdf", sort.f, pset, fields, filter, direction)
  pdf(file=out.f, width=4, height=4)
  print.ggplot(p)
  dev.off()
}

calc.recomb <- function(pset, filter='recomb', test=F) {
  sites <- get.pset.sites(pset, filter=filter, test=test)
  sites <- add.primate.subs(sites, test=test)

  print("Partial correlations...")

  # Compute partial correlations a la 
  # http://mbe.oxfordjournals.org/content/23/5/1068/T2.expansion.html
  cor.df <- subset(sites, select=c('lrt_stat', 'recombM', 'gc', 'nongap_bl'))
  cor.df <- subset(cor.df, !is.na(recombM) & !is.na(gc))
  p.res <- pcor(cor.df, method='spearman')
  rm(cor.df)
  gc(F)
  est <- p.res$estimate
  diag(est) <- 0
  min.v <- min(est, na.rm=T)
  max.v <- max(est, na.rm=T)

  write.csv(est, file=scratch.f(sprintf("recomb_pcor_%s_%s.csv", pset, filter)))
  xt <- xtable(est)
  xt <- color.columns(xt, columns=colnames(xt), limits=c(min.v, max.v))
  out.f <- scratch.f(sprintf("recomb_pcor_%s_%s.txt", pset, filter))
  print.latex(xt, file=out.f)

  quant.f <- function(x, fld) {
    print(fld)
    q.df <- summarize.quantiles(x, fld, primate.subs=T)
    # Keep a subset of primate sub fields.
    q.df <- remove.cols(q.df, c(
      'p.f.cpg', 'h.f.cpg', 'syn.f.cpg',
      'p.f.sw', 'h.f.sw', 'syn.f.sw',
      'syn.f.nsyn',
      'p.f.subs', 'h.f.subs'
    ))
    q.df <- remove.cols(q.df, c('pset', 'label', 'filter'))

    pct.flds <- c('neg.f10', 'pos.f10', 'f.less.0.5', 'f.gt.1.5')
    for (pct.fld in pct.flds) {
      q.df[, pct.fld] <- q.df[, pct.fld] * 100
    }

    print(q.df)
    write.csv(q.df, file=scratch.f(sprintf("recomb_quantiles_%s_%s_%s.csv", fld, pset, filter)), row.names=F)

    xt <- xtable(q.df, digits=3)
    xt <- color.columns(xt,   
      columns=c('med.bl', 'neg.f10', 'pos.f10', 'f.less.0.5', 'f.gt.1.5', 'mean.recombM',
        'mean.recombF', 'mean.gc')
    )
    xt
  }
  out.f <- scratch.f(sprintf("recomb_quantiles_%s_%s.txt", pset, filter))

  print("Quantile summaries...")
  xt <- quant.f(sites, 'gc')
  print.latex(xt, file=out.f)
  cat("\\midrule", file=out.f, append=T)
  xt <- quant.f(sites, 'recombM')
  print.latex(xt, file=out.f, append=T)
  cat("\\midrule", file=out.f, append=T)
  xt <- quant.f(sites, 'recombF')
  print.latex(xt, file=out.f, append=T)

  print("Matrices...")
  summarize.matrix(sites, 'recombM', 'gc', filter=filter, plot=F, test=test)
  gc(F)
  summarize.matrix(sites, 'recombF', 'gc', filter=filter, plot=F, test=test)
}

get.gc.content <- function(width=10000) {
  gc.f <- scratch.f(sprintf("gc_human_%d.Rdata", width))
  if (!file.exists(gc.f)) {
    library(BSgenome.Hsapiens.UCSC.hg19)

    window.width <- width

    chr.names <- seqnames(Hsapiens)[1:25]
    df.list <- list()
    gc.fn <- function(x) {
      af <- alphabetFrequency(x)
      sum(af[2:3]) / sum(af[1:4])
    }
    for (i in 1:length(chr.names)) {
      chr.name <- chr.names[i]
      cur.seq <- Hsapiens[[chr.name]]
      print(chr.name)

      roll.my.own <- T
      if (roll.my.own) {
        cur.len <- length(cur.seq)
        if (cur.len < window.width) {
          window.width <- cur.len / 2
        }
        n.by <- floor(window.width / 10)
        window.starts <- seq(from=1, to=cur.len-window.width, by=n.by)
        window.ends <- seq(from=window.width, to=cur.len, by=n.by)
        window.mids <- window.starts + floor(window.width/2)
        window.mids.l <- window.mids - n.by / 2
        window.mids.h <- window.mids + n.by / 2
        v <- Views(cur.seq, start=window.starts, end=window.ends)
        gcs <- viewApply(v, gc.fn, simplify=T)
        cur.gc.df <- data.frame(
          chr=chr.name,
          start = window.mids.l,
          end = window.mids.h,
          gc=gcs,
          stringsAsFactors=F
        )
        cur.gc.df <- subset(cur.gc.df, !is.nan(gc))
      } else {
        cur.seq <- injectHardMask(cur.seq, '-')
        gc.sum <- letterFrequencyInSlidingView(cur.seq, window.width, c('GC', 'AT'), as.prob=F)
        print(head(gc.sum))
        stop()
        print(summary(gc.sum))
        # gc.sum holds windows starting at i, going to i+window.width-1
        colnames(gc.sum) <- NULL

        half.win <- floor(window.width/2)
        L <- nchar(cur.seq)
        gc.L <- length(gc.sum)

        # we want our GC to be in centered windows.
        # SO, for 1 .. window.width/2 we will use the first value
        # and for nchar(cur.seq) - window.width/2 .. nchar(cur.seq)
        # we will use the last value
        lo.start <- seq(from=1, to=half.win)
        lo.end <- seq(from=2, to=half.win+1)
        lo.gc <- rep(gc.sum[1], times=length(lo.start))
        hi.start <- seq(from=L-half.win+1, to=L-1)
        hi.end <- seq(from=L-half.win+2, to=L)
        hi.gc <- rep(gc.sum[length(gc.sum)], times=length(hi.start))
        mid.start <- seq(from=1 + half.win, to=gc.L + half.win)
        mid.end <- seq(from = 2 + half.win, to=gc.L + half.win + 1)
        mid.gc <- gc.sum

        cur.gc.df <- data.frame(
          chr = chr.name,
          start = c(lo.start, mid.start, hi.start),
          end = c(lo.end, mid.end, hi.end),
          gc = c(lo.gc, mid.gc, hi.gc)
        )

        # Only keep values for every window.width/100 bases.
        every.n <- floor(window.width / 100)
        cur.gc.df <- subset(cur.gc.df, start %% every.n == 0)
        print(nrow(cur.gc.df))
        print(head(cur.gc.df, n=500))
        stop()
      }
      df.list[[i]] <- cur.gc.df
      print(nrow(cur.gc.df))
    }
    gc.df <- rbind.fill(df.list)
    save(gc.df, file=gc.f)
  }
  load(gc.f)
  gc.df
}

add.recomb.gc.to.genes <- function(genes, width=width) {

  genes <- subset(genes, select=c('data_id', 'chr_name', 'chr_start',
  'chr_end', 'chr_strand', 'ref_protein_id'))

  print("Getting recombination rates...")
  recomb.ranges <- get.recomb.ranges(width=width)

  print("Getting GC content...")
  gc.ranges <- get.gc.ranges(width=width)

  print("Creating genes GRanges...")

  # Remove genes w/out coordinates for now... but bring 'em back in after.
  na.genes <- subset(genes, is.na(chr_start))
  has.human <- grepl("ENSP0.*", as.character(genes$ref_protein_id))
  no.human <- genes[!has.human, ]
  genes <- genes[has.human, ]
  genes <- subset(genes, !is.na(chr_start))

  st <- genes$chr_start
  ed <- genes$chr_end
  chr <- genes$chr_name
  strd <- genes$chr_strand

  genes.ranges <- GRanges(
    seqnames = chr,
    ranges = IRanges(st, ed),
    strand = strd,
    as.list(genes)
  )

  print("Adding recomb. GRanges...")
  with.recomb <- add.df(genes.ranges, recomb.ranges,
    keep.fields=c('recombM', 'recombF')
  )
  print("Adding GC Granges...")
  with.gc <- add.df(genes.ranges, gc.ranges,
    keep.fields=c('gc')
  )

  with.gc <- subset(with.gc, select=c('data_id', 'gc'))
  mrgd <- merge(with.recomb, with.gc)

  mrgd <- rbind.fill(mrgd, na.genes, no.human)

  mrgd <- subset(mrgd, select=c('data_id', 'recombM', 'recombF', 'gc'))
  if (width == 1e4) {
    colnames(mrgd) <- c('data_id', 'recombM_10k', 'recombF_10k', 'gc_10k')
  } else if (width == 1e5) {
    colnames(mrgd) <- c('data_id', 'recombM_100k', 'recombF_100k', 'gc_100k')
  } else if (width == 1e6) {
    colnames(mrgd) <- c('data_id', 'recombM_1mb', 'recombF_1mb', 'gc_1mb')
  }
  mrgd
}

add.recomb.gc <- function(sites, test=F, width=width) {
  pset <- sites[1, 'parameter_set_id']

  print("Getting recombination rates...")
  recomb.ranges <- get.recomb.ranges(width=width)

  print("Getting GC content...")
  gc.ranges <- get.gc.ranges(width=width)

  if (!(pset %in% c(1, 6))) {
    print("Loading primate sites for chr coords...")
    # If we're in a species group without human seqs, we have to bring in the human locations
    # from the human sites. Keep fingers crossed that there's enough memory around...!
    primate.sites <- get.pset.sites(1, filter='orig', test=test)
    print(str(primate.sites))
    print(str(sites))
    primate.sites <- subset(primate.sites, !is.na(chr_start),
      select=c('chr_start', 'chr_end', 'chr_name', 'data_id',
      'aln_position')
    )
    print("Merging...")
    sites$chr_start <- NULL
    sites$chr_end <- NULL
    sites$chr_name <- NULL
    sites <- merge(sites, primate.sites, by=c('data_id', 'aln_position'))
    print("Done!")
    rm(primate.sites)
    gc(T)
  }

  sites <- subset(sites, !is.na(chr_start))

  print("Creating sites GRanges...")
  st <- sites$chr_start
  ed <- sites$chr_end
  chr <- sites$chr_name
  sites$chr_start <- NULL
  sites$chr_end <- NULL
  sites$chr_name <- NULL
  sites.ranges <- GRanges(
    seqnames = paste('chr', chr, sep=''),
    ranges = IRanges(st, st+2),
    strand = "*",
    as.list(sites)
  )
  rm(sites)
  
  gc(F)

  print("Adding recomb. GRanges...")
  with.recomb <- add.df(sites.ranges, recomb.ranges,
    keep.fields=c('recombM', 'recombF')
  )
  rm(recomb.ranges)
  print("Adding GC Granges...")
  with.gc <- add.df(sites.ranges, gc.ranges,
    keep.fields=c('gc')
  )
  rm(gc.ranges)
  rm(sites.ranges)

  with.gc <- subset(with.gc, select=c('data_id', 'aln_position', 'gc'))
  print("Merging data frames...")
  mrgd <- merge(with.recomb, with.gc)
  mrgd
}

get.recomb.ranges <- function(width) {
  recomb <- collect.recomb.rates(width=width)
  st <- recomb$start
  ed <- recomb$end
  chr <- recomb$chr
  recomb$start <- NULL
  recomb$end <- NULL
  recomb$chr <- NULL
  GRanges(
    seqnames = chr,
    ranges =
      IRanges(st, ed),
    strand = "*",
    as.list(recomb)
  )
}

get.gc.ranges <- function(width=width) {
  gc.df <- get.gc.content(width=width)
  st <- gc.df$start
  ed <- gc.df$end
  chr <- gc.df$chr
  gc.df$start <- NULL
  gc.df$end <- NULL
  gc.df$chr <- NULL
  GRanges(
    seqnames = chr,
    ranges = IRanges(st, ed),
    gc = gc.df$gc
  )
}

# Merge the GC and recombination data into the sites
add.df <- function(orig.ranges, other.ranges, keep.fields) {
  match.indices <- match(orig.ranges, other.ranges)
  other.data <- as.data.frame(other.ranges)
  rm(other.ranges)
  other.data$mrg.id <- 1:nrow(other.data)
  other.data <- subset(other.data, select=c('mrg.id', keep.fields))
  orig.data <- as.data.frame(orig.ranges, stringsAsFactors=F)
  rm(orig.ranges)
  orig.data$mrg.id <- match.indices
  mrgd <- merge(orig.data, other.data, by='mrg.id', all.x=T, all.y=F)
  mrgd$mrg.id <- NULL
  mrgd
}

add.primate.subs <- function(sites, test=T, keep.cols=c('mut_ws', 'mut_sw', 'mut_sw_nocpg', 'is_cpg', 'mut_nsyn', 'taxon_id')) {
  print("Adding primate subs...")
  p.df <- get.primate.subs(test=test)
  p.df <- subset(p.df, mut_count=1)
  p.df <- subset(p.df, confidence > 0.8)
  p.df$aln_position <- p.df$aln_pos

  p.df$cpg_fwd <- (p.df$nuc_from == 'C' & p.df$nuc_to == 'T') & p.df$mut_cpg == 1
  p.df$cpg_rev <- (p.df$nuc_from == 'G' & p.df$nuc_to == 'A') & p.df$mut_rev_cpg == 1

  p.df$is_cpg <- p.df$cpg_fwd | p.df$cpg_rev
  p.df$mut_sw <- 1 - p.df$mut_ws
  # Remove cpg muts from the ws and sw rates.
  p.df$mut_sw_nocpg <- ifelse(p.df$is_cpg, 0, p.df$mut_sw)

  p.df <- subset(p.df, select=c('data_id', 'aln_position', keep.cols))
  sites <- merge(sites, p.df, by=c('data_id', 'aln_position'), all.x=T)
  sites
}

get.primate.subs <- function(test=T) {
  test.str <- ifelse(test, '_test', '')
  subs.f <- scratch.f(sprintf("primate_subs%s.Rdata", test.str))
  if (!file.exists(subs.f)) {
    print("Collecting subs from table...")
    # Get primate substitution events.
    prim.tx <- all.primates()
    prim.str <- paste('(', paste(prim.tx, collapse=', '), ')', collapse='')
    limit.str <- ifelse(test, 'LIMIT 500000', '')
    cmd <- sprintf("select data_id, aln_pos, taxon_id, mut_nsyn, confidence, mut_cpg, mut_rev_cpg, codon_cpg, mut_ws, aa_from, aa_to, nuc_from, nuc_to from subs where
      aln_pos IS NOT NULL and data_id IS NOT NULL and taxon_id IS NOT NULL
    and taxon_id IN %s %s",
      prim.str, limit.str)
    con <- connect(db())
    subs <- dbGetQuery(con, cmd)
    disconnect(con)
    save(subs, file=subs.f)
  }
  load(subs.f)
  subs
}

filter.default <- function(sites) {
  sites <- filter.branchlength(sites)
  sites <- filter.random(sites)
  sites <- filter.ncod(sites)
  sites <- filter.omega(sites)
  sites
}

test.stringent <- function() {
  sites <- get.pset.sites(1, filter='default', test=T)
  sites <- filter.stringent(sites)
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
  sites <- filter.dups(sites)
  sites <- filter.clusters(sites)
  sites <- filter.ncod.stringent(sites)
  sites <- filter.default(sites)

  sites
}

filter.clusters <- function(sites, return.inverse=F) {
  sites <- filter.default(sites)

  # Collect data_id / aln_positions of alignment segments w/ bad windows.
  baddies.f <- scratch.f("win_baddies.Rdata")
  load(baddies.f)
  baddies <- df

  #con <- connect(db())
  #cmd <- sprintf("select * from win_baddies")
  #baddies <- dbGetQuery(con, cmd)
  #disconnect(con)

  baddies <- subset(baddies, is_leaf == 1)
  print(sprintf("%d bad windows", nrow(baddies)))

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

  #print(quantile(bad.summaries$n.bad.positions, c(0.1, 0.25, 0.5, 0.75, 0.9, 1)))
  #print(quantile(bad.summaries$f.bad.positions, c(0.1, 0.25, 0.5, 0.75, 0.9, 1)))

  #f.bad.threshold <- quantile(bad.summaries$f.bad.positions, 0.75)
  f.bad.threshold <- 0.1

  bad.genes <- subset(bad.summaries, f.bad.positions > f.bad.threshold)
  bad.gene.ids <- bad.genes$data_id
  print(sprintf("Calling bad genes with >%.3f sites within bad windows", f.bad.threshold))
  print(sprintf("Yielded %d bad genes", length(bad.gene.ids)))
  
  rm(bad.summaries)
  rm(bad.genes)
  rm(mgd)

  bad.df <- bad.df[!duplicated(bad.df),]

  # Merge the cluster_sites data frame with the sites data frame...
  # There should be no loss here, so double-check that.  Remove any
  # sites that should be filtered based on high NS subs within
  # 15-codon windows.
  bad.str <- paste(bad.df$data_id, bad.df$aln_position, sep=' ')
  print(sprintf("Yielded %d bad sites", length(bad.str)))
  rm(bad.df)
  sites.str <- paste(sites$data_id, sites$aln_position, sep=' ')

  bad.sites <- sites.str %in% bad.str
  bad.genes <- sites$data_id %in% bad.gene.ids

  print(sprintf("Actually got %d bad gene IDs", sum(bad.genes)))
  print(sprintf("Actually got %d bad site IDs", sum(bad.sites)))
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

filter.dups <- function(sites, return.inverse=F) {
  genes <- get.genes()

  #dup.threshold <- quantile(genes$dup_species_count, 0.75)
  dup.threshold <- 2
  print(sprintf("Duplication species count threshold: > %d", dup.threshold))
  bad.genes <- subset(genes, dup_species_count > dup.threshold)
  bad.ids <- bad.genes$data_id
  print(sprintf("%d genes with too many dups", length(bad.ids)))

  if (return.inverse) {
    subset(sites, data_id %in% bad.ids)
  } else {
    subset(sites, !(data_id %in% bad.ids))
  }
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

filter.constant <- function(sites) {
  subset(sites, is.na(note) | note != 'constant')
}

filter.pfam <- function(sites) {
  sites <- filter.default(sites)
  subset(sites, !is.na(pfam_domain))
}

filter.pfam.stringent <- function(sites) {
  sites <- filter.stringent(sites)
  subset(sites, !is.na(pfam_domain))
}

pset.df <- function(factors=F, prefix='') {
  map.list <- list(
    '1' = 'Primates',
    '4' = 'Atlantogenata',
    '10' = 'HMRD',
    '7' = 'Sparse Glires',
    '9' = 'HQ Mammals',
    '2' = 'Glires',
    '3' = 'Laurasiatheria',
    '8' = 'Sparse Mammals',
    '5' = 'Eutheria',
    '6' = 'Mammals'
  )
  char.list <- list(
    '1' = 'p',
    '4' = 'a',
    '10' = 'f',
    '7' = 'sg',
    '9' = 'h',
    '2' = 'g',
    '3' = 'l',
    '8' = 'sm',
    '5' = 'e',
    '6' = 'm'
  )
  df <- data.frame(
    pset_id = as.numeric(names(map.list)),
    pset = names(map.list),
    label = paste(prefix, unlist(map.list), sep=''),
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

use.type.to.alias <- function(use_type, factors=F) {
  fct <- factor(use_type,
    levels=c('omega', 'ci', 'omega_noconstant', 'ci_noconstant'),
    labels = c(
      '\\omgml',
      '\\ci',
      '\\omgml (Excl. Constant)',
      '\\ci (Excl. Constant)'
    )
  )
  if (factors) {
    fct
  } else {
    as.character(fct)    
  }
}

filter.to.alias <- function(filters, factors=F) {
  fct <- factor(filters,
    levels=c('none', 'default', 'stringent', 'pfam', 'pfam_stringent', 'clusters', 'clusters_inverse', 'dups_inverse'),
    labels = c(
      'None',
      'Relaxed',
      'Conservative',
      'Pfam',
      'Pfam Conservative',
      'Clusters Excluded',
      '(Clusters)',
      '(Paralogs)'
    )
  )
  if (factors) {
    fct
  } else {
    as.character(fct)    
  }
}

pset.to.alias <- function(psets, factors=F, ...) {
  map <- pset.df(factors=factors, ...)
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
    '4' = c(9358, 9361, 9371, 9785, 9813),
    '10' = c(9606, 10090, 10116, 9615),
    '7' = c(10090, 10116, 10020, 43179, 10141),
    '9' = c(9606, 9598, 9544, 10090, 10116, 9615, 9913, 9823, 9796),
    '2' = c(9978, 9986, 10020, 10090, 10116, 10141, 43179),
    '3' = c(9365, 9615, 9646, 9685, 9739, 9796, 9823, 9913, 30538,
  42254, 59463, 132908),
    '8' = c(9258, 9315, 9361, 9785, 9606, 10090, 9615),
    '5' = c(9358, 9361, 9365, 9371, 9478, 9483, 9544, 9593, 9598,
  9601, 9606, 9615, 9646, 9685, 9739, 9785, 9796, 9813, 9823, 9913,
  9978, 9986, 10020, 10090, 10116, 10141, 30538, 30608, 30611, 37347,
  42254, 43179, 59463, 61853, 132908),
    '6' = c(9258, 9315, 9358, 9361, 9365, 9371, 9478, 9483, 9544,
  9593, 9598, 9601, 9606, 9615, 9646, 9685, 9739, 9785, 9796, 9813,
  9823, 9913, 9978, 9986, 10020, 10090, 10116, 10141, 13616, 30538,
  30608, 30611, 37347, 42254, 43179, 59463, 61853, 132908)
  )

  aliases <- llply(taxids.list, function(x) {
    paste(sort(taxid.to.alias(x)), collapse=', ')
  })
  counts <- llply(taxids.list, function(x) {
    length(x)
  })

  genes.df <- get.genes.split()
  pset.mpls <- dlply(genes.df, .(pset), function(x) {
    median(x$mpl, na.rm=T)
  })
  pset.bls <- dlply(genes.df, .(pset), function(x) {
    median(x$bl, na.rm=T)
  })

  out.df <- data.frame(
    'Name' = pset.to.alias(names(taxids.list)),
    'Count' = unlist(counts),
    'Species' = paste("\\tiny{", unlist(aliases), "}", sep=""),
    'slrMPL' = 0,
    'slrTotal' = 0
  )
  out.df[names(pset.mpls), 'slrMPL'] <- unlist(pset.mpls)
  out.df[names(pset.bls), 'slrTotal'] <- unlist(pset.bls)

  xt <- xtable(out.df)
  xt <- color.column(xt, 'slrMPL')
  xt <- color.column(xt, 'slrTotal')
  print.latex(xt, scratch.f("table_species_set_summary.txt"))

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
  df <- subset(df, filter %in% c('none', 'default', 'stringent', 'pfam', 'clusters_inverse', 'dups_inverse'))
  df <- df[order(df$pset.fact, df$flt.fact),]
  df[duplicated(df$pset), 'lbl'] <- ''

  first.split <- nrow(df) / 2
  .write.summary.table(df, prefix='filter_summaries', split.at.indices=c(first.split))

  df <- subset(summaries, filter == 'default')
  df <- df[order(df$pset.fact),]
  df[duplicated(df$filter), 'flt'] <- ''
  .write.summary.table(df, prefix='pset_summaries')

  df <- subset(summaries, filter == 'stringent')
  df <- df[order(df$pset.fact),]
  df[duplicated(df$filter), 'flt'] <- ''
  .write.summary.table(df, prefix='pset_summaries_stringent')
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

  df$n_genes <- as.integer(df$n_genes)
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
    paste('(', R.oo::trim(colb), ')', sep='')
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

remove.cols <- function(x, cols) {
  subset(x, select=setdiff(colnames(x), cols))
}

generic.opts <- function(p) {
    p + opts(
      axis.text.x = theme_text(angle=90, hjust=1, size=8),
      strip.text.y = theme_text(angle=0, hjust=0, size=10),
      strip.background = theme_blank()
    ) + theme_bw()
}