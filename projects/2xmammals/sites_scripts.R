uname  <- Sys.getenv("USER")
if (uname == 'gj1') {
  source("~/src/greg-ensembl/projects/2xmammals/aln.R")
  source("~/src/greg-ensembl/projects/2xmammals/tree.R")
  source("~/src/greg-ensembl/projects/2xmammals/analyze_mammals.R")
  source("~/src/greg-ensembl/projects/2xmammals/analyze_genes.R")
  source("~/src/greg-ensembl/projects/2xmammals/analyze_substitutions.R")
  source("~/src/greg-ensembl/projects/2xmammals/output_data.R")
  source("~/src/greg-ensembl/projects/2xmammals/estimate_trees.R")
  source("~/src/greg-ensembl/projects/2xmammals/calc_recomb.R")
  source("~/src/greg-ensembl/scripts/xtable_utils.R")
} else {
  source("~/lib/greg-ensembl/projects/2xmammals/aln.R")
  source("~/lib/greg-ensembl/projects/2xmammals/tree.R")
  source("~/lib/greg-ensembl/projects/2xmammals/analyze_mammals.R")
  source("~/lib/greg-ensembl/projects/2xmammals/analyze_genes.R")
  source("~/lib/greg-ensembl/projects/2xmammals/analyze_substitutions.R")
  source("~/lib/greg-ensembl/projects/2xmammals/output_data.R")
  source("~/lib/greg-ensembl/projects/2xmammals/estimate_trees.R")
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
      print(as.character(e))
      con <- connect(db())
      out.df <- data.frame(
        tstmp = mysqlEscapeStrings(con, paste(timestamp(), fn.name,
          paste(args, collapse=' '), collapse=' ')),
        fn_name = as.character(fn.name),
        fn_args = mysqlEscapeStrings(con, paste(args, collapse=' ')),
        error = mysqlEscapeStrings(con, as.character(e))
      )
      write.or.update(out.df, 'errors', con, 'tstmp')
      disconnect(con)
    }

    tryCatch(
      do.call(fn.name, as.list(args[2:length(args)])),
      error = error.f
    )
  }
}

collect_genes <- function(pset, filter='default', n.indices=0, test=F, subset.index=NULL) {
  pset <- as.integer(pset)
  if (!is.null(subset.index)) {
    subset.index <- as.integer(subset.index)
    n.indices <- as.integer(n.indices)
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
  process.genes(genes, pset, filter=filter, subset.index=subset.index,
    n.indices=n.indices, test=test)
  print(" done!")
}

write_alns <- function(chr.s, remove.paralogs, mask.clusters, 
  mask.nongaps, jobs.per.chr, job.index) {
  remove.paralogs <- as.logical(remove.paralogs)
  mask.clusters <- as.logical(mask.clusters)
  mask.nongaps <- as.logical(mask.nongaps)
  jobs.per.chr <- as.integer(jobs.per.chr)
  job.index <- as.integer(job.index)

  # Remove all existing functions and objects
  rm(list=ls(pattern='get', envir=.GlobalEnv), envir=.GlobalEnv)

  print("  loading file...")
  out.f <- scratch.f(sprintf("export_%s.Rdata", chr.s))
  load(out.f)

  print("  writing alns...")
  attach(ensembl_alns)

  write.alns(chr.name=chr.s,
    mask.nongaps=mask.nongaps,
    remove.paralogs=remove.paralogs, 
    mask.clusters=mask.clusters,
    jobs.per.chr=jobs.per.chr,
    job.index=job.index
  )
}

write_alns_data <- function(remove.paralogs, mask.clusters, mask.nongaps) {
  write.alns.data(remove.paralogs, mask.clusters, mask.nongaps)
}

combine_concat_alns <- function(remove.paralogs, mask.clusters, mask.nongaps, subdir) {
  combine.concat.alns(remove.paralogs, mask.clusters, mask.nongaps, subdir)
}

concat_alns <- function(remove.paralogs, mask.clusters, mask.nongaps, n.jobs=1, subset.index=1) {
  concat.alns(remove.paralogs, mask.clusters, mask.nongaps, n.jobs, subset.index)
}

estimate_m0_tree <- function(out.f, subdir, pset, aln.in, n.codons, clean=F) {
  estimate.m0.tree(out.f, subdir, pset, aln.in, n.codons, clean)
}

estimate_branch_tree <- function(out.f, subdir, pset, aln.in, n.codons, clean=F) {
  estimate.branch.tree(out.f, subdir, pset, aln.in, n.codons, clean)
}

collect_alns <- function(chr.s, test=F) {
  test <- as.logical(test)
  collect.alns(chr.s, test=test)
}

collect_domains <- function(pset, fltr='default', n.indices=0, test=F, subset.index=NULL) {
  pset <- as.integer(pset)
  if (!is.null(subset.index)) {
    subset.index <- as.integer(subset.index)
    n.indices <- as.integer(n.indices)
  }
  test <- as.logical(test)

  print("Processing...")
  process.domains(pset, fltr=fltr, subset.index=subset.index,
    n.indices=n.indices, test=test)
  print(" done!")
}

plot_pvals <- function(pset) {
  pset <- as.integer(pset)
  plot.pval.example(pset, test=F)
}

collect_clusters <- function(taxon_id, test=F) {
  taxon_id <- as.integer(taxon_id)
  test <- as.logical(test)

  process.cluster(taxon_id, win.size=15, test=test)
}

collect_sites <- function(pset, filter='default', test=F) {
  pset <- as.integer(pset)
  test <- as.logical(test)
  collect.sites(pset, filter, test)
}

summary_table <- function(pset, filter='default', test=F) {
  test <- as.logical(test)
  pset <- as.numeric(pset)
  sites <- get.pset.sites(pset, filter=filter, test=test)
  df <- summarize.sites(sites, filter=filter)
  if (!test) {
    con <- connect(db())
    write.or.update(df, 'summaries', con, 'label')
    disconnect(con)
  } else {
    print(df)
  }
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

empirical_counts <- function(pset, filter, test=F) {
  pset <- as.integer(pset)
  calc.empirical.counts(pset, filter, test=test)
}

cumulative_calc <- function(pset, filter, sort.f, direction, test=F) {
  pset <- as.integer(pset)
  test <- as.logical(test)
  cumulative.calc(pset, filter, sort.f, direction, test)
}

sarah_output <- function(test=F) {
  sarah.output(test=test)
}

eduard_output <- function(test=F) {
  eduard.output(test=test)
}

fit_distr <- function(pset,
  filter,
  distr,
  use,
  test=F
) {

  pset <- as.integer(pset)
  test <- as.logical(test)
  print("Loading sites...")
  sites <- get.pset.sites(pset, filter=filter, test=test)
  sites <- subset(sites, select=c('data_id', 'omega_lower', 'omega_upper', 'omega', 'parameter_set_id', 'note'))

  print("Fitting...")
  for (i in 1:50) {
    print(i)
    con <- connect(db())
    df <- dbGetQuery(con, sprintf("select * from fitdistr where 
      pset=%s and filter='%s' and dist='%s' and use_type='%s' and i=%d", 
      pset, filter, distr, use, i
    ))
    disconnect(con)
    if (nrow(df) == 0) {
      write.t <- ifelse(test, F, T)
      fit.sites(sites, distr=distr, filter=filter, use=use, i=i, write.to.table=write.t)
    }
  }
}

pset_correlation <- function(filter, pset.a, pset.b, test=F) {
  test <- as.logical(test)
  pset.a <- as.integer(pset.a)
  pset.b <- as.integer(pset.b)

  pset.correlation(filter, pset.a, pset.b, test)  
}


parallel_sites <- function(filter, pset.a, pset.b, test=F) {
  pset.a <- as.integer(pset.a)
  pset.b <- as.integer(pset.b)
  test <- as.logical(test)
  parallel.sites(filter, pset.a, pset.b, test)
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

plot_sites_scatters <- function(pset, test) {
  pset <- as.integer(pset)
  test <- as.logical(test)
  plot.sites.scatters(pset, test)
}

plot_qc_hsitograms <- function(pset, filter='default', test=F) {
  pset <- as.integer(pset)
  test <- as.logical(test)  
  plot.qc.histograms(pset, filter, test)
}

nsyn_plot <- function(pset, test=F) {
  pset <- as.integer(pset)
  test <- as.logical(test)
  
  nsyn.plot(pset, test)
}

main()