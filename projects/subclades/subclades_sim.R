uname  <<- Sys.getenv("USER")
lib <- ifelse(uname == 'gj1', 'src', 'lib')
source(sprintf("~/%s/greg-ensembl/projects/2xmammals/aln.R", lib))
source(sprintf("~/%s/greg-ensembl/projects/2xmammals/tree.R", lib))
source(sprintf("~/%s/greg-ensembl/projects/phylosim/tree.plot.R", lib))

work.f <- function(str) {
  lib <- ifelse(uname == 'gj1', 'src', 'lib')
  normalizePath(sprintf("~/%s/greg-ensembl/projects/subclades/%s", lib, str))
}

scratch.f <- function(str) {
  lib <- ifelse(uname == 'gj1', 'src', 'lib')
  path <- sprintf("~/scratch/subclades/%s", str)
  if (file.exists(path)) {
    return(normalizePath(path))
  } else {
    return(path)
  }
}

result.dir <- function(params) {
  str <- paste(params, collapse='/')
  scratch.f(str)
}

result.f <- function(params, f) {
  paste(result.dir(params), '/', f, sep='')
}

get.species.groups <- function() {
  #c('Primates', 'Glires', 'Laurasiatheria', 'Atlantogenata', 'Mammals')
  c('Primates', 'Glires', 'Mammals')
}


get.all.species.groups <- function() {
  c('Primates', 'Glires', 'Laurasiatheria', 'Atlantogenata', 'Mammals')
}

get.pop.sizes <- function() {
  c(2e4, 3e4, 4e4, 6e4, 1e5, 5e5)
}

get.total.bls <- function() {
  c(0.5, 1, 2, 'observed')
}

get.positive.fs <- function() {
  c(0, 0.01, 0.1)
}

bl.for.species <- function(species) {
  tree.f <- work.f(paste(species, '.nh', sep=''))
  tree <- read.tree(tree.f)
  tree <- tree.remove.leaf(tree, tree.node.with.label(tree, 'Platypus'))
  as.numeric(sprintf("%.3f", tree.total.branch.length(tree)))
}


bsub.sim.aln <- function() {
  species.groups <- get.species.groups()
  pop.sizes <- get.pop.sizes()
  total.bls <- get.total.bls()
  positive.fs <- get.positive.fs()
  n.jobs <- 250

  total.jobs <- 0
  for (species in species.groups) {
    for (pop.size in pop.sizes) {
      for (total.bl in total.bls) {
        for (positive.f in positive.fs) {
          xtra <- paste(species, pop.size, total.bl, positive.f, sep=' ')
          total.jobs <- total.jobs + n.jobs
          bsub.function('sim.aln', extra.args=xtra, jobarray=n.jobs, drop.output=T)
        }
      }
    }
  }
  print(total.jobs)
}

test.sim.aln <- function() {
  sim.aln('Primates', 5e4, 1, 0.05)
}

sim.aln <- function(species, pop.size, total.bl, positive.f, rep=1) {
  library(phylosim)

  pop.size <- as.numeric(pop.size)

  if (total.bl != 'observed') {
    total.bl <- as.numeric(total.bl)    
  }
  positive.f <- as.numeric(positive.f)

  params <- list(
    species=species, 
    pop.size=pop.size,
    total.bl=total.bl,
    positive.f=positive.f
  )

  print(unlist(params))

  dir.create(result.dir(params), recursive=T, showWarnings=F)
  aln.f <- result.f(params, paste('aln_', rep, '.fasta', sep=''))
  tree.out.f <- result.f(params, paste('tree_', rep, '.fasta', sep=''))
  slr.out.f <- result.f(params, paste('slr_', rep, '.txt', sep=''))
  sim.out.f <- result.f(params, paste('sim_', rep, '.txt', sep=''))
  results.f <- result.f(params, paste('results_', rep, '.csv', sep=''))

  tree.f <- work.f(paste(species, '.nh', sep=''))
  tree <- read.tree(tree.f)
  tree <- tree.remove.leaf(tree, tree.node.with.label(tree, 'Platypus'))
  cur.total <- tree.total.branch.length(tree)

  if (total.bl == 'observed') {
    # Sample from a realistic distribution of total branch lengths...
    log.mean.sd <- switch(species, 
      'Primates' = c(log(0.82), 0.4),
      'Glires' = c(log(1.88), 0.4),
      'Laurasiatheria' = c(log(2.15), 0.4),
      'Mammals' = c(log(7.48), 0.27)
    )
    bl.smpl <- rlnorm(1, log.mean.sd[1], log.mean.sd[2])
    print(bl.smpl)
    stopifnot(bl.smpl > 0)
    print(sprintf("Scaling tree by %.2f to total bl=%.2f", bl.smpl/cur.total, bl.smpl))
    tree <- tree.scale.to(tree, bl.smpl)
  } else {
    tree <- tree.scale.to(tree, total.bl)
  }
  
  if (file.exists(aln.f)) {
    aln <- aln.read(aln.f)
  } else {
    print("Simulating sequences...")
    xx <- get.w.distribution(100000, pop.size, positive.f=positive.f)
    omegas <- xx$w
    aln <- do.simulation(tree, omegas, sim.out.f)
    aln.write(aln, aln.f)
  }
  sim.df <- read.csv(sim.out.f)

  if (!file.exists(slr.out.f)) {
    print("Running SLR...")
    tree <- tree.remove.node.labels(tree)
    tree.write(tree, tree.out.f)
    run.slr(tree.out.f, aln.f, slr.out.f)
    file.remove(tree.out.f)
  }

  print("Parsing SLR output...")
  slr.df <- parse.slr.output(slr.out.f)

  slr.df$species <- species
  slr.df$pop.size <- pop.size
  slr.df$total.bl <- total.bl
  slr.df$positive.f <- positive.f
  slr.df$lrt_stat <- ifelse(slr.df$omega < 1, -slr.df$lrt_stat, slr.df$lrt_stat)

  print("Merging sim & SLR sites...")
  merged.df <- merge(sim.df, slr.df, by='site')

  print(head(merged.df))

  if (!file.exists(results.f)) {
    print("Writing merged data...")
    write.csv(merged.df, file=results.f, row.names=F)
  }
}

do.simulation <- function(tree, omegas, sim.out.f, seq.length=1000) {
  ### Start simulating the seqs.
  PSIM_FAST <<- TRUE
  seq.length <- seq.length
  kappa <- 2
  process <- GY94(kappa=kappa)

  codon.seq <- CodonSequence(length=seq.length, processes=list(list(process)))
  #sampleStates(codon.seq)
  #print(summary(codon.seq))

  omega.sample <- sample(omegas, size=seq.length, replace=F)
  #print(omega.sample)
  for (i in 1:length(omega.sample)) {
    setParameterAtSites(codon.seq, process, id="omega", value=omega.sample[i], index=i)
  }
  sampleStates(codon.seq)

  out.df <- data.frame(
    site=1:seq.length,
    sim.omega=sprintf("%.4f", omega.sample),
    sim.kappa=kappa
  )
  write.csv(out.df, file=sim.out.f, row.names=F)

  sim <- PhyloSim(root.seq=codon.seq, phylo=tree)
  Simulate(sim)
  aln <- sim$alignment
  aln <- aln.remove.phylosim.internals(aln)
  aln
}

run.slr <- function(tree.f, aln.f, slr.out.f) {
  script.f <- work.f('../2xmammals/run_slr.pl')
  cmd <- sprintf("perl %s --tree=%s --aln=%s --out=%s",
    script.f, tree.f, aln.f, slr.out.f
  )
  system(cmd)
}

parse.slr.output <- function(out.f) {
  input <- readLines(out.f)
  site.lines <- grep('^\\s+(\\d+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)', input)
  rows <- data.frame()
  for (x in site.lines) {
    line <- input[x]
    p <- unlist(strsplit(line, '\\s+'))
    
    type <- grep('(Constant|Synonymous|Single char)', p, value=T)
    type <- ifelse(length(type) > 0, type[1], '')
    note <- grep('\\s+(\\++|-+)\\s+', p, value=T)
    note <- ifelse(length(note) > 0, note[1], '')

    rows <- rbind(rows, data.frame(
      site = as.numeric(p[2]),
      neutral = as.numeric(p[3]),
      optimal = as.numeric(p[4]),
      omega = as.numeric(p[5]),
      lower = as.numeric(p[6]),
      upper = as.numeric(p[7]),
      lrt_stat = as.numeric(p[8]),
      pval = as.numeric(p[9]),
      pval_adj = as.numeric(p[10]),
      qval = as.numeric(p[11]),
      type = type,
      note = note
    ))
  }
  rows
}

test.collect.sims <- function() {
  collect.sims(get.species.groups()[1], get.pop.sizes()[1], get.total.bls()[1])
}

collect.sims <- function(species.group, pop.size, total.bl) {
  library(plyr)

  pop.size <- as.numeric(pop.size)
  if (total.bl != 'observed') {
    total.bl <- as.numeric(total.bl)
  }

  species.dir <- scratch.f(paste(species.group, pop.size, total.bl, sep='/'))
  print(sprintf("Collecting results from %s ...", species.dir))

  csv.files <- Sys.glob(sprintf("%s/*/*.csv", species.dir))

  df.list <- list()
  for (i in 1:length(csv.files)) {
    print(csv.files[i])
    cur.df <- read.csv(csv.files[i], stringsAsFactors=F)
    #out.df <- rbind.fill(out.df, cur.df)
    df.list[[i]] <- cur.df
  }

  print("Rbinding...")
  out.df <- do.call(rbind.fill, df.list)
  print(nrow(out.df))

  combined.results.f <- scratch.f(sprintf("combined_results_%s_%s_%s.csv", species.group, pop.size, total.bl))
  write.csv(out.df, file=combined.results.f, row.names=F)  
  
  print("Ddplying...") 
  ddply.result <- ddply(out.df, c('species', 'pop.size', 'total.bl', 'positive.f'), summary.f)
  print(nrow(ddply.result))
  ddply.results.f <- scratch.f(sprintf("ddply_results_%s_%s_%s.csv", species.group, pop.size, total.bl))
  write.csv(ddply.result, file=ddply.results.f, row.names=F)

  print("Fitting distribution...")
  #fits.results.f <- scratch.f(sprintf("fit_results_%s_%s_%s.csv", species.group, pop.size, total.bl))
  #fits.result <- ddply(out.df, c('species', 'pop.size', 'total.bl', 'positive.f'), fit.distr.f)
  #write.csv(fits.result, file=fits.results.f, row.names=F)  
}

fit.distr.f <- function(xx, distr) {
  library(dfoptim)

  true.mean <- mean(xx$sim.omega)

  sites <- xx[, c('omega', 'lower', 'upper', 'note', 'sim.omega')]
  colnames(sites) <- c('omega', 'omega_lower', 'omega_upper', 'note', 'sim.omega')

  sites$data_id <- 0
  sites$parameter_set_id <- 0

  if (Sys.getenv("USER") == 'greg') {
    source("~/lib/greg-ensembl/projects/2xmammals/analyze_mammals.R")
  } else {
    source("~/src/greg-ensembl/projects/2xmammals/analyze_mammals.R")
  }  

  #pdf(file="~/scratch/test.pdf")
  #p <- ggplot(xx, aes(x=sim.omega))
  #p <- p + geom_histogram()
  #print(p)
  #dev.off()
  
#  sites$omega <- sites$sim.omega
  out.df <- fit.sites(sites, 'beta', use='ci', i=0, boot.sample=F, write.to.table=F)
  out.df$n.sites <- nrow(sites)
  out.df$true_mean <- true.mean
  out.df
}

test.fit.distr.f <- function() {
  #cur.results <- read.csv(scratch.f('run_2/results_Glires_60000_3.csv'))
#  cur.results <- read.csv(scratch.f('run_3/combined_results_Primates_40000_0.823.csv'))
   cur.results <- read.csv(scratch.f('combined_results_Glires_20000_1.721.csv'))
  out.df <- ddply(cur.results, c('species', 'pop.size', 'total.bl', 'positive.f'), fit.distr.f)
  print(out.df)
}

summary.f <- function(xx) {
  library(boot)
  sub.summary.f <- function(x) {
    lrt.thresh <- qchisq(0.95, df=1)

    n.sites <- nrow(x)
    n.pos <- nrow(subset(x, lrt_stat > lrt.thresh))
    n.neg <- nrow(subset(x, lrt_stat < -lrt.thresh))
    f.pos <- n.pos / n.sites

    n.fp <- nrow(subset(x, lrt_stat > lrt.thresh & sim.omega < 1))
    n.fn <- nrow(subset(x, lrt_stat < lrt.thresh & sim.omega > 1))
    n.tp <- nrow(subset(x, lrt_stat > lrt.thresh & sim.omega > 1))
    n.tn <- nrow(subset(x, lrt_stat < lrt.thresh & sim.omega < 1))

    fdr <- n.fp / (n.fp + n.tp + 1)
    fpr <- n.fp / n.sites
    tpr <- n.tp / n.sites
    sens <- n.tp / (n.tp + n.fn)
    spec <- n.tn  / (n.tn + n.fp)

    if (x[1, 'positive.f'] == 0) {
      fdr <- 0
      sens <- 0
    }

    data.frame(
      n.sites = n.sites,
      n.pos = n.pos,
      f.pos = f.pos,
      n.neg = n.neg,
      n.fp = n.fp,
      n.fn = n.fn,
      n.tp = n.tp,
      n.tn = n.tn,
      fdr = fdr,
      fpr = fpr,
      tpr = tpr,
      spec = spec,
      sens = sens
    )
  }

  main.df <- sub.summary.f(xx)

  #boot.cols <- c('sens')
  boot.cols <- c('f.pos', 'n.pos', 'tpr', 'fpr', 'fdr', 'spec', 'sens')
  for (i in 1:length(boot.cols)) {
    boot.col <- boot.cols[i]
    print(boot.col)

    boot.f <- function(D, d) {
      E <- D[d, ]
      df.out <- sub.summary.f(E)
      df.out[, boot.col]
    }
    n.boot.reps <- 100
    #print("booting")
    boot.result <- boot(xx, boot.f, R=n.boot.reps)
    #print("boot.ci-ing")
    ci.res <- boot.ci(boot.result, type='basic', conf=0.95)
    lo.s <- paste(boot.col, '_lo', sep='')
    hi.s <- paste(boot.col, '_hi', sep='')
    #prin("Don!");

    boot.lo <- ci.res$basic[1, 4]
    boot.hi <- ci.res$basic[1, 5]

    main.df[, lo.s] <- boot.lo
    main.df[, hi.s] <- boot.hi
  }

  main.df
}

test.summary.f <- function() {
  cur.results <- Sys.glob(scratch.f("combined_results_*.csv"))
  cur.results <- cur.results[1]
  x <- read.csv(cur.results)
  out.df <- ddply(x, c('species', 'pop.size', 'total.bl', 'positive.f'), summary.f)
  print(out.df)
}

plot.results <- function() {
  results.rdata.f <- scratch.f("all_results.Rdata")
  if (!file.exists(results.rdata.f)) {
    result.files <- Sys.glob(scratch.f("ddply_results_*.csv"))
    df.list <- list()
    for (i in 1:length(result.files)) {
      cur.f <- result.files[i]
      print(cur.f)
      df.list[[i]] <- read.csv(cur.f)
    }
    results <- do.call(rbind.fill, df.list)
    save(results, file=results.rdata.f)
  }
  load(results.rdata.f)
  out.df <- results

  print(head(out.df))

  x <- out.df
  print(str(x$total.bl))


  x$facet_lbl <- paste(x$species, x$total.bl)
  #print(x$facet_lbl)
  x$pop.size <- factor(x$pop.size)

  plot.df <- data.frame()
  plot.cols <- c('f.pos', 'tpr', 'fpr', 'spec', 'sens')
  for (i in 1:length(plot.cols)) {
    cur.col <- plot.cols[i]
    cur.df <- x
    cur.df$value.lbl <- cur.col
    cur.df$value <- cur.df[, cur.col]
    cur.df$value.lo <- cur.df[, paste(cur.col, '_lo', sep='')]
    cur.df$value.hi <- cur.df[, paste(cur.col, '_hi', sep='')]

    plot.df <- rbind(plot.df, cur.df)
  }

  p <- ggplot(plot.df, aes(x=pop.size, y=value, ymin=value.lo, ymax=value.hi, colour=species))
  p <- p + geom_errorbar(position="dodge")
  p <- p + facet_grid(value.lbl ~ positive.f, scales="free_y")
  pdf(file=scratch.f("../plot_fdr.pdf"))
  print(p)
  dev.off()
}

test.plot.results.hist <- function() {

  cur.results <- Sys.glob(scratch.f("combined_results_Mammals_*observed*.csv"))
  x <- data.frame()
  for (csv in cur.results) {
    print(csv)
    x <- rbind.fill(x, read.csv(csv))
  }
  plot.results.hist(x, 6)

  return()

  cur.results <- Sys.glob(scratch.f("combined_results_Primates_*observed*.csv"))
  x <- data.frame()
  for (csv in cur.results) {
    print(csv)
    x <- rbind.fill(x, read.csv(csv))
  }
  plot.results.hist(x, 1)

  cur.results <- Sys.glob(scratch.f("combined_results_Glires_*observed*.csv"))
  x <- data.frame()
  for (csv in cur.results) {
    print(csv)
    x <- rbind.fill(x, read.csv(csv))
  }
  plot.results.hist(x, 2)



}

plot.results.hist <- function(df, pset.id) {
  library(grid)
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

  out.f <- scratch.f(sprintf("pval_example_%s.pdf", pset.id))
  pdf(file=out.f, width=12, height=2.5)
  xx <- subset(df, total.bl == 'observed' & pop.size == 60000 & positive.f == 0.01)
  plot.pval.example(xx)
  dev.off()

  dx <- 0.02
  
  p <- ggplot(df, aes(x=pos.pval))
  p <- p + theme_bw()
  p <- p + geom_histogram(binwidth=dx, colour=NA, fill='black')
  p <- p + scale_x_continuous(limits=c(0, 1))
#  p <- p + scale_y_continuous(limits=c(0, 5000))
  p <- p + facet_grid(pop.size ~ positive.f)

  print("Plotting...")
  pdf(file=scratch.f(sprintf("plot_hist_%s.pdf", pset.id)))
  print(p)
  dev.off()

}

plot.tree.sizes <- function() {
  species.groups <- get.species.groups()
  
  tree.list <- list()
  for (species in species.groups) {
    tree.f <- work.f(paste(species, '.nh', sep=''))
    tree <- read.tree(tree.f)
    tree <- tree.remove.leaf(tree, tree.node.with.label(tree, 'Platypus'))
    total.bl <- tree.total.branch.length(tree)
    tree.lbl <- sprintf("%s, total bl=%.2f", species, total.bl)
    tree.list[[tree.lbl]] <- tree
  }

  pdf(file=scratch.f("../tree_sizes.pdf"), width=7, height=4)
  p <- ggphylo(tree.list, do.plot=F)
#  p <- tree.plot(tree.list, label.size=2, do.plot=F, line.size=1)
  print(p)
  dev.off()

}

plot.pval.example <- function(sites) {
  print(str(sites))
  sites <- subset(sites, select=c('omega', 'pval', 'pos.pval'))

  frmt <- function(x, ...) {
    sprintf("%.1e", x)
  }

  n.cols <- 4
  bw <- 0.025

  vplayout(n.cols, 1)

  sub.neg <- subset(sites, omega < 1)
  p <- ggplot(sub.neg, aes(x=pval, y=..count..))
  p <- p + geom_histogram(binwidth=bw)
  p <- p + scale_x_continuous("Two-tailed (w < 1)")
  p <- p + scale_y_continuous()
  p <- generic.opts(p)
  p <- p + opts(
    axis.text.y = theme_blank(),
    axis.text.x = theme_blank(),
    axis.ticks = theme_blank(),
    axis.title.y = theme_blank(),
    axis.title.x = theme_blank()  
  )
  print(p, vp=subplot(1, 1))

  y.sq <- seq(from=0, to=2e4, by=.5e4)

  sub.pos <- subset(sites, omega > 1)
  p <- ggplot(sub.pos, aes(x=pval, y=..count..))
  p <- p + geom_histogram(binwidth=bw)
  p <- p + scale_x_continuous("Two-tailed (w > 1)")
  p <- p + scale_y_continuous()
#  p <- p + coord_cartesian(ylim=c(0, 2e4))
  p <- generic.opts(p)
  p <- p + opts(
    axis.text.y = theme_blank(),
    axis.text.x = theme_blank(),
    axis.ticks = theme_blank(),
    axis.title.y = theme_blank(),
    axis.title.x = theme_blank()  
  )
  print(p, vp=subplot(2:3, 1))

  y.sq <- seq(from=0, to=4e4, by=1e4)

  p <- ggplot(sites, aes(x=pos.pval, y=..count..))
  p <- p + geom_histogram(binwidth=bw)
  p <- p + scale_x_continuous("One-tailed")
  p <- p + scale_y_continuous()
  p <- p + coord_cartesian(ylim=c(0, 4e3))
  p <- generic.opts(p)
  p <- p + opts(
    axis.text.y = theme_blank(),
    axis.text.x = theme_blank(),
    axis.ticks = theme_blank(),
    axis.title.y = theme_blank(),
    axis.title.x = theme_blank()  
  )
  print(p, vp=subplot(4, 1))
}


test.get.w.distribution <- function() {
  lst <- get.w.distribution(50000, 50000, positive.f=0.01)
  x <- lst$w
  print(summary(x))
  print(sum(x > 1))
}

get.w.distribution <- function(n, Ne, positive.f=0.01, allow.positive=F, dist='nielsen') {
  human.Ne <- 50000

  if (dist == 'nielsen') {
    S <- rnorm(n, mean=-3, sd=2)
  } else if (dist == 'boyko') {
    # boyko 08 normal
    S = rnorm(n, -38.5, 28.6)
    # boyko 08 lnorm
    #S <- -rlnorm(n, 5.02, 5.94)
    S <- pmax(-25636, S)
    S <- pmin(25636, S)
  } else if (dist == 'wilson') {
    S <- c(
      runif(n * .8, -500, -100),
      runif(n * .2, -100, -50),
      runif(n * .05, -50, -10),
      runif(n * .01, -10, -5),
      runif(n * .01, -5, -1),
      runif(n * .03, -1, -0),
      runif(n * .03, 0, 1),
      runif(n * .02, 1, 5),
      runif(n * .02, 5, 10),
      runif(n * .02, 10, 50)
    )
  }

  if (!allow.positive) {
    #S <- pmin(-0.001, S)
    S <- S[S < 0]
  }
  if (positive.f > 0) {
    n.pos <- ceiling(length(S) * positive.f)
    pos.S <- runif(n.pos, min=0.5, max=3)
    pos.indices <- sample(1:length(S), n.pos, replace=F)
    S[pos.indices] <- pos.S
  }

  s <- S / human.Ne / 2
  ws <- s.to.w(s, Ne)
  list(
    w=ws,
    s=s,
    S=S
  )
}

s.to.w <- function(ss, Ne) {
  # w = (2*N*s) / (1 - e^[-2*n*s])
  two.n.s <- Ne * ss
  two.n.s / (-expm1(-two.n.s))
}

bsub.collect.sims <- function() {
  species.groups <- get.species.groups()
  pop.sizes <- get.pop.sizes()
  #total.bls <- get.total.bls()
  total.bls <- c('observed')

  for (species in species.groups) {
    for (pop.size in pop.sizes) {
      for (total.bl in total.bls) {
        xtra <- paste(species, pop.size, total.bl)
        bsub.function('collect.sims', extra.args=xtra)
      }
    }
  }
}


plot.w.distributions <- function() {
  library(plyr)
  library(ggplot2)
  pop.sizes <- get.pop.sizes()
#  pop.sizes <- c(1e1, 1e2, 1e3, 5e3, 1e4, 2e4, 3e4, 4e4, 5e4, 6e4, 8e4)
  n <- 50000
  human.Ne <- 50000

  lst <- get.w.distribution(n, pop.sizes[1], 0, allow.positive=T)
  xx <- data.frame(
    s=lst$s,
    S=lst$S
  )
  print(summary(xx))
  p <- ggplot(xx, aes(x=s))
  p <- p + geom_histogram()
  pdf(file=scratch.f("../s_distributions.pdf"))
  print(p)
  dev.off()

  w.df <- data.frame()
  mark.df <- data.frame()
  for (f.pos in c(0, 0.1)) {
    for (Ne in pop.sizes) {
      if (f.pos == 0) {
        x <- get.w.distribution(n, Ne, f.pos, allow.positive=T)
        ws <- x$w
      } else {
        x <- get.w.distribution(n, Ne, f.pos, allow.positive=F)
        ws <- x$w
      }
      cur.df <- data.frame(
        w=ws,
        Ne=Ne,  
        f.pos=f.pos
      )
      w.df <- rbind.fill(w.df, cur.df)

      mrks <- s.to.w(c(-0.01, -0.001, -0.0001), Ne)

      mark.df <- rbind.fill(mark.df, data.frame(
        orig.s = as.character(c(-0.01, -0.001, -0.0001)),
        w=mrks,
        Ne=Ne,
        f.pos=f.pos
      ))
      print(paste(Ne, nrow(subset(cur.df, w > 1))))
      print(paste(Ne, nrow(subset(cur.df, w < 0.1))))
    }
  }

  w.df$w <- pmin(5, w.df$w)
  w.df$Ne <- factor(w.df$Ne, levels=sort(unique(w.df$Ne)))  

  mark.df$w <- pmin(5, mark.df$w)
  mark.df$Ne <- factor(mark.df$Ne, levels=sort(unique(mark.df$Ne)))  
#  print(mark.df)

  p <- ggplot(w.df, aes(x=w))
  p <- p + theme_bw()
  p <- p + geom_freqpoly(binwidth=0.02)
  p <- p + geom_vline(data=mark.df, aes(xintercept=w, colour=orig.s), linetype='dashed')
  p <- p + scale_y_continuous(limits=c(0, 2000))
  p <- p + facet_grid(Ne ~ f.pos)
  pdf(file=scratch.f("w_distributions.pdf"), width=8, height=5)
  print(p)
  dev.off()
  
  #pdf(file="~/scratch/test2.pdf", width=4, height=5)
  #p <- ggplot(w.df, aes(x=S))
  #p <- p + theme_bw()
  #p <- p + facet_grid(Ne ~ ., scale="free_y")
  #p <- p + geom_freqpoly(binwidth=0.1)
  #print(p)
  #dev.off()

}

archive.results <- function(label) {
  print(sprintf("Archiving results to %s ...", label))

  archive.f <- scratch.f(label)
  dir.create(archive.f)

  old.wd <- getwd()
  setwd(scratch.f(''))

  species.s <- paste(get.all.species.groups(), collapse=' ')
  system(sprintf("mv *.* %s -vf %s", species.s, label))
  setwd(old.wd)
}