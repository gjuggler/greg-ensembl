bsub.estimate.both.trees <- function() {
  subdir <- 'concat_nongaps_20'
  out.f <- output.folder(T, T, T)
  aln.in <- 'all_combined.fasta'

  psets <- c(1, 2, 3, 4)
  reps <- 100
  get.se <- TRUE
  clean <- FALSE
  for (pset in psets) {
      for (n.codons in c(5e5)) {
      xtra <- paste(out.f, subdir, pset, aln.in, n.codons, clean, get.se, sep=' ')

      required.mem <- n.codons / 100000 * 3
      required.mem <- ceiling(required.mem)
      required.mem <- max(6, required.mem)
      required.mem <- min(24, required.mem)
      #bsub.function('estimate_both_trees', mem=required.mem, extra.args=xtra, jobarray=reps, queue='basement')
    }
  }

  psets <- c(6)
  reps <- 100
  get.se <- FALSE
  clean <- TRUE
  for (pset in psets) {
    for (n.codons in c(1e5)) {
      xtra <- paste(out.f, subdir, pset, aln.in, n.codons, clean, get.se, sep=' ')

      required.mem <- n.codons / 100000 * 3
      required.mem <- ceiling(required.mem)
      required.mem <- max(6, required.mem)
      required.mem <- min(18, required.mem)
      bsub.function('estimate_both_trees', mem=required.mem, extra.args=xtra, jobarray=reps, queue='basement')
    }
  }
}

test.estimate.m0 <- function() {
  base.dir <- output.folder(T, T, T)
  subdir <- 'concat_nongaps_20'
  pset <- 1
  #aln.in <- '../00/EIF3A.fasta'
  aln.in <- 'alns_1.fasta'
  n.codons <- 1000
  estimate.m0.tree(base.dir, subdir, pset, aln.in, n.codons, clean=T)
}

test.estimate.branch <- function() {
  base.dir <- output.folder(T, T, T)
  subdir <- 'concat_nongaps_20'
  pset <- 1
  #aln.in <- '../00/EIF3A.fasta'
  aln.in <- 'alns_1.fasta'
  n.codons <- 100000
  estimate.branch.tree(base.dir, subdir, pset, aln.in, n.codons, clean=T)
}

test.estimate.both.trees <- function() {
  base.dir <- output.folder(T, T, T)
  subdir <- 'concat_nongaps_20'
  pset <- 1
  aln.in <- 'alns_1.fasta'
  n.codons <- 500
  clean <- TRUE
  rep <- 1
  estimate.both.trees(base.dir, subdir, pset, aln.in, n.codons, rep, clean=clean, get.se=F)
}

estimate.both.trees <- function(base.dir,
  subdir='concat_full',
  pset=1,
  aln.in='test.fasta',
  n.codons=0,
  rep=0,
  clean=FALSE,
  get.se=FALSE
) {

  estimate.m0.tree(base.dir, subdir, pset, aln.in, n.codons, rep, clean)
  estimate.branch.tree(base.dir, subdir, pset, aln.in, n.codons, rep, clean, get.se)
}

estimate.m0.tree <- function(base.dir,
  subdir='concat_full',
  pset=1,
  aln.in='test.fasta',
  n.codons=0,
  rep=0, 
  clean=FALSE) {

  tree <- get.subset.taxid.tree(pset.id=pset, include.outgroup=T)
  tree <- tree.taxids.to.labels(tree)
  tree$node.label <- gsub('[^[:alnum:]]', '_', tree$node.label)
  
  estimate.tree(tree, 'm0', base.dir, subdir, pset, aln.in, n.codons, rep, clean)
}

plot.subset.trees <- function() {
  load.ggphylo()
  psets <- c(1, 2, 3, 4)

  pdf(file="~/scratch/subset.trees.pdf")

  tree.l <- list()
  for (pset in psets) {
    print(pset)
    tree <- get.subset.taxid.tree(pset.id=pset, include.outgroup=T)
    tree <- tree.taxids.to.labels(tree)
    print(tree$tip.label)
    tree$node.label <- gsub('[^[:alnum:]]', '_', tree$node.label)

    tree <- tree.normalize.branchlengths(tree)

    taxon.lbl <- pset.to.alias(pset)
    tree.l[[taxon.lbl]] <- tree
  }

  ggphylo(tree.l)
  dev.off()
}

load.ggphylo <- function() {
  library(devtools)
  dev_mode(TRUE)
  library(RColorBrewer)
  load_all("~/lib/greg-ensembl/projects/ggphylo", reset=T)
}

estimate.branch.tree <- function(base.dir,
  subdir='concat_full',
  pset=1,
  aln.in='test.fasta',
  n.codons=0,
  rep=0, 
  clean=FALSE,
  get.se=FALSE) {

  library(rphast)

  print("loading M0 tree...")
  m0.tree.f <- get.tree.file('m0', base.dir, subdir, pset, aln.in, n.codons, rep)
  print(m0.tree.f)
  tree <- read.tree(file=m0.tree.f)

  m0.aln.f <- get.tmpaln.file('m0', base.dir, subdir, pset, aln.in, n.codons, rep)
  print(m0.aln.f)
  msa.obj <- read.msa(file=m0.aln.f,
    pointer.only=T,
    tuple.size=3
  )

  model <- 'free'

  estimate.tree(tree, model, base.dir, subdir, pset, aln.in, n.codons, rep, clean, get.se, msa.obj=msa.obj)

  n.codons <- as.integer(n.codons)
  aln.lbl <- gsub('[^[:alnum:]]', '_', aln.in)
  aln.lbl <- gsub('\\.fasta', '', aln.lbl)
  lbl <- paste(pset, '_', model, '_', aln.lbl, '_', n.codons, '_', rep, sep='')

  sub.sub.dir <- paste('est_', pset, '_', n.codons, sep='')
  entire.dir <- paste(base.dir, '/', subdir, '/', sub.sub.dir, sep='')

  nhx.out.f <- paste(entire.dir, '/', 'tree_', lbl, '.nhx', sep='')

  if (file.exists(nhx.out.f)) {
    output.branch.data(nhx.out.f)
  } else {
    print("NHX output not found -- something went wrong??")
    print(nhx.out.f)
  }
  
}

output.branch.data <- function(tree.f) {
  library(devtools)
  dev_mode(TRUE)
  library(RColorBrewer)
  load_all("~/lib/greg-ensembl/projects/ggphylo", reset=T)

  tree <- tree.read.nhx(tree.f)
  #print(str(tree))
  tree <- tree.remove.outgroup(tree)
  #print(str(tree))

  #cur.data <- as.data.frame(tree)
  #cur.data$dN.dS <- as.numeric(cur.data$dN.dS)
  #print(head(cur.data))
  #tree <- tree.load.data(tree, cur.data)

  #print(str(tree))

  tree.pdf.out.f <- gsub('\\.nhx', '_tree.pdf', tree.f)
  dnds.pdf.out.f <- gsub('\\.nhx', '_dnds.pdf', tree.f)
  csv.out.f <- gsub('\\.nhx', '_data.csv', tree.f)

  tree.df <- tree.as.data.frame(tree, order.cladewise=T)

  max.depth <- max(tree.df$Depth)
  #print(tree.df)
  tree.df <- ddply(tree.df, .(id), function(x) {
    x$Label <- paste(x$Label, paste(rep(' ', (x$Depth-1)*2), collapse=''), collapse='')
    x
  })  
  tree.df <- sort.df.by.tree(tree.df, tree)

  #print(str(tree.df))
  tree.df$dN.dS_se <- ifelse(is.na(tree.df$dN.dS_se), 0, tree.df$dN.dS_se)
  tree.df$dN.dS <- as.numeric(tree.df$dN.dS)
  tree.df$dN.dS_se <- as.numeric(tree.df$dN.dS_se)
  tree.df$dN.dS <- ifelse(is.na(tree.df$dN.dS), 0, tree.df$dN.dS)
  tree.df$dN.dS <- ifelse(tree.df$dN.dS > 0.6, 0, tree.df$dN.dS)
  tree.df$low <- tree.df$dN.dS - tree.df$dN.dS_se
  tree.df$low <- pmax(0, tree.df$low)
  tree.df$hi <- tree.df$dN.dS + tree.df$dN.dS_se
  tree.df$hi <- pmin(0.6, tree.df$hi)
  tree.df$Label <- factor(tree.df$Label, levels=tree.df$Label)
  #print(tree.df$Label)

  print(head(tree.df))
  # TODO: Plot the tree w/ branches colored by dnds.
  write.csv(tree.df, file=csv.out.f, row.names=F)
  
  p <- ggplot(tree.df, aes(y=as.numeric(dN.dS), x=Label))
  p <- p + theme_bw()
  p <- p + geom_crossbar(aes(ymin=low, ymax=hi), fill='gray')
  p <- p + opts(
    axis.text.x = theme_text(angle=90, hjust=1)
  )
  pdf(file=dnds.pdf.out.f, height=10, width=6)
  print(p)
  dev.off()

  pdf(file=tree.pdf.out.f, width=8, height=8)
  ggphylo(tree, 
    line.size=2,
#    line.size.by='`dN.dS`',
#    line.size.scale = scale_size_continuous(range=c(1, 3)),
    line.color.by='dN.dS',
    line.color.scale = scale_colour_gradientn(colours=brewer.pal(5, 'YlOrRd')),
    internal.label.angle=15,
    plot.internals=F
  )
  dev.off()

}

test.output.branch.data <- function() {
  #tree.f <- scratch.f('alns_NGPRCW/concat_nongaps_20/tree_6_free_all_combined_fasta_10000.nhx')
  #tree.f <- scratch.f('alns_NGPRCW/concat_nongaps_20/tree_1_free_all_combined_fasta_300000.nhx')
  #tree.f <- scratch.f('alns_NGPRCW/concat_nongaps_20/est_1_100000/tree_1_free_all_combined_fasta_100000_1.nhx')
  tree.f <- scratch.f('alns_NGPRCW/concat_nongaps_20/est_6_100000/tree_6_free_all_combined_fasta_100000_3.nhx')
  #tree.f <- scratch.f('alns_NGPRCW/concat_nongaps_20/est_6_50000/tree_6_free_all_combined_fasta_50000_1.nhx')
  print(tree.f)
  if (file.exists(tree.f)) {
    print("File exists -- loading...")
    output.branch.data(tree.f)
  }
}

test.rphast.msa <- function() {
  require(rphast)
  #ff <- '/nfs/users/nfs_g/gj1/scratch/gj1_2x_63_alt/current/data/alns_NGPRCW/concat_nongaps_20/all_combined.fasta'
  ff <- '/nfs/users/nfs_g/gj1/scratch/gj1_2x_63_alt/current/data/alns_NGPRCW/concat_nongaps_20/alns_1.fasta'
  x <- read.msa(
    file=ff,
    pointer.only=T,
    tuple.size=3
  )

  print(names(x))
  print(ncol.msa(x))
  x <- sub.msa(x, seqs=c('Human', 'Orangutan', 'Chimpanzee', 'Macaque'), keep=T, pointer.only=T)
  strip.gaps.msa(x, strip.mode='all.gaps')
  print(ncol.msa(x))

  #print(sub.msa(x, start.col=1, end.col=50), print.seq=T)
    
}

get.tree.file <- function(
  model,
  base.dir,
  subdir,
  pset,
  aln.in,
  n.codons,
  rep=0
) {
  n.codons <- as.integer(n.codons)

  aln.lbl <- gsub('[^[:alnum:]]', '_', aln.in)
  aln.lbl <- gsub('\\.fasta', '', aln.lbl)
  lbl <- paste(pset, '_', model, '_', aln.lbl, '_', n.codons, '_', rep, sep='')

  sub.sub.dir <- paste('est_', pset, '_', n.codons, sep='')

  entire.dir <- paste(base.dir, '/', subdir, '/', sub.sub.dir, sep='')

  aln.out.f <- paste(entire.dir, '/', 'tree_', lbl, '.nh', sep='')
  aln.out.f
}

get.tmpaln.file <- function(
  model,
  base.dir,
  subdir,
  pset,
  aln.in,
  n.codons,
  rep=0
) {
  n.codons <- as.integer(n.codons)

  aln.lbl <- gsub('[^[:alnum:]]', '_', aln.in)
  aln.lbl <- gsub('\\.fasta', '', aln.lbl)
  lbl <- paste(pset, '_', model, '_', aln.lbl, '_', n.codons, '_', rep, sep='')

  sub.sub.dir <- paste('est_', pset, '_', n.codons, sep='')
  entire.dir <- paste(base.dir, '/', subdir, '/', sub.sub.dir, sep='')

  aln.out.f <- paste(entire.dir, '/', 'tmpaln_', lbl, '.fasta', sep='')
  aln.out.f
}

estimate.tree <- function(tree,
  model,
  base.dir,
  subdir='concat_full',
  pset=1,
  aln.in='test.fasta',
  n.codons=0,
  rep=0,
  clean=FALSE,
  get.se=FALSE,
  msa.obj=NULL
) {

  n.codons <- as.integer(n.codons)
  clean <- as.logical(clean)  
  get.se <- as.logical(get.se)

  library(rphast)
  tree.f <- tempfile()
  aln.f <- tempfile()

  aln.lbl <- gsub('[^[:alnum:]]', '_', aln.in)
  aln.lbl <- gsub('\\.fasta', '', aln.lbl)

  lbl <- paste(pset, '_', model, '_', aln.lbl, '_', n.codons, '_', rep, sep='')

  sub.sub.dir <- paste('est_', pset, '_', n.codons, sep='')
  entire.dir <- paste(base.dir, '/', subdir, '/', sub.sub.dir, sep='')  
  if (!file.exists(entire.dir)) {
    print("Creating directory...")
    print(entire.dir)
    dir.create(entire.dir, recursive=T)
  }

  aln.out.f <- paste(entire.dir, '/', 'tmpaln_', lbl, '.fasta', sep='')
  tree.out.f <- paste(entire.dir, '/', 'tree_', lbl, '.nh', sep='')
  nhx.out.f <- paste(entire.dir, '/', 'tree_', lbl, '.nhx', sep='')
  pdf.f <-      paste(entire.dir, '/', 'tree_', lbl, '.pdf', sep='')
  file.out.f <- paste(entire.dir, '/', 'data_', lbl, '.txt', sep='')

  output.files <- c(tree.out.f, pdf.f, file.out.f)

  skip.long.stuff <- FALSE
  if (all(file.exists(output.files)) && clean == FALSE) {
    skip.long.stuff <- TRUE
  }
  if (!is.null(msa.obj)) {
    skip.long.stuff <- TRUE
    write.msa(msa.obj, file=aln.f)
    write.msa(msa.obj, file=aln.out.f)
  }

  write.tree(tree, file=tree.f)

  print(as.character(tree))

  if (!skip.long.stuff) {  
    # Use rphast for handling large alignments...
    print("loading aln...")
    aln.in.f <- paste(base.dir, subdir, aln.in, sep='/')
    msa <- read.msa(file=aln.in.f, pointer.only=T, tuple.size=3)
    print(sprintf("%d columns", ncol.msa(msa)))
    print(sprintf("%d rows", nrow.msa(msa)))
    print("restrict.to.tree...")
    not.in.tree <- setdiff(names(msa), tree$tip.label)
    #print("not in tree: ")
    #print(not.in.tree)
    msa <- sub.msa(msa, seqs=not.in.tree, keep=F, pointer.only=T)

    if (pset != 6) {
      print("stripping gaps...")
      strip.gaps.msa(msa, strip.mode='all.gaps')
    }

    # Sample the desired number of columns.
    if (n.codons > 0) {
      n <- ncol.msa(msa) / 3
      if (n %% 1 != 0) {
        stop("Codon alignment is not a multiple of 3 nucleotides!")
      }
      if (n.codons > n) {
        print("Trying to sample more codons than the alignment contains -- using the whole alignment...")
        n.codons <- n
        sampled.codons <- 1:n.codons - 1
      } else {
        sampled.codons <- sample(n, size=n.codons, replace=F) - 1
      }
      codon.ind <- sort(c(sampled.codons * 3, sampled.codons * 3 + 1, sampled.codons * 3 + 2))
      codon.ind <- codon.ind + 1
      #print(codon.ind)
      print(sprintf("Sampling %d codons from alignment...", n.codons))
      msa <- msa[, codon.ind, pointer.only=T]
      print(sprintf("Done! %d codons in new alignment", ncol.msa(msa)/3))
    }

    print(sprintf("%d columns", ncol.msa(msa)))
    print(sprintf("%d rows", nrow.msa(msa)))
    #print(msa)

    print("writing aln...")
    write.msa(msa, file=aln.f)
    write.msa(msa, file=aln.out.f)
  }

  if (clean) {
    for (f in output.files) {
      if (file.exists(f)) {
        print(sprintf("Removing output file %s ...", f))
        file.remove(f)
      }
    }
  }

  clean.s <- ifelse(clean, '--clean', '')
  get.se.s <- ifelse(get.se, '--get_se', '')

  print("calculating...")
  script.f <- project.f('estimate_m0_tree.pl')
  cmd.s <- sprintf("perl %s --tree=%s --aln=%s --tree_out=%s --nhx_out=%s --file_out=%s --model=%s %s %s",
    script.f, 
    normalizePath(tree.f),
    normalizePath(aln.f),
    normalizePath(tree.out.f), 
    normalizePath(nhx.out.f),
    normalizePath(file.out.f),
    model,
    clean.s,
    get.se.s
  )
  print(cmd.s)
  system(cmd.s)

  # Plot the tree.
  if (file.exists(tree.out.f)) {
    pdf(file=pdf.f)
    tree <- read.tree(tree.out.f)
    plot(tree)
    dev.off()
  }

  file.remove(tree.f)
  file.remove(aln.f)
}

test.collect.bootstrap.reps <- function() {
  base.dir <- scratch.f('alns_NGPRCW')
  subdir <- 'concat_nongaps_20'
  collect.bootstrap.reps(base.dir, subdir, 6, 100000)
}

collect.bootstrap.reps <- function(base.dir, subdir, pset, n.codons) {  
  library(devtools)
  dev_mode(TRUE)
  library(RColorBrewer)
  load_all("~/lib/greg-ensembl/projects/ggphylo", reset=T)

  n.codons <- sprintf("%d", n.codons)
  sub.sub.dir <- paste('est_', pset, '_', n.codons, sep='')
  entire.dir <- paste(base.dir, '/', subdir, '/', sub.sub.dir, sep='')
  
  data.combined.f <- paste(entire.dir, '/', 'data.combined.Rdata', sep='')
  data.plot.f <- paste(entire.dir, '/', 'data.summary.pdf', sep='')
  data.csv.f <- paste(entire.dir, '/', 'data.summary.csv', sep='')
  data.scatterplot.f <- paste(entire.dir, '/', 'data.bl.dnds.pdf', sep='')
  data.scatterplot2.f <- paste(entire.dir, '/', 'data.bl.dnds_se.pdf', sep='')
  data.treeplot.f <- paste(entire.dir, '/', 'data.tree.pdf', sep='')
  tree.nh.f <- paste(entire.dir, '/', 'data.tree.nh', sep='')
  tree.nhx.f <- paste(entire.dir, '/', 'data.tree.nhx', sep='')

  if (!file.exists(data.combined.f)) {
    print(entire.dir)
    csv.files <- Sys.glob(sprintf("%s/*.csv", entire.dir))

    comb.df <- data.frame()
    for (i in 1:length(csv.files)) {
      print(i)
      x <- read.csv(csv.files[i], stringsAsFactors=F)
      x$rep <- i
      comb.df <- rbind.fill(comb.df, x)
    }

    data <- comb.df
    save(data, file=data.combined.f)
  }

  load(data.combined.f)
  data$label <- data$name

  tree.files <- Sys.glob(sprintf("%s/*_free_*.nh", entire.dir))
  print(paste("Tree #1:",tree.files[1]))
  tree <- tree.read.nhx(tree.files[1])
  print(labels(tree))

  tree <- tree.read.nhx(tree.files[10])
  print(labels(tree))

 summary.df <- ddply(data, .(label), function(x) {
    data.frame(
      n.reps = nrow(x),
      dN.dS_mean = mean(x$dN.dS),
      dN.dS_25 = quantile(x$dN.dS, 0.25),
      dN.dS_75 = quantile(x$dN.dS, 0.75),
      dN.dS_sd = sd(x$dN.dS),
      dS_mean = mean(x$dS),
      dS_sd = sd(x$dS),
      branchlength_mean = mean(x$t),
      label = x[1, 'label']
    )
  })

  summary.df$dnds <- summary.df$dN.dS_mean

  dnds.df <- summary.df[, c('label', 'dnds')]
  tree <- tree.load.data(tree, dnds.df)
  tree <- tree.remove.outgroup(tree)
  tree <- tree.apply.branchlengths(summary.df, tree, column='branchlength_mean')

  nhx.str <- as.character(tree, ignore.tags=F)  
  nh.str <- as.character(tree, ignore.tags=T)
  writeChar(nhx.str, tree.nhx.f)
  writeChar(nh.str, tree.nh.f)

  write.csv(summary.df, file=data.csv.f, row.names=F)

  fctr <- factorize.labels.by.tree(data, tree, indent.depth=T)
  data$label.indent <- fctr

  #data.plot.f <- "~/scratch/data.summary.pdf"
  pdf(file=data.plot.f, width=12, height=8)
  p <- ggplot(data, aes(x=label.indent, y=dN.dS))
  p <- p + theme_bw()
  p <- p + geom_boxplot(size=0.2, alpha=0.7, outlier.size=0)
  p <- p + opts(
    axis.text.x = theme_text(angle=90, hjust=1)
  )
  print(p)
  dev.off()

  pdf(file=data.scatterplot.f)
  p <- ggplot(summary.df, aes(x=dS_mean, y=dN.dS_mean))
  p <- p + theme_bw()
  p <- p + geom_point()
  p <- p + geom_text(aes(label=label), size=2, hjust=0, alpha=0.3)
  print(p)
  dev.off()

  pdf(file=data.scatterplot2.f)
  p <- ggplot(summary.df, aes(x=dS_mean, y=dN.dS_sd))
  p <- p + theme_bw()
  p <- p + geom_point()
  p <- p + geom_text(aes(label=label), size=2, hjust=0, alpha=0.3)
  print(p)
  dev.off()


  pdf(file=data.treeplot.f, width=16, height=8)
  ggphylo(tree,
    line.size=2,
    line.color.by='dnds',
    line.color.scale = scale_colour_gradientn(colours=c('blue', 'red')),
    internal.label.angle=15
  )
  dev.off()

  system(sprintf("cp %s/data*pdf ~/scratch", entire.dir))
}

string.replace <- function(string, map) {
  keys <- names(map)
  for (key in keys) {
    string <- gsub(key, map[[key]], string)
  }
  string
}

package.data <- function(lbl, pset=6, n.codons=100000) {
  base.dir <- scratch.f('alns_NGPRCW')
  subdir <- 'concat_nongaps_20'
  sub.sub.dir <- paste('est_', pset, '_', sprintf("%d",n.codons), sep='')
  entire.dir <- paste(base.dir, '/', subdir, '/', sub.sub.dir, sep='')

  old.wd <- getwd()
  setwd(entire.dir)

  map <- list(
    DATE = format(Sys.time(), "%b %d %Y")
  )

  # Write a README file.
  readme.s <- string.replace("

This package contains genome-wide dN/dS estimates in mammals produced
by Gregory Jordan (greg@ebi.ac.uk) on DATE.

The source data were a concatenated set of alignments of
protein-coding genes incorporating all mammalian species in Ensembl
v63. 

An outline of the alignment processing steps:

1) Ensembl gene annotations and gene trees were used to identify ~16k
genes with largely orthologous patterns in mammals.

2) Genomes with PHRED-like quality scores available were filtered at
Q>25 to remove low-quality sequence data.

2) PRANK (http://code.google.com/p/prank-msa/) was used to re-align
the coding sequences of each set of genes.

3) Paralogous copies were removed and clusters of nonsynonymous
substitutions were filtered out.

4) Alignments for all genes were concatenated, and only columns with
sequence present in at least 20 species were retained.

More complete documentation of the data processing and filtering steps
can be found in Chapters 3 and 4 of Greg Jordan's thesis (sorry!). See
http://github.com/gjuggler/greg-thesis .

To estimate dN/dS for each branch:

1) The phylogenetic tree structure was taken as fixed, following the
structure of recently-published mammalian supertrees.

2) Roughly 100 replicate estimates were produced as follows:

2a) An alignment containing 100k codons was extracted from the complete
concatenated alignment (containing ~2.8M codons).

2b) Branch lengths for the tree were estimated using a PAML M0 model (model=0).

2c) Using those initial branch lengths, branch-specific dN/dS values
were estimated using the PAML free-ratios model (model=1).

3) The resulting set of estimates (1 value per branch x ~100 bootstrap
replicates) is summarized and presented in the attached files:

> data.bl.dnds_se.pdf - A plot of the standard deviation of dN/dS (y
  axis) versus the mean dS value (x axis) for each branch. Shows the
  effect of information content (e.g. more vs. fewer substitutions) on
  the consistency of estimating dN/dS along each branch of the
  tree. Estimates for the shortest branches are expeted to be somewhat
  less reliable.

> data.bl.dnds.pdf - A plot of the mean dN/dS (y axis) versus the mean
  dS value (x axis) for each branch. The general trend towards lower
  mean dN/dS for longer branches suggests a slight bias towards higher
  dN/dS estimates along shorter branches of the tree.

> data.summary.csv - A CSV file summarizing the dN/dS and branch
  length estimates for each branch of the mammalian tree.

> data.summary.pdf - A PDF plotting the data from data.summary.csv. Branches are sorted roughly in 'tree order'.

> data.tree.nh - A Newick-formatted tree of the mammalian species (without Platypus, which was used as an outgroup and subsequently removed).

> data.tree.nhx - A NHX-formatted tree of the mammalian species, as above but with dN/dS values included as NHX annotations.

> data.tree.pdf - A PDF plotting the mammalian tree, with branch lengths and branch colors corresponding to the mean estimates across 100 replicates.
", map)
  
  writeChar(readme.s, "README.txt")

  all.plots <- Sys.glob("data*pdf")
  csv.files <- Sys.glob("data*csv")
  tree.files <- Sys.glob("data*nh*")
  all.files <- c(all.plots, csv.files, tree.files, "README.txt")
  all.files.s <- paste(all.files, collapse=' ')

  out.f <- sprintf("%s.zip", lbl)
  system(sprintf("zip %s %s", out.f, all.files.s))
  system("cp *.zip ~/warehouse")

  setwd(old.wd)
}

archive.results <- function(label) {
  base.dir <- scratch.f('alns_NGPRCW')
  subdir <- 'concat_nongaps_20'
  entire.dir <- paste(base.dir, '/', subdir, '/', label, sep='')

  print(sprintf("Archiving results to %s ...", label))

  archive.f <- entire.dir
  dir.create(archive.f)

  old.wd <- getwd()
  setwd(paste(base.dir, '/', subdir, sep=''))

  system(sprintf("mv est_* -vf %s", label))
  setwd(old.wd)
}