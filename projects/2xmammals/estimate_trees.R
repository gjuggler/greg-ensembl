bsub.estimate.m0 <- function() {
  newicks.df <- get.subset.newicks()

  species.groups <- newicks.df$name
  #species.groups <- c('Primates')

  subdirs <- get.output.subdirs()
  subdirs <- c('concat_nongaps_20')

  out.f <- output.folder(T, T, T)
  for (subdir in subdirs) {
    for (species.group in species.groups) {
      for (n.codons in c(1e4, 1e5, 5e5, 1e6, 2e6)) {
        aln.in <- 'all_combined.fasta'
        #aln.in <- 'alns_1.fasta'
        #aln.in <- '../00/EIF3A.fasta'
        xtra <- paste(out.f, subdir, species.group, aln.in, n.codons, TRUE, sep=' ')
        bsub.function('estimate_m0_tree', mem=16, extra.args=xtra, queue='long')
      }
    }
  }
}

bsub.estimate.branch <- function() {
  newicks.df <- get.subset.newicks()

  species.groups <- newicks.df$name
  species.groups <- c('Primates')

  subdirs <- get.output.subdirs()
  subdirs <- c('concat_nongaps_20')

  out.f <- output.folder(T, T, T)
  for (subdir in subdirs) {
    for (species.group in species.groups) {
#      for (n.codons in c(1e4, 1e5, 5e5, 1e6, 2e6)) {
      for (n.codons in c(1e4, 1e5)) {
        aln.in <- 'all_combined.fasta'
        xtra <- paste(out.f, subdir, species.group, aln.in, n.codons, TRUE, sep=' ')
        bsub.function('estimate_branch_tree', mem=16, extra.args=xtra, queue='long')
      }
      return()
    }
  }
}


test.estimate.m0 <- function() {
  base.dir <- output.folder(T, T, T)
  subdir <- 'concat_nongaps_20'
  species.group <- 'Primates'
  #aln.in <- '../00/EIF3A.fasta'
  aln.in <- 'alns_1.fasta'
  n.codons <- 1000
  estimate.m0.tree(base.dir, subdir, species.group, aln.in, n.codons, clean=T)
}

test.estimate.branch <- function() {
  base.dir <- output.folder(T, T, T)
  subdir <- 'concat_nongaps_20'
  species.group <- 'HMRD'
  #aln.in <- '../00/EIF3A.fasta'
  aln.in <- 'alns_1.fasta'
  n.codons <- 100000
  estimate.branch.tree(base.dir, subdir, species.group, aln.in, n.codons, clean=T)
}

estimate.m0.tree <- function(base.dir,
  subdir='concat_full',
  species.group='HMRD',
  aln.in='test.fasta',
  n.codons=0,
  clean=FALSE) {

  pset <- species.group.to.pset(species.group)
  tree <- get.subset.taxid.tree(pset.id=pset, include.outgroup=T)
  tree <- tree.taxids.to.labels(tree)
  tree$node.label <- NULL

  #print(tree$tip.label)

  #print("loading tree...")
  #tree.csv <- paste(base.dir, 'trees.csv', sep='/')
  #trees <- read.csv(tree.csv, stringsAsFactors=F)
  #cur.row <- subset(trees, name == species.group)
  #tree <- read.tree(text=cur.row$label_newick)
  
  estimate.tree(tree, 'm0', base.dir, subdir, species.group, aln.in, n.codons, clean)
}

estimate.branch.tree <- function(base.dir,
  subdir='concat_full',
  species.group='HMRD',
  aln.in='test.fasta',
  n.codons=0,
  clean=FALSE) {
  library(phylosim)

  print("loading M0 tree...")
  m0.tree.f <- get.tree.file('m0', base.dir, subdir, species.group, aln.in, n.codons)
  print(m0.tree.f)
  tree <- read.tree(file=m0.tree.f)

  estimate.tree(tree, 'free', base.dir, subdir, species.group, aln.in, n.codons, clean)

  nhx.out.f <- paste(base.dir, '/', subdir, '/', 'tree_', species.group, '_', 'free', '.nhx', sep='')
  dnds.plot.f <- paste(base.dir, '/', subdir, '/', 'tree_', species.group, '_', 'free', '_', 'dnds', '.pdf', sep='')

  if (file.exists(nhx.out.f)) {
    tree.s <- readLines(nhx.out.f)
    tree <- read.nhx.tree(tree.s)
    pdf(file=dnds.plot.f)
    x <- PhyloSim()
    x$.phylo <- tree
    p <- plotTree(x, 
#      color.by='dN/dS',
      line.width=2,
      tree.do.plot=F
    )
    p <- p$grob
#    p <- p + scale_colour_gradient()
    print.ggplot(p)
    dev.off()
  }
  
}

test.rphast.msa <- function() {
  library(rphast)
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
  species.group,
  aln.in,
  n.codons
) {
  n.codons <- as.integer(n.codons)

  aln.lbl <- gsub('[^[:alnum:]]', '_', aln.in)
  aln.lbl <- gsub('\\.fasta', '', aln.lbl)
  lbl <- paste(species.group, '_', model, '_', aln.lbl, '_', n.codons, sep='')

  aln.out.f <- paste(base.dir, '/', subdir, '/', 'tree_', lbl, '.nh', sep='')
  aln.out.f
}

estimate.tree <- function(tree,
  model,
  base.dir,
  subdir='concat_full',
  species.group='HMRD',
  aln.in='test.fasta',
  n.codons=0,
  clean=FALSE
) {

  n.codons <- as.integer(n.codons)
  clean <- as.logical(clean)  

  library(rphast)
  tree.f <- tempfile()
  aln.f <- tempfile()

  aln.lbl <- gsub('[^[:alnum:]]', '_', aln.in)
  aln.lbl <- gsub('\\.fasta', '', aln.lbl)

  lbl <- paste(species.group, '_', model, '_', aln.lbl, '_', n.codons, sep='')

  aln.out.f <- paste(base.dir, '/', subdir, '/', 'tmpaln_', lbl, '.fasta', sep='')
  tree.out.f <- paste(base.dir, '/', subdir, '/', 'tree_', lbl, '.nh', sep='')
  nhx.out.f <- paste(base.dir, '/', subdir, '/', 'tree_', lbl, '.nhx', sep='')
  pdf.f <-      paste(base.dir, '/', subdir, '/', 'tree_', lbl, '.pdf', sep='')
  file.out.f <- paste(base.dir, '/', subdir, '/', 'data_', lbl, '.txt', sep='')

  output.files <- c(tree.out.f, nhx.out.f, pdf.f, file.out.f)

  skip.long.stuff <- FALSE
  if (all(file.exists(output.files)) && clean == FALSE) {
    skip.long.stuff <- TRUE
  }

  write.tree(tree, file=tree.f)

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

    if (species.group != 'Mammals') {
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
        file.remove(f)
      }
    }
  }

  clean.s <- ifelse(clean, '--clean', '')

  print("calculating...")
  script.f <- project.f('estimate_m0_tree.pl')
  cmd.s <- sprintf("perl %s --tree=%s --aln=%s --tree_out=%s --nhx_out=%s --file_out=%s --model=%s %s",
    script.f, 
    normalizePath(tree.f), 
    normalizePath(aln.f), 
    normalizePath(tree.out.f), 
    normalizePath(nhx.out.f), 
    normalizePath(file.out.f), 
    model,
    clean.s
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