output.folder <- function(remove.paralogs, mask.clusters, mask.nongaps) {
  remove.paralogs <- as.logical(remove.paralogs)
  mask.clusters <- as.logical(mask.clusters)
  mask.nongaps <- as.logical(mask.nongaps)

  n.s <- ifelse(mask.nongaps, 'NG', 'ng')
  p.s <- ifelse(remove.paralogs, 'PR', 'pr')
  c.s <- ifelse(mask.clusters, 'CW', 'cw')
  out.f <- paste('alns_', n.s, p.s, c.s, sep='')

  scratch.f(out.f)
}

sarah_output <- function(test=F) {
  test <- as.logical(test)

  print("getting sites...")
  sites <- get.pset.sites(6, filter='default', test=test)
  n.sites <- nrow(sites)
  sites <- subset(sites, !is.na(chr_start))
  n.mapped <- nrow(sites)

  print("getting genes...")
  genes <- get.genes()
  genes <- subset(genes, ref_taxon_id == 9606)
  keep.gene.cols <- c(
    'data_id',
    'gene_name',
    'ref_protein_id',
    'm_slr_dnds',
    'chr_strand'
  )
  genes <- subset(genes, select=keep.gene.cols)

  keep.sites.cols <- c(
    'data_id',
    'aln_position',
    'seq_position',
    'chr_name',
    'chr_end',
    'chr_start',
    'lrt_stat',
    'omega',
    'omega_upper',
    'omega_lower',
    'pos.pval',
    'pos.bh'
  )
  sites <- subset(sites, data_id %in% unique(genes$data_id), select=keep.sites.cols)

  print("merging genes & sites...")
  merged <- merge(sites, genes, by=c('data_id'))
  merged <- ren(merged,
    'gene_dnds_mammals'='m_slr_dnds'
  )
  merged$gene_name <- as.factor(merged$gene_name)
  merged$ref_protein_id <- as.factor(merged$ref_protein_id)

  get.seq <- FALSE
  if (get.seq) {
    print("getting seqs...")
    # Get the human seq
    con <- connect(dbname())
    data.ids <- unique(merged$data_id)
    data.ids.s <- paste('(', paste(data.ids, collapse=','), ')', sep='')
    cmd <- sprintf("select data_id, aln_position, `9606` codon from seqs where data_id IN %s order by data_id,aln_position",
      data.ids.s)
    seqs.df <- dbGetQuery(con, cmd)
    disconnect(con)
    merged <- merge(merged, seqs.df, by=c('data_id', 'aln_position'))
  }

  calc.pep <- FALSE
  if (calc.pep) {
    print("calculating peps...")
    library(phylosim)
    assign("PSIM_FAST", TRUE, envir=.GlobalEnv)
    ca <- CodonAlphabet()
    aa.df <- adply(merged$codon, 1, function(x) {
      aa <- translateCodon(ca, x)
      aa <- ifelse(is.null(aa), '-', aa)
      data.frame(
        pep=aa
      )
    })
    merged$pep <- aa.df$pep
  }
  n.merged <- nrow(merged)

  print(sprintf("n: %d mapped: %d merged: %d", n.sites, n.mapped, n.merged))
  #print("writing CSV...")
  #write.csv(merged, file=scratch.f(sprintf("mammal_slr_sarah_%s.csv", Sys.Date())), row.names=F)
  print("writing Rdata...")
  save(merged, file=scratch.f(sprintf("mammal_slr_sarah_%s.Rdata", Sys.Date())))
}

eduard_output <- function(test=F) {
  test <- as.logical(test)

  sites <- get.pset.sites(6, filter='default', test=test)
  n.sites <- nrow(sites)
  n.mapped <- nrow(subset(sites, !is.na(chr_start)))

  # Only keep sites with pos.pval < 0.05, i.e. some evidence for pos-sel
  sites <- subset(sites, pos.pval < 0.05)
  sites <- subset(sites, !is.na(chr_start))
  n.pos <- nrow(sites)
  
  genes <- get.genes()
  genes <- subset(genes, ref_taxon_id == 9606)
  keep.gene.cols <- c(
    'data_id',
    'gene_name',
    'ref_protein_id',
    'm_slr_dnds',
    'chr_strand'
  )
  genes <- subset(genes, select=keep.gene.cols)

  keep.sites.cols <- c(
    'data_id',
    'aln_position',
    'lrt_stat',
    'pos.pval',
    'pos.bh',
    'ncod',
    'chr_name',
    'chr_end',
    'chr_start',
    'omega',
    'omega_upper',
    'omega_lower'
  )
  sites <- subset(sites, data_id %in% unique(genes$data_id), select=keep.sites.cols)

  merged <- merge(sites, genes, by=c('data_id'))
  merged <- ren(merged,
    'gene_dnds_mammals'='m_slr_dnds'
  )

  # Get the human seq
  con <- connect(dbname())
  data.ids <- unique(merged$data_id)
  data.ids.s <- paste('(', paste(data.ids, collapse=','), ')', sep='')
  cmd <- sprintf("select data_id, aln_position, `9606` codon from seqs where data_id IN %s order by data_id,aln_position",
    data.ids.s)
  seqs.df <- dbGetQuery(con, cmd)
  disconnect(con)

  merged <- merge(merged, seqs.df, by=c('data_id', 'aln_position'))

  library(phylosim)
  assign("PSIM_FAST", TRUE, envir=.GlobalEnv)
  ca <- CodonAlphabet()
  aa.df <- adply(merged$codon, 1, function(x) {
    aa <- translateCodon(ca, x)
    aa <- ifelse(is.null(aa), '-', aa)
    data.frame(
      pep=aa
    )
  })
  merged$pep <- aa.df$pep
  n.merged <- nrow(merged)

  print(sprintf("n: %d mapped: %d pos: %d merged: %d", n.sites, n.mapped, n.pos, n.merged))
  
  write.csv(merged, file=scratch.f(sprintf("mammal_pscs_%s.csv", Sys.Date())), row.names=F)
}

bsub.collect.alns <- function(...) {
  genes <- get.genes()
  all.chrs <- unique(genes$chr_name)

  for (chr in all.chrs) {
    xtra <- paste(chr, sep=' ')
    bsub.function('collect_alns', mem=6, extra.args=xtra, ...)
  }
}

collect.alns <- function(chr.s, test=F) {
  genes <- get.genes()
  genes <- subset(genes, chr_name == chr.s)
  
  if (test) {
    genes <- head(genes, n=20)
  }

  id.s <- paste(genes$data_id, collapse=',')

  print(sprintf("%d genes to collect", nrow(genes)))

  print("  getting seqs...")
  con <- connect(db()); 
  cmd <- sprintf("select * from seqs where data_id in (%s)", id.s)
  seqs.df <- dbGetQuery(con, cmd)
  disconnect(con)
  seqs.df <- process.seqs.df(seqs.df)

  # Get the cluster windows and genes data frames.
  print("  getting bad windows...")
  windows.df <- get.all.baddies()
  windows.df <- subset(windows.df, data_id %in% genes$data_id)
  print("  getting genes...")
  genes.df <- get.genes()
  genes.df <- subset(genes.df, data_id %in% genes$data_id)

  # Get the mammalian sites corresponding to these genes, and keep only important
  # columns.
  sites <- get.sitewise()
  sites <- subset(sites, data_id %in% genes$data_id, 
    select=c('data_id', 'aln_position', 'seq_position', 'ncod', 'nongap_bl', 'random')
  )

  print('Creating new env...')
  out.env <- new.env()
  assign("seqs.df", seqs.df, envir=out.env)
  assign("windows.df", windows.df, envir=out.env)
  assign("genes.df", genes.df, envir=out.env)

  assign("get.aln", get.aln, envir=out.env)
  assign("get.aln.df", get.aln.df.export, envir=out.env)
  assign("get.cluster.windows", get.cluster.windows.export, envir=out.env)
  assign("get.genes", get.genes.export, envir=out.env) 
  assign("get.sitewise", get.sitewise.export, envir=out.env)
  assign("sitewise.df", sites, envir=out.env)
  assign("get.species.tree", get.species.tree.export, envir=out.env) 
  assign("species.tree", get.species.tree(), envir=out.env) 

  assign("df.to.aln", df.to.aln, envir=out.env)
  assign("aln.stats", aln.stats, envir=out.env)
  assign("taxids.beneath.node", taxids.beneath.node, envir=out.env)
  assign("node.with.label", node.with.label, envir=out.env) 
  assign("get.taxid.df", get.taxid.df.export, envir=out.env)
  assign("taxid.df", get.taxid.df(), envir=out.env)
  assign("taxid.to.alias2", taxid.to.alias2, envir=out.env)
  assign("remove.branchlengths", remove.branchlengths, envir=out.env)
  assign("depth.to.root", depth.to.root, envir=out.env)
  assign("parent.node", parent.node, envir=out.env)
  assign("leaves.beneath", leaves.beneath, envir=out.env)
  assign("extract.subtree", extract.subtree, envir=out.env)
  assign("sort.aln.by.tree", sort.aln.by.tree, envir=out.env)
  assign("restrict.aln.to.seqs", restrict.aln.to.seqs, envir=out.env)

  assign("genomes", get.taxa.df(), envir=out.env)
  assign("species.subsets", get.subset.newicks(), envir=out.env)
  assign("get.subset.newicks", get.subset.newicks.export, envir=out.env)

  assign("write.aln", write.aln, envir=out.env)
  assign("write.alns", write.alns, envir=out.env)
 
  print("  saving file...")
  out.f <- scratch.f(sprintf("export_%s.Rdata", chr.s))
  ensembl_alns <- out.env
  rm(out.env)
  save(ensembl_alns, file=out.f)
  print("Done!")
}

bsub.write.alns <- function(...) {
  genes <- get.genes()
  all.chrs <- unique(genes$chr_name)
  print(all.chrs)
  
  jobs.per.chr <- 15

  for (mask.nongaps in c(TRUE)) {
    for (remove.paralogs in c(TRUE, FALSE)) {
      for (mask.clusters in c(TRUE, FALSE)) {
        for (chr in all.chrs) {
          xtra <- paste(chr, remove.paralogs, mask.clusters, mask.nongaps, jobs.per.chr, sep=' ')
          bsub.function('write_alns', mem=2, extra.args=xtra, jobarray=jobs.per.chr, ...)
        }
      }
    }
  }
}

write.alns <- function(chr.name, jobs.per.chr=1, job.index=1, ...) {
  genes <- get.genes()

  low.ind <- 1
  if (jobs.per.chr > 1) {
    n.genes <- nrow(genes)
    genes.per.job <- ceiling(n.genes / jobs.per.chr)
    # job.index is 1-based
    low.ind <- (job.index - 1) * genes.per.job + 1
    hi.ind <- min(nrow(genes), low.ind + genes.per.job - 1)
    genes <- genes[low.ind:hi.ind,]
    print(sprintf("%d genes for chromosome, this job gets %d", n.genes, nrow(genes)))
  }

  aln.df <- data.frame()
  for (i in 1:nrow(genes)) {
    cur.gene <- genes[i,]
    gene.name <- cur.gene$gene_name
    ref.taxon <- cur.gene$ref_taxon_id
    if (is.na(gene.name)) {
      gene.name <- cur.gene$ref_gene_id
    }
    #print(sprintf("writing %s (%d of %d)", gene.name, i, nrow(genes)))
    #print(cur.gene)
    data.prefix <- cur.gene$data_prefix

    #print("getting aln...")
    aln <- get.aln(cur.gene$data_id, ...)

    #print("creating folder...")
    folder <- output.folder(...)  
    dir.create(folder)
    if (ref.taxon != 9606) {
      dir.create(file.path(folder, 'no_human'))
      data.dir <- paste('no_human/', data.prefix, sep='')
    } else {
      data.dir <- paste(data.prefix, sep='')
    }
    dir.create(file.path(folder, data.dir))
    aln.f <- paste(file.path(folder, data.dir), '/', gene.name, '.fasta', sep='')
    #print("writing aln...")
    write.aln(aln, aln.f)
    print(sprintf("Wrote %s (%d of %d)", aln.f, i, nrow(genes)))

    # Do some stats on the alignment & append to the CSV file
    x <- aln.stats(aln)
    x$data_id <- cur.gene$data_id
    aln.df <- rbind.fill(aln.df, x)
  }

  dir.create(file.path(folder, 'aln_stats'))
  csv.f <- paste(folder, '/', 'aln_stats/aln_stats_', chr.name, '_', job.index, '.csv', sep='')
  write.csv(aln.df, file=csv.f, row.names=F, quote=F)
}

bsub.write.alns.data <- function() {
  for (mask.nongaps in c(TRUE)) {
    for (remove.paralogs in c(FALSE, TRUE)) {
      for (mask.clusters in c(FALSE, TRUE)) {
        # One job per set of alignments, to write out a CSV summary
        # and README file.
        xtra <- paste(remove.paralogs, mask.clusters, mask.nongaps, sep=' ')
        bsub.function('write_alns_data', mem=2, extra.args=xtra)
      }
    }
  }
}

write.alns.data <- function(remove.paralogs=T, mask.clusters=T, mask.nongaps=T) {
  out.f <- output.folder(remove.paralogs, mask.clusters, mask.nongaps)

  ## Write a CSV of all the genes' info.
  print("getting genes...")
  genes <- get.genes()
  keep.gene.cols <- c(
    'data_id',
    'data_prefix',
    'gene_name',
    'ref_protein_id',
    'ref_gene_description',
    'aln_length',
    'm_slr_dnds',
    'chr_name',
    'chr_strand',
    'chr_start',
    'chr_end'
  )
  genes <- subset(genes, select=keep.gene.cols)
  genes.f <- paste(out.f, '/', 'genes.csv', sep='')

  ## Load up the aln_stats file and merge w/ the genes.
  print("getting aln stats...")
  csv.files <- Sys.glob(paste(out.f, '/aln_stats/*.csv', sep=''))
  stats.df <- data.frame()
  for (f in csv.files) {
    stats.df <- rbind.fill(stats.df, read.csv(f))
  }

  print("merging genes & aln stats...")
  print(str(genes))
  print(str(stats.df))
  genes <- merge(genes, stats.df, by=c('data_id'))

  ## Write out the gene / alignment metadata.
  write.csv(genes, file=genes.f, row.names=F)

  ## Write out the species tree.
  print("getting trees...")
  trees <- get.subset.newicks()
  trees.f <- paste(out.f, '/', 'trees.csv', sep='')
  write.csv(trees, file=trees.f, row.names=F)

  ## Write out the taxid-to-species data.
  taxids <- get.taxid.df()
  taxids.f <- paste(out.f, '/', 'taxon_ids.csv', sep='')
  write.csv(taxids, file=taxids.f, row.names=F)
  
  ## Write out genome build metadata.
  genomes <- get.taxa.df()
  genomes.f <- paste(out.f, '/', 'genomes.csv', sep='')
  write.csv(genomes, file=genomes.f, row.names=F)

  ## Write out a short README file.

  readme.s <- "
{timestamp}
{ensembl.string}
{filter.string}
---

This directory contains alignments of mammalian protein-coding genes,
created by Gregory Jordan at the European Bioinformatics
Institute. For more information, contact Greg at: greg@ebi.ac.uk

These alignments were generated using Ensembl gene annotations and
gene trees from the Compara gene tree pipeline, all sourced from
Ensembl release {ensembl.version}. Only trees with largely orthologous
mammalian relationships (i.e., trees without highly prevalent
mammalian duplications or deletions) were included, yielding {n.genes}
alignments.

Due to the prevalence of misannotation, misassembly, and sequencing
error in low-coverage and non-model organism genomes, additional
alignment filters were designed to mask out suspicious sequence data
or suspicious alignment columns as 'N's. The filters applied to create
the alignments within this directory are indicated above. Each filter
is briefly described below:

  > REMOVE_SHORT: remove short sequences (defined as those with less
    than half of the mean sequence length) from each alignment by
    replacing with all gaps.

  > MASK_NONGAPS: mask out alignment columns where the amount of
    non-gap sequences is below a given threshold. The threshold was
    set to remove roughly 5-10% of alignment columns.

  > REMOVE_PARALOGS: when this filter is *not* applied, paralogs
    (defined as multiple sequences from the same species existing
    within the same largely-orthologous tree) are resolved by keeping
    only the longest or least-divergent paralogous copy. This filter
    removes *all* paralogous copies from each alignment by replacing
    sequences with all gaps.

  > MASK_CLUSTERS: due to sequencing, assembly or annotation errors,
    clusters of nonhomologous sequence were apparent in some
    alignments. This filter masks out regions of the alignment with an
    exceptionally high density of non-synonymous substitutions (in the
    99.9% percentile for a given genome).

Several other files are included with additional information about the
genomes, species and alignments used. These files include:

  > genes.csv - various metadata and calculations
    for each gene alignment.
    
  > trees.csv - Newick-format species tree structures, without branch
    lengths, for the 38 mammals and various species subsets.

  > taxon_ids.csv - NCBI taxon identifiers for all species and
    taxonomic nodes (including internal nodes, e.g. African great
    apes) used to create the alignments.

  > genomes.csv - Genome build information for the Ensembl genomes
    included in the analysis.
---
  "

  remove.short <- TRUE
  filter.s <- sprintf('Filters: %s %s %s',
    ifelse(remove.short, 'REMOVE_SHORT', ''),
    ifelse(mask.nongaps, 'MASK_NONGAPS', ''),
    ifelse(remove.paralogs, 'REMOVE_PARALOGS', ''),
    ifelse(mask.clusters, 'MASK_CLUSTERS', '')
  )

  data.l <- list(
    timestamp = format.Date(Sys.Date(), "%B %d, %Y"),
    ensembl.string = 'Ensembl release 63',
    filter.string = filter.s,
    n.genes = nrow(genes),
    ensembl.version = '63'
  )

  keys <- names(data.l)
  for (i in 1:length(data.l)) {
    key <- keys[i]
    val <- data.l[[i]]
    print(paste(key, val))
    readme.s <- gsub(paste('\\{', key, '\\}', sep=''), val, readme.s)
  }
  #cat(readme.s)
  readme.f <- paste(out.f, '/', 'README.txt', sep='')
  cat(readme.s, file=readme.f)
}

get.output.subdirs <- function() {
  c('concat_full', 'concat_no_blanks', 'concat_nongaps_05',
    'concat_nongaps_10', 'concat_nongaps_20')
}

bsub.concat.alns <- function(n.jobs=50) {
  for (remove.paralogs in c(TRUE, FALSE)) {
    for (mask.clusters in c(TRUE, FALSE)) {
      for (mask.nongaps in c(TRUE)) {
        xtra <- paste(remove.paralogs, mask.clusters, mask.nongaps, n.jobs, sep=' ')
        bsub.function('concat_alns', jobarray=n.jobs, extra.args=xtra)
      }
    }
  }
}

concat.alns <- function(remove.paralogs=T, mask.clusters=T, mask.nongaps=T,
  n.jobs=1, subset.index=1) {
  out.f <- output.folder(remove.paralogs, mask.clusters, mask.nongaps)

  n.jobs <- as.integer(n.jobs)
  subset.index <- as.integer(subset.index)

  aln.files <- Sys.glob(paste(out.f, '/*/*.fasta', sep=''))

  # Remove directories with concatenated files...
  concat.dirs <- grepl('concat', aln.files)
  aln.files <- aln.files[!concat.dirs]

  low.ind <- 1
  if (n.jobs > 1) {
    n.files <- length(aln.files)
    files.per.job <- ceiling(n.files / n.jobs)
    # job.index is 1-based
    low.ind <- (subset.index - 1) * files.per.job + 1
    hi.ind <- min(n.files, low.ind + files.per.job - 1)
    aln.files <- aln.files[low.ind:hi.ind]
    print(sprintf("%d files, this job gets %d", n.files, length(aln.files)))
  }
  #aln.files <- aln.files[1:10]

  print(sprintf("%d alignments to concatenate", length(aln.files)))
  aln <- NULL
  for (aln.f in aln.files) {
    print(aln.f)
    cur.aln <- read.aln(aln.f)

    #print(aln.is.codon(cur.aln))

    if (is.null(aln)) {
      aln <- cur.aln
    } else {
      aln <- aln.concat(aln, cur.aln)
    }
  }
  
  print(aln.length(aln) / 3)

  dir.create(file.path(out.f, 'concat_full'))
  aln.f <- paste(out.f, '/concat_full/alns_', subset.index, '.fasta', sep='')
  aln.write(aln, aln.f)

  aln.noblanks <- aln.remove.allgap.columns(aln)
  dir.create(file.path(out.f, 'concat_no_blanks'))
  aln.f <- paste(out.f, '/concat_no_blanks/alns_', subset.index, '.fasta', sep='')
  aln.write(aln.noblanks, aln.f)
  rm(aln.noblanks)

  aln.nogaps <- aln.remove.ngap.columns(aln, n=5)
  dir.create(file.path(out.f, 'concat_nongaps_05'))
  aln.f <- paste(out.f, '/concat_nongaps_05/alns_', subset.index, '.fasta', sep='')
  aln.write(aln.nogaps, aln.f)
  rm(aln.nogaps)

  aln.nogaps <- aln.remove.ngap.columns(aln, n=10)
  dir.create(file.path(out.f, 'concat_nongaps_10'))
  aln.f <- paste(out.f, '/concat_nongaps_10/alns_', subset.index, '.fasta', sep='')
  aln.write(aln.nogaps, aln.f)
  rm(aln.nogaps)

  aln.nogaps <- aln.remove.ngap.columns(aln, n=20)
  dir.create(file.path(out.f, 'concat_nongaps_20'))
  aln.f <- paste(out.f, '/concat_nongaps_20/alns_', subset.index, '.fasta', sep='')
  aln.write(aln.nogaps, aln.f)
  rm(aln.nogaps)
}

# Should take about 15 minutes for each job...
bsub.combine.concat.alns <- function() {

  subdirs <- get.output.subdirs() 

  for (remove.paralogs in c(TRUE, FALSE)) {
    for (mask.clusters in c(TRUE, FALSE)) {
      for (mask.nongaps in c(TRUE)) {
        for (subdir in subdirs) {
          xtra <- paste(remove.paralogs, mask.clusters, mask.nongaps, subdir, sep=' ')
          bsub.function('combine_concat_alns', mem=12, extra.args=xtra)
        }
      }
    }
  }
}

combine.concat.alns <- function(remove.paralogs=T, mask.clusters=T, mask.nongaps=T, subdir) {
  library(rphast)

  out.f <- output.folder(remove.paralogs, mask.clusters, mask.nongaps)
  out.f <- paste(out.f, '/', subdir, sep='')
  print(out.f)

  aln.files <- Sys.glob(paste(out.f, '/alns*.fasta', sep=''))
  #aln.files <- aln.files[1:2]
  print(sprintf("%d files to combine in %s", length(aln.files), subdir))

  aln.list <- list()
  for (i in 1:length(aln.files)) {
    aln.f <- aln.files[i]
    print(aln.f)
    cur.msa <- read.msa(file=aln.f, pointer.only=T, tuple.size=3)
    #print(cur.msa)
    aln.list[[i]] <- cur.msa
  }

  print("concatenating...")
  all.msa <- concat.msa(aln.list, pointer.only=T)
  aln.f <- paste(out.f, '/', 'all_combined.fasta', sep='')
  print(sprintf("writing %s", aln.f))
  aln.tmp <- tempfile()
  write.msa(all.msa, file=aln.tmp)
  file.copy(aln.tmp, aln.f, overwrite=T)
}


old.combine.concat.alns <- function(remove.paralogs=T, mask.clusters=T, mask.nongaps=T, subdir) {
  out.f <- output.folder(remove.paralogs, mask.clusters, mask.nongaps)
  out.f <- paste(out.f, '/', subdir, sep='')
  print(out.f)

  aln.files <- Sys.glob(paste(out.f, '/', subdir, '/alns*.fasta', sep=''))
  print(sprintf("%d files to combine in %s", length(aln.files), subdir))

  aln <- NULL 
  for (aln.f in aln.files) {
    print(aln.f)
    print("reading...")
    cur.aln <- s.aln.read(aln.f)
    print(Sys.time())

    if (is.null(aln)) {
      aln <- cur.aln
    } else {
      print("concat...")
      aln <- s.aln.concat(aln, cur.aln)
      print(Sys.time())
    }
  }

  aln.f <- paste(out.f, '/', 'all_combined.fasta', sep='')
  print(sprintf("writing %s", aln.f))

  aln.tmp <- tempfile()
  s.aln.write(aln, aln.tmp)
  file.copy(aln.tmp, aln.f, overwrite=T)

  rm(cur.aln)
  rm(aln)
  print("  done!")
}

