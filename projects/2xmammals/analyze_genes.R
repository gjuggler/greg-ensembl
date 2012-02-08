uname  <- Sys.getenv("USER")
if (uname == 'gj1') {
  source("~/src/greg-ensembl/scripts/liftOver.R")
  source("~/src/greg-ensembl/projects/phylosim/PhyloSimPlots.R")
} else {
#  source("~/lib/greg-ensembl/scripts/go_enrichments.R")
  source("~/lib/greg-ensembl/scripts/liftOver.R")
  source("~/lib/greg-ensembl/projects/phylosim/PhyloSimPlots.R")
}

write.aln <- function(x, ff) {
  # Mostly copied from ape, but I removed the space between the Fasta
  # '>' and the sequence ID.
  zz <- file(ff, "w")
  on.exit(close(zz))
  colw <- 120
  colsep <- ''
  indent <- ''
  nbcol <- 1
  N <- dim(x)
  S <- N[2]
  N <- N[1]
  xx <- vector("list", N)
  for (i in 1:N) xx[[i]] <- x[i, ]
  names(xx) <- rownames(x)
  x <- xx
  rm(xx)
  for (i in 1:N) {
    cat(">", names(x)[i], file = zz, sep='')
    cat("\n", file = zz, sep='')
    X <- paste(x[[i]], collapse = "")
    S <- length(x[[i]])
    totalcol <- ceiling(S/colw)
    if (nbcol < 0) 
      nbcol <- totalcol
    nb.lines <- ceiling(totalcol/nbcol)
    SEQ <- character(totalcol)
    for (j in 1:totalcol) {
      SEQ[j] <- substr(X, 1 + (j - 1) * colw, colw + (j - 1) * colw)
    }
    for (k in 1:nb.lines) {
      endsel <- if (k == nb.lines) 
        length(SEQ)
        else nbcol + (k - 1) * nbcol
      cat(indent, file = zz)
      cat(SEQ[(1 + (k - 1) * nbcol):endsel], sep = colsep, file = zz)
      cat("\n", file = zz)
    }
  }
}

save.alns <- function() {
  con <- connect(dbname())
  cmd <- "select * from genes where ref_taxon_id=9606 and chr_name like '%chr%';"
  df <- dbGetQuery(con, cmd)
  disconnect(con)

  print(table(df$chr_name))
  return()

  for (i in 1:nrow(df)) {
    con <- connect(dbname())
    data_id <- df[i, 'data_id']
    cmd <- sprintf("select * from seqs where data_id=%d order by aln_position",
      data_id)
    aln.df <- dbGetQuery(con, cmd)
    disconnect(con)

    print(paste(i, nrow(aln.df)))
  }
}

process.seqs.df <- function(seqs.df) {
  cols <- colnames(seqs.df)
  cols <- setdiff(cols, c('data_id', 'aln_position'))  

  for (i in 1:length(cols)) {
    cc <- as.character(cols[i])
    seqs.df[,cc] <- ifelse(seqs.df[, cc] == '---', NA, seqs.df[, cc])
    seqs.df[,cc] <- as.factor(seqs.df[,cc])
  }
  seqs.df
}

get.subset.newicks.export <- function() {
  species.subsets
}

get.cluster.windows.export <- function(d.id) {
  win.df <- subset(windows.df, data_id==d.id)
  win.df
}

get.genes.export <- function() {
  genes.df
}

get.taxid.df.export <- function() {
  taxid.df
}

get.aln.df.export <- function(d.id) {
  aln.df <- subset(seqs.df, data_id==d.id)
  aln.df
}

is.export <- function() {
  if (exists('ALN_EXPORT')) {
    exprt <- get('ALN_EXPORT')
    if (exprt) {
      return(TRUE)
    }
  }
  return(FALSE)
}

get.aln.df <- function(d.id) {
  if (is.export()) {
    aln.df <- subset(seqs.df, data_id==d.id)
    return(aln.df)
  }
  con <- connect(dbname())
  cmd <- sprintf("select * from seqs where data_id=%d order by aln_position",
    d.id)
  aln.df <- dbGetQuery(con, cmd)
  disconnect(con)
  aln.df
}

get.all.baddies <- function() {
  con <- connect(dbname())
  cmd <- sprintf("select * from win_baddies")
  win.df <- dbGetQuery(con, cmd)
  disconnect(con)
  win.df
}

get.cluster.windows <- function(data_id) {
  con <- connect(dbname())
  cmd <- sprintf("select * from win_baddies where data_id=%d",
    data_id)
  win.df <- dbGetQuery(con, cmd)
  disconnect(con)
  win.df
}

get.aln <- function(data.id,
  aln.pos=NULL,
  remove.paralogs=T,
  mask.clusters=T,
  mask.nongaps=T,
  mask.nongap.threshold=2,
  remove.short.seqs=T,
  keep.species = 'Mammals'
) {
  genes <- get.genes()
  
  cur.gene <- subset(genes, data_id==data.id)
  if (nrow(cur.gene) == 0) {
    # Try finding by gene_name
    cur.gene <- subset(genes, gene_name==data.id)
  }
  if (nrow(cur.gene) == 0) {
    stop(sprintf("Cannot find gene for ID %s!", data.id))
  }

  id <- cur.gene$gene_name

  aln.df <- get.aln.df(cur.gene$data_id)
  taxid.df <- get.taxid.df()
  
  if (!is.null(aln.pos)) {
    aln.df <- subset(aln.df, aln_position %in% aln.pos)
  }

  if (mask.nongaps) {
    #print("mask nongaps...")
    # TODO: use non-gap BL to mask out sites from the alignment.
    #print("get sitewise...")
    sites <- get.sitewise()
    sites <- subset(sites, data_id == data.id)
    #print(sites[1:5,])
    remove.sites <- subset(sites, nongap_bl < mask.nongap.threshold)
    remove.positions <- remove.sites$aln_position

    if (length(remove.positions) > 0) {
      #print(remove.positions)
      species <- colnames(aln.df)
      #print(species)
      species <- setdiff(species, c('data_id', 'aln_position'))  
      #print(species)
      # Mask out bad columns with Ns (or leave as gaps)
      #print("Masking bad cols")
      aln.df[remove.positions, species] <- ifelse(aln.df[remove.positions, species] == '---', '---', 'NNN')
      #print("done!")
      warning(sprintf("Removed %d gappy columns from %s", length(remove.positions), id))
    }
  }

  if (remove.short.seqs) {
    cols <- colnames(aln.df)
    cols <- setdiff(cols, c('data_id', 'aln_position'))  
    seq.df <- subset(aln.df, select=cols)

    len.by.row <- apply(seq.df, 2, function(x) {
      sum(x != '---' & x != 'NNN')
    })
    nonzero.lengths <- len.by.row[len.by.row != 0]
    mean.nonzero <- mean(nonzero.lengths)

    # Remove sequences with less than half the mean length of non-NNN seqs.
    short.genes <- which(len.by.row <= mean.nonzero/2 && len.by.row > 0)
    short.tx <- names(len.by.row)[short.genes]
    #print(short.tx)
    if (length(short.tx) > 0) {
      for (j in 1:length(short.tx)) {
        i2 <- as.character(short.tx[j])
        aln.df[, i2] <- ifelse(aln.df[, i2] == '---', '---', '---')
      }
      warning(sprintf("Removed %d short sequences from %s", length(short.tx), id))
    }
  }

  if (remove.paralogs) {
    paralog.s <- cur.gene$dup_species_list
    if (!is.na(paralog.s)) {
      paralog.tx <- unlist(strsplit(paralog.s, ' '))
      for (j in 1:length(paralog.tx)) {
        i2 <- as.character(paralog.tx[j])
        aln.df[, i2] <- ifelse(aln.df[, i2] == '---', '---', '---')
      }
      warning(sprintf("Removed %d paralogs from %s", length(paralog.tx), id))
    }
  }

  if (mask.clusters) {
    win.df <- get.cluster.windows(cur.gene$data_id)
    if (nrow(win.df) > 0) {
      for (i in 1:nrow(win.df)) {
        x <- win.df[i,]
        cur.tx <- x$taxon_id
        #print("  getting beneath..")
        #print(cur.tx)
        spec.beneath <- taxids.beneath.node(cur.tx)
        #print("  got em!")
        spec.tree <- get.species.tree()
        n.species <- spec.tree$Nnode

        # Only remove windows if the number of nsyn subs is greater than
        # half the window size, and the number of species beneath the bad
        # node is less than half the total tree
        if (length(spec.beneath ) <= n.species/2 && x$n_ns_subs <= x$win_size) {

          tx.lbl <- taxid.to.alias2(cur.tx, taxid.df)
          #print(sprintf("Masking window from %s", tx.lbl))
      
          aln.lo <- x$aln_position
          aln.hi <- x$aln_position + x$win_size

          for (j in 1:length(spec.beneath)) {
            i1 <- aln.lo:aln.hi
            i2 <- as.character(spec.beneath[j])
            # Mask out bad regions with Ns, or leave as gaps
            #print(aln.df[i1, i2])
            aln.df[i1, i2] <- ifelse(aln.df[i1, i2] == '---', '---', 'NNN')
            #print(aln.df[i1, i2])
          }
        }
      }
      warning(sprintf("Removed %d substitution clusters from %s", length(win.df), id))
    }
  }

  aln <- df.to.aln(aln.df)

  # Remove species if necessary.
  if (keep.species != 'all') {
    subset.newicks <- get.subset.newicks()
    cur.ind <- which(subset.newicks$name == keep.species)
    cur.subset <- subset.newicks[cur.ind,]
    taxid.s <- cur.subset$taxid_list
    cur.taxids <- strsplit(taxid.s, ' ')[[1]]
    aln <- restrict.aln.to.seqs(aln, cur.taxids)
  }

  # Sort aln by tree.
  tree <- get.species.tree()
  tree <- remove.branchlengths(tree)
  tree <- extract.subtree(tree, rownames(aln))
  aln <- sort.aln.by.tree(aln, tree)

  # Translate the aln taxon ids into names.
  tx.ids <- rownames(aln)
  species.ids <- taxid.to.alias2(tx.ids, taxid.df)
  rownames(aln) <- species.ids
  aln
}

df.to.aln <- function(aln.df) {
  aln.positions <- aln.df$aln_position
  cols <- colnames(aln.df)
  cols <- setdiff(cols, c('data_id', 'aln_position'))  
  seq.df <- subset(aln.df, select=cols)
  rownames(seq.df) <- aln.positions
  aln <- t(as.matrix(seq.df))
  if (any(is.na(aln))) {
    aln <- ifelse(is.na(aln), '---', aln)
  }
  aln
}

aln.stats <- function(aln) {
  n.seqs <- length(aln[,1])
  n.cols <- length(aln[1,])

  nongap.counts <- apply(aln, 1, function(x) {
    sum(!grepl('-', x))
  })
  nz.counts <- nongap.counts != 0
  n.nonblank.seqs <- sum(nz.counts)
  mean.seq.length <- mean(nongap.counts[nz.counts])

  nuc.counts <- apply(aln, 1, function(x) {
    str <- paste(x, collapse='')
    mtchs <- unlist(gregexpr('[GACT]', str))
    sum(mtchs != -1)
  })
  total.nucs <- sum(nuc.counts)

  nongap.cols <- apply(aln, 2, function(x) {
    sum(!grepl('(---|-)', x))
  })
  n.allgaps <- sum(nongap.cols == 0)
  n.nogaps <- sum(nongap.cols == (n.nonblank.seqs))

  percent.gap.nucs <- (n.seqs * n.cols - total.nucs/3) / (n.seqs * n.cols)
  percent.allgap.columns <- n.allgaps / n.cols
  percent.nogap.columns <- n.nogaps / n.cols

  str <- sprintf('Seqs: %d
Columns: %d
All-gap columns: %d
Mean seq length: %.1f
Total nucs: %d\n',
    n.seqs,
    n.cols,
    n.allgaps,
    mean.seq.length,
    total.nucs
  )
  #cat(str)

  out.df <- data.frame(
    n_species = n.seqs,
    n_nonblank_seqs = n.nonblank.seqs,
    mean_seq_length = mean.seq.length,
    n_columns = n.cols,
    n_allgap_columns = n.allgaps,
    n_nogap_columns = n.nogaps,
    n_nucs = total.nucs,
    percent_allgap_columns = percent.allgap.columns,
    percent_nogap_columns = percent.nogap.columns,
    percent_gap_nucs = percent.gap.nucs
  )
  out.df <- format(out.df, digits=1)
  out.df
}

get.sitewise.export <- function() {
  #print("Getting sitewise export")
  sitewise.df
}

get.sitewise <- function() {
  #print("Getting full sitewise data")
  sites <- get.pset.sites(6, filter='default', test=F)
}

get.species.tree.export <- function() {
  species.tree
}

get.species.tree <- function() {
  f.con <- file("~/src/greg-ensembl/projects/orthologs/compara_63_taxids.nh")
  str <- readLines(con=f.con)
  close(f.con)
  tree <- read.nhx.tree(str)
  tree
}

get.taxa.df <- function() {
  con <- connect.livemirror('ensembl_compara_64')
  cmd <- sprintf('select * from genome_db where name != "ancestral_sequences"')
  tx.df <- dbGetQuery(con, cmd)
  disconnect(con)

  tx.df$alias <- as.character(taxid.to.alias(tx.df$taxon_id))

  tx.df <- subset(tx.df, select=c('taxon_id', 'alias', 'name', 'assembly'))
  tx.df
}

get.subset.newicks <- function() {
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
  pset.ids <- as.integer(names(taxids.list))
  taxid.names = pset.to.alias(names(taxids.list))
  names(taxids.list) <- taxid.names
  
  newick.df <- ldply(taxids.list, function(x) {
    tree <- get.species.tree()
    tree.with.outgroup <- tree

    # Drop tips falling outside the taxid list.
    not.in.set <- setdiff(tree$tip.label, x)
    tree <- drop.tip(tree, not.in.set)
    tree2 <- tree
    tree3 <- tree

    # Include a distant outgroup...
    taxid.df <- get.taxid.df()
    taxid.df <- subset(taxid.df, name == 'Platypus')
    outgroup.taxid <- taxid.df$taxon_id
    xx <- unique(c(x, outgroup.taxid))
    not.in.set <- setdiff(tree.with.outgroup$tip.label, xx)
    tree.with.outgroup <- drop.tip(tree.with.outgroup, not.in.set)
    str.with.outgroup <- write.tree(tree.with.outgroup)

    tree$edge.length <- NULL
    tree$node.label <- NULL
    str <- write.tree(tree)

    tree2$tip.label <- taxid.to.alias(tree2$tip.label, include.internals=T)
    tree2$edge.length <- NULL
    tree2$node.label <- NULL
    str2 <- write.tree(tree2)

    str3 <- write.tree(tree3)

    data.frame(
      taxid_newick=str,
      label_newick=str2,
      taxid_newick_with_internals=str3,
      taxid_newick_with_outgroup=str.with.outgroup,
      taxid_list=paste(x, collapse=' '),
      stringsAsFactors=F
    )
  })
  newick.df$.id <- NULL
  newick.df$name <- taxid.names
  newick.df$pset_id = pset.ids

  newick.df
}

tree.remove.internal.labels <- function(tree) {
  tree$node.label
}

tree.taxids.to.labels <- function(tree, space.to.underscore=T) {
  tree$tip.label <- taxid.to.alias(tree$tip.label, include.internals=T)
  if (any(is.character(tree$node.label))) {
    tree$node.label <- taxid.to.alias(tree$node.label, include.internals=T)
  }

  if (space.to.underscore) {
    tree$tip.label <- gsub(' ', '_', tree$tip.label)
    tree$node.label <- gsub(' ', '_', tree$node.label)
  }

  tree
}

species.group.to.pset <- function(species.group) {
  x <- pset.df()
  x <- subset(x, label == species.group)
  x$pset_id
}

get.subset.taxid.tree <- function(pset.id, include.outgroup=F) {
  subset.df <- get.subset.newicks()
  #print(subset.df)
  cur.row <- subset(subset.df, pset_id==pset.id)
  if (include.outgroup) {
    tree.str <- cur.row$taxid_newick_with_outgroup
  } else {
    tree.str <- cur.row$taxid_newick_with_internals
  }
  tree <- read.tree(text=tree.str)
  tree
}

taxids.beneath.node <- function(taxid) {
  require(ape)

  tree <- get.species.tree()
  #print(str(tree))
  cur.node <- node.with.label(tree, sprintf("%d", taxid))
  #print(cur.node)
  if (cur.node > tree$Nnode) {
    nongap.tree <- extract.clade(tree, cur.node)
    nongap.tree$tip.label
  } else {
    taxid
  }
}

plot.alns <- function(gene.names, ...) {
  for (i in 1:length(gene.names)) {
    cur.gene <- gene.names[i]
    filename <- paste('aln_', cur.gene, '.pdf', sep='')
    plot.aln(gene_name=cur.gene, filename=filename, ...)
  }
}

plot.aln <- function(
  data_id, 
  aln_lo=1, 
  aln_hi=99999, 
  gene_name=NULL,
  filename="test.pdf", 
  taxon_id=9606,
  include.psets = c(1, 2, 3, 4, 5, 6),
  include.subs = F,
  keep.species = 'mammals',
  remove.blank.columns=F,
  flatten.seq=F,
  plot.chars=T,
  plot.protein=F
) {
  library(R.oo)
  library(phylosim)
  source("~/src/greg-ensembl/projects/phylosim/PhyloSimPlots.R")

  if (!is.null(gene_name)) {
    con <- connect(dbname())
    cmd <- sprintf("select * from genes where gene_name='%s'", gene_name)
    df <- dbGetQuery(con, cmd)
    disconnect(con)
    if (nrow(df) == 0) {
      print(sprintf("Gene with name '%s' not found!", gene_name))
      return()
    }
    data_id <- df$data_id
  }

  con <- connect(dbname())
  cmd <- sprintf("select * from seqs where data_id=%d and aln_position between %d and %d order by aln_position",
    data_id, aln_lo, aln_hi)
  df <- dbGetQuery(con, cmd)
  disconnect(con)

  if (include.subs) {
    con <- connect(dbname())
    cmd <- sprintf("select * from subs where data_id=%d and aln_position between %d and %d order by aln_position",
      data_id, aln_lo, aln_hi)
    subs.df <- dbGetQuery(con, cmd)
    disconnect(con)
  }

  con <- connect(dbname())
  cmd <- sprintf("select * from genes where data_id=%d", data_id)
  gene <- dbGetQuery(con, cmd)
  disconnect(con)

  print(sprintf("Gene: %s Align length: %d dN/dS: %.3f", gene$gene_name, gene$aln_length, gene$m_slr_dnds))

  # Get the gene name, MPL, and mammal dN/dS
  gene.name <- gene$gene_name
  gene.mpl <- gene$m_slr_mean_path
  gene.dnds <- gene$m_slr_dnds
  dup.spc <- gene$dup_species_count
  #print(gene)
  gene.xlab <- sprintf("%s Alignment Position (dN/dS=%.2f, MPL=%.2f, Dup. Spec.=%d)", gene.name, gene.dnds, gene.mpl, dup.spc)

  # Get the sites values.
  con <- connect(dbname())
  min.aln.pos <- min(df$aln_position)
  max.aln.pos <- max(df$aln_position)
  cmd <- sprintf("select * from sites where data_id=%d and aln_position between %d and %d", data_id, min.aln.pos, max.aln.pos)
  sites <- dbGetQuery(con, cmd)
  disconnect(con)
  sites <- subset(sites, parameter_set_id %in% include.psets)
  tracks <- NA
  if (nrow(sites) > 0) {
    sites$source <- paste("PNSC", pset.to.alias(sites$parameter_set_id))
    lrt.score <- rep(0, nrow(sites))
    lrt.score <- ifelse(sites$lrt_stat > qchisq(0.95, df=1), 1, lrt.score)
    lrt.score <- ifelse(sites$lrt_stat < -qchisq(0.95, df=1), -1, lrt.score)
    lrt.score <- ifelse(sites$lrt_stat > qchisq(0.99, df=1), 2, lrt.score)
    lrt.score <- ifelse(sites$lrt_stat < -qchisq(0.99, df=1), -2, lrt.score)
    sites$score <- lrt.score
    #sites$score <- sites$lrt_stat

    sites$start <- sites$aln_position - aln_lo + 1
    sites$end <- sites$aln_position - aln_lo + 1

    dnds <- subset(sites, parameter_set_id == 6)
    include.dnds <- F
    if (nrow(dnds) > 0) {
      if (include.dnds) {
        dnds$source <- paste("dN/dS", pset.to.alias(dnds$parameter_set_id))
        dnds$score <- dnds$omega
        sites <- rbind(sites, dnds)
      }

      # Add Pfam annotations.
      pf.sites <- subset(dnds, !is.na(pfam_domain))
      if (nrow(pf.sites) > 0) {
        pfam_id <- as.factor(pf.sites$pfam_domain)
        pf.sites$source <- "Pfam"
        pf.sites$score <- as.integer(pfam_id)
        sites <- rbind(sites, pf.sites)
      }

      # Add splice distances.
      splice.dists <- subset(dnds, !is.na(splice_distance))
      if (nrow(splice.dists) > 0) {
        splice.dists$source <- "Exon Border"
        splice.dists$score <- pmin(2, 1 / (abs(splice.dists$splice_distance)))
        sites <- rbind(sites, splice.dists)
      }

    }

    tracks <- gff.to.tracks(sites)
  }

  aln.positions <- df$aln_position
  cols <- colnames(df)
  cols <- setdiff(cols, c('data_id', 'aln_position'))
  
  seq.df <- subset(df, select=cols)
  rownames(seq.df) <- aln.positions
  aln <- t(as.matrix(seq.df))

  if (keep.species == 'slim_mammals') {
    mamms <- slim.mammals()
  } else if (keep.species == 'mammals') {
    mamms <- mammals()
  } else if (keep.species == 'primates') {
    mamms <- primates()
  }
  mamms <- c(mamms, taxon_id)
  mamms <- unique(mamms)
  aln <- restrict.aln.to.seqs(aln, mamms)

  f.con <- file("~/src/greg-ensembl/projects/orthologs/compara_63_taxids.nh")
  str <- readLines(con=f.con)
  close(f.con)
  tree <- read.nhx.tree(str)
  non.mammals <- setdiff(tree$tip.label, mamms)
  tree <- drop.tip(tree, non.mammals)
  tree <- remove.branchlengths(tree)
  aln <- sort.aln.by.tree(aln, tree)

  aln2 <- aln
  aln2 <- aln.tx(aln2)
  rownames(aln2) <- paste(gene.name, rownames(aln2), sep='_')
  out.f <- paste('aln_', gene.name, '.fasta', sep='')
#  write.seqs(aln2, out.f)
  #write.dna(aln2, file=out.f, format='fasta', colsep='')

#  if (remove.blank.columns) {
#    aln <- remove.blank.columns(aln)
#  }

  tx.ids <- rownames(aln)
  species.ids <- taxid.to.alias(tx.ids, include.internals=T)
  rownames(aln) <- species.ids
  tree$tip.label <- taxid.to.alias(tree$tip.label, include.internals=T)

  # Flatten to human sequence.
  if (flatten.seq) {
    h.gaps <- columns.without.sequence(aln, 'Human')

    xpos <- 1:aln.length(aln)
    xx <- rep(FALSE, aln.length(aln))
    xx[h.gaps] <- TRUE
    xx.cum <- cumsum(xx)
    xx.df <- data.frame(
      pos = xpos,
      offset = xx.cum
    )

    aln <- remove.columns(aln, h.gaps)
    tracks <- llply(tracks, function(x) {
      x <- subset(x, !(pos %in% h.gaps))
      x <- merge(x, xx.df, all.x=T)
      x$pos <- x$pos - x$offset
      x$offset <- NULL
      x
    })
  }

  if (plot.protein) {
    aln <- aln.tx(aln)
  }

  out.w <- .2 * length(aln[1,])
  out.h <- .2 * (length(rownames(aln)) + 3 + length(tracks))

  pdf(file=filename, width=out.w, height=out.h)
  p <- aln.plot(aln, 
    tree=tree,
    color.scheme='auto',
    aln.plot.chars=plot.chars,
    aln.char.text.size=1.5,
    aln.char.alpha = 0.7,
    aln.pos.offset=aln_lo-1,
    aln.plot.xlabel = gene.xlab,
    aln.tracks = tracks,
    tree.labels.in.tree = F,
    tree.xlim.expand = c(0.05, 0.05),
    tree.xlab = '',
    tree.plot.legend = F,
    tree.show.x.axis = T,
    tree.x.axis.color = 'white',
    line.width = 1
  )
  grid.draw(p$grob)
  dev.off()
}

write.seqs <- function(aln, out.f) {
  zz <- file(out.f, "w")  # open an output file connection

  seq.nms <- rownames(aln)
  for (i in 1:length(seq.nms)) {
    cur.row <- aln[i, ]
    cur.row <- paste(cur.row, collapse='')
    #print(paste(cur.row, collapse=''))
    cur.row <- gsub('-', '', cur.row)
    if (nchar(cur.row) == 0) {
      next()
    }
    #print(cur.row)
    cat(paste(">", seq.nms[i], "\n", sep=''), file=zz)
    cat(paste(cur.row, "\n", sep=''), file=zz)
  }
  close(zz)
}

example.alns <- function() {
  plot.aln(1514564, 1, 300, filename="aln_ex_x.pdf", plot.chars=F)
  plot.aln(1640146, 1000, 1500, filename="aln_ex_rhesus.pdf", plot.chars=F)
  # Chimpanzee crap:
  #get.aln(1592585, 700, 750, filename="aln_ex_1.pdf")
  # Gorilla crap:
  plot.aln(gene_name="G3BP1", 1, 100, filename="aln_ex_2.pdf", plot.chars=F)
  # Human crap:
  #get.aln(gene_name="TPM1", filename="aln_ex_3.pdf", flatten.seq=T)
  #get.aln(gene_name="TPM3", filename="aln_ex_3b.pdf", flatten.seq=T)
  # Platypus ain't so bad:
  #get.aln(1898971, 560, 610, filename="aln_ex_4.pdf")  
  # Mouse/rat crap:
  plot.aln(1830947, 370, 450, filename="aln_ex_5.pdf")
}

win.baddies <- function() {

  # Look at the top baddie window from each primate species.
  # Chimp:
  get.aln(1418043, 200, 350, "baddie_chimp.pdf", taxon_id=9598)

  # Human crap, affects primate pos-sel:
  get.aln(1591434, 200, 300, "baddie_human.pdf", taxon_id=9606)  

  # Gorilla crap, but no effect on primate pos-sel:
  get.aln(1653419, 700, 800, "baddie_gorilla_1.pdf", taxon_id=9593)
  get.aln(1833054, 380, 460, "baddie_gorilla_2.pdf", taxon_id=9593)

  # Macaque crap -- looks kind of reasonable...
  get.aln(882461, 100, 200, "baddie_macaque.pdf", taxon_id=9544)

  # Giant panda
  get.aln(822393, 1, 100, "baddie_panda.pdf", taxon_id=9646)
}

get.domains <- function() {
  out.f <- scratch.f("domains.sites.Rdata")
  if (!file.exists(out.f)) {
    con <- connect(dbname())
    cmd <- sprintf("select * from domains_sites")
    domains <- dbGetQuery(con, cmd)
    disconnect(con)
    save(domains, file=out.f)
  }
  load(out.f)
  domains
}

get.domains.merged <- function(filter='default', keep.once=c(), keep.fields=c('hoch_bh'), keep.info=T, do.union=T) {
  domains <- get.domains.sites(filter=filter)

  if (keep.info) {
    always.keep <- c('pfam_domain', 'pfam_name', 'pfam_length', 'pfam_type', 'pfam_descr')
  } else {
    always.keep <- c('pfam_domain')
  }
  keep.fields <- unique(c(keep.fields))

  psets <- sort(unique(domains$pset))
  out.df <- subset(domains, select=unique(c(always.keep, keep.once)))
  out.df <- out.df[!duplicated(out.df$pfam_domain),]

  for (i in 1:length(psets)) {
    #print(i)
    next.sub <- subset(domains, pset == psets[i], select=c('pfam_domain', keep.fields))
    
    colnames(next.sub) <- c('pfam_domain', paste(keep.fields, '.', psets[i], sep=''))
    suffx <- paste('.', psets[i], sep='')
    out.df <- merge(out.df, next.sub, by=c('pfam_domain'), all.x=do.union)
  }
  out.df
}

get.all.merged <- function(filter='default', keep.once=c('gene_name'), keep.fields=c('hoch_bh'), do.union=T) {
  genes <- get.genes.sites(filter=filter)

  keep.fields <- unique(c(keep.fields))

  psets <- sort(unique(genes$pset))
  out.df <- subset(genes, select=unique(c('data_id', keep.once)))
  out.df <- out.df[!duplicated(out.df$data_id),]

  for (i in 1:length(psets)) {
    #print(i)
    next.sub <- subset(genes, pset == psets[i], select=c('data_id', keep.fields))
    
    colnames(next.sub) <- c('data_id', paste(keep.fields, '.', psets[i], sep=''))
    suffx <- paste('.', psets[i], sep='')
    out.df <- merge(out.df, next.sub, by=c('data_id'), all.x=do.union)
  }
  out.df
}

get.domains.sites <- function(filter='default') {
  cache.df <- sprintf("domains.sites.%s", filter)
  if (exists(cache.df, envir=.GlobalEnv)) {
    domains <- get(cache.df, envir=.GlobalEnv)
    return(domains)
  }

  cache.f <- scratch.f(sprintf("domains.sites.%s.Rdata", filter))
  if (!file.exists(cache.f)) {
    print("Domains.sites doesn't exist -- re-fetching...")
    con <- connect(dbname())
    cmd <- sprintf("select * from domain_sites where filter='%s'", filter)
    domains <- dbGetQuery(con, cmd)
    disconnect(con)

    print("Touching up domains...")
    # Do some ddply touch-ups... fix the z-scores, do FDR correction, etc. etc.
    domains <- ddply(domains, .(pset, filter), function(x) {
      x <- subset(x, n_sites > 200)

      # Take the PSG p-values (TPM, empirical, Fisher) and do some FDR control
      x$hoch_bh <- p.adjust(x$hoch_p, method='BH')
      x$emp_005_bh <- p.adjust(x$emp_005, method='BH')
      x$emp_01_bh <- p.adjust(x$emp_01, method='BH')
      x$emp_05_bh <- p.adjust(x$emp_05, method='BH')
      x$emp_10_bh <- p.adjust(x$emp_10, method='BH')

      print(sprintf("%d", x[1, 'pset']))
      x
    })

    save(domains, file=cache.f)
  }
  #print("Loading cache...")
  load(cache.f)

  pos.gene.gaps <- gsub('\\S+', 'X', domains$pos_gene_names)
  pos.gene.gaps <- gsub('[^X]', '', pos.gene.gaps)
  n.pos.genes <- nchar(pos.gene.gaps)
  domains$n.pos.genes <- n.pos.genes

  assign(cache.df, domains, envir=.GlobalEnv)
  domains
}


get.genes.sites <- function(filter='default') {
  cache.df <- sprintf("genes.sites.%s", filter)
  if (exists(cache.df, envir=.GlobalEnv)) {
    genes <- get(cache.df, envir=.GlobalEnv)
    return(genes)
  }

  cache.f <- scratch.f(sprintf("genes.sites.%s.Rdata", filter))
  if (!file.exists(cache.f)) {
    print("Genes.sites doesn't exist -- re-fetching...")
    con <- connect(dbname())
    cmd <- sprintf("select * from genes_sites where filter='%s'", filter)
    genes <- dbGetQuery(con, cmd)
    disconnect(con)

    print("Touching up genes...")
    # Do some ddply touch-ups... fix the z-scores, do FDR correction, etc. etc.
    genes <- ddply(genes, .(parameter_set_id), function(x) {
    
      # Take the PSG p-values (TPM, empirical, Fisher) and do some FDR control
      x$hoch_bh <- p.adjust(x$hoch_p, method='BH')
      x$tpm_05_bh <- p.adjust(x$tpm_05, method='BH')
      x$tpm_10_bh <- p.adjust(x$tpm_10, method='BH')
      x$tpm_20_bh <- p.adjust(x$tpm_20, method='BH')
      x$tpm_50_bh <- p.adjust(x$tpm_50, method='BH')
      x$emp_005_bh <- p.adjust(x$emp_005, method='BH')
      x$emp_01_bh <- p.adjust(x$emp_01, method='BH')
      x$emp_05_bh <- p.adjust(x$emp_05, method='BH')
      x$emp_10_bh <- p.adjust(x$emp_10, method='BH')
      x$fis_bh <- p.adjust(x$fis_p, method='BH')

      print(sprintf("%d", x[1, 'pset']))
      x
    })

    save(genes, file=cache.f)
  }
  #print("Loading cache...")
  load(cache.f)

  # Require at least 50 sites per gene.
  genes <- subset(genes, n_sites > 50)

  assign(cache.df, genes, envir=.GlobalEnv)
  genes
}

filter.nsites <- function(genes) {
  subset(genes, n_sites > 30)
}

plot.scatters <- function() {
  plot.dnds.scatter(1, 2, 'default')
  plot.dnds.scatter(1, 3, 'default')
  plot.dnds.scatter(1, 4, 'default')
  plot.dnds.scatter(1, 6, 'default')
}

plot.dnds.scatter <- function(pset_a, pset_b, filter='default') {
  genes <- get.genes.sites(filter=filter)

  keep.cols <- c('gene_name', 'data_id', 'slr_dnds', 'slr_dnds_q', 'slr_bl', 'slr_dnds_z')
  genes.a <- subset(genes, pset == pset_a, select=keep.cols)
  genes.b <- subset(genes, pset == pset_b, select=keep.cols)
  mrgd <- merge(genes.a, genes.b, by='data_id')
  
  #mrgd$slr_bl.x <- pmin(quantile(mrgd$slr_bl.x, 0.9), mrgd$slr_bl.x)
  #mrgd$slr_bl.x <- pmax(quantile(mrgd$slr_bl.x, 0.1), mrgd$slr_bl.x)

  p <- ggplot(mrgd, aes(x=slr_dnds_z.x, y=slr_dnds_z.y))
  p <- p + geom_abline(colour='black', slope=1, intercept=0, linetype='dashed')
  p <- p + geom_point(size=0.5, alpha=1)
  lmt <- c(-4.5, 2.5)
  p <- p + scale_x_continuous(paste(pset.to.alias(pset_a), "dN/dS Z-Score"), limits=lmt)
  p <- p + scale_y_continuous(paste(pset.to.alias(pset_b), "dN/dS Z-Score"), limits=lmt)
  p <- generic.opts(p)

  out.f <- scratch.f(sprintf("genes_dnds_scatter_%d_%d_%s.pdf", pset_a, pset_b, filter))
  pdf(file=out.f)
  print.ggplot(p)
  dev.off()
}

# Simple: do a histogram of dN/dS (or log dN/dS) for each subgroup.
plot.dnds.hist <- function(filter='default') {
  genes <- get.genes.sites(filter=filter)
  genes$pset <- pset.to.alias(genes$parameter_set_id, factors=T, prefix=' ')

  genes$slr_dnds <- pmin(genes$slr_dnds, 1)
  genes$slr_dnds <- pmax(genes$slr_dnds, 0.001)

  n.bins <- 50
  bw <- diff(range(genes$slr_dnds)) / n.bins

  p <- ggplot(genes, aes(x=slr_dnds))
  p <- p + theme_bw()
  p <- p + geom_histogram(binwidth=bw)
#  p <- p + scale_x_log10()
  p <- p + facet_grid(pset ~ .)
  p <- generic.opts(p)

  out.f <- scratch.f(sprintf("genes_dnds_hist_%s.pdf", filter))
  pdf(file=out.f, width=5, height=8)
  print.ggplot(p)
  dev.off()
}

get.merged.genes <- function(pset_a, pset_b, keep.cols=c('slr_dnds'), filter='default') {
  genes <- get.genes.sites(filter=filter)
  keep.cols <- unique(c('data_id', keep.cols))
  genes.a <- subset(genes, pset == pset_a, select=keep.cols)
  genes.b <- subset(genes, pset == pset_b, select=keep.cols)
  mrgd <- merge(genes.a, genes.b, by='data_id')
  mrgd
}

plot.gene.correlations <- function() {
  gene.correlations(1, 2)
  gene.correlations(1, 3)
  gene.correlations(1, 4)
  gene.correlations(4, 3)
  gene.correlations(3, 2)
  gene.correlations(2, 7)
  gene.correlations(6, 10)
}

gene.correlations <- function(pset_a, pset_b, filter='stringent', plot=T) {
  print(paste(pset_a, pset_b))
  genes <- get.merged.genes(pset_a, pset_b, keep.cols=c('slr_dnds', 'slr_dnds_z'), filter=filter)

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
    xs <- data[indices, 'slr_dnds.x']
    ys <- data[indices, 'slr_dnds.y']
    ws <- data[indices, 'weight']
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


  # Do a PCA without re-centering the data... this forces the 'root' of the best fit line to
  # be at the origin.
  zero.pca <- iwpca.slope.int(genes$slr_dnds.x, genes$slr_dnds.y, center=F)

  genes$weight <- 1 / ((genes$slr_dnds.x + genes$slr_dnds.y)/2 + 0.05)
  wts <- genes$weight
  wt.pca <- iwpca.slope.int(genes$slr_dnds.x, genes$slr_dnds.y, ws=wts, center=F)

  # Do the PCA with re-centering, allowing the line to sit anywhere.
  #nonzero.pca <- iwpca.slope.int(genes$slr_dnds.x, genes$slr_dnds.y, ws=wts, center=T)

  if (plot) {
    dd <- 0.01
    genes$slr_dnds.x <- round_any(genes$slr_dnds.x, dd)
    genes$slr_dnds.y <- round_any(genes$slr_dnds.y, dd)

    p <- ggplot(genes, aes(x=slr_dnds.x, y=slr_dnds.y))
    p <- p + theme_bw()
    p <- p + geom_tile(stat='sum', aes(fill=..prop..))
    p <- p + scale_fill_gradientn("Density", colour=c('white', rgb(0.5, 0.5, 1), 'red'), trans='log')
    p <- p + geom_abline(slope=1, intercept=0, colour='black', linetype='dashed')
    p <- p + geom_abline(slope=zero.pca[1], intercept=0, colour='black')
    lmt <- c(0, 1)
    p <- p + scale_x_continuous(paste(pset.to.alias(pset_a), "dN/dS"), limits=lmt, expand=c(0,0))
    p <- p + scale_y_continuous(paste(pset.to.alias(pset_b), "dN/dS"), limits=lmt, expand=c(0,0))
    p <- generic.opts(p)
    p <- p + opts(
      legend.position='none'
    )

    svg(file=scratch.f(sprintf("genes_cor_%d_%d_%s.svg", pset_a, pset_b, filter)), width=4, height=4)
    print.ggplot(p)
    dev.off()

  }

  zero.pca[1]
}

get.go.annotations <- function(exclude.iea=F) {
  ens.f <- scratch.f("ens_go_63.Rdata")
  if (!file.exists(ens.f)) {
    if (is.ebi()) {
      load.go.csv("~/lib/greg-ensembl/projects/2xmammals/ens_go_63.csv", ens.f)    
    } else {
      load.go.csv("~/src/greg-ensembl/projects/2xmammals/ens_go_63.csv", ens.f)
    }
  }
  load(ens.f)
  if (exclude.iea) {
    go.ens.excl.iea
  } else {
    go.ens
  }
}

bsub.check.data <- function() {
  # TODO: Do some quality checks on the data...
  
}



bsub.genes.enrichment <- function(...) {
  filters <- c('default', 'stringent')

  for (pset in c(6)) {
    for (filter in c('stringent')) {
      for (excl.iea in c(F)) {
        xtra <- sprintf("%d %s %s %s", pset, filter, 'emp_01_bh', as.character(excl.iea))
        bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)

        xtra <- sprintf("%d %s %s %s", pset, filter, 'emp_05_bh', as.character(excl.iea))
        bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)

        xtra <- sprintf("%d %s %s %s", pset, filter, 'emp_10_bh', as.character(excl.iea))
        bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)

        xtra <- sprintf("%d %s %s %s", pset, filter, 'hoch_bh', as.character(excl.iea))
        bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)

        xtra <- sprintf("%d %s %s %s", pset, filter, 'dnds_05', as.character(excl.iea))
        bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)

        xtra <- sprintf("%d %s %s %s", pset, filter, 'dnds_02', as.character(excl.iea))
        bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)        

        xtra <- sprintf("%d %s %s %s", 1, filter, 'core_01', as.character(excl.iea))
        bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)
  
        xtra <- sprintf("%d %s %s %s", 1, filter, 'core_05', as.character(excl.iea))
        bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)
      }
    }
  }

  return()

  for (pset in c(1, 3, 5, 6, 9)) {
    for (filter in filters) {
      for (excl.iea in c(T, F)) {
        xtra <- sprintf("%d %s %s %s", pset, filter, 'emp_01_bh', as.character(excl.iea))
        bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)
      }
    }
  }
  for (pset in c(2, 4, 7, 8, 10)) {
    for (filter in filters) {
      for (excl.iea in c(T, F)) {
        xtra <- sprintf("%d %s %s %s", pset, filter, 'emp_05_bh', as.character(excl.iea))
        bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)
      }
    }
  }

  return()

  # For 5 and 6 test the Hoch BH and top dN/dS
  for (pset in 5:6) {
    for (filter in filters) {
      for (excl.iea in c(T, F)) {
        xtra <- sprintf("%d %s %s %s", pset, filter, 'hoch_bh', as.character(excl.iea))
        bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)
        xtra <- sprintf("%d %s %s %s", pset, filter, 'dnds_05', as.character(excl.iea))
        bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)
        xtra <- sprintf("%d %s %s %s", pset, filter, 'dnds_02', as.character(excl.iea))
        bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)
      }
    }
  }
  
  # Do low-dnds PSGs.
  for (pset in c(1, 3, 6)) {
    for (filter in filters) {
      for (excl.iea in c(T, F)) {
        xtra <- sprintf("%d %s %s %s", pset, filter, 'psg_lowdnds_emp10', as.character(excl.iea))
        bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)
        xtra <- sprintf("%d %s %s %s", pset, filter, 'psg_lowdnds_emp05', as.character(excl.iea))
        bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)
        xtra <- sprintf("%d %s %s %s", pset, filter, 'psg_lowdnds_emp01', as.character(excl.iea))
        bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)
      }
    }
  }

  # Shared / unique PSGs.
  for (filter in filters) {
    for (excl.iea in c(T, F)) {
      xtra <- sprintf("%d %s %s %s", 6, filter, 'all_emp', as.character(excl.iea))
      bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)
      xtra <- sprintf("%d %s %s %s", 6, filter, 'all_hoch', as.character(excl.iea))
      bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)

      xtra <- sprintf("%d %s %s %s", 6, filter, 'm_hoch_noemp', as.character(excl.iea))
      bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)
      xtra <- sprintf("%d %s %s %s", 6, filter, 'm_emp_5_no1', as.character(excl.iea))
      bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)
      xtra <- sprintf("%d %s %s %s", 6, filter, 'm_emp_10_no5', as.character(excl.iea))
      bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)

      xtra <- sprintf("%d %s %s %s", 1, filter, 'unique_m', as.character(excl.iea))
      bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)
      xtra <- sprintf("%d %s %s %s", 1, filter, 'unique_p', as.character(excl.iea))
      bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)
      xtra <- sprintf("%d %s %s %s", 1, filter, 'unique_l', as.character(excl.iea))
      bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)

      xtra <- sprintf("%d %s %s %s", 1, filter, 'core_005', as.character(excl.iea))
      bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)
      xtra <- sprintf("%d %s %s %s", 1, filter, 'core_01', as.character(excl.iea))
      bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)
      xtra <- sprintf("%d %s %s %s", 1, filter, 'core_05', as.character(excl.iea))
      bsub.function('genes_enrichment', mem=4, extra.args=xtra, ...)

    }    
  }

}

bsub.t.test.enrichment <- function(...) {
  method <- 't_test'
  pset <- 1
  for (filter in c('default', 'stringent', 'pfam')) {
    for (excl.iea in c(T, F)) {
      xtra <- sprintf("%d %s %s %s", pset, filter, method, as.character(excl.iea))
      bsub.function('genes_enrichment', queue='normal', mem=6, extra.args=xtra, ...)
    }
  }
  
}

calc.genes.enrichment <- function(pset, filter='default', method='dnds_05', excl.iea=F, other_pset=NA) {
  if (!is.ebi()) {
    source("~/src/greg-ensembl/scripts/go_enrichments.R")
  } else {
    source("~/lib/greg-ensembl/scripts/go_enrichments.R")
  }

  go.ens <- get.go.annotations(exclude.iea=excl.iea)
  scores.df <- get.enrichment.df(pset_id=pset, filter=filter, method=method, other_pset=other_pset)
  scores.df <- subset(scores.df, select=c('score', 'binary.score', 'id', 'name', 'length', 't_test_a', 't_test_b'))

  print(head(scores.df))
  lbl <- sprintf("%s %s %s %s %s", pset, filter, method, excl.iea, other_pset)
  go.res <- do.enrichments(scores.df, lbl, go.ens)

  go.res$filter <- filter
  go.res$method <- method
  go.res$excl.iea <- as.integer(excl.iea)
  go.res$pset <- pset
  go.res$other_pset <- other_pset

  go.res$label <- paste(go.res$label, go.res$go_id, sep=' ')

  con <- connect(dbname())
  i <- 1
  while (i < nrow(go.res)) {
    lo <- i
    hi <- i + 200
    hi <- pmin(hi, nrow(go.res))
    sub.res <- go.res[lo:hi,]
    write.or.update(sub.res, 'go', con, 'label')
    i <- i + 200
  }
  disconnect(con)
}

get.enrichment.df <- function(pset_id=1, method='dnds_5', filter='stringent', other_pset=NA) {
  if (method %in% c('t_test', 'pos_parallel')) {
    genes <- get.merged.genes(pset_id, other_pset,
      keep.cols=c('aln_length', 'ref_protein_id', 'gene_name',
        'slr_dnds', 'slr_dnds_z', 'f_pos_10', 'f_neg_10', 'tpm_05', 'f_pos_01', 'f_pos_gene_bh_10'),
      filter=filter
    )
    genes$gene_name <- genes$gene_name.x
    genes$ref_protein_id <- genes$ref_protein_id.x
    genes$aln_length <- genes$aln_length.x
  } else {
    genes <- get.genes.sites(filter=filter)
    genes <- subset(genes, pset==pset_id)
  }

  # Methods:
  # dnds_05
  # hoch

  if (method == 'dnds_05') {
    genes$score <- genes$slr_dnds_z
    dnds.q <- quantile(genes$score, 0.95)
    genes$binary.score <- ifelse(genes$score >= dnds.q, 1, 0)

  } else if (method == 'dnds_02') {
    genes$score <- genes$slr_dnds_z
    dnds.q <- quantile(genes$score, 0.98)
    genes$binary.score <- ifelse(genes$score >= dnds.q, 1, 0)

  } else if (method %in% c('psg_lowdnds_emp10', 'psg_lowdnds_emp05', 'psg_lowdnds_emp01')) {
    genes$score <- switch(method,
      'psg_lowdnds_emp10' = genes$emp_10_bh,
      'psg_lowdnds_emp05' = genes$emp_05_bh,
      'psg_lowdnds_emp01' = genes$emp_01_bh
    )

    qv <- switch(method,
      'psg_lowdnds_emp10' = 0.75,
      'psg_lowdnds_emp05' = 0.78,
      'psg_lowdnds_emp01' = 0.8
    )
    genes$score <- ifelse(genes$slr_dnds_z > qnorm(qv), 1, genes$score)

    genes$binary.score <- ifelse(genes$score < 0.2, 1, 0)
    genes$score <- -log(pmax(0.00001, genes$score))

  } else if (method %in% c('all_emp', 'all_hoch')) {
    genes <- get.all.merged(filter=filter,
      keep.once=c('gene_name', 'aln_length', 'ref_protein_id'),
      keep.fields=c('hoch_bh', 'emp_01_bh', 'hoch_p', 'emp_01')
    )
    min.p <- rep(Inf, nrow(genes))
    for (i in 1:10) {
      if (method == 'all_emp') {
        emp.s <- sprintf("emp_01_bh.%d", i)
        min.p <- pmin(genes[, emp.s], min.p, na.rm=T)
      } else if (method == 'all_hoch') {
        hoch.s <- sprintf("hoch_p.%d", i)
        min.p <- pmin(genes[, hoch.s], min.p, na.rm=T)
      }
    }
    genes$score <- min.p
    genes$score.bh <- p.adjust(genes$score, method='BH')
    if (method == 'all_emp') {
      genes$binary.score <- ifelse(genes$score.bh < 0.1, 1, 0)
    } else if (method == 'all_hoch') {
      genes$binary.score <- ifelse(genes$score.bh < 0.1, 1, 0)
    }
    genes$score <- -log(pmax(0.0001, genes$score))

  } else if (method == 'hoch_bh') {
    genes$score <- genes$hoch_bh
    genes$binary.score <- ifelse(genes$score < 0.1, 1, 0)
    genes$score <- -log(pmax(0.0001, genes$score))

  } else if (method == 'emp_01_bh') {
    genes$score <- genes$emp_01_bh
    genes$binary.score <- ifelse(genes$score < 0.1, 1, 0)
    genes$score <- -log(pmax(0.0001, genes$score))

  } else if (method == 'emp_05_bh') {
    genes$score <- genes$emp_05_bh
    genes$binary.score <- ifelse(genes$score < 0.1, 1, 0)
    genes$score <- -log(pmax(0.0001, genes$score))

  } else if (method == 'm_hoch_no_emp') {
    genes$score <- genes$hoch_bh
    genes$score <- ifelse(genes$emp_05_bh < 0.1, 1, genes$score)
    genes$binary.score <- ifelse(genes$score < 0.1, 1, 0)
    genes$score <- -log(pmax(0.0001, genes$score))

  } else if (method == 'm_emp_5_no1') {
    genes$score <- genes$emp_05_bh
    genes$score <- ifelse(genes$emp_01_bh < 0.1, 1, genes$score)
    genes$binary.score <- ifelse(genes$score < 0.1, 1, 0)
    genes$score <- -log(pmax(0.0001, genes$score))

  } else if (method == 'm_emp_10_no5') {
    genes$score <- genes$emp_10_bh
    genes$score <- ifelse(genes$emp_05_bh < 0.1, 1, genes$score)
    genes$binary.score <- ifelse(genes$score < 0.1, 1, 0)
    genes$score <- -log(pmax(0.0001, genes$score))

  } else if (method %in% c('unique_p', 'unique_m', 'unique_l')) {
    other.pset <- switch(method,
      unique_p = c(3, 6),
      unique_l = c(1, 6),
      unique_m = c(1, 3)
    )
    this.pset <- switch(method,
      unique_p = c(1),
      unique_l = c(3),
      unique_m = c(6)
    )

    genes <- get.genes.sites(filter=filter)
    genes <- subset(genes, pset == this.pset)

    other.genes <- get.genes.sites(filter=filter)
    other.genes <- subset(other.genes, pset %in% other.pset)

    genes$score <- genes$emp_01_bh
    sig.others <- other.genes$emp_01_bh < 0.1
    sig.other.ids <- other.genes[sig.others, 'data_id']
    print(length(unique(sig.other.ids)))

    genes[genes$data_id %in% sig.other.ids, 'score'] <- 1
    genes$binary.score <- ifelse(genes$score < 0.1, 1, 0)
    genes$score <- -log(pmax(0.0001, genes$score))

  } else if (method %in% c('core_005', 'core_05', 'core_01')) {  

    genes <- get.all.merged(filter=filter,
      keep.once=c('gene_name', 'aln_length', 'ref_protein_id'),
      keep.fields=c(
        'emp_10_bh', 'emp_01_bh', 'emp_05_bh', 'emp_005_bh',
        'emp_01', 'emp_05', 'emp_10', 'emp_005')
    )
    max.p <- rep(-Inf, nrow(genes))
    psets <- c(1, 2, 3)

    for (i in psets) {
      if (method == 'core_005') {
        emp.s <- sprintf("emp_005.%d", i)
      } else if (method == 'core_05') {
        emp.s <- sprintf("emp_05.%d", i)
      } else if (method == 'core_01') {
        emp.s <- sprintf("emp_01.%d", i)
      }
      max.p <- pmax(genes[, emp.s], max.p, na.rm=T)
    }
    genes$score <- max.p
    genes$score <- ifelse(genes$score == -Inf, 1, genes$score)

    genes$score.bh <- p.adjust(genes$score, method='BH')
    genes$binary.score <- ifelse(genes$score.bh < 0.1, 1, 0)
    genes$score <- -log(pmax(0.0001, genes$score))
  
  } else if (method == 't_test') {
    # Use distance from the z-score to the unit line.
    zx <- genes$slr_dnds_z.x
    zy <- genes$slr_dnds_z.y
    z.dist <- (zx - zy) / sqrt(2)
    z.dnds <- (genes$slr_dnds.x + genes$slr_dnds.y) / 2

    z.cum <- pnorm(zy)
    # This alternative finds genes in the more negative z-score range
    # (which are mostly conserved and less interesting...)
    #z.score <- z.dist * (max(z.dnds) - z.dnds)
    z.score <- z.dist * z.cum
#    z.score <- z.dist
    #z.score <- ifelse(zx > qnorm(0.9), z.score, 0)

    genes$score <- z.score
    score.q <- quantile(genes$score, 0.99)

    print(cor(genes$slr_dnds_z.x, genes$slr_dnds_z.y))
    print(cor(genes$f_pos_10.x, genes$f_pos_10.y))

    pdf(file="test.pdf")
    plot(genes$slr_dnds_z.x, genes$slr_dnds_z.y)
    abv <- genes$score > score.q
    points(genes[abv, 'slr_dnds_z.x'], genes[abv, 'slr_dnds_z.y'], col='red') 
    dev.off()

    genes$t_test_a <- zx
    genes$t_test_b <- zy

    genes$binary.score <- ifelse(genes$score > score.q, 1, 0)
  } else {
    stop("Enrichment method not recognized!")
  }

  if (method != 't_test') {
    genes$t_test_a <- NA
    genes$t_test_b <- NA
  }

  genes <- subset(genes, gene_name != 'NA')
  genes$id <- genes$ref_protein_id
  genes$name <- genes$gene_name
  genes$length <- genes$aln_length

  #print(sprintf("%d sig genes", sum(genes$binary.score)))

  xx <- subset(genes, select=c('score', 'binary.score', 'name', 'length'))

  xx <- xx[order(xx$score, decreasing=T),]
  #print(head(xx, n=30))

  genes
}

plot.gene.fractions <- function() {
  genes <- get.genes.sites(filter='default')
  genes <- subset(genes, pset == 6)

  a <- 'f_pos_20'
  b <- 'f_ntr_20'
  c <- 'f_neg_20'

  genes$is_sig <- ifelse(genes$pos_holm_01 > 0, 1, 0)

  p <- ternary.plot(genes, a, b, c, color.fld='slr_dnds', size.fld='is_sig')
  pdf(file=scratch.f("test.pdf"), width=4, height=4)
  print.ggplot(p)
  dev.off()
}

ternary.plot <- function(x, fld.a, fld.b, fld.c, color.fld=NA, size.fld=NA, plot.hist=F) {
  library(ggplot2)

  a <- x[, fld.a]
  b <- x[, fld.b]
  c <- x[, fld.c]

  n <- length(a)

  # a: up
  # b: down-right
  # c: down-left
  yy <- a * sin(pi/3)
  xx <- b + a * cos(pi/3)


  x$xx <- xx
  x$yy <- yy  
  p <- ggplot(x, aes(x=xx, y=yy))

  sz <- 0.3
  if (!is.na(color.fld)) {
    if (!is.na(size.fld)) {
      p <- p + geom_point(aes_string(colour=color.fld, size=size.fld))
    } else {
      p <- p + geom_point(size=sz, aes_string(colour=color.fld))
    }
  } else {
    if (!is.na(size.fld)) {
      p <- p + geom_point(aes_string(size=size.fld))
    } else {
      p <- p + geom_point(size=sz)
    }
  }

  if (!is.na(size.fld)) {
    p <- p + scale_size(to=c(sz / 4, sz * 2))
  }
  if (!is.na(color.fld)) {
    p <- p + scale_colour_gradient()
  }


  plot.labels <- T
  if (plot.labels) {
    p <- p + geom_text(x=cos(pi/3), y=1 * sin(pi/3), label=fld.a, hjust=0.5, vjust=0)
    p <- p + geom_text(x=1, y=0, label=fld.b, hjust=1, vjust=1)
    p <- p + geom_text(x=0, y=0, label=fld.c, hjust=0, vjust=1)
  }

  # Create binned histograms along the edges.
  brks <- seq(from=-1, to=1, length.out=20)
  h <- hist(a-b, breaks=brks, plot=F)
  ab.df <- data.frame(mid=h$mids, count=h$counts/n)
  h <- hist(b-c, breaks=brks, plot=F)
  bc.df <- data.frame(mid=h$mids, count=h$counts/n)
  h <- hist(a-c, breaks=brks, plot=F)
  ac.df <- data.frame(mid=h$mids, count=h$counts/n)

  ab.df$count <- pmax(0.001, ab.df$count)
  bc.df$count <- pmax(0.001, bc.df$count)
  ac.df$count <- pmax(0.001, ac.df$count)
  if (!plot.hist) {
    ab.df$count <- 0.001
    bc.df$count <- 0.001
    ac.df$count <- 0.001
  }

  cf <<- cos(pi/3)
  sf <<- sin(pi/3)
  cff <<- cos(pi/3)
  sff <<- sin(pi/3)
  p <- p + geom_segment(data=ac.df,
    aes(x=.5*cf*(mid+1), y=.5*sf*(mid+1),
        xend=(.5*cf*(mid+1)) - count*sff,
        yend=(.5*sf*(mid+1)) + count*cff),
    colour='black',
    size=2
  )
  p <- p + geom_segment(data=ab.df,
    aes(x=1-(.5*cf*(mid+1)), y=.5*sf*(mid+1),
        xend=1-(.5*cf*(mid+1)) + count*sff,
        yend=(.5*sf*(mid+1)) + count*cff),
    colour='black',
    size=2
  )
  p <- p + geom_segment(data=bc.df,
    aes(y=0, x=.5*(mid+1),
        xend=.5*(mid+1),
        yend=-count*cff),
    colour='black',
    size=2
  )

  if (plot.hist) {
    p <- p + scale_x_continuous(limits=c(-0.5, 1.5))
    p <- p + scale_y_continuous(limits=c(-sin(pi/3)*.5, sin(pi/3)*1.5))
  } else {
    p <- p + scale_x_continuous()
    p <- p + scale_y_continuous()
  }
  p <- p + coord_equal()
  p <- p + theme_bw()
  p <- p + opts(
    legend.position='none'
  )

  return(p)
}


write.psg.table <- function(filter='stringent') {
  genes <- get.genes.sites(filter=filter)
  
  genes$pset.fact <- pset.to.alias(genes$pset, factors=T)
  genes$filter.fact <- filter.to.alias(genes$filter, factors=T)

  print("Summarizing...")
  summary.rows <- ddply(genes, .(pset), function(x) {
    # Number of genes
    n <- nrow(x)

    # Mean dN/dS
    mean.dnds <- mean(x$slr_dnds)
    # dN/dS corresponding to zero z-score (log-mean)
    sorted.zs <- x[order(x$slr_dnds_z),]
    max.neg.zscore.index <- max(which(sorted.zs$slr_dnds_z < 0))
    logmean.dnds <- sorted.zs[max.neg.zscore.index, 'slr_dnds']

    # FWER PSGs
    n.hoch.20 <- nrow(subset(x, hoch_bh < 0.2))
    n.hoch.10 <- nrow(subset(x, hoch_bh < 0.1))
    n.hoch.05 <- nrow(subset(x, hoch_bh < 0.05))

    # TPM PSGs
    n.tpm.50 <- nrow(subset(x, tpm_50_bh < 0.1))
    n.tpm.20 <- nrow(subset(x, tpm_20_bh < 0.1))
    n.tpm.10 <- nrow(subset(x, tpm_10_bh < 0.1))
    n.tpm.05 <- nrow(subset(x, tpm_05_bh < 0.1))
    n.fis <- nrow(subset(x, fis_bh < 0.1))

    # Empirical PSGs
    n.e.005 <- nrow(subset(x, emp_005_bh < 0.1))
    n.e.01 <- nrow(subset(x, emp_01_bh < 0.1))
    n.e.05 <- nrow(subset(x, emp_05_bh < 0.1))
    n.e.10 <- nrow(subset(x, emp_10_bh < 0.1))

    out.df <- data.frame(
      filter = x[1, 'filter.fact'],
      pset = x[1, 'pset.fact'],
      
      n.genes = n,

      mean.dnds = mean.dnds,
      logmean.dnds = logmean.dnds,

#      n.hoch.20 = n.hoch.20,
      n.hoch.10 = n.hoch.10,
#      n.hoch.05 = n.hoch.05,

      n.fis = n.fis,
      n.tpm.50 = n.tpm.50,
      n.tpm.20 = n.tpm.20,
      n.tpm.10 = n.tpm.10,
      n.tpm.05 = n.tpm.05,
      
      n.e.005 = n.e.005,
      n.e.01 = n.e.01,
      n.e.05 = n.e.05,
      n.e.10 = n.e.10
    )
    out.df
  })
  
  summary.rows <- summary.rows[order(summary.rows$pset),]

  summary.rows[duplicated(summary.rows$filter), 'filter'] <- NA

  xt <- xtable(summary.rows)
  xt <- color.columns(xt, columns=c('n.genes', 'mean.dnds', 'logmean.dnds'))

  xt <- color.columns(xt, columns=c('n.hoch.10', 'n.tpm.50', 'n.tpm.20', 'n.tpm.10', 'n.tpm.05', 'n.fis', 'n.e.005', 'n.e.01', 'n.e.05', 'n.e.10'))

  out.f <- scratch.f(sprintf("psg_summary_%s.txt", filter))
  print.latex(xt, file=out.f)
  
  #print(summary.rows)
  
  mamm.row <- subset(summary.rows, pset == 'Mammals')
  print(mamm.row[1, c('n.hoch.10', 'n.tpm.50', 'n.e.01')] / mamm.row$n.genes)
}

pval.histogram <- function() {
  flds <- c('hoch_p', 'fis_p', 'tpm_05', 'tpm_10', 'tpm_20', 'emp_005', 'emp_01', 'emp_05', 'emp_10')
  flds.fact <- pval.to.alias(flds)

  genes <- get.genes.sites(filter='stringent')
  genes <- subset(genes, pset==6)

  out.df <- data.frame()
  for (i in 1:length(flds)) {
    out.df <- rbind(out.df, data.frame(
      fld=flds.fact[i],
      pval=genes[, flds[i]]
    ))
  }

  out.df <- ddply(out.df, .(fld), function(x) {
    ecdf.f <- ecdf(x$pval)
    xx <- seq(from=0, to=1, length.out=400)
    data.frame(
      xx = xx,
      yy = ecdf.f(xx)
    )
  })

  print(head(out.df))

  out.df$fld <- factor(out.df$fld)

  p <- ggplot(out.df, aes(x=xx, y=yy))
  p <- generic.opts(p)
  p <- p + geom_line()
  p <- p + facet_wrap('fld')
  p <- p + coord_equal()
  p <- p + opts(
    strip.background = theme_blank()    
  )
  
  out.f <- scratch.f(sprintf("psg_pval_histogram.pdf"))
  pdf(out.f)
  print.ggplot(p)
  dev.off()

  genes <- get.genes.sites(filter='stringent')
  
  df1 <- genes
  df1$fld <- 'hoch_p'
  df1$pval <- df1$hoch_p

  df2 <- genes
  df2$fld <- 'emp_01'
  df2$pval <- df2$emp_01

  df3 <- genes
  df3$fld <- 'tpm_05'
  df3$pval <- df3$tpm_05

  df4 <- genes
  df4$fld <- 'fis_p'
  df4$pval <- df4$fis_p
  
  pval.df <- rbind(df1, df2, df3, df4)

  out.df <- ddply(pval.df, .(pset, fld), function(x) {
    ecdf.f <- ecdf(x$pval)
    xx <- seq(from=0, to=1, length.out=400)
    data.frame(
      xx = xx,
      yy = ecdf.f(xx)
    )
  })

  out.df$pset <- pset.to.alias(out.df$pset, factors=T)
  out.df$fld <- pval.to.alias(out.df$fld)

  out.df$yy <- pmin(out.df$yy, 0.3)

  out.df <- base::within(out.df, {
    hoch.ind <- pset == 'Hochberg'
    pset <- reorder(pset, yy, function(x) {
      -quantile(x, 0.5)
    })
  })
  

  p <- ggplot(out.df, aes(x=xx, y=yy, colour=pset, linetype=pset))
  p <- generic.opts(p)
  p <- p + geom_line()
  p <- p + scale_colour_discrete("Species Group")
  p <- p + scale_linetype_discrete("Species Group")
  p <- p + scale_x_continuous("P-value")
  p <- p + scale_y_continuous("Cumulative Proportion")
  p <- p + facet_wrap('fld')
  p <- p + coord_equal()
  p <- p + opts(
    strip.background = theme_blank(),
    axis.text.x = theme_text(angle=90, hjust=1)
  )

  out.f <- scratch.f(sprintf("psg_pval_histogram_psets.pdf"))
  pdf(out.f, width=6, height=6)
  print.ggplot(p)
  dev.off()


}

length.bias <- function() {
  genes <- get.genes.sites(filter='stringent')
  genes <- subset(genes, pset==6)

  flds <- pval.fields()

  df.out <- data.frame()
  for (fld in flds) {
    genes$tmp <- genes[, fld]
    sig.genes <- subset(genes, tmp < 0.05)
    nsig.genes <- subset(genes, tmp >= 0.05)

    n.sig <- nrow(sig.genes)

    mean.sig <- mean(sig.genes$n_sites)
    mean.nsig <- mean(nsig.genes$n_sites)

    t.res <- t.test(sig.genes$n_sites, nsig.genes$n_sites)
    t.p <- t.res$p.value

    df.out <- rbind(df.out, data.frame(
      filter = filter.to.alias('stringent'),
      pset = pset.to.alias(6, factors=T),
      pval.thresh = .05,
      fld = fld,
      n.sig = n.sig,
      mean.sig = mean.sig,
      mean.nsig = mean.nsig,
      t.p = t.p
    ))
  }

  print(df.out)

}

pval.fields <- function() {
  c('hoch_p', 'emp_01', 'tpm_05')
}

pval.to.alias <- function(pvals, factors=T) {
  fct <- factor(pvals,
    levels = c('hoch_p', 'fis_p', 'tpm_05', 'tpm_10', 'tpm_20', 'tpm_50',
      'emp_005', 'emp_01', 'emp_05', 'emp_10', 
      'hoch_p.1', 'hoch_p.6', 'emp_10.1', 'emp_10.6', 'emp_01.1', 'emp_01.3', 'emp_01.6',
      'slr_dnds.1', 'slr_dnds.2', 'slr_dnds.3', 'slr_dnds.6'),
    labels = c(
      'Hochberg',
      'Fisher',
      'TPM 0.05',
      'TPM 0.1',
      'TPM 0.2',
      'TPM 0.5',
      'Empirical 0.005',
      'Empirical 0.01',
      'Empirical 0.05',
      'Empirical 0.1',
      'Hochberg (Primates)',
      'Hochberg (Mammals)',
      'Empirical 0.1 (Primates)',
      'Empirical 0.1 (Mammals',
      'Empirical 0.01 (Primates)',
      'Empirical 0.01 (Laurasiatheria)',
      'Empirical 0.01 (Mammals)',
      'M0 dN/dS (Primates)',
      'M0 dN/dS (Glires)',
      'M0 dN/dS (Laurasiatheria)',
      'M0 dN/dS (Mammals)'
    ))
  if (factors) {
    fct
  } else {
    as.character(fct)
  }
}

get.common.psgs <- function() {
  genes <- get.genes.sites(filter='default')
  genes <- subset(genes, pset == 6)
  pubs <- get.published.genes()
  genes <- merge(genes, pubs, all.x=T)
  emp.psgs <- subset(genes, emp_01_bh < 0.1)$gene_name
  k.psgs <- subset(genes, kosiol == 1)$gene_name
  n.psgs <- subset(genes, nielsen == 1)$gene_name

  intersect(emp.psgs, intersect(k.psgs, n.psgs))
}

pval.venn <- function() {
  library(Vennerable)
  genes <- get.genes.sites(filter='stringent')

  # Do the 3 mammalian superorders using emp_05 at FDR<0.1
  psg.p <- subset(genes, pset == 1 & emp_01_bh < 0.1)$gene_name
  psg.g <- subset(genes, pset == 2 & emp_01_bh < 0.1)$gene_name
  psg.l <- subset(genes, pset == 3 & emp_01_bh < 0.1)$gene_name
  psg.a <- subset(genes, pset == 4 & emp_01_bh < 0.1)$gene_name
  psg.e <- subset(genes, pset == 5 & emp_01_bh < 0.1)$gene_name
  psg.m <- subset(genes, pset == 6 & emp_01_bh < 0.1)$gene_name
  psg.sg <- subset(genes, pset == 7 & emp_01_bh < 0.1)$gene_name
  psg.sm <- subset(genes, pset == 8 & emp_01_bh < 0.1)$gene_name
  psg.hq <- subset(genes, pset == 9 & emp_01_bh < 0.1)$gene_name
  psg.hmrd <- subset(genes, pset == 10 & emp_01_bh < 0.1)$gene_name

  psg.m.hoch <- subset(genes, pset == 6 & hoch_bh < 0.1)$gene_name
  psg.m.tpm <- subset(genes, pset == 6 & tpm_05_bh < 0.1)$gene_name
  psg.m.emp10 <- subset(genes, pset == 6 & emp_10_bh < 0.1)$gene_name
  psg.m.emp05 <- subset(genes, pset == 6 & emp_05_bh < 0.1)$gene_name
  psg.m.emp01 <- subset(genes, pset == 6 & emp_01_bh < 0.1)$gene_name

  pan.superorder <- unique(c(psg.p, psg.g, psg.l))
  pan.low <- union.all(psg.g, psg.a, psg.sg, psg.sm, psg.hmrd)
  pan.hi <- union.all(psg.p, psg.l, psg.m, psg.e, psg.hq)

  # Print out some simple counts
  # Get the union of all genes w/ an empirical PSG in some pset.
  genes.all <- unique(genes$gene_name)
  psg.all <- unique(subset(genes, emp_01_bh < 0.1)$gene_name)
  print(sprintf("Union of all Empirical 0.01 PSGs: %d PSG, %d all (%.3f%%)", 
    length(psg.all), length(genes.all), length(psg.all)/length(genes.all)*100
  ))

  psg.all <- unique(subset(genes, emp_01_bh < 0.1 | tpm_05_bh < 0.1 | hoch_bh < 0.1)$gene_name)
  print(sprintf("Union of ALL PSGs: %d PSG, %d all (%.3f%%)", 
    length(psg.all), length(genes.all), length(psg.all)/length(genes.all)*100
  ))

  # Get PSGs from kosiol & nielsen
  genes <- get.genes.sites(filter='stringent')
  genes <- subset(genes, pset == 6)
  emp.psgs <- subset(genes, hoch_bh < 0.1)$data_id

  pubs <- get.published.genes()
  genes <- pubs
  c.psgs <- subset(genes, clark == 1)$data_id
  n.psgs <- subset(genes, nielsen == 1)$data_id
  r.psgs <- subset(genes, rhesus == 1)$data_id 
  k.psgs <- subset(genes, kosiol == 1)$data_id

  w <- Venn(Sets=list(
    Clark = c.psgs,
    Nielsen = n.psgs,
    Rhesus = r.psgs
  ))
  out.f <- scratch.f(sprintf("psg_venn_pubs1.pdf"))
  pdf(file=out.f)
  plotVenn(w)
  dev.off()

  w <- Venn(Sets=list(
    RMGSAC = r.psgs,
    Kosiol = k.psgs,
    Nielsen = n.psgs
  ))
  out.f <- scratch.f(sprintf("psg_venn_pubs2.pdf"))
  pdf(file=out.f)
  plotVenn(w)
  dev.off()

  w <- Venn(Sets=list(
    me = emp.psgs,
    Kosiol = k.psgs,
    Nielsen = n.psgs
  ))
  out.f <- scratch.f(sprintf("psg_venn_pubs3.pdf"))
  pdf(file=out.f)
  plotVenn(w)
  dev.off()

  return()

  w <- Venn(Sets=list(
    emp_10=psg.m.emp10,
    emp_05=psg.m.emp05,
    emp_01=psg.m.emp01
  ))
  out.f <- scratch.f(sprintf("psg_venn_emp_10_5_1.pdf"))
  pdf(file=out.f)
  plotVenn(w)
  dev.off()

  return()
  
  w <- Venn(Sets=list(
    Low=pan.low,
    High=pan.hi
  ))
  out.f <- scratch.f(sprintf("psg_venn_hi_low.pdf"))
  pdf(file=out.f)
  plotVenn(w)
  dev.off()

  w <- Venn(Sets=list(
    Hochberg=psg.m.hoch,
    Empirical=psg.m,
    TPM=psg.m.tpm
  ))
  out.f <- scratch.f(sprintf("psg_venn_methods.pdf"))
  pdf(file=out.f)
  plotVenn(w)
  dev.off()


  w <- Venn(Sets=list(
    HMRD=psg.hmrd,
    'HQ Mammals'=psg.hq,
    'Mammals'=psg.m
  ))
  out.f <- scratch.f(sprintf("psg_venn_mammals.pdf"))
  pdf(file=out.f)
  plotVenn(w)
  dev.off()

  w <- Venn(Sets=list(
    Eutheria=psg.e,
    'Primates+Glires+Laurasiatheria'=pan.superorder
  ))
  out.f <- scratch.f(sprintf("psg_venn_eutheria_allsuperorders.pdf"))
  pdf(file=out.f)
  plotVenn(w)
  dev.off()

  w <- Venn(Sets=list(
    Primates=psg.p,
    Glires=psg.g,
    Laurasiatheria=psg.l
  ))
  out.f <- scratch.f(sprintf("psg_venn_superorders.pdf"))
  pdf(file=out.f)
  plotVenn(w)
  dev.off()

  w <- Venn(Sets=list(
    Mammals=psg.m,
    Primates=psg.p,
    Laurasiatheria=psg.l
  ))
  out.f <- scratch.f(sprintf("psg_venn_superorders2.pdf"))
  pdf(file=out.f)
  plotVenn(w)
  dev.off()


}

pub.goids <- function() {
  x <- get.all.go.ids()

  c.g <- read.table('clark-go.txt')
  r.g <- read.table('rhesus-go.txt')
  k.g <- read.table('kosiol-go.txt')

  x$C <- ifelse(x$go_id %in% c.g$V1, 1, 0)
  x$R <- ifelse(x$go_id %in% r.g$V1, 1, 0)
  x$K <- ifelse(x$go_id %in% k.g$V1, 1, 0)
  x
}

get.all.go.ids <- function() {
  out.f <- scratch.f("go_ids.Rdata")
  if (!file.exists(out.f)) {
    con <- connect(dbname())
    cmd <- sprintf("select distinct(go_id) from go;")
    go.df <- dbGetQuery(con, cmd)
    disconnect(con)  
    go.df
    save(go.df, file=out.f)
  }
  load(out.f)
  go.df
}

write.go.table <- function() {
  source("~/lib/greg-ensembl/scripts/go_enrichments.R")

  top.go <- function(filter, method, pset=NA, excl.iea=0) {
    con <- connect(dbname())
    pset.s <- ifelse(is.na(pset), '', paste('and pset=', pset, sep=''))
    cmd <- sprintf("select method, total_n, total_sig, n, n_sig, n_exp, go_id, descr, fis_pval, topgo_pval, goseq_wall, mwu_pval, mean_length, genes
      from go where filter='%s' and method='%s' %s and excl_iea=%d
#        and n < 500 and n > 2
        order by topgo_pval asc;",
        filter, method, pset.s, excl.iea
    )
    go.df <- dbGetQuery(con, cmd)
    disconnect(con)  
    
    go.df$fis_bh <- go.df$fis_pval
    go.df$mwu_bh <- go.df$msu_pval
    if (nrow(go.df > 0)) {
      go.df$fis_bh <- p.adjust(go.df$fis_pval, method='BH')
      go.df$mwu_bh <- p.adjust(go.df$mwu_pval, method='BH')
    }

    if (nrow(go.df) > 0) {
      go.df <- go.df[order(go.df$n, decreasing=T),]
    }
    go.df
  }

  p.go <- function(x) {
    print(x[1, c('total_n', 'total_sig')])
    go <- subset(x, fis_bh < 0.1 & topgo_pval < 0.25, 
      select=c('go_id', 'n', 'n_sig', 'n_exp', 'descr', 'fis_bh', 
      'topgo_pval', 'goseq_wall', 'genes')
    )
    fis.ids <- go$go_id
    #go$genes <- substr(go$genes, 1, 60)
    max.n <- 100
    if (nrow(go) > max.n) {
      print(head(go, n=max.n))
      print("... Truncated at 30 rows.")
    } else {
      print(go)
    }
    mwu.sig <- subset(x, mwu_bh < 0.1 & !(go_id %in% fis.ids), 
      select=c('go_id', 'descr', 'mwu_pval')
    )
    if (nrow(mwu.sig) > 0) {
      print("Mwu...")
      mwu.sig <- mwu.sig[order(mwu.sig$mwu_pval),]
      print(head(mwu.sig, n=10))
      if (nrow(mwu.sig) > 10) {
        print("... Truncated at 10 rows.")
      }
    }
  }

  # OK, here we go.
  x <- get.all.go.ids()
  f <- 'stringent'
  e <- 0
  enrich.list <- list(
    P = top.go(f, 'emp_01_bh', 1, e),
    g = top.go(f, 'emp_05_bh', 2, e),
    L = top.go(f, 'emp_01_bh', 3, e),
    M = top.go(f, 'emp_01_bh', 6, e),
    m = top.go(f, 'emp_05_bh', 6, e),
    H = top.go(f, 'hoch_bh', 6, e),
    D = top.go(f, 'dnds_02', 6, e),
    i = top.go(f, 'core_05', 1, e)
  )

  nms <- names(enrich.list)
  for (i in 1:length(enrich.list)) {
    chr <- nms[i]
    chr.p <- paste(chr, '.p', sep='')

    cur.df <- enrich.list[[i]]
    cur.df[, chr.p] <- cur.df$fis_pval
    x <- merge(x, cur.df[, c('go_id', chr.p)], by=c('go_id'), all.x=T)

    cur.sub <- subset(cur.df, n_sig >4)
    cur.sub$fis_bh <- p.adjust(cur.sub$fis_pval, method='BH')

    sig.df <- subset(cur.sub, fis_bh < 0.1)
    x[, chr] <- as.integer(x$go_id %in% sig.df$go_id)
  }

  pub.x <- pub.goids()
  x <- merge(x, pub.x, all.x=T)

  # Take the n, n_sig, and genes from the emp_01_bh mammals.
  ms.df <- enrich.list[['M']]
  x <- merge(x, ms.df, by=c('go_id'), all.x=T)

  me.flds <- c('P', 'g', 'L', 'M', 'm', 'H', 'D', 'i')
  x$me.str <- ''
  for (fld in me.flds) {
    x$me.str <- paste(x$me.str, ifelse(x[, fld], fld, '~'), sep='')
  }

  them.flds <- c('C', 'R', 'K')
  x$them.str <- ''
  for (fld in them.flds) {
    x$them.str <- paste(x$them.str, ifelse(x[, fld], fld, '~'), sep='')
  }

  x$me.sig <- x$M + x$m + x$H + x$D + x$i
  x$them.sig <- x$C + x$R + x$K

  x$descr <- getTermsDefinition(x$go_id, 'BP', numChar=70)
  #x$genes <- substr(x$genes, 1, 70)
  regex.res <- regexpr('(\\w*\\W*){1,5}', x$genes, perl=T)
  match.pos <- regex.res
  match.len <- attr(regex.res, 'match.length')
  longer.than <- match.len >= nchar(x$genes) - 1
  x$genes <- substr(x$genes, match.pos, match.len)
  x$genes <- gsub(' (\\w)', ', \\1', x$genes)
  x$genes <- R.oo::trim(x$genes)

  gns <- c('SETX', 'BRCA2', 'HAUS6', 'CEP152', 'CEP250', 'CENPT', 'CENPI', 'CENPQ', 'CENPH')
  for (g in gns) {
    x$genes <- gsub(g, paste('\\\\bf{', g, '}', sep=''), x$genes)
  }

  # Remove some cruft.
  x <- subset(x, !grepl("regulation", descr))
  x <- subset(x, !grepl("activation", descr))
  x <- subset(x, !grepl("acute", descr))
  x <- subset(x, !grepl("other organism", descr))
  x <- subset(x, !grepl("external stimulus", descr))
  x <- subset(x, !grepl("system", descr))

  x$descr <- gsub('tumor necrosis factor', 'TNF', x$descr)
  x$descr <- gsub('pathway', '', x$descr)

  tbl.cols <- c('go_id', 'descr', 'me.str', 'them.str', 'M.p', 'topgo_pval',
    'n', 'n_sig', 'n_exp', 'genes')

  x.top <- subset(x, n < 300 & n > 30 & topgo_pval < 0.1)
  x.top <- x.top[order(x.top$fis_pval),]
  x.top <- head(x.top, n=10)
 
  x.me <- subset(x, them.sig == 0 & (me.sig >= 2 | H == 1) & n < 300 & n > 30 & topgo_pval < 0.1)
  x.me <- subset(x.me, !(go_id %in% x.top$go_id))
  x.me <- x.me[order(x.me$fis_pval),]
  x.me <- head(x.me, n=10)
 
  x.them <- subset(x, them.sig >= 2)
  x.them <- subset(x.them, !(go_id %in% x.top$go_id))
  x.them <- x.them[order(x.them$them.sig, x.them$fis_pval),]

  print.ltx <- function(y, f) {
    y <- subset(y, select=tbl.cols)
    int.cols <- c('n', 'n_sig')
    for (col in int.cols) {
      y[, col] <- as.integer(y[, col])
    }
    print(y)
    xt <- xtable(y)
    xt <- xt.surround(xt, column='me.str', prefix='\\texttt{', suffix='}')  
    xt <- xt.surround(xt, column='them.str', prefix='\\texttt{', suffix='}')  
    #xt <- color.columns(xt, columns=c('me.sig', 'them.sig'))
    below.tg <- which(xt$topgo_pval < 0.05)
    below.fis <- which(xt$M.p < 0.05)
    xt$M.p <- sprintf("%.1e", xt$M.p)
    xt$topgo_pval <- sprintf("%.1e", xt$topgo_pval)
    xt$n_exp <- sprintf("%.1f", xt$n_exp)
    xt <- xt.surround(xt, column='M.p', rows=below.fis, prefix='\\bf{', suffix='}')  
    xt <- xt.surround(xt, column='topgo_pval', rows=below.tg, prefix='\\bf{', suffix='}')  
    xt <- xt.surround(xt, column='genes', prefix='\\tiny{', suffix='}')  
    #print(xt, sanitize.text.function=function(x){x})
    print.latex(xt, file=f)
  }

  out.f <- scratch.f(sprintf("psg_go_top.txt"))
  print.ltx(x.top, f=out.f)
  out.f <- scratch.f(sprintf("psg_go_them.txt"))
  print.ltx(x.them, f=out.f)
  out.f <- scratch.f(sprintf("psg_go_me.txt"))
  print.ltx(x.me, f=out.f)  

  return()

  # Emp_0.01 in Mammals recapitulates much of Kosiol's analysis.
  # - complement activation
  # - innate / adaptive immune response
  # - etc. etc.
  go <- top.go('default', 'emp_01_bh', pset=6, excl.iea=0)
  p.go(go)

  return()

  # Hochberg PSGs in mammals:
  #go <- top.go('stringent', 'hoch_bh', pset=6, excl.iea=0)
  #go.m.hoch <- subset(go, fis_bh < 0.1 & topgo_pval < 0.25)
  #p.go(go.m.hoch)

  # TODO: emp_0.05 but not emp_0.01
  # TODO: emp_0.10 but not emp_0.05

  # Low-dnds PSGs in mammals
  #go <- top.go('stringent', 'psg_lowdnds_emp05', pset=6, excl.iea=0)
  #go.lowdnds <- subset(go, fis_bh < 0.1 & topgo_pval < 0.25)
  #p.go(go.lowdnds)

  f <- 'default'
  
  if (F) {
    # Unique primates / laur's
    go <- top.go(f, 'unique_p', pset=1, excl.iea=0)
    p.go(go)
    go <- top.go(f, 'unique_l', pset=1, excl.iea=0)
    p.go(go)
  }

  if (F) {
    go <- top.go(f, 'psg_lowdnds_emp05', pset=6, excl.iea=0)
    p.go(go)
    go <- top.go(f, 'psg_lowdnds_emp01', pset=6, excl.iea=0)
    p.go(go)
  }

  if (F) {
    go <- top.go(f, 'm_emp_5_no1', pset=6, excl.iea=0)
    p.go(go)
    go <- top.go(f, 'm_emp_10_no5', pset=6, excl.iea=0)
    p.go(go)
  }

  if (F) {
    go <- top.go(f, 'dnds_02', pset=6, excl.iea=0)
    p.go(go)
#    go <- top.go(f, 'hoch_bh', pset=6, excl.iea=0)
#    p.go(go)
  }

  if (F) {
    go <- top.go(f, 'core_05', pset=1, excl.iea=0)
    p.go(go)
  }

  return()

  go <- top.go('stringent', 'dnds_01', pset=6, excl.iea=0)
  go.m.dnds <- subset(go, fis_bh < 0.1 & topgo_pval < 0.25)

  p.go(go.m.dnds)
  p.go(go.m.hoch)
  p.go(go.m.emp)

  go <- top.go('stringent', 'unique_p', pset=1, excl.iea=0)
  go.p <- subset(go, fis_bh < 0.2 & topgo_pval < 0.2)

  return()

  xx <- top.go('default', 'psg_lowdnds_emp01', pset=1, excl.iea=0)
  print(head(xx$descr))
  xx <- top.go('default', 'psg_lowdnds_emp01', pset=3, excl.iea=0)
  print(head(xx$descr))
  xx <- top.go('default', 'psg_lowdnds_emp01', pset=6, excl.iea=0)
  print(head(xx$descr))

  #xx <- top.go('default', 'all_bh', pset=6, excl.iea=1)
  #print(head(xx$descr))

  #xx <- top.go('default', 'core_all', pset=1, excl.iea=0)
  #print(head(xx$descr, n=10))

  #xx <- top.go('default', 'core_all', pset=1, excl.iea=1)
  #print(head(xx$descr, n=10))

  #xx <- top.go('default', 'core_hi', pset=1, excl.iea=0)
  #print(head(xx$descr, n=10))

  xx <- top.go('default', 'unique_l', pset=1, excl.iea=0)
  print(head(xx$descr, n=10))

}

get.nearby.pairs <- function(force.fetch=F) {
  print("Getting nearby pairs")
  out.f <- scratch.f("merged_nearby_pairs.Rdata")
  if (!file.exists(out.f) || force.fetch) {
    genes <- get.genes()
    genes <- add.family.ids(genes)

    close.df <- ddply(genes, .(chr_name), function(x) {
      print(x[1, 'chr_name'])
      summary(x$chr_start)
      x <- x[order(x$chr_start),]
      dd <- as.matrix(dist(x$chr_start))
      diag(dd) <- Inf
      close.ind <- which(dd < 2000000)

      rws <- as.vector(row(dd))
      cls <- as.vector(col(dd))
      close.rws <- rws[close.ind]
      close.cls <- cls[close.ind]
    
      close.a <- x[close.rws, ]
      close.b <- x[close.cls, ]

      min.id <- pmin(close.a$data_id, close.b$data_id)
      max.id <- pmax(close.a$data_id, close.b$data_id)

      ids <- paste(min.id, max.id)

      xx <- data.frame(
        ids = ids,
        id.a = close.a$data_id,
        id.b = close.b$data_id,
        gene.a = close.a$gene_name,
        gene.b = close.b$gene_name,
        fam.a = close.a$ensf,
        fam.b = close.b$ensf,
        stringsAsFactors=F
      )
      xx <- subset(xx, fam.a == fam.b)
      xx
    })

    genes$nearby_sib <- 0
    genes[genes$data_id %in% close.df$id.a, 'nearby_sib'] <- 1

    nearby.genes <- subset(genes, select=c('data_id', 'nearby_sib', 'ensf'))
    save(nearby.genes, file=out.f)
  }
  load(out.f)

  nearby.genes
}

top.domains <- function() {
  filter = 'default'

  pf.df <- 'pfam.sites.df'
  if (!exists(pf.df, envir=.GlobalEnv)) {
    sites <- get.pset.sites(6, 'pfam', test=F)
    assign(pf.df, sites, envir=.GlobalEnv)
  }
  sites <- get(pf.df, envir=.GlobalEnv)
  print(sprintf("Got %d sites", nrow(sites)))
  pos.sites <- subset(sites, !is.na(pfam_domain) & pos.pval < 0.01)
  print(sprintf("%d pos sites", nrow(pos.sites)))
  genes <- get.all.merged(filter='default',
    keep.once=c('gene_name'),
    keep.fields=c('emp_01_bh')
  )

  keep.flds <- c('hoch_bh', 'emp_10_bh', 'emp_01_bh', 'emp_05_bh', 'f_pos_05', 'n_sites', 'n_genes',
    'pos_gene_names', 'all_gene_names', 'n.pos.genes', 'pos_01', 'pos_05', 'f_pos_10', 'f_neg_10', 'f_pos_05',
    'f_neg_05')

  domains.d <- get.domains.merged(filter='default',
    keep.fields = keep.flds
  )
  domains.s <- get.domains.merged(filter='stringent',
    keep.fields=keep.flds
  )
  domains.p <- get.domains.merged(filter='pfam',
    keep.fields=keep.flds
  )
  
  domains <- merge(domains.d, domains.s, by=c('pfam_domain'), all.x=T, suffixes=c('.d', '.s'))
  domains2 <- merge(domains.d, domains.p, by=c('pfam_domain'), all.x=T, suffixes=c('.d', '.p'))

  domains <- merge(domains, domains2)

  domains$pfam_type <- domains$pfam_type.d
  domains$pfam_length <- domains$pfam_length.d
  domains$pfam_name <- domains$pfam_name.d
  domains$pfam_descr <- domains$pfam_descr.d

  sig.ids <- function(filter, method, pset) {
    x <- get.domains.merged(filter=filter, keep.fields=method)
    m.str <- paste(method, '.', pset, sep='')
    x$val <- x[, m.str]
    #print(m.str)
    sig.sub <- subset(x, !is.na(val) & val < 0.1)
    sig.sub$pfam_domain
  }

  get.list <- function(f) {
    list(
      P = sig.ids(f, 'emp_01_bh', 1),
      g = sig.ids(f, 'emp_05_bh', 2),
      L = sig.ids(f, 'emp_01_bh', 3),
      M = sig.ids(f, 'emp_01_bh', 6),
      m = sig.ids(f, 'emp_05_bh', 6),
      H = sig.ids(f, 'hoch_bh', 6)
    )
  }
  
  sig.s <- get.list('stringent')
  sig.d <- get.list('default')
  sig.p <- get.list('pfam')

  domains$s.sig <- 0
  domains$d.sig <- 0
  domains$p.sig <- 0

  sub.g <- domains
  
  nms <- names(sig.s)
  s.s <- rep('', nrow(domains))
  d.s <- rep('', nrow(domains))
  p.s <- rep('', nrow(domains))
  for (i in 1:length(sig.s)) {
    chr <- nms[i]

    cur.ids <- sig.s[[i]]
    sub.g$s.sig <- sub.g$s.sig + ifelse(sub.g$pfam_domain %in% cur.ids, 1, 0)
    sub.g[, paste(chr, '.s', sep='')] <- ifelse(sub.g$pfam_domain %in% cur.ids, 1, 0)
    cur.s <- ifelse(sub.g$pfam_domain %in% cur.ids, chr, '~')
    s.s <- paste(s.s, cur.s, sep='')

    cur.ids <- sig.d[[i]]
    sub.g$d.sig <- sub.g$d.sig + ifelse(sub.g$pfam_domain %in% cur.ids, 1, 0)
    sub.g[, paste(chr, '.d', sep='')] <- ifelse(sub.g$pfam_domain %in% cur.ids, 1, 0)
    cur.s <- ifelse(sub.g$pfam_domain %in% cur.ids, chr, '~')
    d.s <- paste(d.s, cur.s, sep='')

    cur.ids <- sig.p[[i]]
    sub.g$p.sig <- sub.g$p.sig + ifelse(sub.g$pfam_domain %in% cur.ids, 1, 0)
    sub.g[, paste(chr, '.p', sep='')] <- ifelse(sub.g$pfam_domain %in% cur.ids, 1, 0)
    cur.s <- ifelse(sub.g$pfam_domain %in% cur.ids, chr, '~')
    p.s <- paste(p.s, cur.s, sep='')

  }

  out.domains <- function(x, flds, file) {

    x$emp_05_bh.6.p <- sprintf("%.1e", x$emp_05_bh.6.p)

    x$n.pos.genes.6.p <- as.integer(x$n.pos.genes.6.p)
    x$n_genes.6.p <- as.integer(x$n_genes.6.p)
    x$pos_01.6.p <- as.integer(x$pos_01.6.p)
    x$pos_05.6.p <- as.integer(x$pos_05.6.p)
    x$n_sites.6.p <- as.integer(x$n_sites.6.p)

    for (i in 1:nrow(x)) {
      pf.d <- x[i, 'pfam_domain']
      print(pf.d)
      pos.sub <- subset(pos.sites, pfam_domain == pf.d)
      pos.ids <- unique(pos.sub$data_id)
      pos.gns <- subset(genes, data_id %in% pos.ids)
      pos.gns <- pos.gns[order(pos.gns$emp_01_bh.6), ]      

      x[i, 'n.pos.sites'] <- nrow(pos.sub)
      x[i, 'n.pos.genes'] <- length(pos.ids)
      x[i, 'pos.gene.names'] <- paste(pos.gns$gene_name, collapse=' ')
      x[i, 'all.pos.genes'] <- paste(pos.gns$gene_name, collapse=' ')
    }

    print(x[, c('pos.gene.names', 'all_gene_names.6.p')])

    x$genes <- x$pos.gene.names
    regex.res <- regexpr('(\\w*\\W*){1,4}', x$genes, perl=T)
    match.pos <- regex.res
    match.len <- attr(regex.res, 'match.length')
    longer.than <- match.len >= nchar(x$genes) - 1
    x$genes <- substr(x$genes, match.pos, match.len)
    x$genes <- gsub(' (\\w)', ', \\1', x$genes)
    x$genes <- R.oo::trim(x$genes)
    x$pos.gene.names <- x$genes

    
    x <- x[order(x$n.pos.sites, decreasing=T), ]
    x <- subset(x, n.pos.sites >= 10)
    x <- subset(x, n.pos.genes >= 3)

    xt <- xtable(x[, flds])
    xt <- color.columns(xt, columns=c('n.pos.sites', 'n.pos.genes'))
    xt <- xt.surround(xt, column='s.s', prefix='\\texttt{', suffix='}')  
    xt <- xt.surround(xt, column='p.s', prefix='\\texttt{', suffix='}')  
    xt <- xt.surround(xt, column='pos.gene.names', prefix='\\tiny{', suffix='}')  
    xt <- xt.surround(xt, column='pfam_domain', prefix='\\tiny{', suffix='}')  
    
    print(x[, flds])
    print.latex(xt, file=file)
  }

  sub.g$s.s <- s.s
  sub.g$d.s <- d.s
  sub.g$p.s <- p.s

  sub.g$pfam_descr <- gsub("interleukin-8", "IL-8", sub.g$pfam_descr)

  sub.flds <- c('pfam_domain', 'pfam_descr')
  sub.flds <- c(sub.flds, 's.s', 'p.s')

  p.flds <- c(sub.flds, 'n_genes.6.p', 'n_sites.6.p',
    'n.pos.genes', 'n.pos.sites', 'pos.gene.names'
    #, 'all.pos.genes', 'all_gene_names.6.p'
    )

  immune.domains <- c(
    'PF00021',
    'PF00048',
    'PF00059',
    'PF00074',
    'PF00340',
    'PF00354',
    'PF00969',
    'PF01823',
    'PF00993',
    'PF00084',
    'PF00062',
    'PF00047',
    'PF07686',
    'PF02841',
    'PF02886',
    'PF00530'
  )
  protease.domains <- c(
    'PF00656',
    'PF00246',
    'PF07859',
    'PF00089'
  )
  protease.inhibitors <- c(
    'PF00031',
    'PF07678',
    'PF07677',
    'PF00079',
    'PF01835'
  )
  all.cat <- unique(c(immune.domains, protease.domains, protease.inhibitors))
  
  x <- subset(sub.g)
#  x <- subset(x, grepl("lectin", pfam_descr, ignore.case=T))
#  x <- subset(x, grepl("CEN", all_gene_names.6.p, ignore.case=T))
  x <- subset(x, pfam_type == 'Domain')
  x <- subset(x, M.p == 1 & m.p == 1)
  x <- x[order(x$n.pos.genes.6.p, decreasing=T), ]

  out.domains(subset(x, pfam_domain %in% immune.domains), flds=p.flds, file=scratch.f("domains_immune.txt"))
  out.domains(subset(x, pfam_domain %in% protease.domains), flds=p.flds, file=scratch.f("domains_protease.txt"))
  out.domains(subset(x, pfam_domain %in% protease.inhibitors), flds=p.flds, file=scratch.f("domains_inhibitors.txt"))
  out.domains(subset(x, !(pfam_domain %in% all.cat)), flds=p.flds, file=scratch.f("domains_other.txt"))
   
#  p <- ggplot(sub.g, aes(x=(1 - f_neg_05.6.p), y=f_pos_05.6.p))
#  p <- p + theme_bw()
#  p <- p + geom_abline(slope=1/10, linetype='dashed')
#  p <- p + geom_point(size=1, alpha=0.3)
#  out.f <- scratch.f("domain_scatter.pdf")
#  pdf(file=out.f)
#  print.ggplot(p)
#  dev.off()

#  x <- subset(sub.g, f_pos_05.6.p / (1 - f_neg_05.6.p) > 0.1)
#  x <- subset(x, pfam_type == 'Domain')
#  x <- x[order(x$f_pos_05.6.p, decreasing=T),]
#  print(head(x[, p.flds], n=20))

}

top.genes <- function() {
  out.f <- scratch.f("tmp.genes.list.Rdata")
  if (!file.exists(out.f)) {

    filter = 'default'
    genes <- get.all.merged(filter=filter,
      keep.once=c('gene_name', 'aln_length', 'dup_species_count',
        'chr_name', 'chr_start', 'chr_end', 
        'gc_100k', 'aln_length', 'slr_dnds'),
      keep.fields=c('hoch_bh', 'emp_10_bh', 'emp_01_bh', 'emp_05_bh', 'emp_05', 'slr_dnds')
    )
    pubs <- get.published.genes(include.na=F)
    genes <- merge(genes, pubs, all.x=T)
    nearby <- get.nearby.pairs()
    genes <- merge(genes, nearby, all.x=T)

    filter = 'stringent'
    genes.s <- get.all.merged(filter=filter,
      keep.once=c('gc_1mb'),
      keep.fields=c('hoch_bh', 'emp_10_bh', 'emp_01_bh', 'emp_05_bh', 'emp_05', 'slr_dnds')
    )
    genes <- merge(genes, genes.s, by=c('data_id'), all.x=TRUE, suffixes=c('.d', '.s'))

    # Get the max p-value across the three superorders
    psets <- c(1, 2, 3)
    max.p.d <- rep(-Inf, nrow(genes))
    for (pset in psets) {
      cur.p.d <- genes[, paste('emp_05.', pset, '.d', sep='')]
      max.p.d <-  pmax(max.p.d, cur.p.d, na.rm=T)
    }
    max.p.d <- ifelse(max.p.d == -Inf, 1, max.p.d)
    max.p.d.bh <- p.adjust(max.p.d, method='BH')
    genes$max.p.d <- max.p.d.bh

    max.p.s <- rep(-Inf, nrow(genes))
    for (pset in psets) {
      cur.p.s <- genes[, paste('emp_05.', pset, '.s', sep='')]
      max.p.s <-  pmax(max.p.s, cur.p.s, na.rm=T)
    }
    non.nas <- !is.na(genes[, 'gc_1mb'])
    max.p.s <- ifelse(max.p.s[non.nas] == -Inf, 1, max.p.s[non.nas])
    max.p.s.bh <- p.adjust(max.p.s, method='BH')
    genes$max.p.s <- 1
    genes[non.nas, 'max.p.s'] <- max.p.s.bh

    sig.ids <- function(filter, method, pset) {
      x <- get.enrichment.df(pset_id=pset, method=method, filter=filter)
      x[x$binary.score == 1, 'data_id']
    }

    get.list <- function(f) {
      list(
        P = sig.ids(f, 'emp_01_bh', 1),
        g = sig.ids(f, 'emp_05_bh', 2),
        L = sig.ids(f, 'emp_01_bh', 3),
        M = sig.ids(f, 'emp_01_bh', 6),
        m = sig.ids(f, 'emp_05_bh', 6),
        H = sig.ids(f, 'hoch_bh', 6),
        D = sig.ids(f, 'dnds_02', 6),
        i = sig.ids(f, 'core_05', 1),
        C = genes[genes$clark == 1, 'data_id'],
        N = genes[genes$nielsen == 1, 'data_id'],
        R = genes[genes$rhesus == 1, 'data_id'],
        K = genes[genes$kosiol == 1, 'data_id']
      )
    }
    stringent.l <- get.list('stringent')
    default.l <- get.list('default')

    save(stringent.l, default.l, genes, file=out.f)
  }
  load(out.f)

  get.grep <- function(ptrn) {
    sub.g <- subset(genes, grepl(ptrn, gene_name, perl=T))
    
    pub.chars <- c('N', 'R', 'C', 'K')

    s.s <- ''
    d.s <- ''
    p.s <- ''
    sub.g$s.sig <- 0
    sub.g$d.sig <- 0
    sub.g$pub.sig <- 0
    nms <- names(stringent.l)
    # Fill in the pub-based PSGs
    for (i in 1:length(stringent.l)) {
      chr <- nms[i]
      if (chr %in% pub.chars) {
        cur.ids <- stringent.l[[i]]
        sub.g$pub.sig <- sub.g$pub.sig + ifelse(sub.g$data_id %in% cur.ids, 1, 0)
        cur.s <- ifelse(sub.g$data_id %in% cur.ids, chr, ' ')
        sub.g[, chr] <- ifelse(sub.g$data_id %in% cur.ids, 1, 0)
        p.s <- paste(p.s, cur.s, sep='')
      }
    }    
    # Fill in stringent.
    for (i in 1:length(stringent.l)) {
      chr <- nms[i]
      if (chr %in% pub.chars) {
        next()
      }
      cur.ids <- stringent.l[[i]]
      sub.g$s.sig <- sub.g$s.sig + ifelse(sub.g$data_id %in% cur.ids, 1, 0)
      cur.s <- ifelse(sub.g$data_id %in% cur.ids, chr, ' ')
      sub.g[, paste(chr, '.s', sep='')] <- ifelse(sub.g$data_id %in% cur.ids, 1, 0)
      s.s <- paste(s.s, cur.s, sep='')
    }
    # Fill in default.
    nms <- names(default.l)
    for (i in 1:length(default.l)) {
      chr <- nms[i]
      if (chr %in% pub.chars) {
        next()
      }
      cur.ids <- default.l[[i]]
      sub.g$d.sig <- sub.g$d.sig + ifelse(sub.g$data_id %in% cur.ids, 1, 0)
      cur.s <- ifelse(sub.g$data_id %in% cur.ids, chr, ' ')
      sub.g[, paste(chr, '.d', sep='')] <- ifelse(sub.g$data_id %in% cur.ids, 1, 0)
      d.s <- paste(d.s, cur.s, sep='')
    }
    
    sub.g$s.s <- s.s
    sub.g$d.s <- d.s
    sub.g$p.s <- p.s
    sub.g
  }

  seen.ids <- c()
  assign('seen.ids', seen.ids, envir=.GlobalEnv)

  p.genes <- function(x, sig=0, xtra=c('emp_01_bh.6.s'), 
    sort.by=xtra, n=nrow(x), desc=F, skip.seen=T, file=NULL) {
    seen.ids <- get('seen.ids', envir=.GlobalEnv)
    x$chrs <- sprintf("%.1f", x$chr_start / 1000 / 1000)
    x <- subset(x, d.sig >= sig)
    x <- x[order(x[, sort.by], decreasing=desc), ]
    if (skip.seen) {
      x <- subset(x, !(data_id %in% seen.ids))
    }
    print(nrow(x))
    x <- head(x, n=n)
    out.cols <- c('gene_name', 'chr_name', 'chrs', 's.s', 'd.s', 'p.s', 'nearby_sib',
      xtra)
    print(x[, out.cols])
    if (skip.seen) {
      seen.ids <- c(seen.ids, x$data_id)
      assign('seen.ids', seen.ids, envir=.GlobalEnv)
    }
    if (!is.null(file)) {
      x$s.s <- gsub(' ', '~', x$s.s)
      x$d.s <- gsub(' ', '~', x$d.s)
      x$p.s <- gsub(' ', '~', x$p.s)
      x$nearby_sib <- ifelse(x$nearby_sib == 1, '+', '-')
      x$emp_01_bh.6.d <- sprintf("%.1e", x$emp_01_bh.6.d)
      x <- subset(x, select=out.cols)
      xt <- xtable(x)
      xt <- xt.surround(xt, column='s.s', prefix='\\texttt{', suffix='}')  
      xt <- xt.surround(xt, column='d.s', prefix='\\texttt{', suffix='}')  
      xt <- xt.surround(xt, column='p.s', prefix='\\texttt{', suffix='}')  
      print.latex(xt, file=file)
    }
  }

  all <- get.grep(".*")

  # OK, get the mammals top 10 for different methods:

  # Top Emp_01 Mammals genes w/ pub overlap
  top.m <- subset(all, pub.sig >= 2)
  p.genes(top.m, xtra=c('emp_01_bh.6.s'), n=10)

  top.pub <- subset(all, pub.sig >= 2 & s.sig == 0 & d.sig == 0)
  p.genes(top.pub, xtra=c('emp_01_bh.6.d'), n=10)

  # Top "Hochberg" genes
  top.h <- subset(all, H.s == 1)
  p.genes(top.h, xtra=c('hoch_bh.6.s'), n=5)

  # Top "independent" genes with evidence in each mammalian superorder
  top.i <- subset(all, P.s == 1 & g.s == 1 & L.s == 1)
  top.i$max.p <- pmax(top.i$emp_01_bh.1.s, top.i$emp_05_bh.2.s, top.i$emp_01_bh.3.s)
  #p.genes(top.i, xtra=c('max.p'), n=5)

  # Top Primates genes
  top.p <- subset(all, P.s == 1 & g.s == 0 & g.d == 0 & L.s == 0 & L.d == 0)
  p.genes(top.p, xtra=c('emp_01_bh.1.s'), n=5)

  top.l <- subset(all, L.s == 1 & g.s == 0 & g.d == 0 & P.s == 0 & P.d == 0)
  p.genes(top.l, xtra=c('emp_01_bh.3.s'), n=5)

  top.g <- subset(all, g.s == 1 & L.s == 0 & L.d == 0 & P.s == 0 & P.d == 0)
  p.genes(top.g, xtra=c('emp_05_bh.2.s'), n=5)

  # Top dN/dS genes
  #top.d <- subset(all, D.s == 1)
  #p.genes(top.d, xtra=c('slr_dnds.6.s'), n=5, desc=T)

  xy <- get.grep("^CENP.*")
  p.genes(xy, xtra=c('emp_01_bh.6.s'))

  xy <- get.grep("^KLK.*")
  p.genes(xy, xtra=c('emp_01_bh.6.s'))

  xy <- get.grep("^IFIT.*")
  p.genes(xy, xtra=c('emp_01_bh.6.s'))

  #xy <- get.grep("^ADAM.*")
  #p.genes(xy, xtra=c('emp_01_bh.6.s'))
  

  #xy <- get.grep("^SLC.*")
  #p.genes(xy, xtra=c('emp_01_bh.6.s'))

  xx <- subset(genes, !is.na(emp_05_bh.6.d))

  is.nearby <- as.factor(xx$nearby_sib == 1)
  is.pos <- as.factor(xx$emp_05_bh.6.d < 0.1)
  fis.res <- fisher.test(is.nearby, is.pos)
  print(fis.res)
  print(sprintf("%d genes (%d families), %d nearby (%d families), %d pos, %d both",
    nrow(xx),
    length(unique(xx$ensf)),
    sum(xx$nearby_sib == 1),
    length(unique(xx[xx$nearby_sib == 1, 'ensf'])),
    sum(xx$emp_05_bh.6.d < 0.1),
    sum(xx$nearby_sib == 1 & xx$emp_05_bh.6.d < 0.1)
  ))

  genes$ensf <- as.character(genes$ensf)

  fam.tbl <- table(genes$ensf)
  big.fams <- names(fam.tbl[fam.tbl > 3])
  xx <- subset(genes, ensf %in% big.fams)
  xx <- subset(xx, !is.na(emp_05_bh.6.d) & !is.na(nearby_sib))

  summarize.fam <- function(x, do.fis=F) {
    ensf <- ifelse(is.null(x[1, 'ensf']), "ENSF", x[1, 'ensf'])
    x <- subset(x, !is.na(nearby_sib) & !is.na(emp_05_bh.6.d))
    is.nearby <- as.factor(x$nearby_sib == 1)
    is.pos <- as.factor(x$emp_05_bh.6.d < 0.1)
    fis.p <- ''
    fis.est <- ''
    if (length(unique(is.nearby)) > 1 && length(unique(is.pos)) > 1 && do.fis) {
      fis.res <- fisher.test(is.nearby, is.pos, alternative='g')
      fis.p <- fis.res$p.value
      fis.est <- fis.res$estimate
    }
     
    n.pos <- sum(x$emp_05_bh.6.d < 0.1)
    pos.genes <- subset(x, emp_05_bh.6.d < 0.1)
    nrby.genes <- subset(x, nearby_sib == 1)

    pos.nrby <- subset(x, emp_05_bh.6.d < 0.1 & nearby_sib == 1)
    pos.nrby <- pos.nrby[order(pos.nrby$emp_05_bh.6.d),]

    data.frame(
      ensf = ensf,
      n.genes = nrow(x),
      n.nearby = nrow(nrby.genes),
      n.pos = nrow(pos.genes),
      n.both = sum(x$emp_05_bh.6.d < 0.1 & x$nearby_sib == 1),
#      f.nearby = nrow(nrby.genes) / nrow(x),
#      f.pos = nrow(pos.genes) / nrow(x),
#      f.both = sum(x$emp_05_bh.6.d < 0.1 & x$nearby_sib == 1) / nrow(x),
      fis.p = fis.p,
#      fis.est = fis.est,
      pos.genes = paste(head(pos.nrby$gene_name, n=4), collapse=', ')
    )    
  }

  fam.df <- ddply(xx, .(ensf), summarize.fam)
  fam.df <- fam.df[order(fam.df$n.both, decreasing=T),]
  fam.df <- subset(fam.df, n.both >= 3)
  print(fam.df)

  fam.slc <- summarize.fam(get.grep("^SLC.*"), T)
  fam.adam <- summarize.fam(get.grep("^ADAM.*"), T)
  fam.col <- summarize.fam(get.grep("^COL[1234].*"), T)
  fam.tlr <- summarize.fam(get.grep("^TLR.*"), T)
  fam.all <- summarize.fam(get.grep(".*"), T)

  fam.slc$ensf <- "Solute Carrier Family"
  fam.adam$ensf <- "ADAM Family"
  fam.col$ensf <- "Collagen"
  fam.tlr$ensf <- "Toll-like Receptors"
  fam.all$ensf <- "All Genes"

  fam.f <- scratch.f("psg_fams_ensf.txt")
  xt <- xtable(fam.df)
  xt <- xt.surround(xt, column='ensf', prefix='\\scriptsize{', suffix='}')
  xt <- xt.surround(xt, column='pos.genes', prefix='\\tiny{', suffix='}')
  print.latex(xt, file=fam.f)
  
  fam.tbl <- rbind.fill(fam.slc, fam.adam, fam.col, fam.tlr, fam.all)
  fam.tbl <- fam.tbl[order(fam.tbl$fis.p, decreasing=T),]

  fam.f <- scratch.f("psg_fams_custom.txt")
  xt <- xtable(fam.tbl)
  blw.p <- xt$fis.p < 0.05
#  xt <- xt.surround(xt, column='fis.p', rows=which(blw.p), prefix='\\bf{', suffix='}')  
  xt <- xt.surround(xt, column='ensf', prefix='\\scriptsize{', suffix='}')
  xt <- xt.surround(xt, column='pos.genes', prefix='\\tiny{', suffix='}')
  print.latex(xt, file=fam.f)

#  xy <- get.grep(".*")
#  xy <- subset(xy, emp_05_bh.6.d < 0.1 & nearby_sib == 1)
#  p.genes(xy, xtra=c('emp_05_bh.6.d'))
   
#  xy <- get.grep("^CENP.*")
#  p.genes(xy)

#  xy <- get.grep("^SLC.*")
#  p.genes(xy, 2)

  xy <- get.grep("^COL[4].*")
  xy2 <- get.grep("^MMP.*")
  out.f <- scratch.f("psg_col_mmp.txt")
  p.genes(rbind(xy, xy2), 1, xtra=c('emp_01_bh.6.d'), file=out.f, skip.seen=F)

 # xy <- get.grep("^(ACAN|COL|MMP|TIMP|ADAM|SPINK).*")
  xy <- get.grep("^(KLK|SPINK|MASP|MBL|MMP|TIMP).*")
  out.f <- scratch.f("psg_col_mmp.txt")
  p.genes(xy, 0, xtra=c('emp_01_bh.6.d'), file=out.f, skip.seen=F)

  return()

  xy <- get.grep("^TLR.*")
  p.genes(xy, 1)

  xy <- get.grep("^IL.*")
  p.genes(xy, 1)

  xy <- get.grep("^APOL.*")
  p.genes(xy, 1)

  xy <- get.grep("^EIF.*")
  p.genes(xy, 1)

  xy <- get.grep("^TAS.*")
  p.genes(xy, 1)

  xy <- get.grep("^ADAM.*")
  p.genes(xy, 1)

  xy <- get.grep("^SIGLEC.*")
  p.genes(xy, 1)

  xy <- get.grep("^CCL.*")
  p.genes(xy, 1)

  xy <- get.grep("^CCR.*")
  p.genes(xy, 1)

  xy <- get.grep("^(CXC|CX3C).*")
  p.genes(xy, 1)

  xy <- get.grep("^DEFB.*")
  p.genes(xy, 1)

  xy <- get.grep("^HB.*")
  p.genes(xy, 1)

  xy <- get.grep("^CYB.*")
  p.genes(xy, 1)

  xy <- get.grep("^HAV.*")
  p.genes(xy, 0)

  xy <- get.grep("^BRCA.*")
  p.genes(xy, 0)

  return()

#  old.genes <- genes
#  genes <- subset(genes, conversions == 1 | nearby_sib == 1)
#  print(nrow(genes))
  

  #xy <- get.grep(".*")
  #p.genes(xy, 0)

#  genes <- old.genes

#  xy <- get.grep(".*", method='emp_01')
#  xy <- tail(xy, n=30)
#  print(xy[, c('gene_name', 'sig.s', 'm_dnds', 'max.p')])
}


union.all <- function(...) {
  unique(unlist(list(...)))
}

intersect.all <- function(...) {
  l <- list(...)
  cum.intersect <- l[[1]]
  for (i in 2:length(l)) {
    cum.intersect <- intersect(cum.intersect, l[[i]])
  }
  cum.intersect
}

pval.correlations <- function() {
  genes <- get.genes.sites(filter='default')
  genes <- subset(genes, pset==6)

}

pub.to.alias <- function(pubs, factors=T) {
  fct <- factor(pubs,
    levels = c('clark', 'nielsen', 'rhesus', 'kosiol', 'sabeti'),
    labels = c(
      'Clark 2003',
      'Nielsen 2005',
      'RMGSAC 2007',
      'Kosiol 2008',
      'Sabeti 2009'
    ))
  if (factors) {
    fct
  } else {
    as.character(fct)
  }
}

psg.rocs <- function(zoom=F) {
  genes <- get.all.merged(filter='stringent',
    keep.once=c('gene_name'),
    keep.fields=c('emp_005', 'emp_10', 'emp_01', 'emp_01_bh', 'hoch_p', 'slr_dnds')
  )
  pubs <- get.published.genes()
  genes <- merge(genes, pubs, all.x=T)

  flds <- c('slr_dnds.1', 'slr_dnds.6', 'emp_01.6', 'hoch_p.6')
#   flds <- c('emp_01.1', 'emp_01.3', 'emp_01.6')
#  flds <- c('hoch_p.6', 'emp_01.6', 'emp_10.6', 'slr_dnds.6')
  pubs <- c('clark', 'nielsen', 'rhesus', 'kosiol')

  out.df <- data.frame()
  lines.df <- data.frame()
  vlines.df <- data.frame()
  for (pub in pubs) {
    pub.df <- data.frame()
    for (fld in flds) {
      print(paste(pub, fld))
      # Sort by the fld, take a cumsum of truth values.
      cur.vals <- genes
      cur.vals$score <- cur.vals[, fld]
      cur.vals$truth <- as.integer(cur.vals[, pub])
      cur.vals <- subset(cur.vals, select=c('data_id', 'score', 'truth', 'emp_01_bh.6'))
      cur.vals <- subset(cur.vals, !is.na(truth))

      desc <- F
      if (grepl('slr', fld)) {
        desc <- T
      }
      cur.vals <- cur.vals[order(cur.vals$score, decreasing=desc),]
      cur.vals$yy <- cumsum(cur.vals$truth)
      cur.vals$xx <- 1:nrow(cur.vals)
      cur.vals$fld <- fld
      cur.vals$pub <- pub
      pub.df <- rbind(pub.df, cur.vals)

      if (fld == 'emp_01.6') {
        # Find the emp_01 value corresponding to BH < 0.1.
        cur.vals <- cur.vals[order(cur.vals$emp_01_bh.6),]
        max.ind <- max(which(cur.vals$emp_01_bh.6 < 0.1))
        print(max.ind)
        max.raw.pval <- cur.vals[max.ind, 'emp_01.6']
        print(max.raw.pval)
        vlines.df <- rbind(vlines.df, data.frame(
          pub = pub,
          xx = max.ind
        ))
      } 
    }
    max.x <- max(cur.vals$xx)
    max.y <- max(cur.vals$yy)
    lines.df <- rbind(lines.df, data.frame(
      pub = pub,
      slope = max.y/max.x
    ))

    if (zoom) {
      pub.df$xx <- pmin(pub.df$xx, max(pub.df$xx)/10)
      pub.df$yy <- pmin(pub.df$yy, max(pub.df$yy)/2)
    }
    out.df <- rbind(out.df, pub.df)
  }
  print(nrow(out.df))
  print(head(out.df))

  out.df$fld <- pval.to.alias(out.df$fld)
  out.df$fld <- factor(out.df$fld)
  out.df$pub <- pub.to.alias(out.df$pub)

  lines.df$pub <- pub.to.alias(lines.df$pub)
  vlines.df$pub <- pub.to.alias(vlines.df$pub)

  print(unique(out.df$fld))
  print(unique(out.df$pub))

  out.df <- subset(out.df, xx %% 20 == 0)

  p <- ggplot(out.df, aes(x=xx, y=yy, colour=fld, linetype=fld))
  p <- generic.opts(p)
  p <- p + geom_abline(data=lines.df, intercept=0, aes(slope=slope), colour='gray', alpha=1)
  p <- p + geom_vline(data=vlines.df, aes(xintercept=xx), linetype='dashed', colour='gray')
  p <- p + geom_line()
  p <- p + facet_wrap('pub', scales="free")
  p <- p + coord_equal()
  p <- p + opts(
    strip.background = theme_blank()  
  )

  out.f <- scratch.f(sprintf("psg_roc.pdf"))
  pdf(file=out.f, width=8, height=6)
  print.ggplot(p)
  dev.off()
}

get.published.genes <- function(include.na=T, force.fetch=F) {
  out.f <- scratch.f(sprintf("merged_all.Rdata"))

  if (!file.exists(out.f) || force.fetch) {
    print("Merging all pub genes...")
    df.list <- list(
      clark = get.clark.genes(force.fetch=force.fetch),
      nielsen = get.nielsen.genes(force.fetch=force.fetch),
      rhesus = get.rhesus.genes(force.fetch=force.fetch),
      kosiol = get.kosiol.genes(force.fetch=force.fetch),
      conversions = get.conversions()
    )

    genes <- get.genes()
    genes <- subset(genes, select=c('data_id'))

    print(paste("Dup genes: ", sum(duplicated(genes$data_id))))

    pubnames <- names(df.list)
    for (i in 1:length(df.list)) {
      print(pubnames[i])
      pub.genes <- df.list[[i]]
      pub.genes[, pubnames[i]] <- pub.genes[, 'truth']
      print(paste("Dup", pubnames[i], sum(duplicated(pub.genes$data_id))))
      pub.sub <- subset(pub.genes, select=c('data_id', pubnames[i]))
      genes <- merge(genes, pub.sub, by='data_id', all.x=T)
    }
    pub.genes <- genes
    save(pub.genes, file=out.f)
  }
  load(out.f)

  x <- pub.genes
  if (!include.na) {
    x$clark <- ifelse(is.na(x$clark), FALSE, x$clark)
    x$nielsen <- ifelse(is.na(x$nielsen), FALSE, x$nielsen)
    x$rhesus <- ifelse(is.na(x$rhesus), FALSE, x$rhesus)
    x$kosiol <- ifelse(is.na(x$kosiol), FALSE, x$kosiol)
    x$conversions <- ifelse(is.na(x$conversions), FALSE, x$conversions)
  }

  #print("Clark nielsen rhesus kosiol conversions")
  #print(sum(pub.genes$clark, na.rm=T))
  #print(sum(pub.genes$nielsen, na.rm=T))
  #print(sum(pub.genes$rhesus, na.rm=T))
  #print(sum(pub.genes$kosiol, na.rm=T))
  #print(sum(pub.genes$conversions, na.rm=T))

  x
}

get.sabeti.genes <- function() {
  out.f <- scratch.f(sprintf("merged_sabeti.Rdata"))

  if (!file.exists(out.f)) {
    print("Calcultaing Sabeti...")

    sabeti.genes <- read.delim("sabeti-genes.txt", header=F, stringsAsFactors=F, col.names='id')
    ncbi.gene.info <- get.ncbi.gene.info()
    sabeti.ncbi <- merge(sabeti.genes, ncbi.gene.info, by.x='id', by.y='Symbol')
    print(nrow(sabeti.ncbi))

    genes <- get.genes()
    genes$ensg <- genes$ref_gene_id

    mrgd <- merge(genes, sabeti.ncbi, by='ensg', all.x=T)
    print(nrow(mrgd))    

    mrgd$score <- ifelse(is.na(mrgd$id), 0, 1)
    sabeti.genes <- subset(mrgd, select=c('data_id', 'score'))
    
    save(sabeti.genes, file=out.f)
  }
  load(out.f)
  sabeti.genes$truth <- as.logical(sabeti.genes$score)

  sabeti.genes
}

get.clark.genes <- function(force.fetch=F) {
  out.f <- scratch.f(sprintf("merged_clark.Rdata"))

  if (!file.exists(out.f) || force.fetch) {
    #print("Calcultaing Clark...")

    ncbi.gene.info <- get.ncbi.gene.info()
    clark.genes <- read.delim("clark-genes.tsv", header=T, sep="\t")
    clark.genes <- subset(clark.genes, !is.na(LocusLink.ID) & LocusLink.ID != "")
  
    print(sprintf("Clark orig: %d", nrow(clark.genes)))

    clark.ncbi <- merge(clark.genes,ncbi.gene.info,by.x="LocusLink.ID",by.y="GeneID")

    # Get truth values for Clark tests.
    clark.ncbi$pval.c <- clark.ncbi$`Model.2..M2..p.value.chimp`
    clark.ncbi$pval.h <- clark.ncbi$`Model.2..M2..p.value.human`
    clark.ncbi$pval2.h <- clark.ncbi$`P.value.model.1..dN.dS..human..one.sided..test.for.pos.selection.`
    clark.ncbi$pval2.c <- clark.ncbi$`P.value.model.1..dN.dS..chimp..one.sided..test.for.pos.selection.`
    clark.ncbi$ensg_id <- clark.ncbi$ensg
    
    # Get E-SLR scores and merge.
    genes <- get.genes()
    genes$ensg_id <- genes$ref_gene_id

    mrgd <- merge(genes, clark.ncbi, by='ensg_id')
    #print(nrow(mrgd))

    clark.genes <- subset(mrgd, select=c('data_id', 'pval.c', 'pval.h', 'pval2.c', 'pval2.h'))
    save(clark.genes, file=out.f)
  }
  load(out.f)

 #print(nrow(clark.genes))
  clark.genes <- within(clark.genes,
    truth <- pval.c < 0.01 | pval.h < 0.01
  )

  # Order by score decreasing (to keep higher-scored PSGs) then remove
  # duplicates by data_id.
  clark.genes <- clark.genes[order(as.integer(clark.genes$truth), decreasing=T),]
  clark.genes <- clark.genes[!duplicated(clark.genes$data_id),]    

  print(sprintf("Clark: %d mapped, %d psg", nrow(clark.genes), sum(clark.genes$truth)))

  clark.genes
}

get.ncbi.gene.info = function() {
  ncbi.gene.info <- read.delim("Homo_sapiens.gene_info", sep="\t", header=F, stringsAsFactors=F)

  names(ncbi.gene.info) <- c("tax_id", "GeneID", "Symbol", "LocusTag", "Synonyms",
    "dbXrefs", "chromosome", "map.location", "description", "type", "symbol.from.nomenclature", 
    "full.name.from.nomenclature", "nomenclature.status", "other.designations",
    "modification.date"
  )

  ncbi.ensgs <- sapply(strsplit(ncbi.gene.info$dbXrefs,"\\|"), function(x){
    a <- grep("Ensembl", x, value=T)
    a=sub("Ensembl:", "", a)
    a[1]
  })
  ncbi.gene.info$ensg <- ncbi.ensgs
  ncbi.gene.info <- subset(ncbi.gene.info,!is.na(ensg) & ensg != "")
  return(ncbi.gene.info)
}

get.rhesus.genes <- function(force.fetch=F) {
  out.f <- scratch.f(sprintf("merged_rhesus.Rdata"))

  if (!file.exists(out.f) || force.fetch) {
    print("Calculating rhesus...")

    n.df <- readLines("rhesus-genes.txt")
    n.df <- n.df[-1]

    n.split <- strsplit(n.df, ' ')
    #print(head(n.split))

    n.df <- ldply(n.split, function(x) {
      data.frame(
        ncbi_id = x[[2]],
        gene_name = x[[3]]
      )
    })
    print(sprintf("Rhesus orig: %d", nrow(n.df)))
    #print(nrow(n.df))

    genes <- get.genes()
    genes$ensg <- genes$ref_gene_id
    mrg.name <- merge(genes, n.df, by=c('gene_name'))
    #print(nrow(mrg.name))

    ensg.nm <- read.table("ensg_to_nm.txt", header=T, sep="\t")
    colnames(ensg.nm) <- c('ensg', 'ncbi_id')
    n.df.ensg <- merge(n.df, ensg.nm, all.x=T)

    mrg.ensg <- merge(genes, n.df.ensg, by=c('ensg'))
    #print(nrow(mrg.ensg))

    mrg.all <- rbind.fill(mrg.name, mrg.ensg)
    #print(nrow(mrg.all))
    mrg.all <- mrg.all[!duplicated(mrg.all$data_id),]
    #print(nrow(mrg.all))

    genes <- get.genes()
    genes$score <- 0
    genes[genes$data_id %in% mrg.all$data_id, 'score'] <- 1

    rhesus.genes <- subset(genes, select=c('data_id', 'gene_name', 'score'))
    save(rhesus.genes, file=out.f)
  }
  load(out.f)
  
  rhesus.genes$truth <- rhesus.genes$score == 1

  print(sprintf("Rhesus: %d mapped, %d psg", nrow(rhesus.genes), sum(rhesus.genes$truth)))

  rhesus.genes
}

get.nielsen.genes <- function(force.fetch=F) {
  out.f <- scratch.f(sprintf("merged_nielsen.Rdata"))

  if (!file.exists(out.f) || force.fetch) {
    n.df <- read.csv("nielsen-genes.csv", header=T, stringsAsFactors=F)
    colnames(n.df) <- c('id', 'ensg_id', 'refseq_id', 'gene_name', 'length',
      'lnl_ps', 'lnl_null', 'pval', 'restrictive'
    )
    print(sprintf("Nielsen orig: %d", nrow(n.df)))
    #print(head(n.df))
    #n.df <- subset(n.df, !is.na(pval))
    n.df$lrt <- 2 * (n.df$lnl_ps - n.df$lnl_null)
    n.df$lrt <- ifelse(n.df$pval == 1, -n.df$lrt, n.df$lrt)
  
    n.df <- n.df[order(n.df$lrt, decreasing=T),]

    #print(head(subset(n.df, select=c('id', 'ensg_id', 'gene_name',
    #'lrt', 'pval', 'restrictive')), n=50))

    x <- subset(n.df, lrt > 0)
    x$truth <- x$pval < 0.1
    #print(paste("Number of Nielsen true positives: ",(sum(x$truth))))
    
    # Now, try to merge nielsen's genes into ours.
    genes <- get.genes()
    genes$ensg_id <- genes$ref_gene_id

    mrg.name <- merge(genes, n.df, by=c('gene_name'))
    #print(nrow(mrg.name))

    mrg.ens <- merge(genes, n.df, by=c('ensg_id'))
    #print(nrow(mrg.ens))

    # Try using the ensg-to-nm mapping from Ensembl biomart.
    ensg.nm <- read.table("ensg_to_nm_np.txt", header=T, sep="\t", stringsAsFactors=F)
    colnames(ensg.nm) <- c('ensg', 'nm', 'np', 'xp')
    ensg.nm$refseq_id <- ensg.nm$nm
    n.df[n.df$refseq_id == '', 'refseq_id'] <- NA
    nm.mrg <- merge(ensg.nm, n.df, by=c('refseq_id'))
    nm.mrg$ensg_id <- nm.mrg$ensg
    nm.mrg <- merge(genes, nm.mrg, by=c('ensg_id'))    

    mrg.all <- rbind.fill(nm.mrg, mrg.ens, mrg.name)
    #print(nrow(mrg.all))
    mrg.all <- mrg.all[!duplicated(mrg.all$data_id),]
    #print(nrow(mrg.all))

    nielsen.genes <- subset(mrg.all, select=c('data_id', 'lrt', 'pval'))
    #print(head(nielsen.genes))
    save(nielsen.genes, file=out.f)
  }
  load(out.f)

  nielsen.genes$truth <- nielsen.genes$lrt > 1.67

  print(sprintf("Nielsen: %d mapped, %d psg", nrow(nielsen.genes), sum(nielsen.genes$truth, na.rm=T)))

  nielsen.genes
}

get.kosiol.genes <- function(force.fetch=F) {
  out.f <- scratch.f(sprintf("merged_kosiol.Rdata"))

  if (!file.exists(out.f) || force.fetch) {
    print("Calcultaing Kosiol...")

    # Load file.
    k.df <- read.table("kosiol-gene-table.tsv", header=T, sep="\t")
    #print(nrow(k.df))
    #print(nrow(subset(k.df, score >= 700)))

    print(sprintf("Kosiol orig: %d (%d psg)", nrow(k.df), nrow(subset(k.df, score > 400))))

    # First, take chr coords and lift over.
    k.df <- lift.over(k.df, 'hg18ToHg19', start.s='chromStart', end.s='chromEnd', chr.s='chrom')

    #print(nrow(subset(k.df, is.na(chromStart))))

    # Now merge her genes with our genes by position.
    genes <- get.genes()
    genes <- subset(genes, !is.na(chr_start))

    chr_name <- genes$chr_name
    chr_start <- genes$chr_start
    chr_end <- genes$chr_end
    chr_strand <- genes$chr_strand
    genes.r <- GRanges(
      seqnames = chr_name,
      ranges = IRanges(chr_start, chr_end),
      data_id = genes$data_id,
      gene_name = genes$gene_name
    )
    #print(genes.r)

    k.r <- GRanges(
      seqnames = k.df$chrom,
      ranges = IRanges(k.df$chromStart, k.df$chromEnd),
      score = k.df$score
    )
    #print(k.r)

    match.ind <- genes.r %in% k.r
    #print(sum(match.ind))

    olap = findOverlaps(k.r, genes.r)

    k.scores <- values(k.r)[["score"]][queryHits(olap)]
    data_ids <- values(genes.r)[["data_id"]][subjectHits(olap)]
    gene_names <- values(genes.r)[["gene_name"]][subjectHits(olap)]

    kosiol.genes <- data.frame(
      score = k.scores,
      data_id = data_ids,
      gene_name = gene_names
    )

    save(kosiol.genes, file=out.f)
  }
  load(out.f)

  # Order by score decreasing (to keep higher-scored PSGs) then remove
  # duplicates by data_id.
  kosiol.genes <- kosiol.genes[order(kosiol.genes$score, decreasing=T),]
  kosiol.genes <- kosiol.genes[!duplicated(kosiol.genes$data_id),]    

  kosiol.genes$truth <- kosiol.genes$score > 400

  print(sprintf("Kosiol: %d mapped, %d psg", nrow(kosiol.genes), sum(kosiol.genes$truth, na.rm=T)))

  kosiol.genes
}

get.conversions <- function() {
  out.f <- scratch.f("merged_conversions.Rdata")

  if (!file.exists(out.f)) {
    ensg.nm <- read.table("ensg_to_nm_np.txt", header=T, sep="\t", stringsAsFactors=F)
    colnames(ensg.nm) <- c('ensg', 'nm', 'np', 'xp')

    # Read in the gene conversions from
    # Benovoy and Drouin 2009 Genetics
    conv.df <- read.csv("mmc1.csv", header=T, stringsAsFactors=F)
    #print(str(conv.df))
    # Make a longer list of all genes involved in any conversion.
    refseq_ids <- unique(c(conv.df$`Gene1.name`, conv.df$`Gene.2.name`))    
    # Strip off the .1, .2, etc.
    print(head(refseq_ids))
    refseq_ids <- sub("\\.\\d", "", refseq_ids)
    print(head(refseq_ids))

    conv.df <- data.frame(
      refseq_id = refseq_ids,
      stringsAsFactors=F
    )

    print(head(ensg.nm))
    conv.np <- merge(ensg.nm, conv.df, by.x='np', by.y='refseq_id')
    conv.xp <- merge(ensg.nm, conv.df, by.x='xp', by.y='refseq_id')

    all.ensg <- unique(c(conv.np$ensg, conv.xp$ensg))
    ensg.df <- data.frame(
      ref_gene_id=all.ensg,
      gene_conversion=1,
      stringsAsFactors=F
    )
    genes <- get.genes()
    genes <- subset(genes, select=c('data_id', 'gene_name', 'ref_gene_id'))
    genes <- merge(genes, ensg.df, all.x=T)
    
#    print(nrow(subset(genes, !is.na(gene_conversion))))
    genes$gene_conversion <- ifelse(is.na(genes$gene_conversion), 0, 1)

    gene.conversions <- genes
    save(gene.conversions, file=out.f)
  }  
  load(out.f)

  gene.conversions$truth <- gene.conversions$gene_conversion
  gene.conversions$gene_conversion <- NULL

  gene.conversions <- gene.conversions[order(gene.conversions$gene_name),]
  #print(subset(gene.conversions, truth == 1))
  #print(str(gene.conversions))

  gene.conversions
}

add.family.ids <- function(genes) {
  out.f <- scratch.f("merged_families.Rdata")
  if (!file.exists(out.f)) {
    genes <- get.genes()

    ens.df <- read.table('ens_to_fam_id.txt', header=T, sep="\t")
    colnames(ens.df) <- c('ref_gene_id', 'ensf')
    print(head(ens.df))

    genes <- merge(genes, ens.df, all.x=T)

    genes.family <- subset(genes, select=c('data_id', 'ensf'))
    save(genes.family, file=out.f)
  }
  load(out.f)

  genes <- merge(genes, genes.family, by=c('data_id'), all.x=T)
  genes
}

get.go.genes <- function(go_ids) {
  go.f <- scratch.f('godata_6 default hoch_bh FALSE NA.Rdata')
  load(go.f)

  getGenesForTerms(root.godata, go_ids)
}

test.better <- function() {

  print("Neuro")
  svg(file=scratch.f("dnds_neuro.svg"), width=12, height=4)
  vplayout(3, 1)
  go.ids <- c('GO:0051960')
  lbl.genes <- c('CDK5RAP2', 'ASPM', 'XRCC4', 'SPP1', 'MAPT', 'TTC3', 'HAP1', 'XRCC6')
  print.ggplot(test.psg(1, 2, go.ids, lbl.genes), vp=subplot(1,1))
  print.ggplot(test.psg(1, 3, go.ids, lbl.genes), vp=subplot(2,1))
  print.ggplot(test.psg(1, 4, go.ids, lbl.genes), vp=subplot(3,1))
  dev.off()                

  print("Sperm")
  svg(file=scratch.f("dnds_spermatogenesis.svg"), width=12, height=4)
  vplayout(3, 1)
  go.ids <- c('GO:0007283')
  lbl.genes <- c('TXNDC8', 'TAF7', 'ADCY10', 'CATSPER2', 'CYLC1',
    'AKAP14', 'SYCP3', 'SPAG16', 'APOB', 'ODF4', 'FANCG'
  )
  print.ggplot(test.psg(1, 2, go.ids, lbl.genes), vp=subplot(1,1))
  print.ggplot(test.psg(1, 3, go.ids, lbl.genes), vp=subplot(2,1))
  print.ggplot(test.psg(1, 4, go.ids, lbl.genes), vp=subplot(3,1))
  dev.off()                

  print("Chromo")
  svg(file=scratch.f("dnds_chromosome.svg"), width=12, height=4)
  vplayout(3, 1)
  go.ids <- c('GO:0000236', 'GO:0007059')
  lbl.genes <- c('TEX11', 'SGOL2', 'CASC5', 'CDK5RAP2', 'BRCA2',
    'CENPC1', 'BRCA1', 'CENPT', 'CENPE')
  print.ggplot(test.psg(1, 2, go.ids, lbl.genes), vp=subplot(1,1))
  print.ggplot(test.psg(1, 3, go.ids, lbl.genes), vp=subplot(2,1))
  print.ggplot(test.psg(1, 4, go.ids, lbl.genes), vp=subplot(3,1))
  dev.off()

  print("Immune")
  svg(file=scratch.f("dnds_immune.svg"), width=12, height=4)
  vplayout(3, 1)
  go.ids <- c('GO:0045087')
  lbl.genes <- c('DEFB118', 'CLEC4D', 'SAA1', 'PLUNC', 'CLEC4A', 'TLR5', 
    'CXCL16', 'SLAMF7', 'IFNK', 'MBL2')
  print.ggplot(test.psg(1, 2, go.ids, lbl.genes), vp=subplot(1,1))
  print.ggplot(test.psg(1, 3, go.ids, lbl.genes), vp=subplot(2,1))
  print.ggplot(test.psg(1, 4, go.ids, lbl.genes), vp=subplot(3,1))
  dev.off()

}

test.psg <- function(p.a, p.b, go_ids, lbl.genes='all') {
  genes <- get.all.merged(filter='default',
    keep.once=c('gene_name', 'ref_protein_id'),
    keep.fields=c('emp_005', 'emp_10', 'emp_01', 'emp_05', 'emp_01_bh', 
      'emp_05_bh', 'hoch_p', 'slr_dnds',
      'slr_dnds_z')
  )

  genes$slr.a <- genes[, paste('slr_dnds.', p.a, sep='')]
  genes$slr.b <- genes[, paste('slr_dnds.', p.b, sep='')]
  genes$p.a <- genes[, paste('emp_01_bh.', p.a, sep='')]
  genes$p.b <- genes[, paste('emp_01_bh.', p.b, sep='')]

  genes$p.x <- genes$slr.a
  genes$p.y <- genes$slr.b

  p.df <- genes
  p.df$p.x <- round_any(p.df$p.x, 0.02)
  p.df$p.y <- round_any(p.df$p.y, 0.02)

  p <- ggplot(p.df, aes(x=p.x, y=p.y))
  p <- p + theme_bw()
  
  p <- p + geom_point(stat="sum", colour=NA, aes(fill=..prop..), size=2, alpha=1, shape=22)
  p <- p + scale_fill_gradient("Density", low='white', high=gray(0.3), trans='log10')

  p <- p + geom_abline(slope=1, colour='black', linetype='dashed', alpha=0.7)

  xx <- get.go.genes(go_ids)
  sub.inn <- subset(genes, ref_protein_id %in% xx)
  sub.inn2 <- subset(sub.inn, !is.na(p.a) & p.a < 0.1)
  sub.inn3 <- subset(sub.inn, !is.na(p.a) & p.b < 0.1)

  print(sub.inn2[order(sub.inn2$p.a), 'gene_name'])
  print(sub.inn3[order(sub.inn3$p.b), 'gene_name'])

  if (nrow(sub.inn) > 0) {
    p <- p + geom_point(data=sub.inn, size=0.75, colour='black', alpha=0.8)
  }
  
#  shp.a <- switch(p.a,
#    19,
#    19,
#    19
     #  )
#  shp.b <- switch(p.b,
#    19,
#    19,
#    19
#  )
   shp.a <- 19
   shp.b <- 19
  clr.a <- switch(p.a,
    'blue',
    'green',
    'red',
    'orange'
  )
  clr.b <- switch(p.b,
    'blue',
    'green',
    'red',
    'orange'
  )

  sz <- 2.5
  if (nrow(sub.inn2) > 0) {
    p <- p + geom_point(data=sub.inn2, aes(x=p.x+0.005, y=p.y-0.005), 
      size=sz, colour=clr.a, fill=clr.a, shape=shp.a, alpha=0.5
    )
  }
  if (nrow(sub.inn3) > 0) {
    p <- p + geom_point(data=sub.inn3, aes(x=p.x-0.005, y=p.y+0.005),
      size=sz, colour=clr.b, fill=clr.b, shape=shp.b, alpha=0.5
    )
  }

  if (lbl.genes == 'all') {
    sub.txt <- rbind(sub.inn2, sub.inn3)
  } else {
    sub.txt <- subset(sub.inn, gene_name %in% lbl.genes)
  }
  p <- p + geom_text(data=sub.txt, aes(label=gene_name, x=p.x-0.02), hjust=1, vjust=0.5)

  x.s <- sprintf("%s dN/dS", pset.to.alias(p.a))
  y.s <- sprintf("%s dN/dS", pset.to.alias(p.b))
  p <- p + scale_x_continuous(x.s, limits=c(0, 1.2))
  p <- p + scale_y_continuous(y.s, limits=c(0, 1.2))

  descr <- getTermsDefinition(go_ids, 'BP', numChar=40)
  descr <- paste(descr, collapse=' / ')

  p <- p + opts(
   title = descr,
   legend.position = 'none'
  )

  return(p)

#  out.f <- scratch.f("test.png")
#  png(file=out.f, width=1200, height=1200)
#  print.ggplot(p)
#  dev.off()

}

thesis.tables <- function() {

  tab1 <- read.csv(file="gorilla_lrt_results.csv", header=F)
  print(tab1)

  xt <- xtable(tab1)
  xt <- color.columns(xt, columns=c('V2', 'V3', 'V4', 'V5', 'V6', 'V7'))
  print.latex(xt, file=scratch.f("gorilla_lrt_results.txt"), na.string='asdf')

  process.go.t <- function(f) {
    x <- read.csv(file=paste(f, ".csv", sep=''), header=F)

    colnames(x) <- c('name', 'ann', 'sig', 'exp', 'id', 'def',
     'fet', 'topgo', 'goseq', 'length', 'genes')
    x$name <- NULL

    regex.res <- regexpr('(\\w*\\W*){1,5}', x$genes, perl=T)
    match.pos <- regex.res
    match.len <- attr(regex.res, 'match.length')
    longer.than <- match.len >= nchar(x$genes) - 1
    x$genes <- substr(x$genes, match.pos, match.len)
    x$genes <- gsub(' (\\w)', ', \\1', x$genes)
    x$genes <- R.oo::trim(x$genes)
    x$genes <- paste("\\tiny{", x$genes, "}", sep='')

    fet.sig <- x$fet < 0.05
    topgo.sig <- x$topgo < 0.05
    goseq.sig <- x$goseq < 0.05

    x$fet <- sprintf("%.1e", x$fet)
    x$topgo <- sprintf("%.1e", x$topgo)
    x$goseq <- sprintf("%.1e", x$goseq)

    xt <- xtable(x)
    xt <- color.rows(xt, rows=which(!fet.sig), columns='fet', color.bg=F)
    xt <- color.rows(xt, rows=which(!topgo.sig), columns='topgo', color.bg=F)
    xt <- color.rows(xt, rows=which(!goseq.sig), columns='goseq', color.bg=F)
    print.latex(xt, file=scratch.f(paste(f, '.txt', sep='')))
  }

  process.go.t("gorilla_go_1")
  process.go.t("gorilla_go_2")
  process.go.t("gorilla_go_3")
  process.go.t("gorilla_go_4")
  process.go.t("gorilla_go_5")
  process.go.t("gorilla_go_6")
  process.go.t("gorilla_go_7")
  process.go.t("gorilla_go_8")

  # Load the comparison of published studies.
  tab <- read.csv(file="gorilla_studies.csv", header=F)
  #print(str(tab))
  tab$V3 <- NULL # Remove bg species
  tab$V5 <- NULL # Remove aligner
  tab$V8 <- NULL # Remove absolute number of accel. genes
  tab$V12 <- NULL # 

  tab$V11 <- gsub(', ', paste(' \\\\', 'newline ',sep=''), tab$V11)
  tab$V11 <- paste('\\tiny{', tab$V11, '}', sep='')
#  tab$V11 <- paste('\\singlespacing{', tab$V11, '}', sep='')
  #print(tab$V11)

  xt <- xtable(tab)
  print.latex(xt, file=scratch.f("gorilla_studies.txt"))

  tab <- read.csv(file="gorilla_top_genes.csv", header=F, stringsAsFactors=F)

  tab$V1 <- gsub('_', ' ', tab$V1)
  tab$V1 <- gsub('branchsite', 'BS', tab$V1)
  tab$V1 <- gsub('branch', '', tab$V1)
  tab$V1 <- gsub('parallel', "\\\\lrtmin", tab$V1)

  tab$V6 <- as.integer(tab$V6)
  tab$V7 <- as.integer(tab$V7)

  b.sig <- tab$V9 < 0.05
  b.sig.adj <- tab$V10 < 0.1

  bs.sig <- tab$V12 < 0.05
  bs.sig.adj <- tab$V13 < 0.1

  tab[duplicated(tab$V1), 'V1'] <- ''

  i <- 1
  while (i < nrow(tab)) {
    print(nrow(tab))
    if (i >= nrow(tab)) {
      break()
    }
    if (tab[i, 'V1'] != '' && i > 1) {
      tmp <- tab
      tab <- rbind.fill(tmp[1:(i-1),], data.frame(V1=NA))
      tab <- rbind.fill(tab, tmp[(i):nrow(tmp),])
      i <- i + 2
    } else {
      i <- i + 1
    }
    #print(i)
  }

  pv.fields <- c('V9', 'V10', 'V12', 'V13')
  for (fld in pv.fields) {
    tab[, fld] <- sprintf("%.1e", tab[, fld])
    if(any(tab[, fld] == 'NA')) {
      tab[tab[, fld] == 'NA', fld] <- ''
    }
    tab[, fld] <- paste("\\tiny{", tab[, fld], "}", sep='')
  }

  print(head(tab))
  print(str(tab))
  xt <- xtable(tab)
  xt <- xt.surround(xt, column='V1', prefix='\\tiny{', suffix='}', na.string='')
  xt <- xt.surround(xt, column='V2', prefix='\\tiny{\\gene{', suffix='}}', na.string='')
  xt <- color.columns(xt, columns=c('V4', 'V5', 'V8', 'V11'), na.string='') # Color dnds & LRTs
  
  xt <- color.rows(xt, rows=which(!b.sig), columns='V9', color.bg=F, na.string='')
  xt <- color.rows(xt, rows=which(!b.sig), columns='V10', color.bg=F, na.string='')
  xt <- color.rows(xt, rows=which(!bs.sig), columns='V12', color.bg=F, na.string='')
  xt <- color.rows(xt, rows=which(!bs.sig.adj), columns='V13', color.bg=F, na.string='')

  print.latex(xt, file=scratch.f("gorilla_top_genes.txt"), na.string='asdf') 
  
  tab <- read.csv(file="gorilla_dnds.csv", header=F, stringsAsFactors=F)
  print(tab)
  tab$V4 <- sprintf("%.1e",tab$V4)
  tab$V6 <- sprintf("%.1e",tab$V6)
  print(tab)
  xt <- xtable(tab)
  attr(xt, 'digits') <- c(0, 2, 2, 4, 2, 4, 2, 2)
  xt <- color.columns(xt, columns=c('V2', 'V3', 'V5', 'V7'), na.string='')
  print.latex(xt, file=scratch.f("gorilla_dnds.txt"), na.string='')
  

}

ren <- function(`_data`, ...){ 
  # Influenced by transform.data.frame 
  e <- unlist(list(...)) 
  inx <- match(e, names(`_data`)) 
  if(any(is.na(inx))) 
    stop("Couldn't find these columns in the data: ", 
         paste(as.vector(e)[is.na(inx)], collapse=" ") ) 
  names(`_data`)[inx] <- names(e) 
  return(`_data`) 
} 

# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.prettysize <- napply(names, function(x) {
                           capture.output(print(object.size(x), units = "auto")) })
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}

# shorthand
lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}