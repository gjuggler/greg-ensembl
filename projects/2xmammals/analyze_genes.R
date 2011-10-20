source("~/src/greg-ensembl/projects/2xmammals/analyze_mammals.R")
source("~/src/greg-ensembl/scripts/go_enrichments.R")

get.aln <- function(
  data_id, 
  aln_lo=1, 
  aln_hi=99999, 
  filename="test.pdf", 
  taxon_id=9606,
  include.psets = c(1, 2, 3, 6),
  keep.species = 'mammals',
  remove.blank.columns=F
) {
  library(R.oo)
  library(phylosim)
  source("~/src/greg-ensembl/projects/phylosim/PhyloSimPlots.R")

  con <- connect(dbname())
  cmd <- sprintf("select * from seqs where data_id=%d and aln_position between %d and %d order by aln_position",
    data_id, aln_lo, aln_hi)
  df <- dbGetQuery(con, cmd)
  disconnect(con)

  con <- connect(dbname())
  cmd <- sprintf("select * from genes where data_id=%d", data_id)
  gene <- dbGetQuery(con, cmd)
  disconnect(con)

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
    sites$source <- paste("SLR", pset.to.alias(sites$parameter_set_id))
    sites$score <- sites$lrt_stat
    sites$start <- sites$aln_position - aln_lo + 1
    sites$end <- sites$aln_position - aln_lo + 1

    dnds <- subset(sites, parameter_set_id == 6)
    if (nrow(dnds) > 0) {
      dnds$source <- paste("dN/dS", pset.to.alias(dnds$parameter_set_id))
      dnds$score <- dnds$omega
      sites <- rbind(sites, dnds)

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
  }
  mamms <- c(mamms, taxon_id)
  mamms <- unique(mamms)
  aln <- restrict.aln.to.seqs(aln, mamms)
  #aln <- aln.tx(aln)

  f.con <- file("~/src/greg-ensembl/projects/orthologs/compara_63_taxids.nh")
  str <- readLines(con=f.con)
  close(f.con)
  tree <- read.nhx.tree(str)
  non.mammals <- setdiff(tree$tip.label, mamms)
  tree <- drop.tip(tree, non.mammals)
  tree <- remove.branchlengths(tree)
  aln <- sort.aln.by.tree(aln, tree)

#  if (remove.blank.columns) {
#    aln <- remove.blank.columns(aln)
#  }

  tx.ids <- rownames(aln)
  species.ids <- taxid.to.alias(tx.ids, include.internals=T)
  rownames(aln) <- species.ids
  tree$tip.label <- taxid.to.alias(tree$tip.label, include.internals=T)

  out.w <- .2 * length(aln[1,])
  out.h <- .2 * (length(rownames(aln)) + 3 + length(tracks))

  pdf(file=filename, width=out.w, height=out.h)
  p <- aln.plot(aln, 
    tree=tree,
    color.scheme='auto',
    plot.chars=T, 
    aln.plot.chars=T, 
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

example.alns <- function() {
  # Chimpanzee crap:
#  get.aln(1592585, 700, 750, "aln_ex_1.pdf")
  # Gorilla crap:
  get.aln(1863113, 1, 100, "aln_ex_2.pdf")
  # Human crap:
  get.aln(1591434, 210, 260, "aln_ex_3.pdf")
  # Platypus ain't so bad:
  get.aln(1898971, 560, 610, "aln_ex_4.pdf")  
  # Mouse/rat crap:
  get.aln(1830947, 370, 450, "aln_ex_5.pdf")
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

    # Do some ddply touch-ups... fix the z-scores
    genes <- ddply(genes, .(parameter_set_id), function(x) {
      x
    })

    save(genes, file=cache.f)
  }
  #print("Loading cache...")
  load(cache.f)

  # Do a final few touch-ups or filters on the data frame.
  genes$pset <- genes$parameter_set_id

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
  gene.correlations(3, 2)
  gene.correlations(2, 7)
  gene.correlations(6, 10)
}

gene.correlations <- function(pset_a, pset_b, filter='default', plot=T) {
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

    pdf(file=scratch.f(sprintf("genes_cor_%d_%d_%s.svg", pset_a, pset_b, filter)), width=4, height=4)
    print.ggplot(p)
    dev.off()

  }

  zero.pca[1]
}

get.go.annotations <- function(exclude.iea=F) {
  ens.f <- scratch.f("ens_go_63.Rdata")
  if (!file.exists(ens.f)) {
    load.go.csv("~/src/greg-ensembl/projects/2xmammals/ens_go_63.csv", ens.f)
  }
  load(ens.f)
  if (exclude.iea) {
    go.ens.excl.iea
  } else {
    go.ens
  }
}

bsub.genes.enrichment <- function(...) {
  # dnds_05
  # dnds_10
  # tpm
  # f_pos
  # f_pos_bh
  # f_pos_neg 
  methods <- c('dnds_05', 'dnds_10', 'tpm', 'f_pos', 'f_pos_bh', 'pos_neg_2', 'pos_neg_3')
  filters <- c('default', 'stringent', 'pfam')
  for (filter in filters) {
    for (method in methods) {
      for (excl.iea in c(T, F)) {
        #xtra <- sprintf("%s %s %s", filter, method, as.character(excl.iea))
        #bsub.pset.function('genes_enrichment', queue='normal', mem=3, extra.args=xtra, ...)
        xtra <- sprintf("%d %s %s %s", 5, filter, method, as.character(excl.iea))
        bsub.function('genes_enrichment', queue='normal', mem=3, extra.args=xtra, ...)
        return()
      }
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
  go.ens <- get.go.annotations(exclude.iea=excl.iea)
  scores.df <- get.enrichment.df(pset_id=pset, filter=filter, method=method, other_pset=other_pset)
  scores.df <- subset(scores.df, select=c('score', 'binary.score', 'id', 'name', 'length', 't_test_a', 't_test_b'))

  print(head(scores.df))
  lbl <- sprintf("%s %s %s %s %s", pset, filter, method, excl.iea, other_pset)
  go.res <- do.enrichments(scores.df, lbl, go.ens)

  go.res$filter <- filter
  go.res$method <- method
  go.res$excl.iea <- excl.iea
  go.res$pset <- pset
  go.res$other_pset <- other_pset

  go.res$label <- paste(go.res$label, go.res$go_id, sep=' ')

  con <- connect(dbname())
  write.or.update(go.res, 'go', con, 'label')
  disconnect(con)
}

get.enrichment.df <- function(pset_id, filter='default', method='dnds_5', other_pset=NA) {
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
  # dnds_10
  # tpm
  # f_pos
  # f_pos_bh
  # pos_neg_2
  # pos_neg_3

  genes$t_test_a <- NA
  genes$t_test_b <- NA

  if (method == 'dnds_05') {
    genes$score <- genes$slr_dnds
    dnds.q <- quantile(genes$score, 0.95)
    genes$binary.score <- ifelse(genes$score > dnds.q, 1, 0)

  } else if (method == 'dnds_10') {
    genes$score <- genes$slr_dnds
    dnds.q <- quantile(genes$score, 0.90)
    genes$binary.score <- ifelse(genes$score > dnds.q, 1, 0)

  } else if (method == 'tpm') {
    genes$score <- genes$tpm_05
    genes$binary.score <- ifelse(genes$score < 0.05, 1, 0)

  } else if (method == 'f_pos') {
    genes$score <- genes$f_pos_20
    pos.q <- quantile(genes$score, 0.90)
    genes$binary.score <- ifelse(genes$score > pos.q, 1, 0)

  } else if (method == 'f_pos_bh') {
    genes$score <- genes$f_pos_gene_bh_10
    genes$binary.score <- ifelse(genes$score > 0, 1, 0)

  } else if (method == 'pos_neg_a') {
    pos_factor <- 2
    genes$score <- genes$f_pos_10 * pos_factor + genes$f_neg_10
    genes$score <- ifelse(genes$n_pos_05 > 3, genes$score, 0)
    score.q <- quantile(genes$score, 0.95)
    genes$binary.score <- ifelse(genes$score > score.q, 1, 0)

  } else if (method == 'pos_neg_b') {
    med.neg <- median(genes$f_neg_10)
    pos_factor <- 8
    genes$score <- genes$f_pos_gene_bh_10
#    genes$score <- ifelse(genes$f_neg_10 > 0.5, genes$score, 0)
#    genes$score <- genes$f_pos_10 * pos_factor + genes$f_neg_10
#    genes$score <- ifelse(genes$f_pos_10 > 0, genes$score, 0)
#    genes$score <- ifelse(genes$n_pos_gene_bh_10 > 1, genes$score, 0)
    score.q <- quantile(genes$score, 0.95)
    genes$binary.score <- ifelse(genes$score > score.q, 1, 0)

    pdf(file="test.pdf")
    genes$xx <- genes$f_neg_10
    genes$yy <- genes$f_pos_10
    plot(genes$xx, genes$yy, col=gray(0.7))
    abv <- genes$score > score.q
    points(genes[abv, 'xx'], genes[abv, 'yy'], col='red', alpha=0.4) 
    abv <- genes$slr_dnds_z > 1.3
    points(genes[abv, 'xx'], genes[abv, 'yy'], col='green', alpha=0.4) 
    abv <- genes$tpm_05 < 0.1
    points(genes[abv, 'xx'], genes[abv, 'yy'], col='blue', alpha=0.4) 
    bh.q <- quantile(genes$f_pos_gene_bh_05, 0.95)
    dev.off()
    
    pdf(file="test.pdf")
    plot(genes$gc_100k, genes$slr_dnds_z)
    dev.off()

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
  } else if (method == 'pos_parallel') {
    min.pos <- pmin(genes$f_pos_01.x, genes$f_pos_01.y)
    genes$score <- min.pos
    score.q <- quantile(genes$score, 0.99)
    genes$binary.score <- ifelse(genes$score > score.q, 1, 0)
  } else {
    stop("Enrichment method not recognized!")
  }

  genes <- subset(genes, gene_name != 'NA')
  genes$id <- genes$ref_protein_id
  genes$name <- genes$gene_name
  genes$length <- genes$aln_length

  genes <- genes[order(genes$score, decreasing=T),]

  xx <- subset(genes, select=c('score', 'binary.score', 'name', 'length', 't_test_a', 't_test_b'))
  print(head(xx, n=30))

  genes
}

plot.gene.fractions <- function() {
  genes <- get.genes.sites(filter='default')
  genes <- subset(genes, pset == 1)

  a <- 'f_pos_20'
  b <- 'f_ntr_20'
  c <- 'f_neg_20'

  genes$slr_bl <- pmin(genes$slr_bl, 20)

  genes$slr_log <- log(pmax(genes$slr_bl, 0.001))
  p <- ggplot(genes, aes(x=slr_log))
  p <- p + geom_histogram()
#  p <- p + geom_point(size=0.5)
  
  pdf(file="test.pdf")
#  p <- p + facet_grid(pset ~ .)
#  png(file="test.png", height=1400, width=400)
  print.ggplot(p)
  dev.off()

  return()

#  genes$log_p <- ifelse(genes$tpm_05 < 0.05, 1, 0)
   genes$is_sig <- ifelse(genes$n_pos_fwer > 1, 1, 0)

  p <- ternary.plot(genes, a, b, c, color.fld='slr_dnds', size.fld='is_sig')
  pdf(file="test.pdf", width=4, height=4)
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