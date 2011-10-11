source("~/src/greg-ensembl/projects/2xmammals/analyze_mammals.R")

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

  if (remove.blank.columns) {
    aln <- remove.blank.columns(aln)
  }

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