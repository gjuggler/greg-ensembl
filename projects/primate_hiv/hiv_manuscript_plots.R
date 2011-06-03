library(phylosim)
library(RColorBrewer)
library(xtable)
library(plyr)
source("~/src/greg-ensembl/projects/phylosim/PhyloSimPlots.R")
#source("~/src/greg-ensembl/projects/slrsim/slrsim.functions.R")
#source("~/src/greg-ensembl/projects/slrsim/slrsim.plots.R")

run.codeml.analysis <- function() {
  if (!exists('windows')) {
    load(dd.f('windows_all.Rdata'))
    windows <<- data
  }
  x <- windows
  x <- subset(x, peptide_window_width == 9999)

  final.nms <- final.gene.names()
  x <- subset(x, gene_name %in% final.nms)
  x <- x[!duplicated(x$gene_name),]
  print(nrow(x))

  for (i in 1:nrow(x)) {
    row <- x[i,]
    aln.f <- dd.f(paste('data/', row$aln_file, sep=''))
    print(aln.f)
    
    new.aln.f <- codeml.f(paste(row$gene_name, '.fasta', sep=''))
    file.copy(aln.f, new.aln.f)

    cmd <- sprintf("bsub -q normal \"perl %s/hiv_codeml_tests.pl --aln=%s \"",
      wd.f(''), tools::file_path_as_absolute(new.aln.f))
    print(cmd)
    system(cmd)
  }
}

plot.final.genes <- function() {
  if (!exists('windows')) {
    load(dd.f('windows_all.Rdata'))
    windows <<- data
  }
  x <- windows

  x <- subset(x, peptide_window_width == 9999)

  final.nms <- final.gene.names()
  x <- subset(x, gene_name %in% final.nms)
  print(nrow(x))

  for (i in 1:nrow(x)) {
    row <- x[i,]
    aln.f <- dd.f(paste('data/', row$aln_file, sep=''))
    print(aln.f)

    file.copy(aln.f, plot.f(paste(row$gene_name, '.fasta', sep='')))

    aln <- read.aln(aln.f)
    rownames(aln) <- ensp.to.species(rownames(aln))
    aln <- sort.aln.if(aln, sort.order())

    pep.aln <- aln.tx(aln)

    n.pages <- ceiling( ncol(pep.aln) / 125)
    print(n.pages)

    pdf.f <- plot.f(paste(row$gene_name, '_pep.pdf', sep=''))
    pdf(file=pdf.f, width=30, height=3*n.pages)
    sim <- PhyloSim(); sim$.alignment <- pep.aln
    plotAlignment(sim, axis.text.size=12, aln.plot.chars=T, aln.char.text.size=3, num.pages=n.pages)
    dev.off()

    pdf.f <- plot.f(paste(row$gene_name, '_nuc.pdf', sep=''))
    pdf(file=pdf.f, width=30, height=3*n.pages)
    aln[aln == 'N'] <- 'X'
    sim <- PhyloSim(); sim$.alignment <- aln
    plotAlignment(sim, axis.text.size=12, char.text.size=2, aln.plot.chars=F, num.pages=n.pages)
    dev.off()
  }
  
}

final.gene.names <- function() {
  a <- read.csv(wd.f('candidates_final.txt'), header=F, stringsAsFactors=F)
  a$V1
}

codeml.f <- function(f) {
  return(paste("~/src/greg-ensembl/projects/primate_hiv/manuscript_codeml/", f, sep=''))
}

# Plots directory.
plot.f <- function(f) {
  return(paste("~/src/greg-ensembl/projects/primate_hiv/manuscript_plots/", f, sep=''))
}

# Working directory.
wd.f <- function(f) {
  return(paste("~/src/greg-ensembl/projects/primate_hiv/", f, sep=''))
}

# Data directory.
dd.f <- function(f) {
  return(paste("~/scratch/gj1_hiv_full_57_c/2011-02-23_01/", f, sep=''))
}

sort.order <- function() {
  c('Human', 'Chimpanzee', 'Gorilla', 'Orangutan', 'Macaque', 'Marmoset', 'Tarsier',
    'MouseLemur', 'Bushbaby')
}

ensp.to.species <- function(ids) {
  ids <- handle.grep(ids, "ENSP0.*", 'Human')
  ids <- handle.grep(ids, "ENSPTRP0.*", 'Chimpanzee')
  ids <- handle.grep(ids, "ENSGGOP0.*", 'Gorilla')
  ids <- handle.grep(ids, "ENSPPYP0.*", 'Orangutan')
  ids <- handle.grep(ids, "ENSMMUP0.*", 'Macaque')
  ids <- handle.grep(ids, "ENSCJAP0.*", 'Marmoset')
  ids <- handle.grep(ids, "ENSTSYP0.*", 'Tarsier')
  ids <- handle.grep(ids, "ENSMICP0.*", 'MouseLemur')
  ids <- handle.grep(ids, "ENSOGAP0.*", 'Bushbaby')
  ids
}

handle.grep <- function(ids, regex, species) {
  matches <- grepl(regex, ids)

  if (sum(matches) > 1) {
    suffixes <- 1:sum(matches)
    new.ids <- paste(species, '_', suffixes, sep='')
  } else {
    new.ids <- species
  }
  
  ids[matches] <- new.ids
  ids
}