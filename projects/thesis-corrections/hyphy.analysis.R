library(ape)
library(devtools)
load_all('/homes/greg/lib/greg-ensembl/projects/ggphylo')
source('/homes/greg/lib/greg-ensembl/projects/2xmammals/aln.R')

cor.dir <- function() {
  '/homes/greg/lib/greg-ensembl/projects/thesis-corrections/'
}

cor.f <- function(x) {
  paste(cor.dir(), x, sep='')
}

prepare.gorilla.jobs <- function() {
  setwd(cor.dir())

  # Assume gorilla alignments were downloaded to gorilla_alns/00/etc...
  # Assume CSV file with top accelerated gorilla genes is at gorilla_top_genes.csv
  # Assume y.Rdata contains a full data frame with info about the genes.

  # For each gene to be analyzed, create a directory under gorilla_hyphy/00/etc...
  # containing the alignment file, a tree file, and an args.txt

  top.genes <- read.csv(file=cor.f('gorilla_top_genes.csv'), header=F, stringsAsFactors=F)
  load(cor.f('y.Rdata'))
  all.genes <- y
  rm(y)

  colnames(top.genes) <- c(
    'test',
    'gene',
    'length',
    'dnds_mamm',
    'dnds_prim',
    'nsyn',
    'syn',
    'branch_lrt',
    'branch_p',
    'branch_fdr',
    'bs_lrt',
    'bs_p',
    'bs_fdr'
  )

  print(str(top.genes))

  all.top <- unique(top.genes$gene)

  for (cur.gene in all.top) {
    cur.row <- subset(all.genes, gene_name == cur.gene)
    dp <- cur.row$data_prefix
    gn <- cur.row$gene_name
    cur.f <- cor.f(sprintf('gorilla_alns/%s/%s_aln_masked_subs.fasta', dp, gn))
    #print(cur.f)
    if (file.exists(cur.f)) {
      print(paste("Found ", gn))
      
      output.dir <- cor.f(sprintf('gorilla_hyphy/%s/%s/', dp, gn))
      if (!file.exists(output.dir)) {
        dir.create(output.dir, recursive=T)
      }

      # Copy the alignment file over.
      aln.f <- sprintf('%s%s.fasta', output.dir, gn)
      tree.f <- sprintf('%s%s.nh', output.dir, gn)
      args.f <- sprintf('%sargs.txt', output.dir)
      cmd.f <- sprintf('%scmd.sh', output.dir)

      file.copy(cur.f, aln.f, overwrite=T)
      aln <- aln.read(aln.f)
      ens.subs <- list(
        c('ENSP0.*', 'Human'),
        c('ENSPTRP0.*', 'Chimpanzee'),
        c('ENSGGOP.*', 'Gorilla'),
        c('ENSPPYP0.*', 'Orangutan'),
        c('ENSMMUP0.*', 'Macaque'),
        c('ENSCJAP0.*', 'Marmoset'),
        c('ENSTSYP0.*', 'Tarsier'),
        c('ENSMICP0.*', 'MouseLemur'),
        c('ENSOGAP0.*', 'Bushbaby')
      )
      aln <- aln.tx.labels(aln, ens.subs)
      aln.write(aln, aln.f)

      tree <- tree.read('(((((Human,Chimpanzee),Gorilla),Orangutan),Macaque),Marmoset);')
      #print(as.character(tree))
      cat(as.character(tree), file=tree.f)
      
      cat(sprintf("0\n%s\n%s", aln.f, tree.f), file=args.f)

      #cmd.txt <- "bsub -q research-rh6 -n 5 -a openmpi mpirun.lsf -np 5 -mca btl tcp,self /homes/greg/src/hyphy_build/bin/HYPHYMPI BASEPATH=/homes/greg/lib/greg-ensembl/projects/thesis-corrections/hyphy_build/lib/hyphy GABranchFiles/ModelSelectorBranchLocal.bf < args.txt"
      cmd.txt <- sprintf('#!/bin/bash\nbsub -o "/homes/greg/lib/scratch/lsf_logs/gorilla_hyphy_%%J_%%I.txt" -q research-rh6 -n 5 -a openmpi "mpirun.lsf -np 5 -mca btl tcp,self /homes/greg/src/hyphy_build/bin/HYPHYMPI BASEPATH=/homes/greg/lib/greg-ensembl/projects/thesis-corrections/hyphy_build/lib/hyphy GABranchFiles/ModelSelectorBranchLocal.bf < %sargs.txt"', output.dir)
      cat(cmd.txt, file=cmd.f)

      print(cmd.txt)
      sub.jobs <- T
      if (sub.jobs) {
        setwd(output.dir)
        system("module load openmpi-x86_64")
        system("chmod a+x cmd.sh")
        system("./cmd.sh")
        setwd(cor.dir())
      }
    }
    return()
  }
}