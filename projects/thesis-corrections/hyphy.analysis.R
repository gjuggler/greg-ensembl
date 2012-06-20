library(ape)
library(devtools)
library(yaml)
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

  #print(str(top.genes))

  all.top <- unique(top.genes$gene)
  print(length(all.top))

  #all.top <- c('SUPT16H', 'CYB5R3')

  for (cur.gene in all.top) {
    cur.row <- subset(all.genes, gene_name == cur.gene)
    dp <- cur.row$data_prefix
    gn <- cur.row$gene_name
    cur.f <- cor.f(sprintf('gorilla_alns/%s/%s_aln_masked_subs.fasta', dp, gn))

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
      args2.f <- sprintf('%sargs2.txt', output.dir)
      cmd.f <- sprintf('%scmd_calc.sh', output.dir)
      cmd2.f <- sprintf('%scmd_collect.sh', output.dir)
      row.f <- sprintf('%s%s.txt', output.dir, gn)

      cat(as.yaml(cur.row), file=row.f)

      # Translate aln IDs and output.
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

      # Output tree.
      tree <- tree.read('(((((Human,Chimpanzee),Gorilla),Orangutan),Macaque),Marmoset);')
      #print(as.character(tree))
      cat(as.character(tree), file=tree.f)
      
      # Output args.
      cat(sprintf("0\n%s\n%s", aln.f, tree.f), file=args.f)

      # Output cmd script.
      n.nodes <- 17
      cmd.txt <- sprintf('#!/bin/bash\nbsub -o "/homes/greg/scratch/lsf_logs/gorilla_hyphy_%%J_%%I.txt" -q research-rh6 -n %d -a openmpi "%s;%s"',
        n.nodes,
        sprintf('mpirun.lsf -np %d -mca btl tcp,self /homes/greg/src/hyphy_build/bin/HYPHYMPI BASEPATH=/homes/greg/lib/greg-ensembl/projects/thesis-corrections/hyphy_build/lib/hyphy GABranchFiles/ModelSelectorBranchLocal.bf < %sargs.txt', n.nodes, output.dir),
        sprintf('/homes/greg/src/hyphy_build/bin/HYPHYMP BASEPATH=/homes/greg/lib/greg-ensembl/projects/thesis-corrections/hyphy_build/lib/hyphy GABranchFiles/BranchGAResultProcessor.bf < %sargs2.txt', output.dir)
      )
      cat(cmd.txt, file=cmd.f)
      system(sprintf("chmod a+rx %s", cmd.f))

      cmd2.txt <- sprintf('#!/bin/bash\n%s',
        sprintf('/homes/greg/src/hyphy_build/bin/HYPHYMP BASEPATH=/homes/greg/lib/greg-ensembl/projects/thesis-corrections/hyphy_build/lib/hyphy GABranchFiles/BranchGAResultProcessor.bf < %sargs2.txt', output.dir)
      )
      cat(cmd2.txt, file=cmd2.f)
      system(sprintf("chmod a+x %s", cmd2.f))
      
      args2.txt <- sprintf("%s_ga_branch.out\n%s_results", aln.f, aln.f)
      cat(args2.txt, file=args2.f)

      #print(cmd.txt)
      print(output.dir)

      sub.jobs <- T
      if (sub.jobs) {
        setwd(output.dir)
        system("module load openmpi-x86_64")
        system(sprintf("%s", cmd.f))
        setwd(cor.dir())
      }
    }
#    return()
  }
}