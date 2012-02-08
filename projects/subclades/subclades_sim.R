lib <- ifelse(uname == 'gj1', 'src', 'lib')
source(sprintf("~/%s/greg-ensembl/projects/2xmammals/aln.R", lib))
source(sprintf("~/%s/greg-ensembl/projects/2xmammals/tree.R", lib))

work.f <- function(str) {
  lib <- ifelse(uname == 'gj1', 'src', 'lib')
  normalizePath(sprintf("~/%s/greg-ensembl/projects/subclades/%s", lib, str))
}

scratch.f <- function(str) {
  lib <- ifelse(uname == 'gj1', 'src', 'lib')
  normalizePath(sprintf("~/scratch/subclades/%s", dir))
}

bsub.sim.alns <- function() {
  species.groups <- c('Primates', 'Glires', 'Laurasiatheria', 'Atlantogenata')
  pop.sizes <- c(1e4, 5e4, 1e5, 5e5, 1e6)
  total.bls <- c(0.25, 0.5, 1, 2, 3)

  for (species in species.groups) {
    for (pop.size in pop.sizes) {
      for (total.bl in total.bls) {
        xtra <- paste(species, pop.size, total.bl, sep=' ')
        bsub.function('sim.aln', extra.args=xtra)
        return()
      }
    }
  }
}

sim.aln <- function(species, pop.size, total.bl) {
  print(species)
  print(as.numeric(pop.size))
  print(as.numeric(total.bl))

  tree.f <- work.f(paste(species, '.nh', sep=''))
  tree <- read.tree(tree.f)

  print(tree)
}