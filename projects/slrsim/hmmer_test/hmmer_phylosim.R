library(R.oo)
library(ape)
library(ggplot2)
source("~/src/greg-ensembl/projects/phylosim/PhyloSimSource.R")
source("~/src/greg-ensembl/projects/phylosim/PhyloSimPlots.R")

read.hmmer <- function(file) {
  lines <- readLines(file)
  lines.toks <- strsplit(lines, "\\s+")

  seen.hmm <- F
  hmm.alphabet <- NULL
  hmm.states <- NULL
  hmm.start.ind <- -1
  for (i in 1:length(lines)) {
    cur.line <- lines.toks[[i]]
    #print(cur.line[1:3])

    if (cur.line[1] == 'HMM') {
      seen.hmm <- T
      hmm.alphabet <- cur.line[2:length(cur.line)]
    }
    
    if (seen.hmm && cur.line[2] == 'm->m') {
      hmm.states <- cur.line[2:length(cur.line)]
    }

    if (seen.hmm && cur.line[2] == 1) {
      hmm.start.ind <- i
      break;
    }
  }

  hmm.triplets <- lines.toks[hmm.start.ind:length(lines.toks)]

  site.list <- list()
  for (i in seq(from=hmm.start.ind, to=length(lines.toks)-2, by=3)) {
    cur.lines <- lines.toks[i:(i+2)]

    match.line <- cur.lines[[1]]
    ins.line <- cur.lines[[2]]
    trans.line <- cur.lines[[3]]

    cur.site <- match.line[2]
    match.probs <- match.line[3:(length(match.line)-3)]
    ins.probs <- ins.line[2:length(ins.line)]
    trans.probs <- trans.line[2:length(trans.line)]

    match.probs <- conv.probs(match.probs)
    ins.probs <- conv.probs(ins.probs)
    trans.probs <- conv.probs(trans.probs)

    names(match.probs) <- hmm.alphabet
    names(ins.probs) <- hmm.alphabet
    names(trans.probs) <- hmm.states

    cur.list <- list(site=cur.site, match.probs=match.probs, ins.probs=ins.probs, trans.probs=trans.probs)
    site.list[[cur.site]] <- cur.list
  }

  return(site.list)
}

# Converts probabilities from HMM format (-log(x) or "*" where x == 0) to numeric
conv.probs <- function(x) {
  x[x == '*'] <- Inf
  x <- as.numeric(x)
  x <- exp(-x)
  return(x)
}

profile.to.phylosim.seq <- function(profile) {
  n.sites <- length(profile)
  #n.sites <- 10  # Shorten the profile for quicker testing...

  seq <- CodonSequence(length=n.sites)
  gy <- GY94(kappa = 2)
  #ins <- ContinuousInsertor(rate = 1, max.length = 18, dist = expression(which(pzipf(1:18, 18, 1.8) >= runif(1, 0, 1))[1]))
  #del <- ContinuousDeletor(rate = 1, max.length = 18, dist = expression(which(pzipf(1:18, 18, 1.8) >= runif(1, 0, 1))[1]))
  #ins <- DiscreteInsertor(rate = 1, sizes=1:18, probs=dzipf(1:18, 18, 1.8))
  #del <- DiscreteDeletor(rate = 1, sizes=1:18, probs=dzipf(1:18, 18, 1.8))

  # Load up a list of codons and corresponding aa's.
  ca <- CodonAlphabet()
  codons <- as.character(ca)
  aas <- unlist(lapply(as.list(codons), function(x) { translateCodon(ca, x) }))

  ins.rates <- rep(0, n.sites)
  del.rates <- rep(0, n.sites)
  omegas <- rep(1, n.sites)

  sites.list <- list()
  process.list <- list()
  for (i in 1:n.sites) {
    print(paste(i, '/',n.sites))
    cur.site <- profile[[i]]
    
    # Create a separate GY94 process
    # Is it faster to clone or to create a new process?
    #cur.gy <- GY94(kappa = 2)
    cur.gy <- clone(gy)

    # Create sitewise insertors and deletors.
    max.indel.length <- 20
    ins.rate <- cur.site$trans.probs[2]
    del.rate <- cur.site$trans.probs[3]
    ins.to.ins.prob <- cur.site$trans.probs[5]
    ins.probs <- c(1, ins.to.ins.prob ^ (1:(max.indel.length-1)))
    del.to.del.prob <- cur.site$trans.probs[7]
    del.probs <- c(1, del.to.del.prob ^ (1:(max.indel.length-1)))
    
    ins.probs <- ins.probs / sum(ins.probs)
    del.probs <- del.probs / sum(del.probs)

    ins <- BrownianInsertor(
      type="discrete", 
      scale=0.05, 
      sizes=1:max.indel.length, 
      probs=ins.probs,
      rate = ins.rate
    )
    del <- DiscreteDeletor(
      rate = del.rate,
      sizes=1:max.indel.length,
      probs=del.probs
    )

    # Set the codon frequencies to approximate the HMM amino acid probs at this site.
    freqs <- aa.probs.to.codon.freqs(cur.site$match.probs, codons, aas)
    freqs <- freqs / sum(freqs)
    setCodonFreqs(cur.gy, freqs)

    # Set omega to a value scaled by the entropy at the given site (and normalized by the max possible entropy)
    # e.g., a site with completely uniform probs will have omega=1
    site.entropy <- site.entropy(cur.site$match.probs)
    omegas[i] <- site.entropy / log2(20)

    attachProcess(seq, ins, i)
    attachProcess(seq, del, i)
    attachProcess(seq, cur.gy, i)
    setParameterAtSites(seq, cur.gy, id="omega", value=omegas[i], index=i)
  }  

  #ins.rates <- ins.rates * 0.05 / mean(ins.rates) + 0.001
  #del.rates <- del.rates * 0.05 / mean(del.rates) + 0.001
  #ins.rates[ins.rates >= 1] <- 1
  #del.rates[del.rates >= 1] <- 1
  #ins.rates <- ins.rates / max(ins.rates)
  #del.rates <- del.rates / max(del.rates)
  #attachProcess(seq, ins, 1:n.sites)
  #attachProcess(seq, del, 1:n.sites)
  #setInsertionTolerance(seq, ins, ins.rates)
  #setDeletionTolerance(seq, del, del.rates)
  #print(getInsertionTolerance(seq, ins))
  #print(getDeletionTolerance(seq, del))

  sampleStates(seq)
  return(seq)
}

site.entropy <- function(probs) {
  return(-sum(probs * log2(probs)))
}

aa.probs.to.codon.freqs <- function(aa.probs, codons, aas) {
  # Take amino acid probabilities and evenly spread them among consituent codons.
  aa.names <- names(aa.probs)
  codon.freqs <- rep(0, length(codons))
  for (i in 1:length(aa.probs)) {
    cur.aa <- aa.names[i]
    cur.codon.indices <- which(aas == cur.aa)
    codon.freqs[cur.codon.indices] <- aa.probs[i] / length(cur.codon.indices)
  }
  return(codon.freqs)
}

load.or.make.tree <- function(tree.file=NA, phylo=NA, n.leaves=NA) {
  if (!is.na(tree.file)) {
    tree <- read.tree(tree.file)
  } else if (!is.na(phylo)) {
    tree <- phylo
  } else {
    tree <- rcoal(n.leaves, br=runif)
    tree.len <- max.length.to.root(tree)
    tree <- scale.tree.by(tree, 2 / tree.len )  
  }
  print(write.tree(tree))

  tree
}

sim.pfam <- function(pfam.id, ...) {
  tree <- load.or.make.tree(...)

  if (!file.exists(pfam.id)) {
    print("  downloading PFam HMM...")
    system(paste("wget ","http://pfam.sanger.ac.uk/family/hmm/", pfam.id, sep=''))
  }

  orig.aln.f <- paste(pfam.id, '_seed.fasta', sep='')
  if (!file.exists(orig.aln.f)) {
    print("  downloading PFam seed aln...")
    cmd <- sprintf('wget "http://pfam.sanger.ac.uk/family/alignment/download/format?format=fasta&gaps=dashes&acc=%s" -O "%s"', pfam.id, orig.aln.f)
    system(cmd)
  }
  print("  plotting seed aln...")
  seed.sim <- PhyloSim()
  seed.aln <- read.aln(orig.aln.f)
  seed.sim$.alignment <- seed.aln
  pdf(file=paste(pfam.id, "_seed.pdf", sep=''), width=20, height=10)
  plotAlignment(seed.sim, aln.plot.chars=T, plot.ancestors=F)
  dev.off()

  print("  reading HMMER file...")
  profile <- read.hmmer(pfam.id)

  do.sim(pfam.id, profile, tree)
}

sim.hmm <- function(id, aln.file=NA, ...) {
  tree <- load.or.make.tree(...)

  print("  converting aln...")
  stockholm.f <- paste(id, '.aln', sep='')
  cmd <- sprintf("perl fasta2stockholm.pl %s > %s", aln.file, stockholm.f)
  system(cmd)

  hmm.f <- paste(id, '.hmm', sep='')
  if (!file.exists(hmm.f)) {
    print("  calculating HMM...")
    cmd <- sprintf("hmmbuild %s %s", hmm.f, stockholm.f)
    system(cmd)
  }
  print("  reading HMMER file...")
  profile <- read.hmmer(hmm.f)

  print("  plotting orig aln...")
  orig.sim <- PhyloSim()
  orig.aln <- read.aln(aln.file)
  orig.aln <- remove.blank.columns(orig.aln)
  orig.sim$.alignment <- orig.aln
  pdf(file=paste(id, "_orig.pdf", sep=''), width=20, height=10)
  plotAlignment(orig.sim, aln.plot.chars=T, plot.ancestors=F)
  dev.off()

  do.sim(id, profile, tree)  
}

do.sim <- function(id, profile, tree) {
  PSIM_FAST <<- TRUE

  seq.f <- paste(id, '_phylosim_seq.Rdata', sep='')
  if (!file.exists(seq.f)) {
    print("  converting profile to phylosim seq...")
    root.seq <- profile.to.phylosim.seq(profile)
    base::save(root.seq, file=seq.f)
  }
  print("  loading phylosim seq from file...")
  base::load(seq.f)

  sim.f <- paste(id, '_phylosim_sim.Rdata', sep='')
  if (!file.exists(sim.f)) {
    print("  simulating...")
    sim <- PhyloSim()
    setPhylo(sim, tree)

    sampleStates(root.seq)
    print(sprintf("root DNA: %.20s ...", root.seq))
    print(sprintf("root AA : %.20s ...", Translate(root.seq, as.character(root.seq))))

    setRootSeq(sim, root.seq)
    Simulate(sim)
    base::save(sim, file=sim.f)
  }
  print("  loading sim from file...")
  base::load(sim.f)


  rm(PSIM_FAST, envir=.GlobalEnv)

  print("  outputting alns & plots")

  # Save DNA alignment
  aln.f <- paste(id, '_dna.fasta', sep='')
  saveAlignment(sim, file=aln.f, skip.internal=T)

  # Read DNA alignment, tx to pep.
  dna.aln <- read.aln(aln.f)
  pep.aln <- aln.tx(dna.aln)

  # Create a Pep sim to hold the peptide alignment.
  pep.sim <- PhyloSim()
  pep.sim$.alignment <- pep.aln
  pep.aln.f <- paste(id, '_pep.fasta', sep='')
  saveAlignment(pep.sim, file=pep.aln.f, skip.internal=T)

  pdf(file=paste(id, "_dna.pdf", sep=''), width=20, height=10)
  plotAlignment(sim, aln.plot.chars=T, plot.ancestors=F)
  dev.off()

  pdf(file=paste(id, "_pep.pdf", sep=''), width=20, height=10)
  plotAlignment(pep.sim, aln.plot.chars=T, plot.ancestors=F)
  dev.off()
}

gharmonic = function(n, s=1, lognexponent=0) {

    if (!is.Numeric(n, integ = TRUE, posit = TRUE))
        stop("bad input for argument 'n'")
    if (!is.Numeric(lognexponent, allow = 1))
        stop("bad input for argument 'lognexponent'")
    if (length(n) == 1 && length(s) == 1) {
        if (lognexponent != 0) sum(log(1:n)^lognexponent * (1:n)^(-s)) else
            sum((1:n)^(-s))
    } else {
        LEN = max(length(n), length(s))
        n = rep(n, len = LEN)
        ans = s = rep(s, len = LEN)
        if (lognexponent != 0) {
            for(ii in 1:LEN)
                ans[ii] = sum(log(1:n[ii])^lognexponent * (1:n[ii])^(-s[ii]))
        } else
            for(ii in 1:LEN)
                ans[ii] = sum((1:n[ii])^(-s[ii]))
        ans
    }
}

dzipf = function(x, N, s, log = FALSE)
{
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)

    if (!is.Numeric(x))
        stop("bad input for argument 'x'")
    if (!is.Numeric(N, integ = TRUE, posit = TRUE))
        stop("bad input for argument 'N'")
    if (!is.Numeric(s, posit = TRUE))
        stop("bad input for argument 's'")
    nn = max(length(x), length(N), length(s))
    x = rep(x, len = nn); N = rep(N, len = nn); s = rep(s, len = nn);
    ox = !is.finite(x)
    zero = ox | round(x) != x | x < 1 | x > N
    ans = (if (log.arg) log(0) else 0) * x
    if (any(!zero))
        if (log.arg) {
            ans[!zero] = (-s[!zero]) * log(x[!zero]) -
                         log(gharmonic(N[!zero], s[!zero]))
        } else {
            ans[!zero] = x[!zero]^(-s[!zero]) / gharmonic(N[!zero], s[!zero])
        }
    ans
}

is.Numeric <- function(x, allowable.length=Inf, integer.valued=FALSE, positive=FALSE) {
    if (all(is.numeric(x)) && all(is.finite(x)) &&
    (if(is.finite(allowable.length)) length(x)==allowable.length else TRUE) &&
    (if(integer.valued) all(x==round(x)) else TRUE) &&
    (if(positive) all(x>0) else TRUE)) TRUE else FALSE
}