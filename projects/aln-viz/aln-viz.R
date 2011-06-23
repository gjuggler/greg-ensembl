library(ape)
library(phylosim)
library(plyr)


# Simulate an alignment with a given indel rate, root length, and tree root-to-tip length.
sim.aln <- function(indel.rate = 0.03, root.length = 100, tree.length = 1, file.base = paste(indel.rate, tree.length, 
    sep = "_")) {
    PSIM_FAST = TRUE
    assign("PSIM_FAST", TRUE, envir = .GlobalEnv)
    
    # Create the root sequence.
    root.sequence <- Sequence(length = root.length)
    
    # Set up indel distribution:
    max.indel <- 20
    id.dist <- exp(max.indel:1)/sum(exp(max.indel:1))
    
    
    # Coding options.
    #coding.pos <-c(1:root.length)
    #codon.alph <- CodonAlphabet()
    #setAlphabets(root.sequence, list(codon.alph), coding.pos)
    #gy94 <-GY94(kappa=2, omega.default=0.2, scale.nuc=TRUE)
    #del.c <- DiscreteDeletor(rate=indel.rate, sizes=1:6, probs=id.dist )
    #ins.c <- DiscreteInsertor(rate=indel.rate, sizes=1:6, probs=id.dist )
    #ins.c$templateSeq <- CodonSequence(length=1,processes=list(list( gy94, del.c, ins.c ) ))
    #setProcesses(root.sequence, list(list(gy94, del.c, ins.c)), coding.pos)
    #setStates(root.sequence, 'ATG', coding.pos[1]); # Set the state.
    #setRateMultipliers(root.sequence, gy94, 0, coding.pos[1]) # Make the site inv
    #setDeletionTolerance(root.sequence, del.c, 0, coding.pos[1]); # Make the site reject deletions.
    #setInsertionTolerance(root.sequence,ins.c ,0, coding.pos[1]);
    
    # Noncoding options.
    noncoding.pos <- c(1:root.length)
    nuc.alph <- NucleotideAlphabet()
    setAlphabets(root.sequence, list(nuc.alph), noncoding.pos)
    k80 <- K80(rate.params = list(Alpha = 2, Beta = 1), base.freqs = c(2, 1, 2, 1)/4)
    del.nc <- DiscreteDeletor(rate = indel.rate, sizes = 1:max.indel, probs = id.dist)
    ins.nc <- DiscreteInsertor(rate = indel.rate, sizes = 1:max.indel, probs = id.dist)
    ins.nc$templateSeq <- NucleotideSequence(length = 1, processes = list(list(k80, del.nc, ins.nc)))
    setProcesses(root.sequence, list(list(k80, del.nc, ins.nc)), noncoding.pos)
    
    sampleStates(root.sequence)
    print(root.sequence)
    
    # data(bird.families)
    # tree <- bird.families
    # tree <- rcoal(4)
    #tree <- read.tree('~/lib/greg-ensembl/projects/slrsim/trees/2xmammals.nh')
    #tree <- read.tree('~/lib/phylosim/examples/data/ensembl_pax.nh')
    set.seed(1)
    tree <- rtree(80)
    tree <- scale.length.to(tree, tree.length)
    
    aln.file <- paste(file.base, ".fasta", sep = "")
    tree.file <- paste(file.base, ".nh", sep = "")
    write.tree(tree, file = tree.file)
    
    sim <- PhyloSim(phylo = tree, root.seq = root.sequence)
    
    Simulate(sim)
    
    saveAlignment(sim, skip.internal = TRUE, file = aln.file)
    write.tree(sim$phylo, file = tree.file)
}

do.big.sims <- function() {
    counter <- 0
    for (indel in seq(from = 0, to = 0.1, length.out = 3)) {
        for (length in seq(from = 0.1, to = 2, length.out = 3)) {
            r.cmd <- sprintf("\nsource('./aln-viz.R');\nsim.aln(indel.rate=%.3f, root.length=300, tree.length=%.3f);", indel, 
                length)
            dir <- "./tmp/"
            if (!file.exists(dir)) {
                dir.create(dir, recursive = T)
            }
            r.tmp <- paste(dir, counter, ".R", sep = "")
            job.out <- paste(dir, counter, ".out", sep = "")
            file.create(r.tmp)
            
            writeChar(r.cmd, r.tmp)
            
            cmd <- sprintf("bsub -M8000 -R 'rusage[mem=8000]' -o %s 'R-2.12.0 --slave < %s'", job.out, r.tmp)
            #print(cmd)
            system(cmd)
            #return()
            counter <- counter + 1
        }
    }
}

plot.subset <- function() {
    set.seed(1)
    
    fasta.files <- dir("./")
    fasta.files <- fasta.files[grep("\\.(fa|fasta)$", fasta.files)]
    #fasta.files <- fasta.files[4]
    
    fasta.files <- fasta.files[fasta.files == "0.05_1.05.fasta"]
    
    for (file in fasta.files) {
        file.base <- gsub("\\.[^\\.]+$", "", file)
        tree.f <- paste(file.base, ".nh", sep = "")
        aln.f <- file
        for (subset.size in c(64, 2, 8, 32)) {
            for (n.pages in c(1, "auto")) {
                for (tolerance in c(1.1, "full")) {
                  pdf.out <- paste(file.base, "_", subset.size, "_", tolerance, "_", n.pages, ".pdf", sep = "")
                  
                  if (file.exists(pdf.out)) {
                    next()
                  }
                  
                  sim <- PhyloSim()
                  readTree(sim, tree.f)
                  readAlignment(sim, aln.f)
                  
                  names <- sim$phylo$tip.label
                  set.seed(1)
                  sub.names <- sample(names, size = subset.size, replace = F)
                  sim <- tree.aln.subset(sim, sub.names, remove.blank.columns = T)
                  
                  w <- 20
                  h <- 2
                  if (n.pages == "auto") {
                    w <- 7
                    h <- 7
                  }
                  
                  if (tolerance == 1.1) {
                    sim <- remove.gaps(sim, tolerance)
                  }
                  
                  pdf(width = w, height = h, file = pdf.out)
                  print(paste("  ", "plotting", file.base, subset.size, tolerance, n.pages))
                  plot(sim, plot.chars = F, num.pages = n.pages, tree.xlim = c(0, 2))
                  dev.off()
                }
            }
        }
    }
    
    for (tolerance in c(0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2)) {
        print(paste("  tolerance=", tolerance))
        sim <- PhyloSim()
        readTree(sim, "0.1_2.nh")
        readAlignment(sim, "0.1_2.fasta")
        del.sim <- remove.gaps(sim, tolerance)
        
        xlim.file <- paste("tolerance", "_", tolerance, "_", "xlim", ".pdf", sep = "")
        if (!file.exists(xlim.file)) {
            pdf(width = 20, height = 2, file = xlim.file)
            plot(del.sim, plot.chars = F, num.pages = 1, tree.xlim = c(0, 2), aln.xlim = c(1, 1500))
            dev.off()
            
            pdf(width = 20, height = 2, file = paste("tolerance", "_", tolerance, "_", "noxlim", ".pdf", sep = ""))
            plot(del.sim, plot.chars = F, num.pages = 1, tree.xlim = c(0, 2))
            dev.off()
            
            pdf(file = paste("tolerance", "_", tolerance, "_", "full", ".pdf", sep = ""))
            plot(del.sim, plot.chars = F, num.pages = "auto", tree.xlim = c(0, 2))
            dev.off()
        }
        else {
            print(paste("File", xlim.file, "already exists. Skipping this tolerance..."))
        }
        
    }
}

plot.pfam <- function() {
    
    for (tolerance in c(0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5)) {
        file.out <- paste("pfam", "_", tolerance, ".pdf", sep = "")
        if (file.exists(file.out)) {
            next()
        }
        print(paste("  Plotting", file.out))
        
        sim <- PhyloSim()
        readTree(sim, "PF00031_full.nhx")
        readAlignment(sim, "PF00031_full.fasta")
        
        names <- sim$.phylo$tip.label
        set.seed(1)
        sub.names <- sample(names, size = 200, replace = F)
        print("  subset")
        sim <- tree.aln.subset(sim, sub.names, remove.blank.columns = T)
        print("  remove gaps")
        sim <- remove.gaps(sim, tolerance)
        
        print("  plot")
        pdf(width = 20, height = 10, file = file.out)
        plot(sim, plot.chars = F, num.pages = 1)
        dev.off()
        
    }
    
}

# Return an alignment and tree with only the sequences listed in keep.seq
# retained.
tree.aln.subset <- function(sim, keep.seqs, remove.blank.columns = T) {
    tree <- sim$.phylo
    aln <- sim$.alignment
    
    # Remove non-keeper leaves from the tree.
    names <- tree$tip.label
    drop.seqs <- names[!(names %in% keep.seqs)]
    tree <- drop.tip(tree, drop.seqs)
    
    # Remove non-keeper seqs from the alignment.
    aln <- aln[rownames(aln) %in% keep.seqs, ]
    
    # Remove blank columns from the alignment if desired
    if (remove.blank.columns) {
        # Get a vector of column-wise alignment strings.
        column.strings <- apply(aln, 2, function(column) {
            paste(column, collapse = "")
        })
        # List of alignment columns containing only gaps.
        gap.columns <- grep("^-+$", column.strings)
        
        # Remove those columns from the alignment.
        #print(paste(' ',length(gap.columns),' blank columns'))
        if (length(gap.columns) > 0) {
            aln <- aln[, -c(gap.columns)]
        }
    }
    
    sim2 <- PhyloSim()
    sim2$.phylo <- tree
    sim2$.alignment <- aln
    return(sim2)
}


### Some miscellaneous tree utility functions.
branch.length <- function(phylo, node) {
    edge.index <- which(phylo$edge[, 2] == node)
    bl <- phylo$edge.length[edge.index]
    if (length(bl) == 0) {
        bl <- 0
    }
    return(bl)
}

scale.tree.to <- function(tree, new.total.length) {
    cur.total.length <- sum(tree$edge.length)
    tree$edge.length <- (new.total.length/cur.total.length) * tree$edge.length
    return(tree)
}

scale.length.to <- function(tree, new.max.length) {
    cur.length.to.root <- max.length.to.root(tree)
    tree$edge.length <- (new.max.length/cur.length.to.root) * tree$edge.length
    return(tree)
}

node.with.label <- function(tree, label) {
    return(which(tree$tip.label %in% label))
}

# Extracts the length of the branch above the given node. Returns 0 if the node is root.
branch.length <- function(phylo, node) {
    edge.index <- which(phylo$edge[, 2] == node)
    bl <- phylo$edge.length[edge.index]
    if (length(bl) == 0) {
        bl <- 0
    }
    return(bl)
}

# The maximum root-to-tip length in the tree.
max.length.to.root <- function(phylo) {
    max.length <- 0
    for (i in 1:length(phylo$tip.label)) {
        cur.length <- length.to.root(phylo, i)
        max.length <- max(max.length, cur.length)
    }
    return(max.length)
}

# The length from the root to the given node. Can be given either as a node ID or a tip label.
length.to.root <- function(phylo, node) {
    tip.index <- node
    if (is.character(node)) {
        tip.index <- which(phylo$tip.label == node)
    }
    
    cur.node.b <- tip.index
    
    p.edges <- phylo$edge
    p.lengths <- phylo$edge.length
    
    length <- 0
    while (length(which(p.edges[, 2] == cur.node.b)) > 0) {
        cur.edge.index <- which(p.edges[, 2] == cur.node.b)
        cur.edge.length <- p.lengths[cur.edge.index]
        length <- length + cur.edge.length
        cur.node.a <- p.edges[cur.edge.index, 1]
        cur.node.b <- cur.node.a
    }
    return(length)
}

 
