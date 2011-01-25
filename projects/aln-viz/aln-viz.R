library(ape)
library(phylosim)
library(plyr)

score.aln.columns <- function(tree,aln) {
  aln.length <- length(aln[1,])

  score.df <- data.frame()
  for (i in 1:aln.length) {
    aln.column <- aln[,i]

    nongap.seqs <- names(aln.column[aln.column != '-'])
    gap.seqs <- names(aln.column[aln.column == '-'])

    # Get the non-gap branch length.
    if (length(nongap.seqs) == 1) {
      nongap.node <- node.with.label(tree,nongap.seqs[1])
      nongap.bl <- branch.length(tree,nongap.node)
    } else {
      nongap.tree <- drop.tip(tree,gap.seqs)
      nongap.bl <- sum(nongap.tree$edge.length)
    }

    nongap.str <- paste(nongap.seqs,collapse=';')
    cur.df <- data.frame(
      pos=i,
      score=nongap.bl,
      nongap.str=nongap.str,
      stringsAsFactors=F
    )
    score.df <- rbind(score.df,cur.df)
  }
  score.df <- score.df[order(score.df$score,score.df$pos),]
  return(score.df)
}

remove.gaps <- function(sim,tolerance) {
  aln <- sim$alignment
  tree <- sim$phylo
  col.scores <- score.aln.columns(tree,aln)

  # Store the deletion markers in a separate data frame.
  deletion.df <- NULL
  
  # Get the mean sequence length
  lengths <- apply(aln,1,function(x) {
    x <- x[x != '-']
    return(stringLength(paste(x,collapse='')))
  })
  mean.seq.length <- mean(lengths)

  repeat {
    aln.length <- length(aln[1,])  
    ratio <- aln.length / mean.seq.length
    if (ratio < tolerance) {
      break;
    }

    # Take the next site from the sorted scores
    lowest.scores <- col.scores[1,]
    col.scores <- col.scores[-c(1),]
    cur.pos <- lowest.scores$pos # Current position of lowest-scoring column.
    cur.score <- lowest.scores$score # The current column score.
    cur.str <- lowest.scores$nongap.str # The current column's nongap pattern.

    # Grab the entire 'current chunk' of alignment which has the same
    # score and non-gap string.
    repeat {
      first.pos <- col.scores[1,'pos']
      first.score <- col.scores[1,'score']
      first.str <- col.scores[1,'nongap.str']
      if (cur.score == max(col.scores$score)) {
        lowest.scores <- NULL
        cur.ratio <- length(aln[1,]) / mean.seq.length
        print(sprintf("Nothing left to remove at ratio %.2f!",cur.ratio))
        break;
      }
      if (first.pos == cur.pos + 1 && first.score == cur.score && first.str == cur.str) {
        cur.pos <- col.scores[1,]$pos
        lowest.scores <- rbind(lowest.scores,col.scores[1,])
        col.scores <- col.scores[-1,]
      } else {
        #print("Done!")
        break;
      }
    }

    if (is.null(lowest.scores)) {
      break;
    }
    
    # remove.us should be a contiguous vector of integers,
    # representing the set of columns to remove.
    remove.us <- lowest.scores$pos

    if (any(diff(remove.us) > 1)) {
      print("ERROR: Removing a non-consecutive set of columns!")
    }
    
    #print(paste("Removing at ",paste(remove.us[1])))

    # Go through columns from right to left, making sure to update
    # the new positions of columns on the right side of the splice.
    rev.pos <- rev(remove.us)
    for (i in 1:length(rev.pos)) {
      cur.pos <- rev.pos[i]
      aln <- splice.column(aln, cur.pos)

      # Update new positions of column scores
      above <- which(col.scores$pos > cur.pos)
      col.scores[above,]$pos <- col.scores[above,]$pos - 1

      # Update new positions of deletion locations.
      if (!is.null(deletion.df) && nrow(deletion.df) > 0) {
        above <- which(deletion.df$pos > cur.pos)
        deletion.df[above,]$pos <- deletion.df[above,]$pos - 1
      }
    }


    # Add a single entry to the data frame of deletions, to be used
    # by the plot function to indicate deletion points.
    cur.deletion <- data.frame(
      pos = cur.pos, # The first position AFTER the deletion splice.
      length = length(rev.pos), # The length of the block removed.
      nongap.str = as.character(cur.str), # nongap sequence IDs for this deletion
      stringsAsFactors=F
    )
    #print(cur.deletion)
    deletion.df <- rbind(deletion.df,cur.deletion)
    
  }

  sim.temp <- PhyloSim();
  setAlignment(sim.temp,aln)
  sim.temp$phylo <- tree

  plot(sim.temp,plot.chars=F,indels=deletion.df)
}

splice.column <- function(aln,pos) {
  return(aln[,-pos])
}

dummy.aln <- function(indel.rate=0.03,root.length=100,tree.length=1) {
  PSIM_FAST = TRUE;
  assign('PSIM_FAST',TRUE,envir=.GlobalEnv)

  # Create the root sequence.
  root.sequence<-Sequence(length=root.length)

  # Set up indel distribution:
  max.indel <- 20
  id.dist<-exp(max.indel:1)/sum(exp(max.indel:1))

  
  # Coding options.
  #coding.pos      <-c(1:root.length)
  #codon.alph  <- CodonAlphabet()
  #setAlphabets(root.sequence, list(codon.alph), coding.pos)
  #gy94    <-GY94(kappa=2, omega.default=0.2, scale.nuc=TRUE)
  #del.c    <- DiscreteDeletor(rate=indel.rate, sizes=1:6, probs=id.dist  )
  #ins.c   <- DiscreteInsertor(rate=indel.rate, sizes=1:6, probs=id.dist )
  #ins.c$templateSeq   <- CodonSequence(length=1,processes=list(list( gy94, del.c, ins.c ) ))
  #setProcesses(root.sequence, list(list(gy94, del.c, ins.c)), coding.pos)
  #setStates(root.sequence, "ATG", coding.pos[1]);            # Set the state.
  #setRateMultipliers(root.sequence, gy94, 0, coding.pos[1])  # Make the site inv
  #setDeletionTolerance(root.sequence, del.c, 0, coding.pos[1]);        # Make the site reject deletions.
  #setInsertionTolerance(root.sequence,ins.c ,0, coding.pos[1]);

  # Noncoding options.
  noncoding.pos   <- c(1:root.length)
  nuc.alph    <- NucleotideAlphabet();
  setAlphabets(root.sequence, list(nuc.alph), noncoding.pos)
  k80     <-K80(rate.params=list("Alpha"=2,"Beta"=1), base.freqs=c(2,1,2,1)/4)
  del.nc   <- DiscreteDeletor(rate=indel.rate, sizes=1:max.indel, probs=id.dist )
  ins.nc  <- DiscreteInsertor(rate=indel.rate, sizes=1:max.indel, probs=id.dist )
  ins.nc$templateSeq  <- NucleotideSequence(length=1,processes=list(list( k80, del.nc, ins.nc ) ))
  setProcesses(root.sequence, list(list(k80, del.nc, ins.nc)), noncoding.pos)

  sampleStates(root.sequence)
  print(root.sequence)

  data(bird.families)
  tree <- bird.families
#  tree <- read.tree("~/work/phylosim/examples/data/slr_bigtree.nh")
  tree <- scale.length.to(tree,tree.length)

  sim <- PhyloSim(
    phylo=tree,
    root.seq=root.sequence
  )

  Simulate(sim)
  saveAlignment(sim,skip.internal=TRUE,file="asdf.fasta")
  write.tree(sim$phylo,file="asdf.nh")
  readAlignment(sim,file="asdf.fasta")
  return(sim)
}

do.big.sim <- function() {
  sim <- dummy.aln(indel.rate=0.03, root.length=200, tree.length = 1)
  assign("sim",sim,envir=.GlobalEnv)
}

tree.aln.subset <- function(tree,aln,keep.seqs,remove.blank.columns=T) {
  drop.seqs <- tree[!(tree$tip.label %in% keep.seqs)]$tip.label
  print(drop.seqs)
  sub.tree <- drop.tip(tree,gap.seqs)

  
  
}



### Some miscellaneous tree utility functions.
branch.length <- function(phylo,node) {
  edge.index <- which(phylo$edge[,2]==node)
  bl <- phylo$edge.length[edge.index]
  if (length(bl)==0) {
    bl <- 0
  }
  return(bl)
}

scale.tree.to <- function(tree,new.total.length) {
  cur.total.length <- sum(tree$edge.length)
  tree$edge.length <- (new.total.length / cur.total.length) * tree$edge.length
  return(tree)
}

scale.length.to <- function(tree,new.max.length) {
  cur.length.to.root <- max.length.to.root(tree)
  tree$edge.length <- (new.max.length / cur.length.to.root) * tree$edge.length
  return(tree)
}

node.with.label <- function(tree,label) {
  return(which(tree$tip.label %in% label))
}

    # Extracts the length of the branch above the given node. Returns 0 if the node is root.
    branch.length <- function(phylo,node) {
      edge.index <- which(phylo$edge[,2]==node)
      bl <- phylo$edge.length[edge.index]
      if (length(bl)==0) {
        bl <- 0
      }
      return(bl)
    }

    # The maximum root-to-tip length in the tree.
    max.length.to.root <- function(phylo) {
      max.length <- 0
      for (i in 1:length(phylo$tip.label)) {
        cur.length <- length.to.root(phylo,i)
        max.length <- max(max.length,cur.length)
      }
      return(max.length)
    }

    # The length from the root to the given node. Can be given either as a node ID or a tip label.
    length.to.root <- function(phylo,node) {
      tip.index <- node
      if (is.character(node)) {
        tip.index <- which(phylo$tip.label==node)
      }
      
      cur.node.b <- tip.index

      p.edges <- phylo$edge
      p.lengths <- phylo$edge.length
      
      length <- 0
      while(length(which(p.edges[,2]==cur.node.b)) > 0) {
        cur.edge.index <- which(p.edges[,2]==cur.node.b)
        cur.edge.length <- p.lengths[cur.edge.index]
        length <- length + cur.edge.length
        cur.node.a <- p.edges[cur.edge.index,1]
        cur.node.b <- cur.node.a # Move up to the next edge
      }
      return(length)
    }
