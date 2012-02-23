require(ape)

# Extracts the NHX annotations from a tree and returns a list with the annotations and the
# tree string (with NHX stuff removed).
tree.read.nhx <- function(str) {
  if (file.exists(str)) {
    # We're reading a file -- load the file into a string.
    str <- readLines(str)
  }

  nhx.matches <- gregexpr("[,\\)\\( ]?([^,\\)\\(]+)?(:?\\d*\\.?\\d*)?\\[.*?\\]", str)
  matches <- nhx.matches[[1]]
  match.pos <- as.numeric(matches)
  match.len <- attr(matches, 'match.length')

  if (match.len[1] == -1) {
    return(tree.read(text=str))
  }

  match.pos <- match.pos + 1
  match.len <- match.len - 1

  #print(match.pos)
  nhx.strings <- substring(str, match.pos, match.pos+match.len-1)
  #print(nhx.strings)

  labels <- gsub("\\[&&NHX.*", "", nhx.strings)
  labels.no.bl <- gsub(":.*", "", labels)

  #print(labels.no.bl)

  # Go through and splice the stripped-down labels back into the string.
  for (i in 1:length(match.pos)) {
    new.label <- gsub(labels.no.bl[i], paste('zzz',i,'zzz',sep=''), labels[i])
    str <- paste(substr(str, 0, match.pos[i]-1 ), new.label, substr(str, match.pos[i] + match.len[i], nchar(str)), sep='')
    match.pos <- match.pos - match.len[i] + nchar(new.label)
  }
  
  # Parse the Phylo object from the cleaned-up string.
  #print(str)
  tree <- tree.read(text=str)

  #print(str(tree))

  # Create a list of NHX annotations keyed by the node ID.
  tag.values <- gsub(".*\\[&&NHX:(.*)\\]", "\\1", nhx.strings)
  tagval.list <- strsplit(tag.values, ":")
  names(tagval.list) <- labels.no.bl

  library(plyr)
  map.list <- llply(tagval.list, function(x) {
    list.out <- list()
    cur.list <- strsplit(x, "=")
    if (length(cur.list) > 0) {
      for (i in 1:length(cur.list)) {
        vec <- cur.list[[i]]
        list.out[vec[1]] = vec[2]
      }
    }
    return(list.out)
  })

  # Replace the labels with the true labels.
  tree$.tags <- list()
  for (i in 1:(tree$Nnode+length(tree$tip.label))) {
    tree$.tags[[i]] <- list()
  }

  for (i in 1:length(match.pos)) {
    cur.node <- node.with.label(tree, paste('zzz', i, 'zzz', sep=''))
    
    leaf <- tree.is.leaf(tree, cur.node)
    real.node.name <- names(map.list)[i]
   if (leaf) {
      tree$tip.label[cur.node] <- real.node.name
    } else {
      tree$node.label[cur.node-length(tree$tip.label)] <- real.node.name
    }
    tree$.tags[[cur.node]] <- map.list[[i]]
  }

  #tree$.tags <- map.list
  return(tree)
}

tree.write <- function(phylo, f) {
  write.tree(phylo, f)
}

# Tries to get a tag from 
tree.get.tag <- function(phylo, node, tag) {
  tag <- phylo$.tags[[node]][[tag]]
  if(is.null(tag)) {return('')}
  return(tag)
}

tree.get.tags <- function(phylo, node) {
  tags <- phylo$.tags[[node]]
  if (is.null(tags)) {
    return(list())
  } else {
    return(tags)
  }
}

tree.foreach <- function(phylo, fn, leaves=T, nodes=T) {
  n.leaves <- length(phylo$tip.label)
  n.internals <- phylo$Nnode

  indices <- c()
  if (leaves) {
    indices <- c(indices, 1:n.leaves)
  }
  if (nodes) {
    indices <- c(indices, (n.leaves+1):(n.leaves+n.internals))
  }
  
  for (node in indices) {
    do.call(fn, as.list(c(phylo, node)))
  }
}

tree.as.data.frame <- function(tree, order.cladewise=T) {
  tags <- c()

  tree.foreach(tree, function(phylo, node) {
    cur.tags <- tree.get.tags(phylo, node)
    tags <<- c(tags, names(cur.tags))
  })

  tree.df <- data.frame(stringsAsFactors=F)
  tree.foreach(tree, function(phylo, node) {
    cur.tags <- tree.get.tags(phylo, node)
    cur.tags[['Label']] <- tree.label.for.node(phylo, node)
    cur.tags[['Depth']] <- tree.leaves.beneath(phylo, node)
    cur.tags[['id']] <- node
    tree.df <<- rbind.fill(tree.df, as.data.frame(cur.tags, stringsAsFactors=F))
  })

  if (order.cladewise) {
    tree.df <- tree.df[tree.order.nodes(tree, include.internals=T),]
  }
  tree.df
}

## read.tree.R (2010-09-27)

##   Read Tree Files in Parenthetic Format

## Copyright 2002-2010 Emmanuel Paradis, Daniel Lawson and Klaus Schliep

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

tree.build <- function(tp)
{
    add.internal <- function() {
        edge[j, 1] <<- current.node
        edge[j, 2] <<- current.node <<- node <<- node + 1L
        index[node] <<- j # set index
        j <<- j + 1L
    }
    add.terminal <- function() {
        edge[j, 1] <<- current.node
        edge[j, 2] <<- tip
        index[tip] <<- j # set index
        X <- unlist(strsplit(tpc[k], ":"))
        tip.label[tip] <<- X[1]
        edge.length[j] <<- as.numeric(X[2])
        k <<- k + 1L
        tip <<- tip + 1L
        j <<- j + 1L
    }
    go.down <- function() {
        l <- index[current.node]
        X <- unlist(strsplit(tpc[k], ":"))
        node.label[current.node - nb.tip] <<- X[1]
        edge.length[l] <<- as.numeric(X[2])
        k <<- k + 1L
        current.node <<- edge[l, 1]
    }
    if (!length(grep(",", tp))) {
        obj <- list(edge = matrix(c(2L, 1L), 1, 2))
        tp <- unlist(strsplit(tp, "[\\(\\):;]"))
        obj$edge.length <- as.numeric(tp[3])
        obj$Nnode <- 1L
        obj$tip.label <- tp[2]
        if (tp[4] != "") obj$node.label <- tp[4]
        class(obj) <- "phylo"
        return(obj)
    }

    tpc <- unlist(strsplit(tp, "[\\(\\),;]"))
    tpc <- tpc[nzchar(tpc)]
    ## the following 2 lines are (slightly) faster than using gsub()
    tsp <- unlist(strsplit(tp, NULL))
    skeleton <- tsp[tsp %in% c("(", ")", ",", ";")]
    nsk <- length(skeleton)
    nb.node <- sum(skeleton == ")")
    nb.tip <- sum(skeleton == ",") + 1
    ## We will assume there is an edge at the root;
    ## if so, it will be removed and put into a vector
    nb.edge <- nb.node + nb.tip
    node.label <- character(nb.node)
    tip.label <- character(nb.tip)

    edge.length <- numeric(nb.edge)
    edge <- matrix(0L, nb.edge, 2)
    current.node <- node <- as.integer(nb.tip + 1) # node number
    edge[nb.edge, 2] <- node
    index <- numeric(nb.edge + 1) # hash index to avoid which
    index[node] <- nb.edge

    ## j: index of the line number of edge
    ## k: index of the line number of tpc
    ## tip: tip number
    j <- k <- tip <- 1L

    for (i in 2:nsk) {
        if (skeleton[i] == "(") add.internal() # add an internal branch (on top)
        if (skeleton[i] == ",") {
            if (skeleton[i - 1] != ")") add.terminal() # add a terminal branch
        }
        if (skeleton[i] == ")") {
            if (skeleton[i - 1] != ")") { # add a terminal branch and go down one level
                add.terminal()
                go.down()
            }
            if (skeleton[i - 1] == ")") go.down() # go down one level
        }
    }

    edge <- edge[-nb.edge, ]
    obj <- list(edge = edge, Nnode = nb.node, tip.label = tip.label)
    root.edge <- edge.length[nb.edge]
    edge.length <- edge.length[-nb.edge]
    if (!all(is.na(edge.length))) # added 2005-08-18
        obj$edge.length <- edge.length
    if (is.na(node.label[1])) node.label[1] <- ""
    if (any(nzchar(node.label))) obj$node.label <- node.label
    if (!is.na(root.edge)) obj$root.edge <- root.edge
    class(obj) <- "phylo"
    obj
}

tree.read <- function(file = "", text = NULL, tree.names = NULL, skip = 0, remove.whitespace=F,
    comment.char = "#", keep.multi = FALSE, ...)
{
    unname <- function(treetext) {
        nc <- nchar(treetext)
	tstart <- 1
	while (substr(treetext, tstart, tstart) != "(" && tstart <= nc)
            tstart <- tstart + 1
	if (tstart > 1)
            return(c(substr(treetext, 1, tstart - 1),
                     substr(treetext, tstart, nc)))
	return(c("", treetext))
    }
    if (!is.null(text)) {
        if (!is.character(text))
          stop("argument `text' must be of mode character")
        tree <- text
    } else {
        tree <- scan(file = file, what = "", sep = "\n", quiet = TRUE,
                     skip = skip, comment.char = comment.char, ...)
    }
    ## Suggestion from Eric Durand and Nicolas Bortolussi (added 2005-08-17):
    if (identical(tree, character(0))) {
        warning("empty character string.")
        return(NULL)
    }
    
    if (remove.whitespace) {
      tree <- gsub("[ \t]", "", tree)
    }
    tree <- unlist(strsplit(tree, NULL))
    y <- which(tree == ";")
    Ntree <- length(y)
    x <- c(1, y[-Ntree] + 1)
    ## Suggestion from Olivier François (added 2006-07-15):
    if (is.na(y[1])) return(NULL)
    STRING <- character(Ntree)
    for (i in 1:Ntree)
        STRING[i] <- paste(tree[x[i]:y[i]], sep = "", collapse = "")

    tmp <- unlist(lapply(STRING, unname))
    tmpnames <- tmp[c(TRUE, FALSE)]
    STRING <- tmp[c(FALSE, TRUE)]
    if (is.null(tree.names) && any(nzchar(tmpnames)))
        tree.names <- tmpnames

    colon <- grep(":", STRING)
    if (!length(colon)) {
        obj <- lapply(STRING, clado.build)
    } else if (length(colon) == Ntree) {
        obj <- lapply(STRING, tree.build)
    } else {
        obj <- vector("list", Ntree)
        obj[colon] <- lapply(STRING[colon], tree.build)
        nocolon <- (1:Ntree)[!1:Ntree %in% colon]
        obj[nocolon] <- lapply(STRING[nocolon], clado.build)
    }
    for (i in 1:Ntree) {
        ## Check here that the root edge is not incorrectly represented
        ## in the object of class "phylo" by simply checking that there
        ## is a bifurcation at the root
        ROOT <- length(obj[[i]]$tip.label) + 1
        if(sum(obj[[i]]$edge[, 1] == ROOT) == 1 && dim(obj[[i]]$edge)[1] > 1)
            stop(paste("There is apparently two root edges in your file: cannot read tree file.\n  Reading Newick file aborted at tree no.", i, sep = ""))
    }
    if (Ntree == 1 && !keep.multi) obj <- obj[[1]] else {
        if (!is.null(tree.names)) names(obj) <- tree.names
        class(obj) <- "multiPhylo"
    }
    obj
}

tree.remove.node.labels <- function(phylo) {
  phylo$node.label <- NULL
  phylo
}

tree.remove.branchlengths <- function(phylo, push.to.tips=F) {
  n.leaves <- length(phylo$tip.label)
  n.nodes <- length(phylo$tip.label)+phylo$Nnode

  max.depth <- 0
  for (i in 1:n.nodes) {
    depth <- tree.depth.to.root(phylo, i)
    max.depth <- max(depth, max.depth)
  }
  max.depth <- max.depth + 1

  for (i in 1:n.nodes) {
    cur.depth <- tree.depth.to.root(phylo,i)
    parent.node <- tree.parent.node(phylo, i)
    edge.index <- which(phylo$edge[,2]==i)

    is.leaf <- i <= n.leaves

    if (is.leaf) {
      cur.count <- 1
    } else {
      cur.count <- tree.leaves.beneath(phylo, i)
    }

    if (parent.node > -1) {
      parent.count <- tree.leaves.beneath(phylo, parent.node)
      if (push.to.tips) {
        count.diff <- parent.count - cur.count
        #print(paste(i, count.diff))
        phylo$edge.length[edge.index] <- count.diff
      } else {
        if (is.leaf) {
          # Branch length should equal diff. between depth to root and max depth.
          cur.d <- tree.depth.to.root(phylo, i)
          phylo$edge.length[edge.index] <- max.depth - cur.d
        } else {
          phylo$edge.length[edge.index] <- 1
        }
      }
    } else {
      phylo$edge.length[edge.index] <- 0
    }
  }

  phylo
}

tree.leaves.beneath <- function(phylo, node) {
  if (is.leaf(phylo, node)) {
    return(1)
  }
  cld <- tree.extract.clade(phylo, node)
  length(cld$tip.label)
}

# The length from the root to the given node. Can be given either as a node ID or a tip label.
tree.depth.to.root <- function(phylo,node) {
  tip.index <- node
  if (is.character(node)) {
    tip.index <- which(phylo$tip.label==node)
  }
  cur.node.b <- tip.index
  p.edges <- phylo$edge

  length <- 0
  while(length(which(p.edges[,2]==cur.node.b)) > 0) {
    cur.edge.index <- which(p.edges[,2]==cur.node.b)
    cur.edge.length <- 1
    length <- length + cur.edge.length
    cur.node.a <- p.edges[cur.edge.index,1]
    cur.node.b <- cur.node.a # Move up to the next edge
  }
  return(length)
}

# Finds the node with a given label.
tree.node.with.label <- function(tree,label) {
  all.labels <- c(tree$tip.label,tree$node.label)
  return(which(all.labels %in% label))
}

tree.label.for.node <- function(tree, node) {
  if (node <= length(tree$tip.label)) {
    return(tree$tip.label[node])
  } else if (node <= (tree$Nnode + length(tree$tip.label))) {
    node.label.index <- node - length(tree$tip.label)
    return(tree$node.label[node.label.index])
  }
}

# Extracts the length of the branch above the given node. Returns 0 if the node is root.
tree.branch.length <- function(phylo,node) {
  edge.index <- which(phylo$edge[,2]==node)
  bl <- phylo$edge.length[edge.index]
  if (length(bl)==0 || is.na(bl)) {
    bl <- 0
  }
  return(bl)
}

tree.scale.by <- function(phylo, factor) {
  phylo$edge.length <- phylo$edge.length * factor
  return(phylo)
}

tree.scale.to <- function(phylo, total.length) {
  cur.total <- tree.total.branch.length(phylo)
  scale.factor <- total.length / cur.total
  tree.scale.by(phylo, scale.factor)
}

tree.total.branch.length <- function(phylo) {
  sum(phylo$edge.length)
}

# The maximum root-to-tip length in the tree.
tree.max.length.to.root <- function(phylo) {
  max.length <- max(tree.lengths.to.root(phylo))
  if (is.na(max.length)) {
    max.depth <- 0
    for (i in 1:length(phylo$tip.label)) {
      cur.depth <- tree.depth.to.root(phylo,i)
      max.depth <- max(max.depth, cur.depth)
    }
    return(max.depth)
  }
  return(max.length)
}

tree.mean.path.length <- function(phylo) {
  mean(tree.lengths.to.root(phylo))
}

tree.lengths.to.root <- function(phylo) {
  lengths <- c()
  if (length(phylo$tip.label) == 0) {
    return(NA)
  }
  for (i in 1:length(phylo$tip.label)) {
    lengths[i] <- tree.length.to.root(phylo,i)
  }
  lengths
}

# The length from the root to the given node. Can be given either as a node ID or a tip label.
tree.length.to.root <- function(phylo,node) {
  tip.index <- node
  if (is.character(node)) {
    tip.index <- which(phylo$tip.label==node)
  }
  cur.node.b <- tip.index

  p.edges <- phylo$edge
  p.lengths <- phylo$edge.length

  if(is.null(p.lengths)) {
    p.lengths <- rep(1, length(p.edges[,1]))
  }

  length <- 0
  while(length(which(p.edges[,2]==cur.node.b)) > 0) {
    cur.edge.index <- which(p.edges[,2]==cur.node.b)
    cur.edge.length <- p.lengths[cur.edge.index]
    if (length(cur.edge.length) == 0 || is.na(cur.edge.length)) {
      cur.edge.length <- 0
    }
    length <- length + cur.edge.length
    cur.node.a <- p.edges[cur.edge.index,1]
    cur.node.b <- cur.node.a # Move up to the next edge
  }
  return(length)
}

tree.remove.leaf <- function(tree, x) {
  tree <- drop.tip(tree, x)
  tree
}

tree.extract.subtree <- function(tree, x) {
  not.in.set <- setdiff(tree$tip.label, x)
  tree <- drop.tip(tree, not.in.set)
  tree
}

tree.is.leaf <- function(phylo,node) {
  return(node <= length(phylo$tip.label))
}

 # Extracts a list of child node IDs for the given node. Returns (-1,-1) if the node is a leaf.
tree.child.nodes <- function(phylo,node) {
  edge.indices <- which(phylo$edge[,1]==node)
  nodes <- phylo$edge[edge.indices,2]
  if (length(nodes)==0) {
    nodes <- list(c(-1,-1))
  } else {
    nodes <- list(nodes)
  }
  return(list(nodes))
}

# Extracts the parent node ID for the given node. Returns -1 if the node is root.
tree.parent.node <- function(phylo,node) {
  edge.index <- which(phylo$edge[,2]==node)
  node <- phylo$edge[edge.index,1]
  if (length(node)==0) {
    node <- -1
  }
  return(node)
}

tree.order.nodes <- function(phylo, include.internals=T) {
  phylo <- reorder(phylo, order="cladewise");

  df.list <- phylo.layout.df(phylo,
    layout.ancestors=T
  )
  
  nodes <- df.list$nodes

  if (!include.internals) {
    nodes <- subset(nodes, is.leaf==FALSE)
  }

  nodes[order(nodes$y), 'node']
}

sort.df.by.tree <- function(tree.df, tree) {
  phylo.order <- tree.order.nodes(tree)

  df.order <- match(phylo.order, tree.df$id)
  tree.df[df.order, ]
}

is.tree <- function(x) {
  if (!is.null(x$edge) && x$Nnode > 0) {
    TRUE
  } else {
    FALSE
  }
}