get.numbers <- function(sites) {
  if (!exists('gene.types')) {
    gene.types <- classify.genes(sites)
    assign('gene.types',gene.types,envir=.GlobalEnv)
  }

  # 'Neutral' type genes are ones with < 0.5 neg.f
  neg.thresh <- 0.5
  # 'positive' types are ones with > 0.01 pos.f
  pos.thresh <- 0.01

  neutral.genes <- subset(gene.types,neg.f < neg.thresh & pos.f >= pos.thresh)
  positive.genes <- subset(gene.types,neg.f >= neg.thresh & pos.f >= pos.thresh)
  purifying.genes <- subset(gene.types,(neg.f >= neg.thresh & pos.f < pos.thresh) | 
                                       (neg.f < neg.thresh & pos.f < pos.thresh) )

  df <- data.frame(
    neutral.n = nrow(neutral.genes),
    purifying.n = nrow(purifying.genes),
    positive.n = nrow(positive.genes)
  )
  print(df)
}

classify.genes <- function(sites) {
  sub.above.one = subset(sites, omega > 1)
  pvals = sub.above.one$pval
  pvals.adj = p.adjust(pvals,method="fdr")
  sub.above.one[,'pval.adj'] = pvals.adj

  classify.gene <- function(gene.sites) {
    mean.dnds <- mean(gene.sites$omega)
    pos.sites <- subset(gene.sites,omega > 1 & pval < 0.01)
    neg.sites <- subset(gene.sites,omega < 1 & pval < 0.01)
    neutral.sites <- subset(gene.sites,pval > 0.01)
    
    n <- nrow(gene.sites)

    n.genes <- length(unique(gene.sites$node_id))
    n.domains <- length(unique(gene.sites$domain))

    asdf <- data.frame(
      n = n,
      n.genes = n.genes,
      n.domains = n.domains,
      pos.n = nrow(pos.sites),
      neg.n = nrow(neg.sites),
      neutral.n = nrow(neutral.sites),
      pos.f = nrow(pos.sites)/n,
      neg.f = nrow(neg.sites)/n,
      neutral.f = nrow(neutral.sites)/n,
      dnds = mean.dnds
    )
    return(asdf)
  }
  
  library(plyr)
  gene.types <- ddply(sites,.(node_id),classify.gene)

  domain.sites <- subset(sites,!is.na(domain))
  domain.types <- ddply(domain.sites,.(domain),classify.gene)

  assign('domain.types',domain.types,envir=.GlobalEnv)
  assign('gene.types',gene.types,envir=.GlobalEnv)
}

