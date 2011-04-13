get.numbers <- function(sites) {
  if (!exists('gene.types')) {
    classify.genes(sites)
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

summarize.sites <- function(sites) {
  # BH-correction is applied separately to positive and negative p-values
  sites[,'pval'] <- 1 - pchisq(sites[,'lrt_stat'],1)
  sites[sites$omega < 1,'pval.adj'] <- p.adjust(sites[sites$omega < 1,'pval'],'fdr')
  sites[sites$omega >= 1,'pval.adj'] <- p.adjust(sites[sites$omega >= 1,'pval'],'fdr')

  sub.above.one <- subset(sites, omega > 1)

  pos.sites <- subset(sub.above.one, pval.adj < 0.05)
  pos.domains <- subset(pos.sites, !is.na(pfam_domain))
  
  n.total.domains <- nrow(unique(sites[,c('node_id','pfam_domain')]))
  n.total.domain.types <- length(unique(sites$pfam_domain))
  n.total.genes <- length(unique(sites$node_id))

  n.pos.domains <- nrow(unique(pos.domains[,c('node_id','pfam_domain')]))
  n.pos.domain.types <- length(unique(pos.domains$pfam_domain))
  n.pos.genes <- length(unique(pos.sites$node_id))

  n.pos.sites <- nrow(pos.sites)
  n.sites <- nrow(sites)  
  f.pos.sites <- n.pos.sites / n.sites

  pos.sites <- subset(sites, omega > 1 & pval < 0.01)
  neg.sites <- subset(sites, omega < 1 & pval < 0.01)
  neutral.sites <- subset(sites, pval > 0.01)

  `f.less.0.5` <- nrow(subset(sites, omega < 0.5)) / nrow(sites)
  `f.less.1` <- nrow(subset(sites, omega < 1)) / nrow(sites)
  `f.gt.1` <- nrow(subset(sites, omega > 1)) / nrow(sites)
  `f.gt.1.5` <- nrow(subset(sites, omega > 1.5)) / nrow(sites)

  return(data.frame(
    n.sites = n.sites,
    'f.less.0.5' = `f.less.0.5`,
    'f.less.1' = `f.less.1`,
    'f.gt.1' = `f.gt.1`,
    'f.gt.1.5' = `f.gt.1.5`,

    n.total.domains = n.total.domains,
    n.total.domain.types = n.total.domain.types,
    n.total.genes = n.total.genes,
    
    n.pos.domains = n.pos.domains,
    n.pos.domain.types = n.pos.domain.types,
    n.pos.genes = n.pos.genes,

    pos.n = n.pos.sites,
    pos.f = nrow(pos.sites)/n.sites,
    neg.f = nrow(neg.sites)/n.sites,
    neutral.f = nrow(neutral.sites)/n.sites
   ))
}

classify.genes <- function(sites) {
  # BH-correction is applied separately to positive and negative p-values
  sites[sites$omega < 1,'pval.adj'] <- p.adjust(sites[sites$omega < 1,'pval'],'fdr')
  sites[sites$omega >= 1,'pval.adj'] <- p.adjust(sites[sites$omega >= 1,'pval'],'fdr')

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
  library(doBy)
  gene.types <- ddply(sites,.(node_id),classify.gene)

  domain.sites <- subset(sites,!is.na(domain))
  domain.types <- ddply(domain.sites,.(domain),classify.gene)
  
  # Restrict to well-covered domains: > 5 genes, > 500 sites
  domain.types <- subset(domain.types, n.genes > 5 & n > 500)
  domain.types <- orderBy(~-pos.f,data=domain.types)

  assign('domain.types',domain.types,envir=.GlobalEnv)
  assign('gene.types',gene.types,envir=.GlobalEnv)
}

summarize.sites.sets <- function() {
  collect.df <- data.frame()
  for (set in c(1,2,3,4)) {
    file.name <- paste('sites_',set,'.Rdata',sep='')
    print(file.name)
    load(file.name)
    print("Summarizing...")
    cur.df <- summarize.sites(sites)
    collect.df <- rbind(collect.df,cur.df)
  }
  write.csv(collect.df,file='sites_summaries.csv')
}

classify.genes.sets <- function() {

  for (set in c(1)) {
    file.name <- paste('sites_',set,'.Rdata',sep='')
    print(file.name)
    load(file.name)
    print("Classifying...")
    classify.genes(sites)
    print("Saving...")
    out.file <- paste('pos_domains_',set,'.csv',sep='')
    write.csv(pos.domains,file=out.file)
  }
}
