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
                                       (neg.f < neg.thresh & pos.f < pos.thresh) 
  )

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

classify.genes <- function(sites, genes) {
  # BH-correction is applied separately to positive and negative p-values
  #sites[sites$omega < 1,'pval.adj'] <- p.adjust(sites[sites$omega < 1,'pval'],'fdr')
  #sites[sites$omega >= 1,'pval.adj'] <- p.adjust(sites[sites$omega >= 1,'pval'],'fdr')

  sites$lrt_stat <- sites$lrt_stat * sign(sites$omega - 0.999)  

  q.t <- qchisq(0.95, df=1)

  neg.stats <- sites$lrt_stat[sites$lrt_stat < 0]
  neg.q <- quantile(neg.stats, probs=c(0, 0.33, 0.67, 1))
  print(neg.q)

  pos.stats <- sites$lrt_stat[sites$lrt_stat > 0]
  pos.q <- quantile(pos.stats, probs=c(0, 0.33, 0.67, 1))
  print(pos.q)

  classify.gene <- function(gene.sites) {
    t <- 3

    # Classify into pos, neg, and neutral.    
    n <- nrow(gene.sites)
    pos.sites <- sum(gene.sites$lrt_stat > q.t)
    neg.sites <- sum(gene.sites$lrt_stat < -q.t)
    neutral.sites <- sum(gene.sites$lrt_stat >= -q.t & gene.sites$lrt_stat < q.t)
    
    # Classify negative sites into strong, moderate, weak
    n.neg <- sum(gene.sites$lrt_stat < 0) + 1
    neg.1 <- sum(gene.sites$lrt_stat > neg.q[1] & gene.sites$lrt_stat <= neg.q[2])
    neg.2 <- sum(gene.sites$lrt_stat > neg.q[2] & gene.sites$lrt_stat <= neg.q[3])
    neg.3 <- sum(gene.sites$lrt_stat > neg.q[3] & gene.sites$lrt_stat <= neg.q[4])

    # Classify positive sites into strong, moderate, weak
    n.pos <- sum(gene.sites$lrt_stat > 0) + 1
    pos.1 <- sum(gene.sites$lrt_stat > pos.q[1] & gene.sites$lrt_stat <= pos.q[2])
    pos.2 <- sum(gene.sites$lrt_stat > pos.q[2] & gene.sites$lrt_stat <= pos.q[3])
    pos.3 <- sum(gene.sites$lrt_stat > pos.q[3] & gene.sites$lrt_stat <= pos.q[4])

    asdf <- data.frame(
      n = n,
      #n.genes = n.genes,
      #n.domains = n.domains,
      pos = pos.sites/n,
      neg = neg.sites/n,
      neutral = neutral.sites/n,
      neg.1 = neg.1/n.neg,
      neg.2 = neg.2/n.neg,
      neg.3 = neg.3/n.neg,
      pos.1 = pos.1/n.pos,
      pos.2 = pos.2/n.pos,
      pos.3 = pos.3/n.pos
    )
    print(asdf)
    return(asdf)
  }
  
  library(plyr)
  library(doBy)
  node.ids <- unique(sites$node_id)
  gene.types <- ddply(sites,.(node_id),classify.gene)

  gene.types <- merge(genes, gene.types, by=c('node_id'))

  assign('gene.types', gene.types, envir=.GlobalEnv)
  save(gene.types, file="gene.types.Rdata")

  return()

  domain.sites <- subset(sites,!is.na(domain))
  domain.types <- ddply(domain.sites,.(domain),classify.gene)

  # Restrict to well-covered domains: > 5 genes, > 500 sites
  domain.types <- subset(domain.types, n.genes > 5 & n > 500)
  domain.types <- orderBy(~-pos.f,data=domain.types)

  assign('domain.types',domain.types,envir=.GlobalEnv)
  assign('gene.types',gene.types,envir=.GlobalEnv)

  good.domains <- subset(domain.types, n.genes > 5 & n > 500)
  sorted.domains <- good.domains

  sorted.domains <- orderBy(~ -pos.f, data=sorted.domains)
  assign('pos.domains',sorted.domains,envir=.GlobalEnv)
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

ternary.plot <- function(a, b, c, clr=NA) {
  library(ggplot2)
  n <- length(a)

  # a: up
  # b: down-right
  # c: down-left
  y <- a * sin(pi/3)
  x <- b + a * cos(pi/3)

  # Create binned histograms along the edges.
  brks <- seq(from=-1, to=1, length.out=20)
  h <- hist(a-b, breaks=brks, plot=F)
  ab.df <- data.frame(mid=h$mids, count=h$counts/n)
  h <- hist(b-c, breaks=brks, plot=F)
  bc.df <- data.frame(mid=h$mids, count=h$counts/n)
  h <- hist(a-c, breaks=brks, plot=F)
  ac.df <- data.frame(mid=h$mids, count=h$counts/n)

  clr[is.na(clr)] <- median(clr, na.rm=T)
  md <<- median(clr, na.rm=T)
  print(summary(clr))

  d <- data.frame(xx=x, yy=y, a=a, b=b, c=c, clr=clr)
  p <- ggplot(d, aes(x=xx, y=yy))
  p <- p + theme_bw()

  if (!is.na(clr[1])) {
    p <- p + geom_point(size=0.5, aes(colour=clr))
    p <- p + scale_colour_gradient()
  } else {
    p <- p + geom_point(size=0.5)
  }

  cf <<- cos(pi/3)
  sf <<- sin(pi/3)
  cff <<- cos(pi/3)
  sff <<- sin(pi/3)
  p <- p + geom_segment(data=ac.df,
    aes(x=.5*cf*(mid+1), y=.5*sf*(mid+1),
        xend=(.5*cf*(mid+1)) - count*sff,
        yend=(.5*sf*(mid+1)) + count*cff),
    colour='black',
    size=2
  )
  p <- p + geom_segment(data=ab.df,
    aes(x=1-(.5*cf*(mid+1)), y=.5*sf*(mid+1),
        xend=1-(.5*cf*(mid+1)) + count*sff,
        yend=(.5*sf*(mid+1)) + count*cff),
    colour='black',
    size=2
  )
  p <- p + geom_segment(data=bc.df,
    aes(y=0, x=.5*(mid+1),
        xend=.5*(mid+1),
        yend=-count*cff),
    colour='black',
    size=2
  )

  p <- p + scale_x_continuous(limits=c(-0.3, 1.3))
  p <- p + scale_y_continuous(limits=c(-sin(pi/3)*.3, sin(pi/3)*1.3))

  p <- p + opts(
    legend.position='none'
  )

  pdf(file="ternary.pdf")
  print(p)
  dev.off()

}