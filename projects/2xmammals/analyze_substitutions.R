bsub.nsyn.plot <- function() {
  bsub.pset.function('nsyn_plot', mem=14, extra.args='FALSE', psets=c(1, 2, 3, 4))
}

plot.taxon.id.trees <- function() {
  library(phylosim)
  psets <- pset.df()$pset_id
  print(psets)

  pdf(file=scratch.f("pset_trees.pdf"))
  for (pset in psets) {
    for (translate.labels in c(F, T)) {
      nm <- pset.to.alias(pset)
      tree <- get.subset.taxid.tree(pset)
      if (translate.labels) {
        tree <- tree.taxids.to.labels(tree)
      }
      psim <- PhyloSim()
      psim$.phylo <- tree
      print(paste("Plotting", nm))
      p <- plotTree(psim, 
        tree.do.plot = F,
        line.color = 'gray',
        tree.xlim.expand = c(0.5, 0.5),
        include.internal.labels=T,
        internal.size = 0.9,
        internal.angle = -20,
        internal.color = 'black'
      )
      p <- p$grob
      p <- p + opts(title=nm)
      print.ggplot(p)
      #grid.newpage()
    }
  }
  dev.off()
}

test.nsyn.plot <- function() {
  nsyn.plot(pset.id=1, test=F, clean=F)
}

nsyn.plot <- function(pset.id=1, test=F, clean=F) {
  alias <- pset.to.alias(pset.id, factors=F)
  print(alias)

  tree <- get.subset.taxid.tree(pset.id)
  all.taxids <- tree.all.labels(tree)

  print(all.taxids)  
  taxids.str <- paste('(', paste(all.taxids, collapse=', '), ')', collapse='')
  print(taxids.str)

  analysis.dir <- scratch.f('nsyn_plot')
  dir.create(analysis.dir, showWarnings=F)

  subs.out.f <- file.path(analysis.dir, paste('subs_', alias, '.Rdata', sep=''))
  print(subs.out.f)
  if (!file.exists(subs.out.f) || clean) {
    # Get the subs from MySQL.
    print("Getting subs from mysql...")
    limit.str <- ifelse(test, 'LIMIT 500000', '')
    cmd <- sprintf("select data_id, aln_pos, taxon_id, mut_nsyn, confidence, mut_cpg, mut_rev_cpg, codon_cpg, mut_ws, 
      aa_from, aa_to, nuc_from, nuc_to, codon_from, codon_to from subs where
      aln_pos IS NOT NULL and data_id IS NOT NULL and taxon_id IS NOT NULL
      and taxon_id IN %s %s",
      taxids.str, limit.str)
    con <- connect(db())
    subs <- dbGetQuery(con, cmd)
    disconnect(con)
    print(nrow(subs))
    print("Saving subs file...")
    save(subs, file=subs.out.f)
  } else {
    load(subs.out.f)
  }
  subs$aln_position <- subs$aln_pos
  subs <- subset(subs, confidence > 0.9)

  sites.out.f <- file.path(analysis.dir, paste('sites_', alias, '.Rdata', sep=''))
  sites.pos.f <- file.path(analysis.dir, paste('pos_sites_', alias, '.Rdata', sep=''))
  if (!file.exists(sites.pos.f) || clean) {
    if (!file.exists(sites.out.f) || clean) {
      print("loading sites...")
      sites <- get.pset.sites(pset.id, filter='default', test=F)
      if (test) {
        sites <- sites[1:500000,]
      }
      print(nrow(sites))
      print("Saving sites file...")
      save(sites, file=sites.out.f)
    } else {    
      print(sites.out.f)
      load(sites.out.f)
    }
    pos.sites <- subset(sites, pos.pval < 0.05)
      print("Saving positive sites file...")
      save(pos.sites, file=sites.pos.f)
  } else {
    print(sites.pos.f)
    load(sites.pos.f) 
  }
  print(sprintf("Positive sites: %d", 
    nrow(pos.sites)
  ))

  #pos.sites <- pos.sites[1:100,]

  print("merging...")
  subs.sites <- merge(pos.sites, subs, by=c('data_id', 'aln_position'), all.x=F)
  subs.sites$taxon_alias <- taxid.to.alias(subs.sites$taxon_id, include.internals=T)
  subs.sites$site_id <- paste(subs.sites$data_id, subs.sites$aln_position, sep='.')

  print(sprintf("Site IDs: %d", length(unique(subs.sites$site_id))))

  #print(sort(table(subs.sites$aa_from)))
  #print(sort(table(subs.sites$aa_to)))
  nsyn.sites <- subset(subs.sites, mut_nsyn==1)
  nsyn.nocpg <- subset(subs.sites, mut_nsyn==1 & mut_cpg==0 & mut_rev_cpg==0)
  df.list <- list(
    all.nsyn=nsyn.sites,
    no.cpg=nsyn.nocpg
  )
  
  for (lbl in c('all.nsyn', 'no.cpg')) {
    xx <- df.list[[lbl]]
    print(sprintf("%s (subs=%d, sites=%d)", lbl, nrow(xx), length(unique(xx$site_id))))
    print(tail(sort(table(paste(xx$aa_from, xx$aa_to, sep='->'))), n=10))
    print(tail(sort(table(paste(xx$codon_from, xx$codon_to, sep='->'))), n=10))
    print(tail(sort(table(paste(xx$nuc_from, xx$nuc_to, sep='->'))), n=10))
  }

  print("ddplying...")
  out.df <- ddply(subs.sites, .(taxon_alias), function(x) {
    n.nsyn <- nrow(subset(x, mut_nsyn == 0))
    n.syn <- nrow(subset(x, mut_nsyn == 1))
    nsyn.adj <- n.nsyn / (n.syn + 1)

    data.frame(
      taxon_id = x[1, 'taxon_id'],
      nsyn = n.nsyn,
      syn = n.syn,
      nsyn.adj = nsyn.adj
    )
  })
  print(out.df)

  out.df <- ddply(subs.sites, .(site_id), function(x) {
    x <- subset(x, mut_nsyn == 1)
    nsyn.taxa <- paste(sort(unique(x$taxon_alias)), collapse=', ')
    n.nsyn.taxa <- length(unique(x$taxon_id))
    data.frame(
      nsyn.taxa = nsyn.taxa,
      n.nsyn.taxa = n.nsyn.taxa,
      stringsAsFactors=F
    )
  })
  print(tail(sort(table(out.df$nsyn.taxa)), n=20))
  print(table(out.df$n.nsyn.taxa))

  xx <- subset(subs.sites, mut_nsyn == 1)
  unique.positions <- as.factor(xx$site_id)
  nsyn.sub.counts <- table(unique.positions)
  print(table(nsyn.sub.counts))

  #results <- out.df
  #sites.out.f <- file.path(analysis.dir, paste('results_', alias, '.Rdata', sep=''))
  #save(results, file=sites.out.f)

}