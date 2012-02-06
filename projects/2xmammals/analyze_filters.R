library(plyr)
library(ggplot2)
uname  <- Sys.getenv("USER")
if (uname == 'gj1') {
  source("~/src/greg-ensembl/scripts/mysql_functions.R")
  source("~/src/greg-ensembl/scripts/xtable_utils.R")
} else {
  source("~/lib/greg-ensembl/scripts/mysql_functions.R")
  source("~/lib/greg-ensembl/scripts/xtable_utils.R")
}

dbname <- function() {
   'gj1_2x_63_alt'
}

plot.phred.filtered <- function() {
  con <- connect(dbname())
  df <- dbGetQuery(con, 'select * from qual_filt')
  dbDisconnect(con)

  df <- ddply(df, .(taxon_id), function(x) {
    total.nucs <- sum(x$length_nucs)
    filtered.nucs <- sum(x$filtered_nucs)
    data.frame(
      total = total.nucs,
      filtered = filtered.nucs
    )
  })
    
  df <- subset(df, filtered > 0)
  df$name <- taxid.to.alias(df$taxon_id)
  
  df <- base::within(df, {
    name <- reorder(name, filtered / total)
  })

  p <- ggplot(df, aes(x=name, y=filtered / total))
  p <- p + theme_bw()
  p <- p + geom_bar()
  p <- p + scale_x_discrete("Species")
  p <- p + scale_y_continuous("Fraction of Nucleotides Filtered")
  p <- p + opts(
    axis.text.x = theme_text(angle=90, hjust=1)
  )
  
  pdf(file="filtered_qual_bars.pdf", width=6, height=4)
  print(p)
  dev.off()
}

plot.paralog.histogram <- function() {
  con <- connect(dbname())
  df <- dbGetQuery(con, 'select * from paralogs')
  dbDisconnect(con)

  # Show an overall distribution of the length ratio. Combine all removed + kept genes.
  # Get the unique rows for the kept gene value.
  df$value <- df$removed_len_ratio
  df$distval <- df$removed_dist
  df$typ <- paste('Discarded', ' (n=', nrow(df), ')', sep='')
  print(nrow(df))
  df2 <- df
  df2 <- df2[!duplicated(df2[, c('kept_id')]),]
  df2$value <- df2$kept_len_ratio
  df2$distval <- df2$kept_dist
  df2$typ <- paste('Kept', ' (n=', nrow(df2), ')', sep='')

  df.comb <- rbind(df, df2)
  df.comb$typ <- paste('All', ' (n=', nrow(df.comb), ')', sep='')

  df.split <- rbind(df, df2, df.comb)
  df.split$value <- pmin(2, df.split$value)

  p <- ggplot(df.split, aes(x=value))
  p <- p + theme_bw()
  p <- p + geom_histogram(binwidth=0.025)
  p <- p + scale_x_continuous("Relative Gene Length")
  p <- p + scale_y_continuous("Count")
  p <- p + facet_grid(typ ~ ., scales="free_y")
  p <- p + opts(
    strip.text.y = theme_text(angle=0, hjust=0.5)
  )

  svg(file="filtered_paralogs_hist.svg", height=4, width=6)
  print(p)
  dev.off()

  short.paralogs <- df$removed_len_ratio < 0.5 & df$kept_len_ratio > 0.5
  equal.distance <- !short.paralogs & df$kept_dist == df$removed_dist
  r.a <- 'Higher Divergence'
  r.b <- 'Length < 0.5'
  r.c <- 'Shorter Length'
  df[, 'reason'] <- r.a
  df[short.paralogs, 'reason'] <- r.b
  df[equal.distance, 'reason'] <- r.c

  n.dist <- nrow(subset(df, reason == r.a))
  n.length.short <- nrow(subset(df, reason == r.b))
  n.length.divergence <- nrow(subset(df, reason == r.c))
  df[df$reason == r.a, 'reason.str'] <- 1
  str.1 <- paste(r.a, '\n(n=',n.dist, ')', sep='')
  df[df$reason == r.b, 'reason.str'] <- 3
  str.3 <- paste(r.b, '\n(n=',n.length.short, ')', sep='')
  df[df$reason == r.c, 'reason.str'] <- 2
  str.2 <- paste(r.c, '\n(n=',n.length.divergence, ')', sep='')

  df$reason.str <- factor(df$reason.str, levels=c(1,2,3), labels=c(str.1, str.2, str.3))

  df$kept_dist <- pmin(df$kept_dist, 2)
  df$removed_dist <- pmin(df$removed_dist, 2)

  xyz <- df[df$reason == r.b,]
  xyz <- subset(xyz, kept_dist > 0 & removed_dist > 0)
  lm.res <- lm(removed_dist ~ kept_dist + 0, data=xyz)
  print(summary(lm.res))
  lm.res <- lm(kept_dist ~ removed_dist + 0, data=xyz)
  print(summary(lm.res))

  r.f <- function(x, a, b) {
    x$kept_dist <- round_any(x$kept_dist, a * 1)
    x$removed_dist <- round_any(x$removed_dist, b * 1)
    x
  }
  df <- r.f(df, 0.05, 0.05)
  df <- subset(df, kept_dist >= 0 & removed_dist >= 0)
  p <- ggplot(df, aes(x=kept_dist, y=removed_dist))
  p <- p + theme_bw()

  p <- p + geom_point(stat='sum', aes(colour=..prop..), size=.8)
  p <- p + scale_colour_gradientn("Density", colour=c('white', rgb(0.5, 0.5, 1), 'red'), trans='log10')

  p <- p + scale_x_continuous("Kept Divergence (subst. / site)", limits=c(0, 2))
  p <- p + scale_y_continuous("Discarded Divergence (subst. / site)", limits=c(0, 2))
  p <- p + facet_grid(. ~ reason.str)
  p <- p + coord_equal()
  p <- p + opts(
    strip.text.y = theme_text(angle=0, hjust=0.5),
    axis.text.x = theme_text(angle=90, hjust=1),
    strip.background = theme_blank()
  )

  svg(file="filtered_paralogs_scatter.svg", width=7, height=4)
  print(p)
  dev.off()
}

subs.sim <- function(
  win=15,
  n = 10 * 1000,
  bl = seq(from=0.01, to=1, length.out=10)
  ) {
  # Given dS branch length l, simulate substitution counts within N-codon windows.
  sub.df <- data.frame()

  # Create a vector of branch lengths sampled from the given distribution
  # with each branch length repeated 'win' times.
  bl.sample <- sample(bl, size=n, replace=T)
  bl.rep <- sapply(bl.sample, function(x) rep(x, win))
  bl.rep <- as.vector(bl.rep)
  #bl.rep <- bl.rep + rnorm(n=length(bl.rep), mean=0, sd=bl.rep/2)
  bl.rep <- pmax(bl.rep, 0.005)

  #f.inv <- 0.1
  #inv.sites <- sample.int(length(bl.rep), length(bl.rep) * f.inv, replace=F)
  #bl.rep[inv.sites] <- 0

  #syn.subs <- rpois(n * win, lambda=bl.rep)
  syn.subs <- rnbinom(n*win, size=.1, mu=bl.rep * .3)
  syn.subs <- ifelse(syn.subs >= 1, 1, 0)

  nsyn.subs <- rnbinom(n*win, size=.1, mu=bl.rep * .7)
  nsyn.subs <- ifelse(nsyn.subs >= 1, 1, 0)

  syn.mat <- matrix(syn.subs, nrow=n, ncol=win, byrow=T)
  syn.counts <- apply(syn.mat, 1, sum)
  nsyn.mat <- matrix(nsyn.subs, nrow=n, ncol=win, byrow=T)
  nsyn.counts <- apply(nsyn.mat, 1, sum)

  #f.inv.windows <- 0.1
  #inv.windows <- sample.int(length(syn.counts), length(syn.counts) * f.inv.windows, replace=F)
  #syn.counts[inv.windows] <- 0
  #nsyn.counts[inv.windows] <- 0

  print(table(syn.counts))
  print(table(nsyn.counts))

  out.df <- data.frame(
    taxon_id = 9999,
    is_leaf = 1,
    n_s_subs = syn.counts,
    n_ns_subs = nsyn.counts,
    n_nongap_codons = 15,
    orig_node_bl = bl.sample,
    node_bl_se = 0,
    data_id = 0,
    aln_position = 0
  )
  out.df  
}

pois.sim <- function(x, win.size) {
  q.val <- quantile(x, 0.99)
  zz <- x[x <= q.val]

  method <- 'simple'
  if (method == 'complex') {
    # Unpack the window counts into a series of 1s and 0s
    subs <- zz
    non.subs <- win.size - zz
    xx <- c(rep(1, times=sum(subs)), rep(0, times=sum(non.subs)))

    print(table(xx))
    ft.nb <- fitdist(xx, 'nbinom')
    sim.vals <- rnbinom(length(xx), size=ft.nb$estimate[1], mu=ft.nb$estimate[2])

    syn.mat <- matrix(sim.vals, nrow=length(zz), ncol=win.size, byrow=T)
    syn.counts <- apply(syn.mat, 1, sum)
  }

  if (method == 'simple') {
    f.zero <- sum(x == 0) / length(x)
    lambda <- -log(f.zero)
    syn.counts <- rpois(length(x), lambda=lambda)
  }
  return(syn.counts)
}

plot.subs.hist <- function(win.size=15, clean=F, include.sim=F, include.pois=F) {
  df.lbl <- paste('sub.windows', win.size, 'df', sep='.')
  if (!exists(df.lbl, envir=.GlobalEnv) || clean) {
    con <- connect(dbname())
    cmd <- sprintf("
      select
        data_id, aln_position, is_leaf, taxon_id, orig_node_bl, node_bl_se, 
          n_ns_subs, n_s_subs, n_nongap_codons
        from windows_%s
        where taxon_id in (9606, 9604, 9598, 9593, 9443, 9258, 42254, 10020, 9483, 61853, 10090, 39107, 314145)
    ", win.size)
    df <- dbGetQuery(con, cmd)
    dbDisconnect(con)
    assign(df.lbl, df, envir=.GlobalEnv)
  }
  cur.df <- get(df.lbl, envir=.GlobalEnv)
  cur.df$sim <- FALSE

  if (include.sim) {
    sim.lbl <- 'sim.df'
    if (!exists(sim.lbl, envir=.GlobalEnv)) {
      print("  simulating...")
      ### Simulate some sites using a poisson process.
      sim.df <- ddply(cur.df, .(taxon_id), function(x) {
        print(x[1, 'taxon_id'])
        bls <- x$orig_node_bl
        bls <- pmin(bls, quantile(bls, 0.9))
        bls <- pmax(bls, quantile(bls, 0.1))
        simulated.df <- subs.sim(n=nrow(x), bl=bls)
        simulated.df$taxon_id <- x[1, 'taxon_id']
        simulated.df$sim <- TRUE  
        simulated.df
      })
      assign(sim.lbl, sim.df, envir=.GlobalEnv)
      print("  done!")
    }
    sim.df <- get(sim.lbl, envir=.GlobalEnv)
    cur.df <- rbind(cur.df, sim.df)  
  }

  unq.tx <- unique(cur.df$taxon_id)
  cur.df$taxon_id <- factor(cur.df$taxon_id, levels=unq.tx, labels=taxid.to.alias(unq.tx, include.internals=T))
  cur.df <- base::within(cur.df, {
    taxon_id <- reorder(taxon_id, orig_node_bl, mean)
  })

  tx.lvls <- dlply(cur.df, .(taxon_id), function(x) {
    cur.tx <- as.character(x[1, 'taxon_id'])
    sprintf(" %s (%.3f)", cur.tx, mean(x$orig_node_bl))
  })  
  cur.df$taxon_id <- factor(cur.df$taxon_id, labels=unlist(tx.lvls))

  if (include.pois) {
    sim.lbl <- paste('sim.windows', win.size, 'df', sep='.')
    clean <- TRUE
    if (!exists(sim.lbl, envir=.GlobalEnv) || clean) {
      print("  simulating...")
      ### Simulate some sites using a poisson process.
      sim.df <- ddply(cur.df, .(taxon_id), function(x) {
        sim.vals <- pois.sim(x$n_s_subs, win.size)
        simulated.df <- data.frame(
          value = sim.vals,
          taxon_id = x[1, 'taxon_id'],
          fld = 'Simulated Synonymous',
          yval = 0
        )
        simulated.df
      })
      assign(sim.lbl, sim.df, envir=.GlobalEnv)
      print("  done!")
    }
    sim.df <- get(sim.lbl, envir=.GlobalEnv)
  }

  med.df <- ddply(cur.df, .(taxon_id), function(x) {
    keep.count <- floor(nrow(x) * 0.001)
    top.s <- x[order(x$n_s_subs, decreasing=T),]
    top.s <- top.s[keep.count, ]
    top.n.s <- top.s$n_s_subs

    top.ns <- x[order(x$n_ns_subs, decreasing=T),]
    top.ns <- top.ns[keep.count, ]
    top.n.ns <- top.ns$n_ns_subs

#    n.ns.cor <- cor(x$n_ns_subs, x$n_s_subs, method='s')

     ns.ovr.s <- function(y) {sum(y$n_ns_subs) / sum(y$n_s_subs)}
     ns.s <- ns.ovr.s(x)

    data.frame(
      is_leaf = x[1, 'is_leaf'],
      n_win = nrow(x),
      n_genes = length(unique(x$data_id)),
      bl_mean = mean(x$orig_node_bl),
      bl_lo = quantile(x$orig_node_bl, 0.25),
      bl_hi = quantile(x$orig_node_bl, 0.75),
      avg_per_codon = mean(x$n_ns_subs + x$n_s_subs) / win.size,
      ns_ovr_s = ns.s,
      ns_sum = sum(x$n_ns_subs),
      s_sum = sum(x$n_s_subs),
      ns_top = top.n.ns,
      s_top = top.n.s
    )
  })

  med.df <- med.df[order(med.df$is_leaf, med.df$bl_mean),]
  print(med.df)

  xt <- xtable(med.df)
  xt <- color.columns(xt, c('ns_top', 's_top'))
  print.latex(xt, "sub_hist_summary.txt")

  ns.df <- cur.df
  ns.df$fld <- 'Nonsynonymous'
  ns.df$value <- ns.df$n_ns_subs

  s.df <- cur.df
  s.df$fld <- 'Synonymous'
  s.df$value <- s.df$n_s_subs

  comb.df <- rbind(ns.df, s.df)
  comb.df <- subset(comb.df, sim==FALSE)

  # Create an extra dataframe for using med.df and the 's_top' and 'ns_top' values
  # to plot triangle marks in the histogram.
  med.df$yval <- 10^3.5
  mark.ns <- med.df
  mark.ns$value <- mark.ns$ns_top - 0.25
  mark.ns$fld <- 'Nonsynonymous'
  mark.s <- med.df
  mark.s$yval <- 10^3
  mark.s$value <- mark.s$s_top - 0.25
  mark.s$fld <- 'Synonymous'
  marks <- rbind(mark.ns, mark.s)
  marks <- subset(marks, select=c('yval', 'value', 'fld', 'taxon_id'))

  print(str(marks))

  comb.df$yval <- 0
  comb.df <- subset(comb.df, select=c('yval', 'value', 'fld', 'taxon_id'))
  
  if (include.pois) {
    print(str(sim.df))
    comb.df <- rbind(comb.df, sim.df)
  }

  print(str(comb.df))

  p <- ggplot(comb.df, aes(x=value, y=..count..+1, fill=fld))
  p <- p + theme_bw()
  p <- p + geom_histogram(binwidth=1, position="dodge")
  
  p <- p + geom_point(data=marks, aes(y=yval, x=value), shape=25, colour=NA, size=2)

  #p <- p + scale_colour_discrete("Substitution Type")
  p <- p + scale_fill_discrete("Substitution Type")
  p <- p + scale_y_log10("Count (Log Scale)")
  p <- p + scale_x_continuous(sprintf("Substitutions (Non-Overlapping %d-Codon Windows)", win.size), limits=c(-1, win.size+3))
  p <- p + facet_grid(taxon_id ~ ., scales="free_y")
  p <- p + opts(
    strip.text.y = theme_text(angle=0, hjust=0, size=10),
    strip.background = theme_blank()
  )
  
  pdf(file=paste("subs_hist_combined_", win.size, ".pdf", sep=''), width=10, height=10)
  print.ggplot(p)
  dev.off()

}

plot.subs.histogram <- function() {
  if (!exists('subs.df', envir=.GlobalEnv)) {
    con <- connect(dbname())
    subs.df <- dbGetQuery(con, '
      select
        taxon_id, orig_node_bl, n_ns_subs, n_s_subs, n_nongap_codons 
        from sub_windows
        where taxon_id in (9606, 9598, 9593, 13616, 30611, 61853)
        limit 300000
    ')
    dbDisconnect(con)
    assign('subs.df', subs.df, envir=.GlobalEnv)
  }
  df <- get('subs.df', envir=.GlobalEnv)
  print("  calculating plots...")

  cum.p <- function(xx, field, rng) {
    print(paste("  ",field,sep=''))
    xx$value <- xx[, field]
    if (rng[2] == 30) {
      bins <- 1
    } else {
      bins <- 0.25
    }
    p <- ggplot(xx, aes(x=value, y=..count..+1))
    p <- p + theme_bw()
    p <- p + geom_histogram(binwidth=bins)
    p <- p + scale_x_continuous(limits=rng)
    p <- p + scale_y_log10()
    p <- p + facet_grid(typ ~ ., scales="free_y")
    p <- p + opts(
      title=field,
      strip.text.y = theme_text(angle=0, hjust=0, size=8),
      axis.text.y = theme_blank()      
    )
    p
  }

  df <- subset(df, n_ns_subs >= 3)

  df$typ <- 'Other Terminal'
  df[is.na(df$taxon_id), 'typ'] <- 'Ancestral Node'
  df[is.na(df$taxon_id), 'taxon_id'] <- 0

  df[df$typ == 'Ancestral Node' & df$orig_node_bl < 0.3 & df$orig_node_bl > 0.1, 'typ'] <- 'Ancestral Short'
  df[df$typ == 'Ancestral Node' & df$orig_node_bl > 0.6, 'typ'] <- 'Ancestral Long'

  df[df$taxon_id > 0, 'typ'] <- 'Other Terminal Node'
  taxid.lbls <- taxid.to.alias(df$taxon_id)
  df[!is.na(taxid.lbls), 'typ'] <- taxid.lbls[!is.na(taxid.lbls)]
  df$typ <- factor(df$typ)

  print(levels(df$typ))

  df$node_bl <- df$orig_node_bl
  df$node_bl <- pmax(0.02, df$orig_node_bl)
  df <- base::within(df, {
    typ <- reorder(typ, orig_node_bl, FUN=median)
  })

  typ.lvls <- dlply(df, .(typ), function(x) {
    cur.typ <- as.character(x[1, 'typ'])
    sprintf("%s (%.3f)", cur.typ, median(x$orig_node_bl))
  })  
  df$typ <- factor(df$typ, labels=unlist(typ.lvls))

  df <- ddply(df, .(typ), function(x) {
    med.bl <- median(x$orig_node_bl)
    x$median_bl <- med.bl
    x
  })
  df$s_per_codon <- df$n_s_subs / df$n_nongap_codons
  df$ns_per_codon <- df$n_ns_subs / df$n_nongap_codons
  df$s_per_codon_bl <- df$n_s_subs / df$n_nongap_codons / df$node_bl
  df$ns_per_codon_bl <- df$n_ns_subs / df$n_nongap_codons / df$node_bl

#  s.p <- cum.p(df, 'n_s_subs', c(0,30))
  s.bl.p <- cum.p(df, 's_per_codon_bl', c(0, 15))
#  ns.p <- cum.p(df, 'n_ns_subs', c(0,30))
  ns.bl.p <- cum.p(df, 'ns_per_codon_bl', c(0,15))

  pdf(file="filtered_subs_hist.pdf", width=8, height=10)
  vplayout(2, 1)
  print(s.bl.p, vp=subplot(1,1))
  print(ns.bl.p, vp=subplot(2,1))
  dev.off()
}

clustered.subs.sim <- function() {

  n <- 100000
  win <- 15
  mtx <- matrix(nrow=n, ncol=win)

  for (i in 1:win) {
    mtx[1:n, i] <- rbinom(n, 1, prob=0.01)
  }

  sums <- apply(mtx, 1, sum)
  table(sums)
}

slim.mammals <- function() {
  c(9986, 9365, 9258, 9315, 9361, 9785, 9606, 9598, 9593, 10090, 10116, 9785,
  9478, 9615, 30611, 9371, 9358, 9796, 9913)
}

mammals <- function() {
  c(9258, 9315, 9358, 9361, 9365, 9371, 9478,
  9483, 9544, 9593, 9598, 9601, 9606, 9615, 9646, 9685, 9739, 9785,
  9796, 9813, 9823, 9913, 9978, 9986, 10020, 10090, 10116, 10141,
  13616, 30538, 30608, 30611, 37347, 42254, 43179, 59463, 61853,
  132908)
}

primates <- function() {
  c(9478, 9483, 9544, 9593, 9598, 9601, 9606, 30608, 30611, 61853)
}

all.primates <- function() {
  c(primates(), 207598, 9604, 314295, 9526, 314293, 376913, 376911, 9443, 37347)
}

eutheria <- function() {
  c(9358, 9361, 9365, 9371, 9478, 9483, 9544, 9593, 9598, 9601, 9606,
  9615, 9646, 9685, 9739, 9785, 9796, 9813, 9823, 9913, 9978, 9986,
  10020, 10090, 10116, 10141, 30538, 30608, 30611, 37347, 42254,
  43179, 59463, 61853, 132908)
}

taxid.to.alias2 <- function(taxids, map.df) {
  match.inds <- match(taxids, map.df$taxon_id)
  map.df$name[match.inds]
}

get.taxid.df <- function() {
  spec.tree <- get.species.tree()
  taxids <- c(spec.tree$tip.label, spec.tree$node.label)

  name_class <- 'ensembl alias'

  unique.ids <- unique(taxids)
  unique.ids <- unique.ids[unique.ids != '']
  taxid_list <- paste(unique.ids, collapse=', ')
  cmd <- sprintf('select n.taxon_id, n.name from ncbi_taxa_name n where
    n.taxon_id IN (%s) and n.name_class like "%%%s%%" order by n.name',
      taxid_list, name_class)
  con <- connect.livemirror('ensembl_compara_64')
  df <- dbGetQuery(con, cmd)
  dbDisconnect(con)

  df  
}

taxid.to.alias <- function(taxids, binomials=F, include.internals=F, save.df=T) {
  df.key <- paste('taxid', as.character(binomials), as.character(include.internals), sep='.')
  if (!exists(df.key, envir=.GlobalEnv) || include.internals || !save.df) {
    con <- connect.livemirror('ensembl_compara_64')
    name_class <- "ensembl alias"
    if (binomials) {
      name_class <- 'scientific'
    }

    if (include.internals) {
      unique.ids <- unique(taxids)
      unique.ids <- unique.ids[unique.ids != '']
      taxid_list <- paste(unique.ids, collapse=', ')
      cmd <- sprintf('select n.taxon_id, n.name from ncbi_taxa_name n where
        n.taxon_id IN (%s) and n.name_class like "%%%s%%" order by n.name',
          taxid_list, name_class)
    } else {
      cmd <- sprintf('select n.taxon_id, n.name from genome_db g, ncbi_taxa_name n
        where n.taxon_id=g.taxon_id and n.name_class like "%%%s%%"
        order by n.name', name_class)
    }
    df <- dbGetQuery(con, cmd)
    dbDisconnect(con)
    if (any(df$name == 'Lesser hedgehog tenrec')) {
      df[df$name == 'Lesser hedgehog tenrec', 'name'] <- 'Tenrec'
    }
    if (save.df) {
      assign(df.key, df, envir=.GlobalEnv)
    }
  }
  if (save.df) {
    df <- get(df.key, envir=.GlobalEnv)
  }

  #print(df)

  if (missing(taxids)) {
    print("  Missing some taxids!")
    print(df)
  }

  match.inds <- match(taxids, df$taxon_id)
  df$name[match.inds]
}

subplot <- function(x, y) viewport(layout.pos.col=x, layout.pos.row=y)
vplayout <- function(x, y) {
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(y,x)))
}
