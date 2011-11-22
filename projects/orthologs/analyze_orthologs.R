library(plyr)
library(ggplot2)
library(xtable)
uname  <- Sys.getenv("USER")
if (uname == 'gj1') {
  source("~/src/greg-ensembl/scripts/mysql_functions.R")
  source("~/src/greg-ensembl/scripts/xtable_utils.R")
} else {
  source("~/lib/greg-ensembl/scripts/mysql_functions.R")
  source("~/lib/greg-ensembl/scripts/xtable_utils.R")
}

do.all <- function() {
  plot.hist()
  plot.dups()
  plot.eutheria.dups()
  print.some.numbers()
  summarize.nodesets()
}

print.some.numbers <- function() {
  con <- connect('gj1_orthologs_63');
  df <- dbGetQuery(con, 'select * from trees')
  dbDisconnect(con)
  
  root.trees <- subset(df, method_id == "Ensembl Roots")
  print(nrow(root.trees))
  print(sum(root.trees$human_count))
  small.trees <- subset(root.trees, leaf_count <= 15)
  large.trees <- subset(root.trees, leaf_count > 15)
  print(sum(small.trees$human_count) / sum(root.trees$human_count))

  fish <- subset(df, method_id == "Fish")
  fish.dups <- subset(fish, human_count == 0  & species_count <= 5 & leaf_count <= 40)
  print(nrow(fish.dups))

  print(nrow(fish.dups) / (19991 * 2))

  hum <- subset(df, method_id == "Human Orthologs")
  mouse.zero <- nrow(subset(hum, mouse_count == 0)) / nrow(hum)
  mouse.one <- nrow(subset(hum, mouse_count == 1)) / nrow(hum)
  mouse.two <- nrow(subset(hum, mouse_count >= 2)) / nrow(hum)

  print(paste(mouse.zero, mouse.one, mouse.two))
}

summarize.set <- function(x) {
    method.id <- x[1, 'method_id']

    # Sort trees from big to small.
    sorted.counts <- x[order(x$leaf_count, decreasing=TRUE), 'leaf_count']
    sorted.sum <- cumsum(sorted.counts)
    # Find the first tree where the cumsum is >50%
    half.index <- min(which(sorted.sum > sum(sorted.counts) / 2))
    n.fifty <- sorted.counts[half.index]

    n.total <- sum(x$leaf_count)
    n.total <- sprintf("%.1e", n.total)

    median.n <- median(x$leaf_count)
    mean.n <- mean(x$leaf_count)
    min.n <- min(x$leaf_count)
    max.n <- max(x$leaf_count)

    n.zero.h <- nrow(subset(x, human_count == 0))
    n.one.h <- nrow(subset(x, human_count == 1))
    n.two.h <- nrow(subset(x, human_count >= 2))

    total.h <- sum(x$human_count)

    mean.mpl <- mean(x$tree_mpl)
    median.mpl <- median(x$tree_mpl)
    sd.mpl <- sd(x$tree_mpl)

    mean.species <- mean(x$species_count)
    tbl <- table(x$species_count)
    max.i <- which.max(tbl)
    mode.species <- dimnames(tbl)[[1]][max.i]
    median.species <- median(x$species_count)

    n <- nrow(x)
    data.frame(
      'Category' = '',
      method_id = method.id,
      'Count' = n,
      'Size' = paste(as.integer(median.n), ' (', min.n, ' / ', max.n, ')', sep=''),
      'Sequence Count' = n.total,
      'N50' = n.fifty,
      
      'Zero' = n.zero.h / n,
      'One' = n.one.h / n,
      'TwoPlus' = n.two.h / n,
      'Total' = total.h,
      
      'MPL' = median.mpl,
      'Species' = as.integer(median.species),
      stringsAsFactors=F
    )
}


summarize.nodesets <- function() {
  con <- connect('gj1_orthologs_63');
  df <- dbGetQuery(con, 'select * from trees')
  dbDisconnect(con)

  methds <- factor(nodeset.levels(), levels=nodeset.levels())
  df$method_id <- factor(df$method_id, levels=nodeset.levels())

  df.summary <- ddply(df, .(method_id), summarize.set)

  # Get the comparison to the Optic database.
  source("compare_to_optic.R")
  optic.summary <- get.optic.summary()
  df.summary <- rbind(df.summary, optic.summary)
  df.summary$method_id <- as.character(df.summary$method_id)

  df.summary[df.summary$method_id == 'Ensembl Roots', 'method_id'] <- 'Ensembl Trees'
  df.summary[df.summary$method_id == 'Optic', 'method_id'] <- 'Optic Trees'
  df.summary[df.summary$method_id == 'MammalSubgroupsPlusOutgroup', 'method_id'] <- '\\footnotesize{MammalSubgroupsPlusOutgroup}'

  # Label ingroups, outgroups, subgroups, orthologs, Compara Trees, Optic Trees
  df.summary[df.summary$method_id == 'Primates', 'Category'] <- 'Ingroups'
  df.summary[df.summary$method_id == 'Eutheria', 'Category'] <- 'Outgroups'
  df.summary[df.summary$method_id == 'MammalSubgroups', 'Category'] <- 'Subgroups'
  df.summary[df.summary$method_id == 'Human Orthologs', 'Category'] <- 'Orthologs'
  df.summary[df.summary$method_id == 'Ensembl Trees', 'Category'] <- 'Default Trees'

  print(df.summary)

  xt <- xtable(df.summary)
  cc.f <- function(x, column) {
    color.column(x, column, low='white', high=rgb(0.3, 0.3, 1))
  }
  cc.all <- function(x) {
    x <- cc.f(x, 'N50')
    x <- cc.f(x, 'Zero')
    x <- cc.f(x, 'One')
    x <- cc.f(x, 'TwoPlus')
    x <- cc.f(x, 'MPL')
    x <- cc.f(x, 'Total')
    x <- cc.f(x, 'Species')
    x <- cc.f(x, 'Count')
  }
  xt <- cc.all(xt)
  print.latex(xt, "ortholog_summary.txt")
  assign("xt", xt, envir=.GlobalEnv)

  root.trees <- subset(df, method_id == "Ensembl Roots")
  small.trees <- subset(root.trees, leaf_count <= 15)
  large.trees <- subset(root.trees, leaf_count > 15)
  a <- summarize.set(root.trees)
  a$method_id <- "Ensembl Roots"
  b <- summarize.set(small.trees)
  b$method_id <- "($<=$15)"
  d <- summarize.set(large.trees)
  d$method_id <- "($>$15)"
  summaries <- rbind(a,b,d)
  summaries$Category <- NULL
  xt <- xtable(summaries)
  print(summaries)
  print(xt)
  print.latex(xt, "ortholog_roots_summary.txt")
}

print.latex <- function(xt, filename) {
  print(xt, 
    file=filename,
    only.contents=T,
    include.colnames=F,
    include.rownames=F,
    sanitize.text.function=function(x){x},
    hline.after=NULL
  )
}

root.nodes.hist <- function() {
  con <- connect('gj1_orthologs_63');
  df <- dbGetQuery(con, 'select * from trees where method_id="Ensembl Roots"')
  dbDisconnect(con)

  lines.df <- data.frame(
    leaf_count = c(1, 2, 3, 4, 5, 6, 7) * 48,
    clr = rgb(0.2, 0.2, 0.8)
  )

  df <- df[order(df$leaf_count),]
  df$cum <- cumsum(df$leaf_count) / sum(df$leaf_count) * 1200
  df$humcum <- cumsum(df$human_count) / sum(df$human_count) * 1200

  print(nrow(subset(df, leaf_count <= 15)))
  df <- subset(df, leaf_count > 15)
  print(nrow(df))

  df <- df[order(df$species_count),]
  df$cumspec <- cumsum(rep(1, times=nrow(df))) / nrow(df) * 1200


  p <- ggplot(df, aes(x=leaf_count))
  p <- p + theme_bw()
  p <- p + geom_vline(data=lines.df, aes(xintercept=leaf_count, colour=clr), linetype='dashed')
  p <- p + scale_colour_identity()
  p <- p + geom_histogram(binwidth=5)
  p <- p + geom_line(aes(y=cum), colour='gray')
#  p <- p + geom_line(aes(y=humcum), colour='gray', linetype='dashed')

  p <- p + scale_x_continuous("Sequence Count")
  p <- p + scale_y_continuous("Count", limits=c(0, 1200))

  pdf(file="ensembl_roots_hist.pdf", height=3, width=6)
  print(p)
  dev.off()


  # Plot the Optic roots histogram, for comparison.
  lines.df <- data.frame(
    leaf_count = c(1, 2) * 8,
    clr = rgb(0.2, 0.2, 0.8)
  )
  source("compare_to_optic.R")
  df <- get.optic.data()
  print(str(df))

  top.y <- 5000

  df <- df[order(df$leaf_count),]
  df$cum <- cumsum(df$leaf_count) / sum(df$leaf_count) * top.y

  print(df[df$group_id == 214,])

  print(nrow(subset(df, leaf_count >= 20)))
  df$leaf_count <- pmin(df$leaf_count, 20)

  p <- ggplot(df, aes(x=leaf_count))
  p <- p + theme_bw()
  p <- p + geom_vline(data=lines.df, aes(xintercept=leaf_count, colour=clr), linetype='dashed')
  p <- p + scale_colour_identity()
  p <- p + geom_histogram(binwidth=1)
  p <- p + geom_line(aes(y=cum), colour='gray')

  p <- p + scale_x_continuous("Sequence Count")
  p <- p + scale_y_continuous("Count", limits=c(0, top.y))

  pdf(file="optic_roots_hist.pdf", height=3, width=6)
  print(p)
  dev.off()
}

fish.nodes.hist <- function() {
  con <- connect('gj1_orthologs_63');
  df <- dbGetQuery(con, 'select * from trees where method_id="Fish"')
  dbDisconnect(con)

  lines.df <- data.frame(
    leaf_count = c(1, 2, 3, 4, 5, 6, 7) * 48,
    clr = rgb(0.2, 0.2, 0.8)
  )

  big.y <- 4100

  df <- df[order(df$leaf_count),]
  df$cum <- cumsum(df$leaf_count) / sum(df$leaf_count) * big.y
  df$humcum <- cumsum(df$human_count) / sum(df$human_count) * big.y


  df$leaf_count <- pmin(df$leaf_count, 100)

  p <- ggplot(df, aes(x=leaf_count))
  p <- p + theme_bw()
  p <- p + geom_vline(data=lines.df, aes(xintercept=leaf_count, colour=clr), linetype='dashed')
  p <- p + scale_colour_identity()
  p <- p + geom_histogram(binwidth=2)
  p <- p + geom_line(aes(y=cum), colour='gray')
#  p <- p + geom_line(aes(y=humcum), colour='gray', linetype='dashed')

  p <- p + scale_x_continuous("Sequence Count", limits=c(0, 100))
  p <- p + scale_y_continuous("Count", limits=c(0, big.y))

  pdf(file="ensembl_fish_hist.pdf", height=3, width=6)
  print(p)
  dev.off()
}

euth.nodes.hist <- function() {
  con <- connect('gj1_orthologs_63');
  df <- dbGetQuery(con, 'select * from trees where method_id="Eutheria"')
  dbDisconnect(con)

  lines.df <- data.frame(
    leaf_count = c(1, 2, 3, 4, 5, 6, 7) * 48,
    clr = rgb(0.2, 0.2, 0.8)
  )

  big.y <- 2100

  df <- df[order(df$leaf_count),]
  df$cum <- cumsum(df$leaf_count) / sum(df$leaf_count) * big.y
  df$humcum <- cumsum(df$human_count) / sum(df$human_count) * big.y

  df$leaf_count <- pmin(df$leaf_count, 100)

  p <- ggplot(df, aes(x=leaf_count))
  p <- p + theme_bw()
  p <- p + geom_vline(data=lines.df, aes(xintercept=leaf_count, colour=clr), linetype='dashed')
  p <- p + scale_colour_identity()
  p <- p + geom_histogram(binwidth=2)
  p <- p + geom_line(aes(y=cum), colour='gray')
#  p <- p + geom_line(aes(y=humcum), colour='gray', linetype='dashed')

  p <- p + scale_x_continuous("Sequence Count", limits=c(0, 100))
  p <- p + scale_y_continuous("Count", limits=c(0, big.y))

  pdf(file="ensembl_euth_hist.pdf", height=3, width=6)
  print(p)
  dev.off()
}


plot.hist <- function() {
  #Fungi/Metazoa group  53
  #Vertebrata  48
  #Amniota  42
  #Mammalia  38
  #Eutheria  35
  #Primates  10
  #Glires  7
  #Laurasiatheria  12
  #Clupeocephala  5
  #Sauria  4
  methds <- factor(nodeset.levels(), levels=nodeset.levels())
  counts <- nodeset.sizes()
  clrs <- nodeset.colors()
  lines.df <- data.frame(
    method_id = methds,
    leaf_count = counts,
    clr = clrs
  )
  print(lines.df)

  con <- connect('gj1_orthologs_63');
  df <- dbGetQuery(con, 'select * from trees')
  dbDisconnect(con)

  df$method_id <- factor(df$method_id, levels=nodeset.levels())

  df$leaf_count <- pmin(300, df$leaf_count)

  p <- ggplot(df, aes(x=leaf_count, group=method_id))
  p <- p + theme_bw()
  p <- p + geom_histogram(binwidth=2)
  p <- p + geom_vline(data=lines.df, aes(xintercept=leaf_count, colour=clr), linetype='dashed')
  p <- p + scale_colour_identity()

  p <- p + facet_grid(method_id ~ .)
  p <- p + scale_x_continuous("Sequence Count")
  p <- p + scale_y_continuous("Count")

  p <- p + opts(
    strip.text.y = theme_text(angle=0, hjust=0, size=14)
  )

  pdf(file="ortholog_size_histogram.pdf", width=10, height=14)
  print(p)
  dev.off()
}

nodeset.colors <- function() {
  c(rep('red', times=9),
    rep('black', times=7)
  )
}

nodeset.sizes <- function() {
  c(10, # primates
    7, # glires
    12, # laur
    4, # sauria
    5, # fish
    35, # eutheria
    42, # amniotes
    48, # vertebrates
    53, # fungi/metazoa
    38, # mammal subgroups
    38, # mammal subgr + outgroup
    53, # human orths
    53, # mouse orths
    53, # zebr orths
    53, # dros orths
    53
    )
}

nodeset.levels <- function() {
  c('Primates', 'Glires', 'Laurasiatheria',
    'Sauria', 'Fish', 'Eutheria', 'Amniotes',
    'Vertebrata', 'Fungi/Metazoa group',
    'MammalSubgroups', 'MammalSubgroupsPlusOutgroup',
    'Human Orthologs', 'Mouse Orthologs', 'Zebrafish Orthologs', 'Drosophila Orthologs',
    'Ensembl Roots'
  )
}

ingroup.levels <- function() {
  c('Primates', 'Glires', 'Laurasiatheria',
  'Sauria', 'Fish')
}

outgroup.levels <- function() {
  c('Eutheria', 'Amniotes', 'Vertebrata', 'Fungi/Metazoa group')
}

ortho.levels <- function() {
  c('Human Orthologs', 'Mouse Orthologs', 'Zebrafish Orthologs', 'Drosophila Orthologs')
}

root.levels <- function() {
  c('Ensembl Roots')         
}

hi.coverage.species <- function() {
  rev(c('yeast', 'c_elegans', 'drosophila', 'c_intestinalis', 'c_savignyi',
    'zebrafish', 'fugu', 'tetraodon', 'medaka', 'stickleback',
    'xenopus', 'lizard', 'zebrafinch', 'chicken', 'turkey',
    'platypus', 'opossum', 'elephant', 'horse', 'pig', 'cow',
    'microbat', 'dog', 'panda', 'rabbit', 'guinea_pig', 'mouse', 'rat',
    'marmoset', 'rhesus', 'gibbon', 'orang', 'gorilla', 'chimp', 'human'))
}

low.coverage.species <- function() {
  rev(c('wallaby', 'sloth', 'armadillo', 'tenrec', 'hyrax', 'hedgehog', 'shrew',
  'alpaca', 'dolphin', 'megabat', 'cat', 'pika', 'squirrel', 'kangaroo_rat',
  #'tree_shrew', 
  'mouse_lemur', 'bushbaby', 'tarsier'))
}

all.species <- function() {
  c('human', 'chimp', 'gorilla', 'orang', 'gibbon', 'rhesus', 'marmoset',
  'tarsier', 'bushbaby', 'mouse_lemur', 
  # 'tree_shrew',
  'guinea_pig', 'squirrel', 'kangaroo_rat', 'mouse', 'rat', 'pika', 'rabbit',
  'hedgehog', 'shrew', 'horse', 'pig', 'alpaca', 'cow', 'dolphin', 'megabat', 'microbat', 'cat', 'dog', 'panda',
  'sloth', 'armadillo', 'tenrec', 'hyrax', 'elephant',
  'opossum', 'wallaby', 'platypus',
  'lizard', 'zebrafinch', 'chicken', 'turkey',
  'xenopus',
  'zebrafish', 'fugu', 'tetraodon', 'medaka', 'stickleback',
  'c_intestinalis', 'c_savignyi',
  'drosophila', 'c_elegans', 'yeast')
}

get.plot.df <- function() {
  if (!exists('plot.df', envir=.GlobalEnv)) {
    con <- connect('gj1_orthologs_63');
    df <- dbGetQuery(con, 'select * from trees')
    dbDisconnect(con)

    all.species <- c(hi.coverage.species(), low.coverage.species())
    all.coverage <- c(rep("High", length(hi.coverage.species())), rep("Low", length(low.coverage.species())))

    plot.df <- data.frame()
    for (i in 1:length(all.species)) {
      species <- all.species[i]
      print(species)
      coverage <- all.coverage[i]
      fld <- paste(species, '_count', sep='')
      
      count.field <- df[, fld]
      count.field[count.field == 2] <- 2
      count.field[count.field >= 3] <- 3

      tree.size <- df[, 'leaf_count']

      plot.df <- rbind(plot.df, data.frame(
        id=df$tree_id,
        method=df$method_id,
        species=species,
        coverage=coverage,
        count=count.field,
        leaf_count=tree.size
      ))
    }
    plot.df$method <- factor(plot.df$method, levels=nodeset.levels())
    plot.df$coverage <- factor(plot.df$coverage, levels=c("High", "Low"))
    plot.df$count <- factor(plot.df$count, levels=c(0, 1, 2, 3), labels=c('0', '1', '2', '3+'))
    save(plot.df, file="plot.df.Rdata")
    assign('plot.df', plot.df, envir=.GlobalEnv)
  }
  plot.df <- get('plot.df', envir=.GlobalEnv)
  plot.df
}

plot.dups <- function() {
  df <- get.plot.df()

#  plot.df <- subset(plot.df, method %in% c('Eutheria', 'Amniotes', 'Vertebrates'))
  plot.df <- subset(plot.df, coverage == "High")
  plot.df <- subset(plot.df, !(method %in% c('Primates', 'Glires', 'Laurasiatheria', 'MammalSubgroupsPlusOutgroup', 'MammalSubgroups')))

  plot.df$species <- factor(plot.df$species, levels=all.species())
  plot.df$count <- factor(plot.df$count, labels=c('0', '1', '2', '3+'))

  print(head(plot.df))
  p <- ggplot(plot.df, aes(x=species, fill=count))
  p <- p + theme_bw()
  p <- p + geom_bar()
  p <- p + facet_grid(method ~ ., scales="free")
  p <- p + scale_fill_discrete("Gene copy count")
  p <- p + scale_x_discrete('Species')
  p <- p + scale_y_continuous('Count')
  p <- p + opts(
    axis.text.x = theme_text(angle=90, hjust=1, size=9),
    strip.text.y = theme_text(angle=0, hjust=0, size=9)
  )

  pdf(file="ortholog_stacked_bar.pdf", width=8, height=8)
  print(p)
  dev.off()
}

plot.eutheria.dups <- function() {
  plot.df <- get.plot.df()
  plot.df <- subset(plot.df, method == 'Eutheria')

  plot.df$species <- factor(plot.df$species, levels=all.species())
  plot.df$n.dups <- factor(plot.df$count, labels=c('0', '1', '2', '3+'))
  plot.df$count <- NULL

  print(head(plot.df))
  p <- ggplot(plot.df, aes(x=species, fill=coverage))
  p <- p + theme_bw()
  p <- p + geom_bar(colour='black')

  p <- p + geom_text(aes(
    label = ifelse(..count.. > min(..count..) + base::diff(range(..count..)) / 2, paste(..count.., " "), paste(" ",..count..)),
    hjust = ifelse(..count.. > min(..count..) + base::diff(range(..count..)) / 2, 1, 0),
    colour = ifelse(..count.. > min(..count..) + base::diff(range(..count..)) / 2, 'black', 'black'),
    y=..count..),
    stat="bin", angle=90, size=2.5)
  p <- p + scale_fill_manual("Coverage", values=c(rgb(0.2, 0.2, 0.6), gray(0.8)))
  p <- p + scale_colour_identity("Coverage")
  p <- p + facet_grid(n.dups ~ ., scales="free")
  p <- p + scale_x_discrete('Species')
  p <- p + scale_y_continuous('Count')
  p <- p + opts(
    axis.text.x = theme_text(angle=90, hjust=1, size=12),
    strip.text.y = theme_text(angle=0, hjust=0.5, size=12)
  )

  pdf(file="ortholog_euth_dups.pdf", width=11, height=6)
  print(p)
  dev.off()  
}

plot.root.dups <- function() {
  plot.df <- get.plot.df()
  plot.df <- subset(plot.df, method == 'Ensembl Roots')

  plot.df$species <- factor(plot.df$species, levels=all.species())
  plot.df$count <- factor(plot.df$count, labels=c('0', '1', '2', '3+'))

  print(head(plot.df))
  p <- ggplot(plot.df, aes(x=species, fill=coverage))
  p <- p + theme_bw()
  p <- p + geom_bar(colour='black')
  p <- p + scale_fill_manual("Coverage", values=c(rgb(0.2, 0.2, 0.6), gray(0.8)))
  p <- p + scale_colour_identity("Coverage")
  p <- p + facet_grid(count ~ ., scales='free')
  p <- p + scale_x_discrete('Species')
  p <- p + scale_y_continuous('Count')
  p <- p + opts(
    axis.text.x = theme_text(angle=90, hjust=1, size=12),
    strip.text.y = theme_text(angle=0, hjust=0.5, size=12)
  )

  pdf(file="ortholog_root_dups.pdf", width=11, height=6)
  print(p)
  dev.off()  
  
}

plot.root.cdf <- function() {
  con <- connect('gj1_orthologs_63');
  df <- dbGetQuery(con, 'select * from trees where method_id="Ensembl Roots"')
  dbDisconnect(con)

  df <- df[order(df$leaf_count, decreasing=F),]
  df$cum <- cumsum(df$leaf_count) / sum(df$leaf_count)

  p <- ggplot(df, aes(x=leaf_count, y=cum))
  p <- p + theme_bw()
  p <- p + geom_line()
  p <- p + scale_x_continuous('Tree Sequence Count')#, limits=c(rev(range(df$leaf_count))))
  p <- p + scale_y_continuous('Cumulative Sequence Count')

  pdf(file="ortholog_root_cum.pdf")
  print(p)
  dev.off()  
}