source("~/src/greg-ensembl/projects/2xmammals/manuscript_numbers.R")
source("~/src/greg-ensembl/scripts/collect_sitewise.R")
library(ggplot2)
library(plyr)
library(GenomicRanges)

sites.correlations <- function() {
  if (!exists("sites")) {
    load("~/scratch/sites_60_g.Rdata")  
    sites$has.pf <- !is.na(sites$pfam_domain)
    sites <- subset(sites, !is.na(splice_distance))
  }

  dx <- 20
  mult <- 3
  sites$trunc_splice <- round(sites$splice_distance / mult)
  sites$trunc_splice <- sites$trunc_splice * mult
  sites[sites$trunc_splice < -dx, 'trunc_splice'] <- -dx
  sites[sites$trunc_splice > dx, 'trunc_splice'] <- dx

  summ.f <- function(x) {
    data.frame(
      n=nrow(x),
      median=median(x$omega),
      mean=mean(x$omega),
      mean_lrt=mean(x$signed_lrt),
      pos_neg=nrow(subset(x, signed_lrt > 1)) / nrow(x)
    )
  }

  pf.sites <- subset(sites, !is.na(pfam_domain))
  
  res <- ddply(pf.sites, .(trunc_splice), summ.f)
  print(res)
  
  res <- ddply(pf.sites, .(exon_position), summ.f)
  print(res)
  
  res <- ddply(sites, .(has.pf), summ.f)
  print(res)
}

recomb.rates <- function() {
  if(!exists("sites")) {
    load("/nfs/users/nfs_g/gj1/scratch/gj1_2x_57/2011-01-21_01/sites_1.Rdata")
    assign("sites", sites, envir=.GlobalEnv)
  }

  dir = "~/src/greg-ensembl/projects/2xmammals"
  if (!exists("recomb")) {
    recomb <- read.table(paste(dir,"/recombRate.txt",sep=''), header=F)
    names(recomb) <- c(
      'chr', 'start', 'end',
      'name',
      'decodeAvg', 'decodeF', 'decodeM',
      'marshfieldAvg', 'marshfieldF', 'marshfieldM',
      'genethonAvg', 'genethonF', 'genethonM'
    )  
    assign("recomb",recomb, envir=.GlobalEnv)
  }

  if (!exists("gc.content")) {
    library(BSgenome.Hsapiens.UCSC.hg19)
    
    window.width <- 10000

    chr.names <- seqnames(Hsapiens)[1:25]
    gc.content <- data.frame()
    for (i in 1:length(chr.names)) {
      chr.name <- chr.names[i]
      cur.seq <- Hsapiens[[chr.name]]
      print(chr.name)
      cur.len <- length(cur.seq)
      
      window.starts <- seq(from=1, to=cur.len-window.width, by=window.width)
      window.ends <- seq(from=window.width, to=cur.len, by=window.width)

      v <- Views(cur.seq, start=window.starts, end=window.ends)

      gc.f <- function(x) {
        af <- alphabetFrequency(x)
        sum(af[2:3])/sum(af)
      }
      gcs <- viewApply(v, gc.f, simplify=T)

      cur.gc.df <- data.frame(
        chr=chr.name,
        start=window.starts,
        end=window.ends,
        gc=gcs
      )
      gc.content <- rbind(gc.content, cur.gc.df)
    }

    assign("gc.content", gc.content, envir=.GlobalEnv)    
  }

  if (!exists("comb.df")) {
    recomb.ranges <-
      GRanges(seqnames = recomb$chr,
        ranges =
          IRanges(recomb$start, recomb$end),
        strand = "*",
        recomb_m = recomb$decodeM,
        recomb_f = recomb$decodeF,
        recomb_avg = recomb$decodeAvg
    )

    sites <- sign.lrt(sites)
    pos.sites <- sites
  #  pos.sites <- subset(pos.sites, signed_lrt > 0)
    pos.sites <- subset(pos.sites, !is.na(chr_start))
#    pos.sites <- subset(pos.sites, ncod > 10)
#    pos.sites <- subset(pos.sites, nongap_bl > 1)  

    pos.ranges <-
      GRanges(
        seqnames = pos.sites$chr_name,
        ranges = IRanges(pos.sites$chr_start, pos.sites$chr_start+2),
        strand = "*",
        signed_lrt = pos.sites$signed_lrt,
        omega = pos.sites$omega,
        lrt_stat = pos.sites$lrt_stat,
        pfam_domain = pos.sites$domain,
        node_id = pos.sites$node_id
    )

    match.indices <- match(pos.ranges, recomb.ranges)
    pos.ranges <- pos.ranges[which(!is.na(match.indices))]
    match.indices <- match.indices[!is.na(match.indices)]
    comb.df <- cbind(
      values(pos.ranges),
      values(recomb.ranges[match.indices])
    )
    assign("pos.ranges", pos.ranges, envir=.GlobalEnv)

    comb.df <- as.data.frame(comb.df)
    comb.df$chr_start <- start(pos.ranges)
    comb.df$chr_name <- as.data.frame(pos.ranges)$seqnames

    pos.ranges <-
      GRanges(
        seqnames = comb.df$chr_name,
        ranges = IRanges(comb.df$chr_start, comb.df$chr_start+2),
        strand = "*",
        signed_lrt = comb.df$signed_lrt,
        omega = comb.df$omega,
        lrt_stat = comb.df$lrt_stat,
        pfam_domain = comb.df$pfam_domain,
        node_id = comb.df$node_id,
        recomb_m = comb.df$recomb_m,
        recomb_f = comb.df$recomb_f,
        recomb_avg = comb.df$recomb_avg
    )

    gc.ranges <-
      GRanges(seqnames = gc.content$chr,
        ranges =
          IRanges(gc.content$start, gc.content$end),
        strand = "*",
        gc = gc.content$gc
    )

    match.indices <- match(pos.ranges, gc.ranges)
    pos.ranges <- pos.ranges[which(!is.na(match.indices))]
    match.indices <- match.indices[!is.na(match.indices)]
    comb.df <- cbind(
      values(pos.ranges),
      values(gc.ranges[match.indices])
    )

    comb.df <- as.data.frame(comb.df)
    comb.df$chr_start <- start(pos.ranges)
    comb.df$chr_name <- as.data.frame(pos.ranges)$seqnames

    assign("comb.df", comb.df, envir=.GlobalEnv)
  }

  x <- comb.df

  higc <- subset(x, gc >= quantile(gc, 0.9) & recomb_m > 8.73534)
  print(summ.f(higc))

  v.pos.sites <- subset(x, signed_lrt > 10)
  pos.sites <- subset(x, signed_lrt > 2.7)
  neg.sites <- subset(x, signed_lrt < -2.7)
  v.neg.sites <- subset(x, signed_lrt < -10)


  pos.hi.gc <- subset(x, signed_lrt >= quantile(signed_lrt, 0.99) & gc >= quantile(gc, 0.75))
  hi.gc <- subset(x, gc >= quantile(gc, 0.75))
  print(nrow(pos.hi.gc))

  res <- t.test(pos.hi.gc$recomb_m, hi.gc$recomb_m)
  print(res)
  res <- t.test(pos.hi.gc$recomb_f, hi.gc$recomb_f)
  print(res)

  res <- t.test(pos.sites$recomb_m, x$recomb_m)
  print(res)
  res <- t.test(pos.sites$recomb_f, x$recomb_f)
  print(res)
  res <- t.test(pos.sites$gc, x$gc)
  print(res)

#  res <- wilcox.test(pos.sites$recomb_f, x$recomb_f, alternative="greater")
#  print(res)
#  res <- wilcox.test(pos.sites$recomb_f, x$recomb_f, alternative="less")
#  print(res)

#  p <- ggplot(comb.df, aes(x=recomb_m, y=signed_lrt))
#  p <- p + stat_bin2d(bins=50)
#  pdf(file="~/src/greg-ensembl/projects/2xmammals/recomb.pdf")
#  print(p)
#  dev.off()

  print("sites")

  summarize.quantiles <- function(df, var, quantiles) {
    print(paste(var,paste(quantiles, collapse=' / '), collapse=': '))
    df$var <- df[, var]
    cur.list <- list()
    quantiles <- sort(quantiles)
    q.values <- quantile(df$var, probs=quantiles, names=F)
    lower.quantiles <- c()
    upper.quantiles <- c()
    for (i in 1:length(quantiles)) {
      quantile.lbl.hi <- quantiles[i]
      quantile.lbl.lo <- 0
      quantile.value <- q.values[i]
      quantile.lower <- min(df$var) - 0.01
      if (i > 1) {
        quantile.lower <- q.values[i-1]
        quantile.lbl.lo <- quantiles[i-1]
      }
      lower.quantiles <- c(lower.quantiles, quantile.lower)
      upper.quantiles <- c(upper.quantiles, quantile.value)
      cur.df <- subset(df, var > quantile.lower & var <= quantile.value)
      print(paste(var, quantile.lbl.lo, quantile.lbl.hi, nrow(cur.df)))
      cur.id <- paste(quantile.lbl.lo, '-', quantile.lbl.hi, sep='')
      cur.list[[cur.id]] <- cur.df
    }
    tbl <- ldply(cur.list, summ.f)
    tbl$lo <- lower.quantiles
    tbl$hi <- upper.quantiles
    tbl$variable <- var
    return(tbl)
  }

  res <- summarize.quantiles(x, 'signed_lrt', c(0.01, 0.25, 0.5, 0.75, 0.99, 1.0))
  write.csv(res, file="~/src/greg-ensembl/projects/2xmammals/signed_lrt.csv", row.names=F)
  res <- summarize.quantiles(x, 'recomb_m', c(0.01, 0.25, 0.5, 0.75, 0.99, 1.0))
  write.csv(res, file="~/src/greg-ensembl/projects/2xmammals/recomb_m.csv", row.names=F)
  res <- summarize.quantiles(x, 'gc', c(0.01, 0.25, 0.5, 0.75, 0.99, 1.0))
  write.csv(res, file="~/src/greg-ensembl/projects/2xmammals/gc.csv", row.names=F)

  pos.sites <- subset(x, signed_lrt > 2.7)
  res <- summarize.quantiles(pos.sites, 'recomb_m', c(0.01, 0.25, 0.5, 0.75, 0.99, 1.0))
  write.csv(res, file="~/src/greg-ensembl/projects/2xmammals/pos_sites_recomb_m.csv", row.names=F)
  res <- summarize.quantiles(pos.sites, 'signed_lrt', c(0.01, 0.25, 0.5, 0.75, 0.99, 1.0))
  write.csv(res, file="~/src/greg-ensembl/projects/2xmammals/pos_sites_signed_lrt.csv", row.names=F)

  gc_hi <- subset(x, gc >= quantile(gc, 0.9))
  res.lo <- summarize.quantiles(gc_hi, 'recomb_f', c(0.01, 0.25, 0.5, 0.75, 0.99, 1.0))
  res.hi <- summarize.quantiles(gc_hi, 'recomb_m', c(0.01, 0.25, 0.5, 0.75, 0.99, 1.0))
  res.lo$gc_type <- 'Low GC'
  res.hi$gc_type <- 'High GC'
  res.combined <- rbind(res.lo, res.hi)
  write.csv(res.combined, file="~/src/greg-ensembl/projects/2xmammals/gc_recomb_m.csv", row.names=F)

  pos.hi.gc <- subset(x, signed_lrt >= quantile(signed_lrt, 0.99) & gc >= quantile(gc, 0.8))
  res.m <- summarize.quantiles(pos.hi.gc, 'recomb_m', c(0.95, 0.99, 1.0))
  res.f <- summarize.quantiles(pos.hi.gc, 'recomb_f', c(0.95, 0.99, 1.0))
  res.m$recomb_type <- 'Male Recomb'
  res.f$recomb_type <- 'Female Recomb'
  write.csv(rbind(res.m, res.f), file="~/src/greg-ensembl/projects/2xmammals/hi_gc_pos_sites.csv", row.names=F)

  return()
  
  # Test for recombination correlation in sites with low and high GC content.
  gc.tests <- list(
    all.sites = x,
    `gc_lo_recomb_lo` = subset(x, gc <= quantile(gc, 0.25) & recomb_m <= quantile(recomb_m, 0.25)),
    `gc_lo_recomb_hi` = subset(x, gc <= quantile(gc, 0.25) & recomb_m > quantile(recomb_m, 0.75)),
    `gc_hi_recomb_lo` = subset(x, gc > quantile(gc, 0.75) & recomb_m <= quantile(recomb_m, 0.25)),
    `gc_hi_recomb_hi` = subset(x, gc > quantile(gc, 0.75) & recomb_m > quantile(recomb_m, 0.75))
  )
  gc.test.tbl <- ldply(gc.tests, summ.f)
  print(head(gc.test.tbl))
  write.csv(gc.test.tbl, file="~/src/greg-ensembl/projects/2xmammals/recomb_gc_recomb.csv", row.names=F)

  return()

  # http://www.yilab.gatech.edu/pcor.R


  mult <- 1
  max.r <- 6

  print(summ.f(x))

  print("Recomb F")
  x$var <- x$recomb_f
  x[x$var > max.r,]$var <- max.r
  x$var <- round(x$var * mult) / mult
  res <- ddply(x, .(var), summ.f)
  print(res)  

  print("Recomb M")
  x$var <- x$recomb_m
  x[x$var > max.r,]$var <- max.r
  x$var <- round(x$var * mult) / mult
  res <- ddply(x, .(var), summ.f)
  print(res)

  print("Signed LRT")
  mult <- 0.5
  max.r <- 25
  min.r <- -15
  x$var <- x$signed_lrt
  x[x$var > max.r,]$var <- max.r
  x[x$var < min.r,]$var <- min.r
  x$var <- round(x$var * mult) / mult
  res <- ddply(x, .(var), summ.f)
  print(res)
  
  #up <- subset(x, signed_lrt > 0)
  #down <- subset(x, signed_lrt < 0)
  #print(cor.test(up$recomb_avg, up$signed_lrt))
  #print(cor.test(down$recomb_avg, down$signed_lrt))

  # Test for a difference in the distribution of recombination rates for pos vs all sites

  

}

summ.f <- function(x) {
  lrt.thresh <- 2.7
  n <- nrow(x)
  data.frame(
    n=nrow(x),
    mean_w=mean(x$omega),
    mean_lrt=mean(x$signed_lrt),
    median_lrt=median(x$signed_lrt),
    f_blw_p3 = nrow(subset(x, omega < 0.3)) / n,
    f_abv_one = nrow(subset(x, omega > 1)) / n,
    f_abv_two = nrow(subset(x, omega > 2)) / n,
    recomb_m = mean(x$recomb_m),
    recomb_f = mean(x$recomb_f),
    recomb_avg = mean(x$recomb_avg),
    pos_f=(nrow(subset(x, signed_lrt > lrt.thresh)) / nrow(x)),
    neg_f=(nrow(subset(x, signed_lrt < -lrt.thresh)) / nrow(x)),
    neut_f=(nrow(subset(x, signed_lrt > -lrt.thresh & signed_lrt < lrt.thresh)) / nrow(x)),
    gc = mean(x$gc)
  )
}

