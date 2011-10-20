source("~/src/greg-ensembl/projects/2xmammals/analyze_mammals.R")

get.chrs <- function() {
  male.r <- scratch.f("rec_male.txt")
  map.df <- read.table(male.r, header=T, stringsAsFactors=F)
  unique(map.df$chr)
}

collect.recomb.rates <- function(width=10000) {
  recomb.f <- scratch.f(sprintf("recomb_%d.Rdata", width))
  if (!file.exists(recomb.f)) {
    male.r <- .collect.recomb.rates(sex='male', width=width)
    male.r$recombM <- male.r$recombRate
    female.r <- .collect.recomb.rates(sex='female', width=width)
    female.r$recombF <- female.r$recombRate
    male.r$recombRate <- NULL
    female.r$recombRate <- NULL
    
    recomb.df <- merge(male.r, female.r, by=c('chr', 'start', 'end'))

    # WAIT -- need to lift over from hg18 to hg19!!!
    print("Lifting over to hg19...")
    print(head(recomb.df))
    source("~/src/greg-ensembl/scripts/liftOver.R")
    recomb.df <- lift.over(recomb.df, 'hg18ToHg19', start.s='start', end.s='end', chr.s='chr')    
    print("Done!")
    print(head(recomb.df))

    save(recomb.df, file=recomb.f)
  }
  load(recomb.f)

  #print("Fraction zero from our calcs:")
  #print(nrow(subset(recomb.df, recombM == 0)) / nrow(recomb.df))

  #print("Fraction zero from Kong 2010 (10kb):")
  #male.r <- scratch.f("rec_male.txt")
  #map.df <- read.table(male.r, header=T, stringsAsFactors=F)
  #print(nrow(subset(map.df, stdrate == 0)) / nrow(map.df))

  recomb.df
}

.collect.recomb.rates <- function(sex='male', width=10000) {
  con <- connect(db())
  cmd <- sprintf("select * from recomb where sex='%s' and width=%d;", sex, width)
  x <- dbGetQuery(con, cmd)
  disconnect(con)

  # Get standardized rates: the mean across all windows with seqbin=1 should be 1.
  good.sites <- subset(x, seqbin==1)
  mean.rate <- mean(good.sites$recombRate)
  x$recombRate <- x$recombRate / mean.rate
  
  # Test against chr22
  # Kong 2010:
  #chr22 20687698  1 2.723375
  #chr22 20697698  1 4.483294
  #chr22 20707698  1 4.483294
  #chr22 20717698  1 4.466250
  #chr22 20727698  1 35.532378
  #chr22 20737698  1 22.766245
  # this approach at 10kb:
  #chr22 20687698  1 2.70319028
  #chr22 20697698  1 4.44998717
  #chr22 20707698  1 4.44998717
  #chr22 20717698  1 4.43306720
  #chr22 20727698  1 35.26840887
  #chr22 20737698  1 22.59713811
  # ... not bad, eh?

  x$start <- x$pos - floor(width/2)
  x$end <- x$pos + floor(width/2) - 1

  x <- subset(x, select=c('chr', 'start', 'end', 'recombRate'))
  x
}

bsub.calc.recomb.rate <- function() {
  chrs <- get.chrs()
  print(chrs)
  for (sex in c('male', 'female')) {
    for (width in c(1e4, 1e5, 1e6)) {
      for (chr in chrs) {
        xtra.s <- sprintf("%s %d %s", sex, width, chr)
        fn <- 'recomb_rate_calc'
        bsub.function(fn, extra.args=xtra.s, mem=6)
      }
    }
  }
}

calc.recomb.rate <- function(sex='male', win.width=10000, chr_name='chr1') {
  female.f <- scratch.f("gmap_female.txt")
  male.f <- scratch.f("gmap_male.txt")
  female.r <- scratch.f("rec_female.txt")
  male.r <- scratch.f("rec_male.txt")
  if (!file.exists(female.f)) {
    system(sprintf("wget http://www.decode.com/addendum/female.gmap -O %s", female.f))
  }
  if (!file.exists(male.f)) {
    system(sprintf("wget http://www.decode.com/addendum/male.gmap -O %s", male.f))
  }
  if (!file.exists(female.r)) {
    system(sprintf("wget http://www.decode.com/addendum/female.rmap -O %s", female.r))
  }
  if (!file.exists(male.r)) {
    system(sprintf("wget http://www.decode.com/addendum/male.rmap -O %s", male.r))
  }

  if (sex == 'male') {
    f <- male.f
  } else {
    f <- female.f
  }

  map.df <- read.table(f, header=T, stringsAsFactors=F)
  print(nrow(map.df))
  map.df <- subset(map.df, chr == chr_name)
  print(nrow(map.df))

  recomb.df <- process.map(map.df, win.width=win.width)
  recomb.df$sex <- sex
  recomb.df$width <- win.width
  recomb.df$label <- paste(recomb.df$sex, recomb.df$width, recomb.df$pos, sep=' ')
  print(nrow(recomb.df))

  con <- connect(db())
  write.or.update(recomb.df, 'recomb', con, 'label')
  disconnect(con)
}

process.map <- function(map.df, win.width=10000, skip.mb=5) {
  w <- win.width

  # Remember -- this is all done on NCBI36 / hg18!!!
  library(BSgenome.Hsapiens.UCSC.hg18)

  rec.fn <- function(x) {
    x <- x[order(x$pos),]
    min.pos <- min(x$pos)
    max.pos <- max(x$pos)
    
    chr <- x[1, 'chr']
    chr.seq <- Hsapiens[[chr]]

    # Start the windows at 5mb after the first SNP
    min.w <- min.pos + skip.mb*1e6
    max.w <- max.pos - skip.mb*1e6
    window.starts <- seq(from=min.w, to=max.w-w, by=w)
    window.ends <- window.starts + w
    window.mids <- window.starts + floor(w/2)

    interp.fn <- function(lo, hi) {
      # Interpolate out to the marker below the window region.
      cM <- 0

      cur.seq <- subseq(chr.seq, start=lo, end=hi)
      # seqbin is 0 if there are missing seq's in the assembly.
      seqbin <- ifelse(maskedwidth(cur.seq) > 0, 0, 1)

      prev.i <- max(which(x$pos < lo))
      next.i <- min(which(x$pos >= hi))
      within.i <- which(x$pos >= lo & x$pos < hi)
      if (!any(within.i)) {
        # No markers within the window, so we interpolate from prev and next.
        snp.phys.d <- x[next.i, 'pos'] - x[prev.i, 'pos']
        snp.gen.d <- x[next.i, 'cM']
        mid.phys.f <- (hi - lo) / snp.phys.d
        cM <- cM + snp.gen.d * mid.phys.f
      } else {
        if (length(within.i) == 1) {
          # One marker within the window -- set it to first and last, so we get fractional
          # genetic distance from each fragment to the left & right.
          first.i <- within.i
          last.i <- within.i
        } else {
          # Multiple markers within the window. Add the sum of the 'within' distance.
          # Make sure not to use the cM from the first marker, since this is taken
          # care of with the 'low.phys.d' below.
          first.i <- min(within.i)
          last.i <- max(within.i)
          cM <- cM + sum(x[(first.i+1):last.i, 'cM'], na.rm=T)
        }
        low.phys.d <- x[first.i, 'pos'] - x[prev.i, 'pos']
        low.phys.f <- (x[first.i, 'pos'] - lo) / low.phys.d
        low.gen.d <- x[first.i, 'cM'] * low.phys.f
        cM <- cM + low.gen.d  

        hi.phys.d <- x[next.i, 'pos'] - x[last.i, 'pos']
        hi.phys.f <- (hi - x[last.i, 'pos']) / hi.phys.d
        hi.gen.d <- x[next.i, 'cM'] * hi.phys.f
        cM <- cM + hi.gen.d  
      }
      mB <- (hi - lo) / 1e6
      c(cM / mB, seqbin)
    }

    lo <- window.starts
    hi <- window.ends

    r.r <- c()
    seqbin <- c()
    pos <- c()
    for (i in 1:length(window.starts)) {
      retval <- interp.fn(lo[i], hi[i])
      r.r[i] <- retval[1]
      seqbin[i] <- retval[2]
      midp <- (lo[i] + hi[i]) / 2
      pos[i] <- midp
    }
    data.frame(
      pos = pos,
      seqbin = seqbin,
      recombRate = r.r
    )
  }

  rec.df <- ddply(map.df, .(chr), rec.fn)
  rec.df
}