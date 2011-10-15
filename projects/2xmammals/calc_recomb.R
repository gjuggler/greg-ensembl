source("~/src/greg-ensembl/projects/2xmammals/analyze_mammals.R")

get.recomb.rate <- function(sex='male', win.width=10000) {
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

  # First, try to recreate the deCODE map.
  rec.f <- scratch.f(sprintf("rec_%s_%d.Rdata", sex, win.width))
  if (!file.exists(rec.f)) {
    map.df <- read.table(f, header=T, stringsAsFactors=F)
    recomb.df <- process.map(map.df, win.width=win.width)
    save(recomb.df, file=rec.f)
  } else {
    load(rec.f)
  }
  recomb.df
}

process.map <- function(f, win.width=10000, skip.mb=5) {
  w <- win.width

  library(BSgenome.Hsapiens.UCSC.hg18)

  rec.df <- ddply(map.df, .(chr), function(x) {
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
      has.snp <- 0
      has.ns <- 0

      cur.seq <- subseq(chr.seq, start=lo, end=hi)
      print(maskedwidth(cur.seq))

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
        has.snp <- 1
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
      cM / mB
    }

    lo <- window.starts
    hi <- window.ends

    r.r <- c()
    pos <- c()
    for (i in 1:length(window.starts)) {
      recomb.rate <- interp.fn(lo[i], hi[i])
      r.r[i] <- recomb.rate
      midp <- (lo[i] + hi[i]) / 2
      pos[i] <- midp
    }
    data.frame(
      recombRate = r.r,
      pos = pos
    )
  })

  # Find the mean recombination value.

  rec.df
}