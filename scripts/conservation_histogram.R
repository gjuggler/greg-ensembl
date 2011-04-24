# This function takes a list of background rates (nonvar.sites), variant information
# (var.sites), and the name of the column we want to plot (hist.column) and creates
# a histogram on top with a factored "rug" plot of variants below.
do.histogram <- function(
  nonvar.sites, 
  var.sites, 
  hist.column,
  types.column,
  types.label="Variant",
  rug.track.size=1/3, 
  do.plot=F
) {

  # Store the chosen histogram column as a new 'rate' column.
  nonvar.sites$rate <- nonvar.sites[, hist.column]
  var.sites$rate <- var.sites[, hist.column]

  print(summary(nonvar.sites$rate))
  print(summary(var.sites$rate))

  nonvar.sites$variant.type <- nonvar.sites[, types.column]
  var.sites$variant.type <- var.sites[, types.column]

  # Use the hist() function to generate breakpoints.
  h <- hist(nonvar.sites$rate, n=30, plot=F)
  dx <- diff(h$breaks)
  dx <- dx[1]/2
  max.count <- max(h$counts)
  # Create a new data frame storing the histogram info.
  df <- data.frame(
    count=h$counts,
    rate=as.numeric(h$mids),
    dx=as.numeric(dx)
  )

  # Sort the variants
  types <- sort(unique(var.sites$variant.type))
  tracks.below <- length(types)
  # Size the rug plot relative to the histogram height.
  t.h <- max.count * rug.track.size

  # Calculate the y-coord range and axis labels.
  y.lo <- -(t.h/5 + t.h * (length(types)))
  y.hi <- max(df$count) * 1.5
  brk.hi <- y.hi
  if (brk.hi > 100) {
    brk.hi <- y.hi - (y.hi %% 100)
  }
  y.brk <- c(y.lo, seq(from=0, to=brk.hi, length.out=5))
  y.pct <- sprintf("%.0f / %.0f%% ",y.brk, (y.brk / brk.hi) * 100)
  y.labels <- c('', y.pct[2:length(y.pct)])

  all.ecdf <- ecdf(nonvar.sites$rate)
  nonvar.sites$cumsum <- all.ecdf(nonvar.sites$rate) * brk.hi

  # Plot the main histogram.
  p <- ggplot(data=df, aes(x=rate))
  p <- p + geom_rect(data=df, aes(xmin=rate-dx, xmax=rate+dx, ymin=0, ymax=count))

  # Overlay the ecdf.
  p <- p + geom_line(data=nonvar.sites, color='gray', aes(y=cumsum))

  # Add a geom_segment rug plot for each variant type
  for (i in 1:length(types)) {
    typ <- types[i]
    y_hi <- -t.h/5 - t.h*(i-1)
    y_lo <- -t.h/5 - t.h*(i)
    sub.vars <- subset(var.sites, variant.type==typ)
    sub.vars$y_hi <- y_hi
    sub.vars$y_lo <- y_lo

    # Order the factor so we have the same order as the PAML plot.
    sub.vars <- sub.vars[order(sub.vars$rate),]
    p <- p + geom_segment(data=sub.vars, alpha=0.7, aes(y=y_hi, yend=y_lo, xend=rate, x=rate, colour=variant.type, size=grantham_score))
  }

  name = ""
  title = ""
  if (hist.column == "paml") {
    title = "Evolutionary rate for nucleotides in NBEAL2"
    name = "Nucleotide evolutionary rate (PAML)"
  } else if (hist.column == "slr_m") {
    title = "Evolutionary constraint for codons in NBEAL2"
    name = "Codon evolutionary constraint (SLR)"
  } else if (hist.column == "slr_dnds") {
    title = "Evolutionary constraint (dN/dS) for codons in NBEAL2"
    name = "dN/dS (SLR)"
  } else {
    title = "Selection for codons"
    name = "Selection (SLR)"
  }

  p <- p + scale_x_continuous(name=name)
  p <- p + scale_y_continuous(name="Bin count / cumulative %", breaks=y.brk, labels=y.labels)
  p <- p + scale_colour_discrete(name=types.label)
  #p <- p + scale_linetype_discrete(name="Variant effect")
  p <- p + opts(title=title)
  if (do.plot) {
    print(p)
  }
  return(p)
}

# Formats numeric columns of a data frame to a specified precision.
format.numeric.df <- function(x,digits=3) {  
  nums <- unlist(lapply(x,is.double))
  print(nums)
  for (col in names(x)) {
    if (nums[col] == TRUE) {
      x[,col] <- formatC(x[,col],digits=digits,format='fg')
    }
  }
  return(x)
}