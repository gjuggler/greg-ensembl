
## grid/lattice version of smoothScatter
## courtesy of D. Sarkar

# Require the KernSmooth library.
library(KernSmooth)

panel.smoothScatter <-
    function (x, y = NULL,
              nbin = 64,
              cuts = 255,
              bandwidth,
              colramp = colorRampPalette(c("white", brewer.pal(9, "Blues"))),
              nrpoints = 100,
              transformation = function(x) x^0.25,
              pch = ".",
              cex = 1, col="black",
              ...,
              subscripts)

    ## subscripts is ignored (and recomputed), but the formal argument
    ## necessary to catch another version passed in as ...

    ## If time permits, replacing the call to panel.levelplot by
    ## panel.rect might make the code more transparent (which would
    ## also make the subscripts thing unnecessary).  The only drawback
    ## I can think of is that we will no longer get contourplots
    ## instead of levelplots, but I don't think this is the right
    ## place for that anyway.

{
    x <- as.numeric(x)
    y <- as.numeric(y)
    if (!is.numeric(nrpoints) | (nrpoints < 0) | (length(nrpoints) != 1))
        stop("'nrpoints' should be numeric scalar with value >= 0.")
    xy <- xy.coords(x, y)
    x <- cbind(xy$x, xy$y)[!(is.na(xy$x) | is.na(xy$y)), ]
    map <- .smoothScatterCalcDensity(x, nbin, bandwidth)
    xm <- map$x1
    ym <- map$x2
    dens <- map$fhat
    dens <- array(transformation(dens), dim = dim(dens))
    panel.levelplot(x = rep(xm, length(ym)),
                    y = rep(ym, each = length(xm)),
                    z = as.numeric(dens),
                    subscripts = TRUE,
                    at = seq(from = 0, to = 1.01 * max(dens), length = cuts + 2),
                    col.regions = colramp(cuts + 1),
                    ...)
    if (nrpoints != 0)
    {
        stopifnot(length(xm) == nrow(dens), length(ym) == ncol(dens))
        ixm <- round((x[, 1] - xm[1])/(xm[length(xm)] - xm[1]) *
                     (length(xm) - 1))
        iym <- round((x[, 2] - ym[1])/(ym[length(ym)] - ym[1]) *
                     (length(ym) - 1))
        idens <- dens[1 + iym * length(xm) + ixm]
        nrpoints <- min(nrow(x), ceiling(nrpoints))
        sel <- order(idens, decreasing = FALSE)[1:nrpoints]
        panel.points(x[sel, 1:2], pch = pch, cex = cex, col = col)
    }
    ## panel.abline(h=0, col="#fe0020") # why was this here?
}


.smoothScatterCalcDensity <- function(x, nbin, bandwidth) {
  
  if (length(nbin) == 1)
    nbin <- c(nbin, nbin)
  if (!is.numeric(nbin) || (length(nbin)!=2))
    stop("'nbin' must be numeric of length 1 or 2")

  if (missing(bandwidth)) {
    bandwidth <- diff(apply(x, 2, quantile, probs=c(0.05, 0.95), na.rm=TRUE)) / 25
  } else {
    if(!is.numeric(bandwidth))
      stop("'bandwidth' must be numeric")
  }

  ## create density map
  rv <- bkde2D(x, gridsize=nbin, bandwidth=bandwidth)
  rv$bandwidth <- bandwidth
  return(rv)
}

smoothScatter <- function(x, y=NULL, 
                          nbin=128,
                          bandwidth,
                          colramp=colorRampPalette(c("white", brewer.pal(9, "Blues"))),
                          nrpoints=100,
                          transformation=function(x) x^.25,
                          xlab=NULL, ylab=NULL, postPlotHook=box,
                          pch=".", cex=1,
                          xlim, ylim, col="black",
                          xaxs=par("xaxs"), yaxs=par("yaxs"), ...) {
  
  if (!is.numeric(nrpoints) | (nrpoints<0) | (length(nrpoints)!=1) )
    stop("'nrpoints' should be numeric scalar with value >= 0.")

  ## similar as in plot.default
  xlabel <- if (!missing(x)) 
    deparse(substitute(x))
  ylabel <- if (!missing(y)) 
    deparse(substitute(y))
  xy <- xy.coords(x, y, xlabel, ylabel)
  xlab <- if (is.null(xlab)) 
    xy$xlab
  else xlab
  ylab <- if (is.null(ylab)) 
    xy$ylab
  else ylab


  ## eliminate NA
  x <- cbind(xy$x, xy$y)[!(is.na(xy$x)|is.na(xy$y)), ]

  ## xlim and ylim
  if(!missing(xlim)) {
    stopifnot(is.numeric(xlim), length(xlim)==2, !any(is.na(xlim)))
    x <- x[ (x[,1]>=xlim[1]) & (x[,1]<=xlim[2]), ]
  } else {
    xlim <- range(x[,1], na.rm=TRUE)
  }
  if(!missing(ylim)) {
    stopifnot(is.numeric(ylim), length(ylim)==2, !any(is.na(ylim)))
    x <- x[ (x[,2]>=ylim[1]) & (x[,2]<=ylim[2]), ]
  } else {
    ylim <- range(x[,2], na.rm=TRUE)
  }
  
  ## create density map
  map  <- .smoothScatterCalcDensity(x, nbin, bandwidth)
  xm   <- map$x1
  ym   <- map$x2
  dens <- map$fhat
  dens <- array(transformation(dens), dim=dim(dens))
  
  ## plot color image
  image(xm, ym, z=dens, col=colramp(256),
        xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, xaxs=xaxs, yaxs=yaxs, ...)
  if(!is.null(postPlotHook)) postPlotHook()
  
  ## plot selection of dots
  if (nrpoints!=0){
    ## we assume that map$x1 and map$x2 go linearly from
    ## their first to their last value in nbin steps
    stopifnot(length(xm)==nrow(dens), length(ym)==ncol(dens))
    ixm <- round((x[,1]-xm[1])/(xm[length(xm)]-xm[1])*(length(xm)-1))
    iym <- round((x[,2]-ym[1])/(ym[length(ym)]-ym[1])*(length(ym)-1))
    idens <- dens[1 + iym*length(xm) + ixm]
    nrpoints <- min(nrow(x), ceiling(nrpoints))
    sel <- order(idens, decreasing=FALSE)[1:nrpoints]
    points(x[sel,1:2], pch=pch, cex=cex, col=col)
  }
}



densCols <- function(x, y=NULL,
                     nbin=128,
                     bandwidth,
                     colramp=colorRampPalette(brewer.pal(9, "Blues")[-(1:3)])) {

  ## similar as in plot.default
  xy <- xy.coords(x, y)

  ## deal with NA
  select <- !(is.na(xy$x)|is.na(xy$y))
  x <- cbind(xy$x, xy$y)[select, ]
  
  ## create density map
  map  <- .smoothScatterCalcDensity(x, nbin, bandwidth)

  ## bin x-values
  xbin <- cut(x[,1], map$x1-(diff(range(map$x1))/(length(map$x1)-1)/2),
              labels=FALSE)

  ## bin y-values
  ybin <- cut(x[,2], map$x2-(diff(range(map$x2))/(length(map$x2)-1)/2),
              labels=FALSE)
    
  dens <- map$fhat[xbin + nrow(map$fhat) * (ybin-1)]
  dens[is.na(dens)]<- 0

  ## transform densities to colors
  colpal <- cut(dens, length(dens), labels=FALSE)
  cols   <- rep(as.character(NA), nrow(x))
  cols[select] <- colramp(length(dens))[colpal]
    
  return(cols)
}
