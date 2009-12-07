### slr.tools.R -- SLR plot tools for R  ###
###  written by Greg Jordan, fall 2008   ###
###  email: greg@ebi.ac.uk               ###

read.slr = function(file) {
  return(read.table(file,header=T,sep="\t"))
}

plot.slr.types = function(slr,xlim=NULL) {
  minr=1
  maxr=max(slr$site)

  ylim=c(-.2,.2)
  if (is.null(xlim)) xlim = c(1,maxr)

  pos_nocorr = paste("positive",c(1,2),sep="")
  pos_corr = paste("positive",c(3,4),sep="")

  pn = slr[slr$type %in% pos_nocorr,]
  pc = slr[slr$type %in% pos_corr,]
  plot.new()
  plot.window(xlim=xlim,ylim=ylim)
  rect(xleft=pn$site-.5,ybottom=rep(0,nrow(pn)),xright=pn$site+.5,ytop=rep(.1,nrow(pn)),col=rgb(1,.5,.5,.5),border=NA)
  rect(xleft=pc$site-.5,ybottom=rep(.1,nrow(pn)),xright=pc$site+.5,ytop=rep(.2,nrow(pc)),col=rgb(1,0,0,1),border=NA)

  syn = slr[slr$note == "synonymous",]
  const = slr[slr$note == "constant",]
  rect(xleft=syn$site-.5,ybottom=rep(-.1,nrow(syn)),xright=syn$site+.5,ytop=rep(0,nrow(syn)),col='gray',border=NA)
  rect(xleft=const$site-.5,ybottom=rep(-.2,nrow(const)),xright=const$site+.5,ytop=rep(-.1,nrow(const)),col='blue',border=NA)
}

plot.slr.blocks = function(slr,
  overlay=FALSE,
  xlim=NULL,
  ylim=c(0,5),
  block.xlim = c(0,1),
  colorByType = TRUE,
  color.lrt = TRUE,
  log=F,
  ...
) {
  if (is.character(slr)) {
    print("Reading SLR file...");
    slr = read.table(slr,header=TRUE,sep="\t")
  }

  names(slr)
  print(slr[1,])
  minr = 1
  maxr = nrow(slr)

  log_lo = 0.001

  # Make sure the axis type is correct.
  par(xaxs='i')

  if (!overlay) {
    if(is.null(xlim)) xlim = c(minr,maxr)
    if (log) {
      plot(NULL,xlim=xlim,ylim=c(log_lo,100),xlab="Site",ylab="Omega",log='y',axes=F,...)
    } else {
      plot(NULL,xlim=xlim,ylim=c(0,5),xlab="Site",ylab="Omega",log='',axes=F,...)
    }
  }

  if (nrow(slr) == 0) {return();};

  # Rectangles for error bars.
  #mids = minr:maxr
  mids = slr$site
  vals = slr$omega
  los = slr$lower
  his = slr$upper
  if (log) {
    los[los < log_lo] = log_lo
    his[his < log_lo] = log_lo
    vals[vals < log_lo] = log_lo
  }
  # use the block.xlim to re-position our x values.
  lo_xs = mids - 0.5 + (block.xlim[1])
  hi_xs = mids - 0.5 + (block.xlim[2])
  
  colors = rep(gray(0.7),maxr)

  #print(slr$upper)

  if (colorByType) {
    colors[slr$type == "positive1"] = rgb(.6,.5,.5)
    colors[slr$type == "positive2"] = rgb(.8,.4,.4)
    colors[slr$type == "positive3"] = rgb(.9,.3,.3)
    colors[slr$type == "positive4"] = rgb(1,0,0)
    colors[slr$type == "negative1"] = rgb(.4,.4,.8)
    colors[slr$type == "negative2"] = rgb(.5,.5,.6)
    colors[slr$type == "negative3"] = rgb(.3,.3,.9)
    colors[slr$type == "negative4"] = rgb(0,0,1)
  } else {
    colors[slr$lower > 1] = rgb(1,0,0)
    colors[slr$upper < 1] = rgb(0,0,1)
  }

  if (color.lrt) {
    cols = slr$lrt_stat
    my.cols = c()
    max = 20;
    cols[cols > max] = max;
    cols[cols < -max] = -max;
    abv = cols >= 0;
    blw = cols < 0;
    my.cols[abv] = hsv(0,cols[abv]/max,.9-cols[abv]/max/2);
    my.cols[blw] = hsv(2/3,-cols[blw]/max,.9+cols[blw]/max/2);
    colors = my.cols;
  }

  rect(xleft=lo_xs,ybottom=los,xright=hi_xs,ytop=his,
       col=colors, border=NA
  )

  lines(c(0,max(mids)),c(1,1),col='black',lwd=2);

  barHeight = 0.05
  rect(xleft=lo_xs,ybottom=vals,xright=hi_xs,ytop=vals,
       col='black',border='black'
  )
}


plot.slr = function(slr,
  overlay=FALSE
)
{
  if (is.character(slr)) {
    slr = read.table(slr,header=TRUE)
  }

  names(slr)
  print(slr[1,])
  minr = 1
  maxr = nrow(slr)

  par(mar=c(2.5,5,1,0.5))

  if (!overlay) {
    plot(NULL,xlim=c(minr,maxr),ylim=c(0,5),xlab="Site", ylab="Omega")
  }

  # Omega values.
  points(minr:maxr, slr[minr:maxr,]$omega, pch=19, xlim=c(minr,maxr), ylim=c(0,5),
)

  # Error bars
  plot.errbars = function(x,y1,y2) {
    npt = length(x)
    for ( i in 1:npt) {
      lcol = 'black'
      if (y2[i]<1){ lcol = 'blue'}
      if (y1[i]>1){ lcol = 'red'}
      lines( c(x[i],x[i]) , c(y1[i],y2[i]) , col=lcol, lwd=1.5)
    }
  }
  plot.errbars(minr:maxr, slr[minr:maxr,]$lower, slr[minr:maxr,]$upper)

  # "rolling average"
  lines(minr:maxr,smooth(slr[minr:maxr,]$omega),col='green')

  # Neutral line.
  abline(h=1,col='orange')
}

color.slr = function(slr,
  constrained=rgb(0,0,1,1),
  positive=rgb(1,0,0,1),
  normal=rgb(1,1,1,0.8)
){

  hue.map = function(x,hueLo,hueHi,valLo,valHi,...) {
    new.val = map.range(x,hueLo,hueHi,valLo,valHi)
    return(hcl(h=new.val,...))
  }
  
  # Fill all sites with the "dull" color.
  slr.colors = rep(normal,max(slr$site)-min(slr$site))

  # Sites confidently above 1 are bright red.
  slr.colors[slr$site] = ifelse(slr$lower > 1,positive,slr.colors)

  # Map the rest of sites 
  slr.colors[slr$site] = ifelse(slr$omega < 1, hue.map(slr$omega,240,360,0,1,l=50,c=90),slr.colors)

  return(slr.colors)
}

map.range = function(x,newLo,newHi,origLo=NULL,origHi=NULL,within.bounds=TRUE) {
  if (is.null(origLo)) origLo = min(x)
  if (is.null(origHi)) origHi = max(x)
  if (within.bounds) {
    x[x < origLo] = origLo
    x[x > origHi] = origHi
  }
  x = x - origLo
  x = x * (newHi-newLo)/(origHi-origLo)
  x = x + newLo
  return(x)
}
