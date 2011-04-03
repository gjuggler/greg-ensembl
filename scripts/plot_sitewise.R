constrain = function(data,lo,hi) {
  data = pmin(data,hi)
  data = pmax(data,lo)
  return(data)
}

plot.cum = function(df, new.plot=T, x.lim=NULL, y.lim=NULL, ...) {
  require(doBy)

  df$mid <- (df$right-df$left)/2.0
  df$signed_lrt = df$lrt_stat
  df[df$omega < 1,]$signed_lrt = -df[df$omega < 1,]$signed_lrt
  df <- orderBy(~signed_lrt,data=df)
  df$index <- (1:nrow(df))/nrow(df)
  df$index_lo <- (1:nrow(df))/nrow(df) - 1/nrow(df)
  print(df[1,])
  if (is.null(x.lim)) {
    x.lim = c(min(df$left),max(df$right))
  }
  if (is.null(y.lim)) {
    y.lim = c(0,1)
  }

  if (new.plot) {
    plot.new()
    plot(NULL,xlim=x.lim,ylim=y.lim,xlab='omega',ylab='cumulative density')
  }

  X = 6
  Y = 1/3
  lrt = df$signed_lrt
  lrt.blue = -constrain(lrt,-X,0)
  lrt.red = constrain(lrt,0,X)
  lrt.blue = (lrt.blue/X) ^ (Y)
  lrt.red = (lrt.red/X) ^ (Y)

  #W = 0.75
  #Z = df$lrt_stat < 6.6
  #lrt.blue[Z] = lrt.blue[Z]*W
  #lrt.red[Z] = lrt.red[Z]*W

  #V = 0.75
  #Z = df$lrt_stat < 3.8
  #lrt.blue[Z] = lrt.blue[Z]*V
  #lrt.red[Z] = lrt.red[Z]*V

  sig.lines <- function(thresh,...) {
    a <- subset(df,signed_lrt < thresh)
    abline(h=a[nrow(a),]$index,col='gray',...)
  }
  sig.lines(-6.6)
  sig.lines(-3.8,lty='dashed')
  sig.lines(3.8,lty='dashed')
  sig.lines(6.6)

  rect(xleft=df$left,xright=df$right,ybottom=df$index_lo,ytop=df$index,
    col=rgb(lrt.red,0,lrt.blue,.5),border=NA,...)
}

#plot.cum(df,x.lim=c(0,5))

plot.fn = function(ft,col) {
  a <- ft$estimate[1]
  b <- ft$estimate[2]
  xx <- seq(from=0,to=10,length.out=200)
  str <- paste("p",fn,"(xx,a,b)",sep="")
  yy <- eval(parse(text=str))
  lines(xx,yy,col=col,lwd=1.5)
}

#fit.fn(df,'lnorm','orange')
#fit.fn(df,'gamma','green')
#fit.fn(df,'weibull','purple')

ggplot.gene <- function(df,field) {
  df <- sign.lrt(df)
  df$x_pos <- df[,field]
#  max_val <- max(abs(df$signed_lrt))
  max_val <- 20
  values = c(-max_val,-5,5,max_val)

  # Get colors.
  cols <- df$signed_lrt
  my.cols <- c()
  ###
  max <- 15
  offset <- 0.9
  rate <- 0.5
  ###
  cols[cols > max] <- max
  cols[cols < -max] <- -max
  abv <- cols >= 0
  blw <- cols < 0
  my.cols[abv] <- hsv(0,cols[abv]/max,offset-cols[abv]/max * rate)
  my.cols[blw] <- hsv(2/3,-cols[blw]/max,offset+cols[blw]/max * rate)
  colors <- my.cols
  df$colors <- colors
  # /end Get colors.

  p <- ggplot(df)
  p <- p + geom_rect(aes(xmin=x_pos-0.5,xmax=x_pos+0.5,ymin=omega_lower+0.01,ymax=omega_upper+0.01,fill=colors))
  p <- p + scale_x_continuous(name="Alignment Position")
  p <- p + scale_fill_identity()
  p <- p + geom_hline(yintercept=1,colour='black',linetype='dashed',size=1,alpha=0.7)
  p <- p + scale_y_log10(name="dN/dS (log scale)")
  return(p)
}
