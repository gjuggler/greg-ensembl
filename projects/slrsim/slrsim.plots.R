library(ape)

subplot <- function(x, y) viewport(layout.pos.col=x, layout.pos.row=y)
vplayout <- function(x, y) {
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(y,x)))
}

plot.indel.distribution <- function(indel_label) {
  library(VGAM)
  toks <- unlist(strsplit(indel_label,' '))
  if (toks[1] == 'NB') {
    # for NB .5 1, a=1 and b=1 - .5
    distr <- dnbinom
    a <- toks[3]
    b <- 1-as.numeric(toks[2])
  } else {
    # for POW 4 30, a=30 and b=4
    distr <- dzipf
    a <- toks[3]
    b <- toks[2]
  }
  a <- as.numeric(a)
  b <- as.numeric(b)  

  xs <- seq(from=1,to=20)
  ys <- distr(xs,a,b)
  df <- data.frame(x=xs,y=ys)
  p <- ggplot(df,aes(x=x,y=y)) + geom_bar(stat='identity')
  return(p)
}

plot.roc = function(
  data,
  plot.x = 'tn',
  plot.y = 'tp',
  color.by = 'slrsim_label',
  plot.at.threshold = NA,
  alpha = 1,
  plot = T
) {
  data$roc.x = data[, plot.x]
  data$roc.y = data[, plot.y]
  data$roc.color = data[, color.by]

  require(ggplot2)
  require(grid)
  #print(paste("Plot.roc row count: ",nrow(data)))

  p <- ggplot(data,aes(x=roc.x,y=roc.y,colour=roc.color))
  p <- p + geom_line(alpha=alpha)
  p <- p + xlab(plot.x)
  p <- p + ylab(plot.y)

  if (any(!is.na(plot.at.threshold))) {
    for (lbl in unique(data$slrsim_label)) {
      data.sub <- subset(data,slrsim_label==lbl)

      for (threshold in plot.at.threshold) {
        t <- threshold
        row.index <- min(which(data.sub$score <= t))
        row <- data.sub[row.index,]
        #print(row)
        p <- p + geom_point(data=row,colour=I('black'),shape=I(1)) + scale_shape(solid=FALSE)
      }
    }
  }

  if (plot) {
    print(p)
  }
  return(p)
}
