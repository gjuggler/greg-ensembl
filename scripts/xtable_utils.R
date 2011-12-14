library(xtable)

# color.column - color the columns of a table for LaTeX output.
#
# Takes an xtable object and a column name as input, and returns the
# xtable with the chosen column stringified and prepended with color
# values using the \cellcolor command from the xcolor package. For an
# easy copy-paste into a LaTeX document, use xtable's print function
# with the following parameter to preserve the \cellcolor command:
# "sanitize.text.function=function(x){x}". I also like to add
# "include.rownames=F". Here's a worked example:
#
# data(iris)
# xt <- xtable(iris[1:10,])
# xt.clr <- color.column(xt, 'Sepal.Length', low=rgb(1,1,1), high='red')
# print(xt.clr, include.rownames=F, sanitize.text.function=function(x){x})
# # Prints a table with rows like: \cellcolor[HTML]{FF8E8E}5.10 & 3.50 & 1.40 & 0.20 & setosa \\
#
# Within your LaTeX preamble, be sure to add the xcolor package with a
# 'table' argument and you're good to go:
#
# \usepackage[table]{xcolor}
#
color.column <- function(xt, column, 
  limits=range(xt[, column], na.rm=T), 
  low='white', 
  color.bg=T,
  high=rgb(0.6, 0.6, 1),
  skip.coloring=F,
  log = F,
  ...
) {
  ## Make sure the colorspace library is loaded to convert rgb to hex.
  library(colorspace)
  color.f <- colorRamp(c(low, high))
  vals <- xt[, column]
  if (log) {
    min.val <- min(vals)
    if (min.val <= 0) {
      vals <- vals - min.val + 0.01
    }
    vals <- log(vals)
    limits <- range(vals, na.rm=T)
  }
  if (diff(limits) > 0) {
    ## Rescale the values from (limits[1], limits[2]) to (0, 1)
    vals <- (vals - limits[1])/diff(limits) * diff(c(0,1)) + 0
    ## Clip values outside the desired range.
    vals <- ifelse(!is.finite(vals) | vals < 1, vals, 1)
    vals <- ifelse(!is.finite(vals) | vals > 0, vals, 0)
  } else {
    # All values are the same -- put them all at zero.
    vals <- 0
  }

  clrs <- color.f(vals)
  clrs[is.na(clrs)] <- color.f(0)
  clr.string <- hex(RGB(clrs/255))
  clr.string <- substring(clr.string, 2) # Remove the hash prefix.

  if (color.bg) {
    clr.command <- paste('\\cellcolor[HTML]{', clr.string, '}', sep='')
  } else {
    clr.command <- paste('\\color[HTML]{', clr.string, '}', sep='')
  }

  if (skip.coloring) {
    clr.command <- ''
  }

  xt <- xt.surround(xt, column=column, prefix=clr.command, ...)
  xt
}

xt.surround <- function(xt, column, rows=c(1:nrow(xt)), prefix='', suffix='', na.string='-') {
  vals <- xt[, column]
  indx <- which(colnames(xt) == column)
  # Grab the digits and display values set by xtable
  digt <- attr(xt, 'digits')[indx+1]
  disply <- attr(xt, 'display')[indx+1]

  if (any(is.na(vals))) {
    vals[is.na(vals)] <- na.string
  }

  formatted.vals <- R.oo::trim(formatC(vals, digits=digt, format=disply))
  formatted.column <- paste(prefix, formatted.vals[rows], suffix, sep='')
  vals <- formatted.vals
  vals[rows] <- formatted.column

  xt[, column] <- vals
  xt
}

color.columns <- function(xt, columns,  ...) {
  for (i in 1:length(columns)) {
    #print(columns[i])
    xt <- color.column(xt, columns[i], ...)
  }
  xt
}

bold.t <- function(xt, column, t=0.05, dir='below', ...) {
  vals <- xt[, column]
  if (dir == 'below') {
    rws <- which(vals < t)
  } else {
    rws <- which(vals > t)
  }

  xt.surround(xt, column=column, rows=rws, prefix='\\bf{', suffix='}')
}

color.rows <- function(xt, rows=c(1:nrow(xt)), columns=colnames(xt), color.bg=T, color=gray(0.5), ...) {
  clr.string <- substring(color, 2) # Remove the hash prefix.
  if (color.bg) {
    clr.command <- paste('\\cellcolor[HTML]{', clr.string, '} ', sep='')
  } else {
    clr.command <- paste('\\color[HTML]{', clr.string, '} ', sep='')
  }

  for (cl in columns) {
    xt <- xt.surround(xt, column=cl, prefix=clr.command, rows=rows, ...)
  }
  xt
}

xt.tiny <- function(xt, column) {
  xt.surround(xt, column, prefix='\\tiny{', suffix='}')
}
xt.small <- function(xt, column) {
  xt.surround(xt, column, prefix='\\small{', suffix='}')
}

print.latex <- function(xt, filename, ...) {
  print(xt, 
    file=filename,
    only.contents=T,
    include.colnames=F,
    include.rownames=F,
    sanitize.text.function=function(x){x},
    hline.after=NULL,
    ...
  )
  print(sprintf("Wrote LaTeX to '%s'", filename))
}
