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
color.column <- function(xt, column, limits=range(xt[, column], na.rm=T), low='white', high=rgb(0.3, 0.3, 1)) {
  ## Make sure the colorspace library is loaded to convert rgb to hex.
  library(colorspace)
  color.f <- colorRamp(c(low, high))
  vals <- xt[, column]
  ## Rescale the values from (limits[1], limits[2]) to (0, 1)
  vals <- (vals - limits[1])/diff(limits) * diff(c(0,1)) + 0
  ## Clip values outside the desired range.
  vals <- ifelse(!is.finite(vals) | vals %inside% c(0,1), vals, NA)

  clrs <- color.f(vals)
  clr.string <- hex(RGB(clrs/255))
  clr.string <- substring(clr.string, 2) # Remove the hash prefix.

  clr.command <- paste('\\cellcolor[HTML]{', clr.string, '} ', sep='')
  xt <- xt.surround(xt, column, prefix=clr.command)
  xt
}

xt.surround <- function(xt, column, prefix='', suffix='') {
  vals <- xt[, column]
  indx <- which(colnames(xt) == column)
  # Grab the digits and display values set by xtable
  digt <- attr(xt, 'digits')[indx+1]
  disply <- attr(xt, 'display')[indx+1]
  print(digt)
  print(disply)
  formatted.vals <- formatC(xt[, column], digits=digt, format=disply)
  formatted.column <- paste(prefix, trim(formatted.vals), suffix, sep='')
  print(formatted.column)
  xt[, column] <- formatted.column
  xt
}

xt.tiny <- function(xt, column) {
  xt.surround(xt, column, prefix='\\tiny{', suffix='}')
}
xt.small <- function(xt, column) {
  xt.surround(xt, column, prefix='\\small{', suffix='}')
}
