g.plotdistcens = function (censdata, distr, para, leftNA = -Inf, rightNA = Inf, 
x.lim=NULL, y.lim=NULL,
    ...) 
{
    def.par <- par(no.readonly = TRUE)
    if (missing(censdata) || !(is.vector(censdata$left) & is.vector(censdata$right) & 
        length(censdata[, 1]) > 1)) 
        stop("datacens must be a dataframe with two columns named left \n            and right and more than one line")
    if ((missing(distr) & !missing(para)) || (missing(distr) & 
        !missing(para))) 
        stop("distr and para must defined")
    if (is.finite(leftNA) & any(is.na(censdata$left))) 
        censdata[is.na(censdata$left), ]$left <- leftNA
    if (is.finite(rightNA) & any(is.na(censdata$right))) 
        censdata[is.na(censdata$right), ]$right <- rightNA
    lcens <- censdata[is.na(censdata$left), ]$right
    if (any(is.na(lcens))) 
        stop("An observation cannot be both right and left censored, coded with two NA values")
    rcens <- censdata[is.na(censdata$right), ]$left
    noricens <- censdata[!is.na(censdata$left) & !is.na(censdata$right), 
        ]
    midnoricens <- (noricens$left + noricens$right)/2
    ordmid <- order(midnoricens)
    ordlcens <- order(lcens)
    ordrcens <- order(rcens)
    nlcens <- length(lcens)
    nrcens <- length(rcens)
    nnoricens <- length(noricens$left)
    n <- length(censdata$left)
    xminright <- min(censdata[!is.na(censdata$right), ]$right)
    xminleft <- min(censdata[!is.na(censdata$left), ]$left)
    xmin = min(xminright, xminleft)
    xmaxright <- max(censdata[!is.na(censdata$right), ]$right)
    xmaxleft <- max(censdata[!is.na(censdata$left), ]$left)
    xmax = max(xmaxright, xmaxleft)
    xrange <- xmax - xmin
    xmin <- xmin - 0.3 * xrange
    xmax <- xmax + 0.3 * xrange
    xlim <- c(xmin, xmax)
    if (!is.null(x.lim)) {
      xlim <- x.lim
    }

    plot(c(0, 0), c(0, 0), type = "n", xlim = xlim, ylim = c(0, 
        1), xlab = "censored data", ylab = "CDF", main = "Cumulative distribution")
    plotlcens <- function(i) {
        y <- i/n
        lines(c(xmin, lcens[ordlcens[i]]), c(y, y))
    }
    if (nlcens >= 1) 
        toto <- sapply(1:nlcens, plotlcens)
    plotnoricens <- function(i) {
        y <- (i + nlcens)/n
        if (noricens[ordmid[i], ]$left != noricens[ordmid[i], 
            ]$right) 
            lines(c(noricens[ordmid[i], ]$left, noricens[ordmid[i], 
                ]$right), c(y, y))
        else points(noricens[ordmid[i], ]$left, y, pch = 4)
    }
    if (nnoricens >= 1) 
        toto <- sapply(1:nnoricens, plotnoricens)
    plotrcens <- function(i) {
        y <- (i + nlcens + nnoricens)/n
        lines(c(rcens[ordrcens[i]], xmax), c(y, y))
    }
    if (nrcens >= 1) 
        toto <- sapply(1:nrcens, plotrcens)
    if (!missing(distr)) {
        if (!is.character(distr)) 
            distname <- substring(as.character(match.call()$distr), 
                2)
        else distname <- distr
        if (!is.list(para)) 
            stop("'para' must be a named list")
        ddistname <- paste("d", distname, sep = "")
        if (!exists(ddistname, mode = "function")) 
            stop(paste("The ", ddistname, " function must be defined"))
        pdistname <- paste("p", distname, sep = "")
        if (!exists(pdistname, mode = "function")) 
            stop(paste("The ", pdistname, " function must be defined"))
        densfun <- get(ddistname, mode = "function")
        nm <- names(para)
        f <- formals(densfun)
        args <- names(f)
        m <- match(nm, args)
        if (any(is.na(m))) 
            stop(paste("'para' specifies names which are not arguments to ", 
                ddistname))
        s <- seq(xmin, xmax, by = (xmax - xmin)/100)
        theop <- do.call(pdistname, c(list(q = s), as.list(para)))
        lines(s, theop, lty = 1)
    }
    #par(def.par)
}
