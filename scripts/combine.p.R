# combine.p --- A short and ugly function to combine probabilities
# from independent tests of significance using either the Fisher,
# Winer or Stouffer's methods
# see:
# Fisher (1958, section 21.1 --- pp 99--101)
# Sokal & Rohlf (1995, section 18.1, pp 794--797)
# Wolf (1986, pp 18--23)
#
# To do:
# 0) fix winer & stouffer to handle signed cases, this is broken!
# 1) add bulletproofing for y stuff, as it stands now, winer
#    will accept all sorts of crappy input
# 2) Beautify using switch on method, or not...
#
# Pete Hurd June 20 2002

combine.p <- function(x, y = NULL, method="fisher") {
  x <- as.vector(na.omit(x))
  if (any(x <= 0) || any(is.na(x)) || any(x>1))
    stop("all entries of x must be in the range (0,1]")
    
  if( method=="fisher"){
    STATISTIC <- 2*(sum(-1*log(x)))
    PARAMETER <- 2*length(x)
    PVAL <- pchisq(STATISTIC, PARAMETER, lower = FALSE)
    METHOD <- "Fisher's combined probability"
    names(STATISTIC) <- "Chi Square"
    names(PARAMETER) <- "X-squared df"
    names(PVAL) <- "p.value"
  }
  else if (method=="winer"){
    t <- qt(x,y,lower.tail=FALSE)
    STATISTIC <- sum(t) / (sqrt(sum(dv/(dv-2))))
    PARAMETER <- NA
    PVAL <-pnorm(STATISTIC,lower.tail=FALSE)
    METHOD <- "Winer's combined probability"
    names(STATISTIC) <- "Z score"
    names(PARAMETER) <- "No df"
    names(PVAL) <- "p.value"    
  }
  else if (method=="stouffer"){
    STATISTIC <- sum(qnorm(x,lower.tail=FALSE)) / sqrt(length(x))
    PARAMETER <- NA
    PVAL <-pnorm(STATISTIC,lower.tail=FALSE)
    METHOD <- "Stouffer's combined probability"
    names(STATISTIC) <- "Z score"
    names(PARAMETER) <- "No df"
    names(PVAL) <- "p.value"
  }
  
  structure(list(statistic=STATISTIC, parameter=PARAMETER,
                 p.value=PVAL,method=METHOD), class="htest")
}
  

