#############################################################################
#   Copyright (c) 2009 Marie Laure Delignette-Muller, Regis Pouillot, Jean-Baptiste Denis                                                                                                  
#                                                                                                                                                                        
#   This program is free software; you can redistribute it and/or modify                                               
#   it under the terms of the GNU General Public License as published by                                         
#   the Free Software Foundation; either version 2 of the License, or                                                   
#   (at your option) any later version.                                                                                                            
#                                                                                                                                                                         
#   This program is distributed in the hope that it will be useful,                                                             
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                                          
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                 
#   GNU General Public License for more details.                                                                                    
#                                                                                                                                                                         
#   You should have received a copy of the GNU General Public License                                           
#   along with this program; if not, write to the                                                                                           
#   Free Software Foundation, Inc.,                                                                                                              
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                                                             
#                                                                                                                                                                         
#############################################################################
### maximum likelihood estimation for censored or non-censored data
###
###         R functions
### 
### many ideas are taken from the fitdistr function of the MASS package and 
### the mle function of the stat package.

mledist <- function (data, distr, start=NULL, fix.arg=NULL, optim.method="default",
    lower=-Inf, upper=Inf, custom.optim=NULL, ...)
    # data may correspond to a vector for non censored data or to
    # a dataframe of two columns named left and right for censored data 
{
    if (!is.character(distr)) 
#        distname <- substring(as.character(match.call()$distr), 2)
    stop("distr must be a character string naming a distribution")
    else 
        distname <- distr
    ddistname <- paste("d",distname,sep="")
    
    if (!exists(ddistname, mode="function"))
        stop(paste("The ", ddistname, " function must be defined"))
    if (distname == "unif")
        stop("Maximum likelihood estimation is not available for the uniform distribution")

    if (!is.null(fix.arg) & is.null(start))
        stop("Starting values must be defined when some distribution parameters are fixed")
    
    if (is.vector(data)) {
        cens <- FALSE
        if (!(is.numeric(data) & length(data)>1)) 
            stop("data must be a numeric vector of length greater than 1 for non censored data
            or a dataframe with two columns named left and right and more than one line for censored data")
    }
    else {
        cens <- TRUE
        censdata <- data
        if (!(is.vector(censdata$left) & is.vector(censdata$right) & length(censdata[,1])>1))
        stop("data must be a numeric vector of length greater than 1 for non censored data
        or a dataframe with two columns named left and right and more than one line for censored data")
        pdistname<-paste("p",distname,sep="")
        if (!exists(pdistname,mode="function"))
            stop(paste("The ",pdistname," function must be defined to apply maximum likelihood to censored data"))

    }
    
    if (cens) {
        # Definition of datasets lcens (left censored)=vector, rcens (right censored)= vector,
        #   icens (interval censored) = dataframe with left and right 
        # and ncens (not censored) = vector
        lcens<-censdata[is.na(censdata$left),]$right
        if (any(is.na(lcens)) )
            stop("An observation cannot be both right and left censored, coded with two NA values")
        rcens<-censdata[is.na(censdata$right),]$left
        ncens<-censdata[censdata$left==censdata$right & !is.na(censdata$left) & 
            !is.na(censdata$right),]$left
        icens<-censdata[censdata$left!=censdata$right & !is.na(censdata$left) & 
            !is.na(censdata$right),]
        # Definition of a data set for calculation of starting values
        data<-c(rcens,lcens,ncens,(icens$left+icens$right)/2)
    }
    
    # MLE fit 
    # definition of starting values if not previously defined
    if (is.null(start)) {
        if (distname == "norm") {
            n <- length(data)
            sd0 <- sqrt((n - 1)/n) * sd(data)
            mx <- mean(data)
            start <- list(mean=mx, sd=sd0)
        }
        if (distname == "lnorm") {
            if (any(data <= 0)) 
                stop("values must be positive to fit a lognormal distribution")
            n <- length(data)
            ldata <- log(data)
            sd0 <- sqrt((n - 1)/n) * sd(ldata)
            ml <- mean(ldata)
            start <- list(meanlog=ml, sdlog=sd0)
        }
        if (distname == "pois") {
            start <- list(lambda=mean(data))
        }
        if (distname == "exp") {
            start <- list(rate=1/mean(data))
        }
        if (distname == "gamma") {
            n <- length(data)
            m <- mean(data)
            v <- (n - 1)/n*var(data)
            start <- list(shape=m^2/v,rate=m/v)
        }
        if (distname == "nbinom") {
            n <- length(data)
            m <- mean(data)
            v <- (n - 1)/n*var(data)
            size <- if (v > m) m^2/(v - m)
                else 100
            start <- list(size = size, mu = m) 
        }
        if (distname == "geom" ) {
            m <- mean(data)
            prob <- if (m>0) 1/(1+m)
                    else 1
            start <- list(prob=prob)        
        }
        if (distname == "beta") {
            if (any(data < 0) | any(data > 1)) 
                stop("values must be in [0-1] to fit a beta distribution")
            n <- length(data)
            m <- mean(data)
            v <- (n - 1)/n*var(data)
            aux <- m*(1-m)/v - 1
            start <- list(shape1=m*aux,shape2=(1-m)*aux)
        }
        if (distname == "weibull") {
            m <- mean(log(data))
            v <- var(log(data))
            shape <- 1.2/sqrt(v)
            scale <- exp(m + 0.572/shape)
            start <- list(shape = shape, scale = scale)
        }
        if (distname == "logis") {
            n <- length(data)
            m <- mean(data)
            v <- (n - 1)/n*var(data)
            start <- list(location=m,scale=sqrt(3*v)/pi)
        }
        if (distname == "cauchy") {
            start <- list(location=median(data),scale=IQR(data)/2)
        }
        if (!is.list(start)) 
            stop("'start' must be defined as a named list for this distribution") 
   } # end of the definition of starting values 
   
   ############# MLE fit using optim or custom.optim ##########
    vstart <- unlist(start)
    vfix.arg <- unlist(fix.arg)
    # check of the names of the arguments of the density function
    argddistname <- names(formals(ddistname))   
    m <- match(names(start), argddistname)
    mfix <- match(names(vfix.arg), argddistname)
    if (any(is.na(m)))
        stop("'start' must specify names which are arguments to 'distr'")
    if (any(is.na(mfix)))
        stop("'fix.arg' must specify names which are arguments to 'distr'")
    # check that some parameters are not both in fix.arg and start
    minter <- match(names(start), names(fix.arg))
    if (any(!is.na(minter)))
        stop("a distribution parameter cannot be specified both in 'start' and 'fix.arg'")

    # definition of the function to minimize : - log likelihood
    # for non censored data
    if (!cens) {
        # the argument names are:
        # - par for parameters (like in optim function)
        # - fix.arg for optional fixed parameters
        # - obs for observations (previously dat but conflicts with genoud data.type.int argument)
        # - ddistnam for distribution name
        if ("log" %in% argddistname){
            fnobj <- function(par, fix.arg, obs, ddistnam){
                -sum(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg), log=TRUE) ) )
            }
        }
        else{
        fnobj <- function(par, fix.arg, obs, ddistnam) {
            -sum(log(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) )
            }
        }
    }
    else {# if !cens
        argpdistname<-names(formals(pdistname))
        if (("log" %in% argddistname) & ("log.p" %in% argpdistname))
            fnobjcens <- function(par,fix.arg,rcens,lcens,icens,ncens,ddistnam,pdistnam)
                -sum(do.call(ddistnam,c(list(x=ncens),as.list(par),as.list(fix.arg),list(log=TRUE)))) -
                sum(do.call(pdistnam,c(list(q=lcens),as.list(par),as.list(fix.arg),list(log=TRUE)))) -
                sum(do.call(pdistnam,c(list(q=rcens),as.list(par),as.list(fix.arg),list(lower.tail=FALSE),list(log=TRUE)))) -
                sum(log(do.call(pdistnam,c(list(q=icens$right),as.list(par),as.list(fix.arg))) - # without log=TRUE here
                do.call(pdistnam,c(list(q=icens$left),as.list(par),as.list(fix.arg))) )) # without log=TRUE here
        else
            fnobjcens <- function(par,fix.arg, rcens,lcens,icens,ncens,ddistnam,pdistnam) {
                out.val <- -sum(log(do.call(ddistnam,c(list(x=ncens),as.list(par),as.list(fix.arg))))) -
                sum(log(do.call(pdistnam,c(list(q=lcens),as.list(par),as.list(fix.arg))))) -
                sum(log(1-do.call(pdistnam,c(list(q=rcens),as.list(par),as.list(fix.arg))))) -
                sum(log( 
                  do.call(pdistnam,c(list(q=icens$right),as.list(par),as.list(fix.arg))) - 
                  do.call(pdistnam,c(list(q=icens$left),as.list(par),as.list(fix.arg)))
                ))
                print(out.val)
                out.val
                }
    }
    # Choice of the optimization method    
    if (optim.method == "default")
    {
        if(is.infinite(lower) && is.infinite(upper))
        { 
            if (length(vstart) > 1) 
                meth <- "Nelder-Mead"
            else 
                meth <- "BFGS"
        }else
            meth <- "L-BFGS-B"
    }else
        meth <- optim.method
        
    # Try to minimize the minus (log-)likelihood using the base R optim function
    if(is.null(custom.optim))
    {

    #print(lower)
    #print(upper)
    #print(meth)
    #print(fix.arg)
    #print(vstart)

        if (!cens) {
            opt <- optim(par=vstart, fn=fnobj, fix.arg=fix.arg, obs=data, ddistnam=ddistname,
            hessian=TRUE, method=meth, lower=lower, upper=upper, ...)
       } else {
            opt<-optim(par=vstart,fn=fnobjcens, fix.arg=fix.arg, rcens=rcens,lcens=lcens,icens=icens,ncens=ncens,
            ddistnam=ddistname,pdistnam=pdistname,hessian=TRUE,
            method=meth,lower=lower,upper=upper,...)
        }

        return(list(estimate = opt$par, convergence = opt$convergence, loglik = -opt$value, 
                    hessian = opt$hessian, optim.function="optim"))  
        
    }
    else # Try to minimize the minus (log-)likelihood using a user-supplied optim function 
    {
        if (!cens) {
          opt <- custom.optim(par=vstart, fn=fnobj, fix.arg=fix.arg, obs=data, ddistnam=ddistname, ...)
        } else {
            if (is.finite(lower)) {
              opt <- custom.optim(par=vstart,
              fn=fnobjcens,fix.arg=fix.arg,
              rcens=rcens,lcens=lcens,icens=icens,ncens=ncens,
              ddistnam=ddistname,pdistnam=pdistname,lower=lower,upper=upper,...)
            } else {
              opt <- custom.optim(par=vstart,
              fn=fnobjcens,fix.arg=fix.arg,
              rcens=rcens,lcens=lcens,icens=icens,ncens=ncens,
              ddistnam=ddistname,pdistnam=pdistname,...)
            }
        }
        return(list(estimate = opt$par, convergence = opt$convergence, loglik = -opt$value, 
                    hessian = NA, optim.function=custom.optim))
    }   
        
     
}

## old function with previous name for censored data
mledistcens<-function (censdata, distr, start=NULL,optim.method="default",lower=-Inf,upper=Inf)
{
    stop("The function \"mledistcens\" is no more used. Now the same function \"mledist\" must be used for censored and non censored data.")
}
