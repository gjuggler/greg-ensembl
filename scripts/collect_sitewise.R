if (exists('drv')) {
  lapply(dbListConnections(drv),dbDisconnect)
} else {
  library(RMySQL)
  drv <- dbDriver('MySQL')
}

if (!exists('dbname')) {
  dbname = 'gj1_eslr_57'
}
#print(paste("dbname is: ",dbname))
connect.db = function(dbname) {
  con <- dbConnect(drv, host='ens-research', port=3306, user='ensro', password='', dbname=dbname)
  return(con)
}
con = connect.db(dbname)

get.vector = function(con,query,columns='all') {
  res = dbSendQuery(con,query)
  if (columns == 'all') {
    data = fetch(res,n=-1)  # Important: the n=-1 causes the entire resultset to be fetched, not just the first 500.
  } else {
    data = fetch(res,n=-1)[[columns]]
  }
  dbClearResult(res)
  return(data)
}

factorize = function(data,columns=names(data)) {
  data[columns] = lapply(data[columns],as.factor)
  return(data)
}

# Get the parameter sets.
get.psets = function(db=dbname) {
  query = sprintf(
    paste('SELECT n.parameter_set_id AS id,sn.parameter_value AS parameter_set_shortname, n.parameter_value AS name FROM %s.parameter_set n, %s.parameter_set sn',
      'WHERE n.parameter_set_id=sn.parameter_set_id AND n.parameter_name="name" AND sn.parameter_name="parameter_set_shortname";'),
    db,db)
  return(get.vector(con,query,columns='all'))
}

get.genes = function(parameter.set.id=1,db=dbname) {
  query = sprintf("SELECT * FROM %s.stats_genes where parameter_set_id=%s",db,parameter.set.id);
  data = get.vector(con,query,columns='all')

  return(data)
}

get.genes.list = function(db=dbname) {
  query = sprintf("SELECT * FROM %s.stats_genes",db);
  data = get.vector(con,query,columns='all')
  return(data)
}

get.genes.merged = function(db=dbname,exclude.cols=NULL,pset.cols=NULL) {
  pset.cols <- c(pset.cols,
    'slr.dnds','hyphy.dnds'
#'hyphy_dnds_lo','hyphy_dnds_hi',
#    'positive_1','positive_2',
#    'tree_length', 'tree_max_path'
  )
  pset.cols <- unique(pset.cols)

  if (is.null(exclude.cols)) {
    exclude.cols <- c('tree_newick',pset.cols)
  }

  param.sets = get.psets(db=db)
  print(param.sets)
  
  genes = get.genes(1,db=db)
  genes <- within(genes,rm(list=exclude.cols))
  
  for (i in 1:nrow(param.sets)) {
    pset = param.sets[i,]$id
    print(param.sets[i,]$parameter_set_shortname)
    cur.genes = get.genes(pset,db=db)
    print(str(cur.genes))
    create.name = function(ext) {paste(param.sets[pset,]$parameter_set_shortname,ext,sep='_')}

    genes.subset = data.frame(
      a=cur.genes$data_id # Need to make sure the data_id is here to merge back with the main genes vector.
    )

    for (i in 1:length(pset.cols)) {
      print(pset.cols[i])
      genes.subset[,paste('col.',i,sep="")] <- cur.genes[,pset.cols[i]]
    }

    col.names = c('data_id')
    for (i in 1:length(pset.cols)) {
      col.names = c(col.names,create.name(pset.cols[i]))
    }

    colnames(genes.subset) = col.names
    genes = merge(genes,genes.subset,all.x=T)
  }
  return(genes);
}

get.positive.sites = function(parameter.set.id=1,db=dbname,limit="") {
  query = sprintf("SELECT * from %s.stats_sites WHERE parameter_set_id=%s AND type IN ('positive1','positive2','positive3','positive4') %s",db,parameter.set.id,limit)
  sites = get.vector(con,query,columns='all')
  return(sites)  
}

get.sites = function(parameter.set.id=1,db=dbname,limit="") {
  query = sprintf("SELECT * from %s.stats_sites WHERE parameter_set_id=%s %s",db,parameter.set.id,limit)
  sites = get.vector(con,query,columns='all')
  return(sites)
}

get.sites.slim = function(fields.list=NA,parameter.set.id=NA,db=dbname,limit.string="") {
  if (is.na(fields.list)) {
    fields.list <- c('node_id','parameter_set_id','filter_value','domain','omega','note','type')
  }
  fields.string <- paste(fields.list,collapse=",")
  
  param.set.string <- ""
  if (!is.na(parameter.set.id)) {
    param.set.string <- paste("WHERE parameter_set_id =",parameter.set.id)
  } 

  query <- sprintf("SELECT %s from %s.stats_sites %s %s",fields.string,db,param.set.string,limit.string)
  sites <- get.vector(con,query,columns='all')
  return(sites)
}

get.sites.for.parameter.set = function(dir='.',parameter.set.id,filter.fn=NULL) {
  load(paste(dir,"/sites_",parameter.set.id,".Rdata",sep=""))
  if (!is.null(filter.fn)) {
    sites <- filter.fn(sites)
  }

  return(sites)
}

sites.summary = function(db=dbname, sites.dir='.', filter.fn=NULL) {
  param.sets = get.psets(db=db)

  df = data.frame()
#  for (i in c(5)) {
  for (i in 1:nrow(param.sets)) {
    pset <- param.sets[i,]$id
    name <- param.sets[i,]$parameter_set_shortname
    print(pset)
    print(name)

    sub.sites <- get.sites.for.parameter.set(dir=sites.dir,parameter.set.id=pset,filter.fn=filter.fn)
    
    n <- nrow(sub.sites)
    if (n == 0) {
      next
    }

    # Collect FDR threshold values from simulated datasets.
    for (thresh in c(0.01,0.05,0.1)) {
      #print(thresh)
      x <- subset(fdr.thresholds, parameter_set_name==name & fdr.threshold==thresh)
      #print(x)
      var.name <- paste('fdr.thresh.',thresh,sep="")
      if (nrow(x)>0) {
        assign(var.name,x[,'score'])
      } else {
        assign(var.name,1000)
      }
    }

    # Collect from indel-simulated datasets.
    for (thresh in c(.01,.05,.1)) {
      x <- subset(fdr.indel.thresholds,parameter_set_name==name & fdr.threshold==thresh)
      var.name <- paste('fdr.indel.thresh.',thresh,sep="")
      if (nrow(x)>0) {
        assign(var.name,x$score)
      } else {
        assign(var.name,1000)
      }
    }
    
    # Use p.adjust to do B-H multiple testing correction.
    sub.above.one = subset(sub.sites, omega > 1)
    pvals = sub.above.one$pval
    pvals.adj = p.adjust(pvals,method="fdr")
    sub.above.one[,'pval.adj'] = pvals.adj

    fdr.pos.sites <- subset(sub.sites, lrt_stat >= `fdr.thresh.0.05` & omega > 1)
    fdr.strong.pos.sites <- subset(sub.sites, lrt_stat >= `fdr.thresh.0.01` & omega > 1)
    n.pos.sites <- nrow(fdr.pos.sites)
    n.strong.pos.sites <- nrow(fdr.strong.pos.sites)
    n.neg.sites <- nrow(subset(sub.sites, type %in% c('negative1','negative2','negative3','negative4')))
    n.strong.neg.sites <- nrow(subset(sub.sites, type %in% c('negative2','negative3','negative4')))

    # Define positive domains and genes using the simulated FDR thresholds.
    pos.domains <- subset(fdr.strong.pos.sites,!is.na(domain))
    x <- pos.domains[,c('node_id','domain')]
    y <- pos.domains[,c('domain')]
    n.pos.domains <- nrow(unique(x))
    n.pos.domain.types <- length(unique(y))

    x <- fdr.strong.pos.sites[,c('node_id')]
    n.pos.genes <- length(unique(x))
    x <- sub.sites[,c('node_id')]
    n.total.genes <- length(unique(x))

    summary.df <- data.frame(
      parameter_set_id = pset,
      parameter_set_name = name,
      total = nrow(sub.sites),
      constant = nrow(subset(sub.sites,note=='constant')),
      synonymous = nrow(subset(sub.sites,note=='synonymous')),
      `omega_lt_p1` = nrow(subset(sub.sites, omega < 0.1)) / n,
      `omega_lt_p2` = nrow(subset(sub.sites, omega < 0.2)) / n,
      `omega_lt_p5` = nrow(subset(sub.sites, omega < 0.5)) / n,
      `omega_lt_one` = nrow(subset(sub.sites, omega < 1)) / n,
      `omega_gt_one` = nrow(subset(sub.sites, omega > 1)) / n,
      `omega_gt_onep5` = nrow(subset(sub.sites, omega > 1.5)) / n,
      `n.positive.domains` = n.pos.domains,
      `n.positive.domain.types` = n.pos.domain.types,
      `n.positive.genes` = n.pos.genes,
      `n.total.genes` = n.total.genes,
      `n.pval.01` = nrow(subset(sub.sites, pval <= 0.01 & omega > 1)),
      `n.pval.05` = nrow(subset(sub.sites, pval <= 0.05 & omega > 1)),
      `n.pval.1` = nrow(subset(sub.sites, pval <= 0.1 & omega > 1)),
      `n.adj.01` = nrow(subset(sub.above.one, pval.adj <= 0.01)),
      `n.adj.05` = nrow(subset(sub.above.one, pval.adj <= 0.05)),
      `n.adj.1` = nrow(subset(sub.above.one, pval.adj <= 0.1)),
      `n.fdr.01` = nrow(subset(sub.sites, lrt_stat >= `fdr.thresh.0.01` & omega > 1)),
      `n.fdr.05` = nrow(subset(sub.sites, lrt_stat >= `fdr.thresh.0.05` & omega > 1)),
      `n.fdr.1` = nrow(subset(sub.sites, lrt_stat >= `fdr.thresh.0.1` & omega > 1)),
      `n.fdr.indel.01` = nrow(subset(sub.sites, lrt_stat >= `fdr.indel.thresh.0.01` & omega > 1)),
      `n.fdr.indel.05` = nrow(subset(sub.sites, lrt_stat >= `fdr.indel.thresh.0.05` & omega > 1)),
      `n.fdr.indel.1` = nrow(subset(sub.sites, lrt_stat >= `fdr.indel.thresh.0.1` & omega > 1)),
      `fdr.thresh.01` = `fdr.thresh.0.01`,
      `fdr.thresh.05` = `fdr.thresh.0.05`,
      `fdr.thresh.1` = `fdr.thresh.0.1`,
      `fdr.indel.thresh.01` = `fdr.indel.thresh.0.01`,
      `fdr.indel.thresh.05` = `fdr.indel.thresh.0.05`,
      `fdr.indel.thresh.1` = `fdr.indel.thresh.0.1`
    )
    df <- rbind(df,summary.df)
    rm(sub.sites)
    gc()
  }

  print(df)
  return(df)
}


sites.fit.distributions = function(db=dbname, sites.dir='.', filter.fn=NULL, fit.type='values') {
  library(fitdistrplus)
  param.sets = get.psets(db=db)
  df = data.frame()
  for (i in c(1)) {
#  for (i in 1:nrow(param.sets)) {
    pset <- param.sets[i,]$id
    name <- param.sets[i,]$parameter_set_shortname
    print(name)
    sites <- get.sites.for.parameter.set(dir=sites.dir,parameter.set.id=pset,filter.fn=filter.fn)

    sites.adj <- sites
    sites.adj = sites.adj[sample(1:nrow(sites.adj),size=100000),]
    sites.adj[sites.adj$omega == 0,]$omega <- 0.0001
    omegas.adj <- sites.adj$omega
    sites.adj[,'left'] = sites.adj[,'omega_lower']    
    sites.adj[,'right'] = sites.adj[,'omega_upper']

    set.seed(0)
    omegas.below.one <- omegas.adj[omegas.adj < 1]
    param.min = sqrt(.Machine$double.eps)

    for (distr in c('weibull','beta','exp','gamma','lnorm')) {
      print(paste("Fit",distr,"..."))
      start <- NULL
      values <- omegas.adj
      data <- sites.adj
      lower <- -Inf
      if (distr == 'beta') {
        #data <- subset(data,right < 1);values <- omegas.below.one;
	start <- list(shape1=1,shape2=1); lower<-c(param.min,param.min)}
      if (distr == 'weibull') {start <- list(shape=0.5,scale=0.5); lower <- c(param.min,param.min)}
      if (distr == 'gamma') {lower=c(param.min,param.min)}
#      if (distr == 'lognormal') {start <- list(-5,4);}
      do.fit = function() {
        fit = NA
        if (fit.type == 'ci') {
	  # Use 'fitdistcens' from the Fitdistrplus package.
  	  library(fitdistrplus)
          if (!is.null(start)) {
            fit <- fitdistcens(data,distr,start=start)
          } else {
            fit <- fitdistcens(data,distr)
          }
        } else {
          if (!is.null(start)) {
            fit <- fitdist(values,distr,start=start,lower=lower)
          } else {
            fit <- fitdist(values,distr,lower=lower)
          }
	  if (!is.na(fit) && is.null(fit$aic)) {
  	    fit$aic = AIC(fit)
          }
        }
        return(fit)
      }
      error.fit = function(e) {
        #print(e)
        return(NA)
      }
      fit = tryCatch(do.fit(),error=error.fit)
      #fit = do.fit()

      fit.str = NA
      aic = NA
      if (!is.na(fit)) {
        aic = fit$aic
	fit.str = as.character(fit)[1]
      }

      fit.df = data.frame(
        parameter_set_id = pset,
        parameter_set_name = name,
        distr = distr,
	aic = aic,
	params = fit.str
      )    
      df <- rbind(df,fit.df)
      print(df)
    }
  }

  print(df)
  return(df)
}

sites.fit.intervals = function(db=dbname, sites.dir='.',filter.fn=NULL) {
  library(fitdistrplus)
  param.sets = get.psets(db=db)
  df = data.frame()
  for (i in 1:nrow(param.sets)) {
    pset <- param.sets[i,]$id
    name <- param.sets[i,]$parameter_set_shortname
    print(name)
    sites <- get.sites.for.parameter.set(dir=sites.dir,parameter.set.id=pset,filter.fn=filter.fn)
    sites <- sites[sample(100000),]
    sites$left <- sites$omega_lower
    sites$right <- sites$omega_upper

    for (distr in c('weibull','exp','gamma','lnorm')) {
      print(paste("Fit",distr,"..."))
 
      error.fit = function(e) {return(NA)}
      fit = tryCatch(fit.fn(sites,distr),error=error.fit)
      fit.str = NA
      aic = NA
      if (!is.na(fit)) {
        aic = fit$aic
	fit.str = as.character(fit)[1]
      }
      fit.df = data.frame(
        parameter_set_id = pset,
        parameter_set_name = name,
        distr = distr,
	aic = aic,
	params = fit.str
      )    
      df <- rbind(df,fit.df)
      print(df)
    }
  }
  return(df)
}

fit.fn = function(df,fn) {
  require(fitdistrplus)
  ft <- fitdistcens(df,fn,lower=-Inf)
  return(ft)
}

sites.plot.distributions = function(db=dbname, sites.dir='.', filter.fn=NULL) {
  library(MASS)
  param.sets = get.psets(db=db)
#  param.sets = subset(param.sets,id==1)
  for (i in 1:nrow(param.sets)) {
    pset <- param.sets[i,]$id
    name <- param.sets[i,]$parameter_set_shortname
    print(name)
    sites <- get.sites.for.parameter.set(dir=sites.dir,parameter.set.id=pset)
    #sites <- get.sites.for.parameter.set(dir=sites.dir,parameter.set.id=0) # Test set.
    print(nrow(sites))
    if (!is.null(filter.fn)) {
      sites <- filter.fn(sites)
    }
    omega.distribution.plot(sites)
  }
}

omega.distribution.plot = function(sites) {
# Plot site-wise big histogram.
par(
  mfrow=c(1,1)
  ,mar=c(4,17,1,2) # bltr
  ,mgp=c(9,7.5,6.5) # title, labels, line
  ,cex.axis=1
  ,mex=1
)
plot(NULL,xlab=expression(omega),ylab="",main="",
  xlim=c(0,1),ylim=c(0,1.5),axes=F)
par(mgp=c(16,14,13.5))
#axis(side=2,at=c(0,2,4,5,6.4),labels=c(0,2,4,5,105))
axis(side=2)
par(mgp=c(2,1,0))
axis(side=1,at=c(0,0.5,1))
title(xlab=expression(omega))
noZ = sites[sites$omega > 0,]
dens.lines("omega",data=noZ,n=200,orig_count=nrow(sites),lwd=1, col=rgb(1,1,1,0.8),plot.lines=F,plot.poly=F,plot.hist=T,fill=T)
dx = diff(c(0,1.2))/200
# Put a bar to the left of the histogram for the sites at zero.
z = sites[sites$omega==0,]
n_noZ = nrow(noZ)
n_sites = nrow(sites)
n_z = nrow(z)
n_const = nrow(sites[!is.na(sites$note) & sites$note=="constant",])
n_syn =   nrow(sites[!is.na(sites$note) & sites$note=="synonymous",])
print(paste(n_const,n_syn))
bar_height = 1.5
bar_width = (n_syn / n_sites) / bar_height
offset = -0.018
y_off = 0
par(xpd=NA)
# Constant on far left.
rect(offset-bar_width*2,y_off,offset-bar_width,y_off+bar_height*n_const/n_syn,col=rgb(1,1,1,1),density=NA,border=NA)
rect(offset-bar_width*2,y_off,offset-bar_width,y_off+bar_height*n_const/n_syn,col=rgb(.5,.5,.5,,1),density=10,angle=45,border='black')
# Synonymous on near left.
rect(offset-bar_width,y_off,offset,y_off+bar_height*1,col=rgb(1,1,1,1),density=NA,border=NA)
rect(offset-bar_width,y_off,offset,y_off+bar_height*1,col=rgb(.5,.5,.5,1),density=10,angle=-45,border='black')
#lines(x=c(offset-bar_width*2,0),y=c(y_off,0),lwd=1,col='black',lty="dashed")
#lines(x=c(offset,dx),y=c(y_off,0),lwd=1,col='black',lty="dashed")


# Legend.
l.col = c("red",rgb(.7,.7,.7,.8),"blue",gray(0.5),gray(0.5))
l.dens = c(NA,NA,NA,10,10)
l.ang = c(NA,NA,NA,45,-45)
l.text = c(expression(omega[lower]),expression(omega[ML]),expression(omega[upper]),
  "constant","synonymous")
legend(x=0.48,y=6-.75,xjust=1,legend=l.text,angle=l.ang,density=l.dens,fill=l.col,box.lwd=1)
m = expression("Site-wise"~omega~estimates)
text(x=0.48,y=6,m,cex=1.7,adj=1)
#xo = 0
#polygon(x=c(xo-0.015,xo-0.015,xo+0.025,xo+0.025),y=c(5.3,5.5,6,5.8),col='white',border=NA,lty=0)
#lines(x=c(xo-0.015,xo+0.025),y=c(5.3,5.8),lwd=1,lty=1);
#lines(x=c(xo-0.015,xo+0.025),y=c(5.5,6),lwd=1,lty=1);
#xo = -bar_width*2-.06
#polygon(x=c(xo-0.015,xo-0.015,xo+0.025,xo+0.025),y=c(5.3,5.5,6,5.8),col='white',border=NA,lty=0)
#lines(x=c(xo-0.015,xo+0.025),y=c(5.3,5.8),lwd=1,lty=1);
#lines(x=c(xo-0.015,xo+0.025),y=c(5.5,6),lwd=1,lty=1);

# Cumulative plot.
cum.sum = function(data,len=200,xlim=c(0,2)) {
  data = data[!is.na(data)]
  data = data[data >= xlim[1] & data <= xlim[2]]
  breaks=seq(from=xlim[1],to=xlim[2],length.out=len)
  bin = cut(data,breaks=breaks,include.lowest=T)
  est <- tabulate(bin, length(levels(bin)))
  est <- cumsum(est) / length(data)
  return(est)
}
print(str(sites))
cum.lo = cum.sum(sites$omega_lower)
cum.ml = cum.sum(sites$omega)
cum.hi = cum.sum(sites$omega_upper)
par(fig=c(0.6,0.9,0.48,0.95)
  ,new=TRUE
  ,mar=c(0,0,0,0)
  ,mgp=c(3,1,0)
  ,xpd=T)
plot.window(xlim=c(0,2),ylim=c(0,1))
axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"))
axis(side=1,at=c(0,1,2))
box()
xs = c(seq(from=0,to=2,length.out=length(cum.lo)))
abline(v=1,lty='dashed',col='gray')
lines(xs,c(cum.ml),lwd=2,col='black')
lines(xs,c(cum.lo),lwd=2,col='red')
lines(xs,c(cum.hi),lwd=2,col='blue')
###  
}

dens.lines = function(field,data,
  plot.poly=F,plot.lines=F,plot.hist=F,fill=F,log=F,
  mode="hist", n=50, orig_count=length(data[[field]]),
  border=NA,col='gray',line.col='black',lwd=2,
  xlim=NULL,ylim=NULL,
  adj=1,...) {
  usr = par("usr")
  x = data[[field]]
  x = x[!is.na(x)]
  if (is.null(xlim)) {xlim=c(usr[1],usr[2])}
  if (is.null(ylim)) {ylim=c(usr[3],usr[4])}
#  orig_count = length(x)
  if (log) x = log(x)
  if (mode=="dens") {
    d = density(x,from=0,to=xlim[2],adj=adj)
} else {
    # Largely from truehist (MASS pkg)
    x = x[x>xlim[1] & x<xlim[2]]
    
    breaks = seq(from=min(x),to=xlim[2],length.out=n)
    dx = diff(breaks)[1]
    bin = cut(x,breaks,include.lowest=T)
    est <- tabulate(bin, length(levels(bin)))
    est <- est/(dx * orig_count)
    n_breaks = length(breaks)
    max_y = ylim[2]
    d = list(x=breaks[-n_breaks]+dx/2,y=est)

    print(paste(field,sum(est*dx)))
    print(paste(field,"Dens at zero:",est[1]))
    if (plot.hist) {
      if (fill) {
        rect(breaks[-n], 0, breaks[-1], pmin(max_y,est), col=col,border=NA,...)
      }
      lines(breaks[-n],pmin(max_y,est), col=line.col,type='s',lwd=lwd)
  }
}
  poly_x = c(d$x)
  poly_y = c(d$y)
  if (plot.poly) {
    polygon(poly_x,poly_y,border=1,col=col,...)
}
  if (plot.lines) {
    lines(poly_x,poly_y,lwd=lwd,col=line.col,...)
  }
}

