if (exists('drv')) {
  lapply(dbListConnections(drv),dbDisconnect)
} else {
  library(RMySQL)
  drv <- dbDriver('MySQL')
}
connect.db = function(dbName) {
  con <- dbConnect(drv, host='ens-research', port=3306, user='ensro', password='', dbname='gj1_eslr_57')
  return(con)
}
con = connect.db("gj1_eslr_57")

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
get.psets = function(db='gj1_eslr_57') {
  query = sprintf(
    paste('SELECT n.parameter_set_id AS id,sn.parameter_value AS parameter_set_shortname, n.parameter_value AS name FROM %s.parameter_set n, %s.parameter_set sn',
      'WHERE n.parameter_set_id=sn.parameter_set_id AND n.parameter_name="name" AND sn.parameter_name="parameter_set_shortname";'),
    db,db)
  return(get.vector(con,query,columns='all'))
}

get.genes = function(parameter.set.id=1,db="gj1_eslr_57") {
  query = sprintf("SELECT * FROM %s.stats_genes where parameter_set_id=%s",db,parameter.set.id);
  data = get.vector(con,query,columns='all')

  return(data)
}

get.genes.list = function(db="gj1_eslr_57") {
  query = sprintf("SELECT * FROM %s.stats_genes",db);
  data = get.vector(con,query,columns='all')
  return(data)
}

get.genes.merged = function(db="gj1_eslr_57",exclude.cols=c()) {
  #all.node.ids = get.vector(con,sprintf('SELECT DISTINCT(node_id) FROM %s.stats_genes',db))
  
  param.sets = get.psets(db=db)
  print(param.sets)
  
  genes = get.genes(1,db=db)
  genes <- within(genes,rm(list=exclude.cols))
  
  for (i in 1:nrow(param.sets)) {
    pset = param.sets[i,]$id
    print(param.sets[i,]$parameter_set_shortname)
    cur.genes = get.genes(pset,db=db)

    create.name = function(ext) {paste(param.sets[pset,]$parameter_set_shortname,ext,sep='')}
    
    col.slr.dnds = create.name('.slr.dnds')
    col.hyphy.dnds = create.name('.hyphy.dnds')
    col.hyphy.dnds.lo = create.name('.hyphy.dnds.lo')
    col.hyphy.dnds.hi = create.name('.hyphy.dnds.hi')
    col.psc.count = create.name('.psc.count')
    col.weak.psc.count = create.name('.weak.psc.count')

    genes.subset = data.frame(
      a=cur.genes$node_id,
      b=cur.genes$slr_dnds,
      c=cur.genes$hyphy_dnds,
      d=cur.genes$hyphy_dnds_lo,
      d=cur.genes$hyphy_dnds_hi,
      e=cur.genes$psc_count,
      f=cur.genes$weak_psc_count
      )
    colnames(genes.subset) = c(
              'node_id',
              col.slr.dnds,
              col.hyphy.dnds,
              col.hyphy.dnds.lo,
              col.hyphy.dnds.hi,
              col.psc.count,
              col.weak.psc.count
              )
    
    genes = merge(genes,genes.subset,all.x=T)
  }
  return(genes);
}

get.positive.sites = function(parameter.set.id=1,db="gj1_eslr_57",limit="") {
  query = sprintf("SELECT * from %s.stats_sites WHERE parameter_set_id=%s AND type IN ('positive1','positive2','positive3','positive4') %s",db,parameter.set.id,limit)
  sites = get.vector(con,query,columns='all')
  return(sites)  
}

get.sites = function(parameter.set.id=1,db="gj1_eslr_57",limit="") {
  query = sprintf("SELECT * from %s.stats_sites WHERE parameter_set_id=%s %s",db,parameter.set.id,limit)
  sites = get.vector(con,query,columns='all')
  return(sites)
}

get.sites.slim = function(fields.list=NA,parameter.set.id=NA,db="gj1_eslr_57",limit.string="") {
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

get.sites.for.parameter.set = function(dir='.',parameter.set.id) {
  load(paste(dir,"/sites_",parameter.set.id,".Rdata",sep=""))
  return(sites)
}

sites.summary = function(db="gj1_eslr_57", sites.dir='.', filter.fn=NULL) {
  param.sets = get.psets(db=db)

  df = data.frame()
  for (i in c(1)) {
#  for (i in 1:nrow(param.sets)) {
    pset <- param.sets[i,]$id
    name <- param.sets[i,]$parameter_set_shortname
    print(pset)
    print(name)

    sub.sites <- get.sites.for.parameter.set(dir=sites.dir,parameter.set.id=pset)

    if (!is.null(filter.fn)) {
      sub.sites <- filter.fn(sub.sites)
    }
    
    n <- nrow(sub.sites)
    if (n == 0) {
      next
    }
    pos.sites <- subset(sub.sites, type %in% c('positive1','positive2','positive3','positive4'))
    strong.pos.sites <- subset(sub.sites, type %in% c('positive2','positive3','positive4'))
    n.pos.sites <- nrow(pos.sites)
    n.strong.pos.sites <- nrow(strong.pos.sites)
    n.neg.sites <- nrow(subset(sub.sites, type %in% c('negative1','negative2','negative3','negative4')))
    n.strong.neg.sites <- nrow(subset(sub.sites, type %in% c('negative2','negative3','negative4')))
    pos.domains <- subset(strong.pos.sites,!is.na(domain))
    x <- pos.domains[,c('node_id','domain')]
    y <- pos.domains[,c('domain')]
    n.pos.domains <- nrow(unique(x))
    n.pos.domain.types <- length(unique(y))

    x <- strong.pos.sites[,c('node_id')]
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
      `positive` = n.pos.sites / n,
      `negative` = n.neg.sites / n,
      `strong.positive` = n.strong.pos.sites / n,
      `strong.negative` = n.strong.neg.sites / n,
      `n.positive` = n.pos.sites,
      `n.negative` = n.neg.sites,
      `n.strong.positive` = n.strong.pos.sites,
      `n.strong.negative` = n.strong.neg.sites,
      `n.positive.domains` = n.pos.domains,
      `n.positive.domain.types` = n.pos.domain.types,
      `n.positive.genes` = n.pos.genes,
      `n.total.genes` = n.total.genes
    )
    df <- rbind(df,summary.df)
    rm(pos.sites,pos.domains,sub.sites)
    gc()
  }

  print(df)
  return(df)
}

sites.fit.distributions = function(db="gj1_eslr_57", sites.dir='.', filter.fn=NULL) {
  library(MASS)
  param.sets = get.psets(db=db)
  df = data.frame()
  for (i in 1:nrow(param.sets)) {
    pset <- param.sets[i,]$id
    name <- param.sets[i,]$parameter_set_shortname
    print(name)
    sites <- get.sites.for.parameter.set(dir=sites.dir,parameter.set.id=pset)
    if (!is.null(filter.fn)) {
      sites <- filter.fn(sites)
    }
    
    omegas <- sites$omega
    sites.adj <- sites
    sites.adj[sites.adj$omega == 0,]$omega <- 0.0001
    omegas.adj <- sites.adj$omega

    set.seed(0)
    omegas.adj <- sample(omegas.adj, size=1000000)
    omegas.below.one <- omegas.adj[omegas.adj < 1]
    param.min = sqrt(.Machine$double.eps)

    for (distr in c('weibull','beta','exponential','gamma','lognormal')) {
      print(paste("Fit",distr,"..."))
      start = NULL
      values = omegas.adj
      lower = -Inf
      if (distr == 'beta') {values = omegas.below.one; start = list(shape1=1,shape2=1); lower=c(param.min,param.min)}
      if (distr == 'weibull') {start = list(shape=0.5,scale=0.5); lower=c(param.min,param.min)}
      if (distr == 'gamma') {lower=c(param.min,param.min)}
      do.fit = function() {
        if (!is.null(start)) {
          fit <- fitdistr(values,distr,start=start,lower=lower)
        } else {
          fit <- fitdistr(values,distr,lower=lower)
        }
        return(fit)
      }
      error.fit = function(e) {
        #print(e)
        return(NA)
      }
      fit = tryCatch(do.fit(),error=error.fit)

      fit.str = NA
      aic = NA
      if (!is.na(fit)) {
        aic <- AIC(fit)
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
    }
  }

  print(df)
  return(df)
}

sites.plot.distributions = function(db="gj1_eslr_57", sites.dir='.', filter.fn=NULL) {
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
