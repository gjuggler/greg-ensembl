if (exists('drv')) {
  lapply(dbListConnections(drv),dbDisconnect)
} else {
  library(RMySQL)
  drv <- dbDriver('MySQL')
}

if (!exists('dbname')) {
  dbname = 'slrsim_anisimova'
  host = 'mysql-greg.ebi.ac.uk'
  port=4134
  user='slrsim'
  password='slrsim'
  userpass='slrsim:slrsim'
} else if (exists('dbname') && dbname == 'gj1_slrsim') {
  host = 'ens-research'
  port=3306
  user='ensro'
  password=''
  userpass='ensro'
} else {
  host = 'mysql-greg.ebi.ac.uk'
  port=4134
  user='slrsim'
  password='slrsim'
  userpass='slrsim:slrsim'
}
con <- dbConnect(drv, host=host, port=port, user=user, password=password, dbname=dbname)
dbURL = paste("mysql://",userpass,"@",host,":",port,"/",dbname,sep="")
print(paste("Connected to:",user,"@",host,":",port,"/",dbname))
print(paste("[",dbURL,"]"))
#url = paste("mysql://slrsim:slrsim@mysql-greg.ebi.ac.uk:4134/",dbname,sep="")

get.vector = function(con,query,columns=1) {
  res = dbSendQuery(con,query)
  if (columns == 'all') {
    data = fetch(res,n=-1)  # Important: the n=-1 causes the entire resultset to be fetched, not just the first 500.
  } else {
    data = fetch(res,n=-1)[[columns]]
  }
  dbClearResult(res)
  return(data)
}

# Grabs from the database the true and inferred omegas for the given reference sequence.
get.all.data = function(sites.cols=NULL,genes.cols=NULL,where=NULL) {
  if (is.null(sites.cols)) {
    sites.cols = "s.aln_dnds,s.aln_lrt,s.seq_position,s.true_dnds,s.true_ncod,s.true_type";
  }
  if (is.vector(sites.cols)) {
    sites.cols = paste(sites.cols,collapse=",")
  }
  if (is.null(genes.cols)) {
    genes.cols = "g.data_id,g.slrsim_label";
  }
  if (is.vector(genes.cols)) {
    genes.cols = paste(genes.cols,collapse=",")
  }

  if (is.null(where)) {
    where.statement = ''
  } else {
    where.statement = where
  }

  for (necessary in c(
    'g.experiment_name','g.slrsim_label',
    'g.data_id','g.node_id',
    'g.phylosim_insertrate','g.phylosim_insertmodel'
    )) {
    if (!necessary %in% genes.cols) {
      genes.cols = paste(genes.cols,necessary,sep=',')
    }
  }

  for (necessary in c(
    's.aln_dnds','s.aln_type','s.aln_note','s.aln_lrt',
    's.true_dnds','s.true_type','s.seq_position'
    )) {
    if (!necessary %in% sites.cols) {
      sites.cols = paste(sites.cols,necessary,sep=',')
    }
  }

  print(genes.cols)
  print(sites.cols)
  query = sprintf("select %s,%s from stats_genes g JOIN stats_sites s ON g.data_id=s.data_id %s",genes.cols,sites.cols,where.statement)
  all = get.vector(con,query,columns='all')
 
  return(all)
}

get.all.by.experiment = function() {
  query = "SELECT distinct(experiment_name) FROM stats_genes"
  experiments = get.vector(con,query)

  data.list = list()
  i=1
  for (experiment in experiments) {
    print(experiment)
    cur.data = get.all.data(where=paste('WHERE experiment_name="',experiment,'"',sep=''))
    data.list[[i]] = cur.data
    i = i + 1
  }

  return(data.list)
}

get.all.nodes = function() {
  query = sprintf("SELECT * FROM stats_sites GROUP BY node_id")
  data = get.vector(con,query,columns='all')
  return(data)
}

get.all.first.reps = function() {
  nodes = get.all.nodes()
  first.reps = subset(nodes,slrsim_rep == 1);
  return(first.reps)
}

dump.protein = function(node_id,base=node_id,tree_file=NA,aln_file=NA,sw_file=NA,params_file=NA,sw_table=NA) {
  if (is.na(tree_file)) {
    tree_file = paste(base,".nh",sep="")
  }
  if (is.na(aln_file)) {
    aln_file = paste(base,".fa",sep="")
  }
  if (is.na(sw_file)) {
    sw_file = paste(base,".csv",sep="")
  }
  if (is.na(params_file)) {
    params_file = paste(base,".txt",sep="")
  }
  if (is.na(sw_table)) {
    sw_table = "stats_sites"
  }

  url = "mysql://slrsim:slrsim@mysql-greg.ebi.ac.uk:4134/slrsim_anisimova"
  system(paste("perl ~/lib/greg-ensembl/scripts/tree_dump.pl",
    " --url=",url,
    " --id=",node_id,
    " --tree=",tree_file,
    " --aln=",aln_file,
    " --sw=",sw_file,
    " --sw_table=",sw_table,
    " --params=",params_file,
    sep="")
  )  
}

dump.nodes = function(nodes,base) {
  for (i in 1:nrow(nodes)) {
    node = nodes[i,]
    dump.protein(node$node_id,base=paste(base,node$node_id,sep=""))
  }
}

dump.all.nodes = function(base) {
  nodes = get.all.nodes()
  dump.nodes(nodes,base)
}

get.test.data = function() {
  query = sprintf("SELECT * FROM stats_slrsim LIMIT 500000");
  data = get.vector(con,query,columns='all')
  return(data)
}

is.paml = function(df) {
  if (!any(colnames(df) %in% c('sitewise_action'))) {
    return(FALSE)
  }
  if (grepl("paml",df[1,'sitewise_action'],ignore.case=T)) {
    return(TRUE)
  } else {
    return(FALSE)
  } 
}

df.stats = function(df,
    thresh=3.8,
    paml_thresh=0.95,type='all') {

  aln_thresh = 1
  if (is.paml(df)) {
    thresh = paml_thresh
    aln_thresh = -1
  }

    # Collect stats for SLR-type runs.
    pos_pos = nrow(subset(df,true_type=="positive1" & aln_lrt>thresh & aln_dnds>aln_thresh))
    neg_pos = nrow(subset(df,true_type!="positive1" & aln_lrt>thresh & aln_dnds>aln_thresh))
    neg_neg = nrow(subset(df,true_type!="positive1" & !(aln_lrt>thresh & aln_dnds>aln_thresh)))
    pos_neg = nrow(subset(df,true_type=="positive1" & !(aln_lrt>thresh & aln_dnds>aln_thresh)))

    pos_all = nrow(subset(df,true_type=="positive1"))
    neg_all = nrow(subset(df,true_type!="positive1"))
    all_pos = nrow(subset(df,aln_lrt>thresh & aln_dnds>aln_thresh))
    all_neg = nrow(subset(df,!(aln_lrt>thresh & aln_dnds>aln_thresh)))

  all = nrow(df)

  # see http://en.wikipedia.org/wiki/Receiver_operating_characteristic
  sens = pos_pos/pos_all  # also: tpr, power, hit rate, recall (anisimova 2002 as power)
  spec = neg_neg/neg_all  # also: tnr
  ppv = pos_pos/all_pos   # also: precision, (anisimova 2002 has it as accuracy)
  npv = neg_neg/all_neg
  fdr = neg_pos/all_pos
  fpr = neg_pos/neg_all   # also: false alarm rate, fallout
  dlr = sens/(1-spec)     # diagnostic likelihood ratio, see http://www.rapid-diagnostics.org/accuracy.htm#posdlr

  pos_true = pos_all/all # proportion of true positives
  pos_inf = all_pos/all # proportion of inferred positives

  true_pos_count = pos_pos
  true_neg_count = neg_neg

  return(list(
    sens=sens,
    spec=spec,
#    ppv=ppv,
#    npv=npv,
    fdr=fdr,
    fpr=fpr,
#    dlr=dlr,
    pos_true=pos_true,
    pos_inf=pos_inf,
    true_pos_count = true_pos_count,
    true_neg_count = true_neg_count
  ))
}

df.swfwer = function(df,
    thresh=3.8,
    paml_thresh=0.95) {

  aln_thresh = 1
  if (is.paml(df)) {
    thresh = paml_thresh
    aln_thresh = -1
  }

  df.fp = function(df,thresh) {
    fps = nrow(subset(df,true_type!="positive1" & aln_lrt>thresh & aln_dnds>aln_thresh))
    return(fps)
  }

  error_by_protein = by(df,df$node_id,df.fp,thresh=thresh)
  max_errors = max(error_by_protein)
  min_errors = min(error_by_protein)
  return(list(swfwer=sum(error_by_protein >= 1)/length(error_by_protein),
              max=max_errors,min=min_errors))
}

df.gene.detection.rate = function(df) {
  # Note: the detection rate is based on PAML's LRT test or SLR's multiple testing-corrected estimates.
  # (This means that a different SLR threshold will *not* change these calculated detection rates.)
  df.detected.positive = function(df) {
    if (is.paml(df)) {
      return(nrow(subset(df, paml_lrt > 5.99)) >= 1)
    } else {
      return(nrow(subset(df, aln_type == 'positive3' | aln_type == 'positive4')) >= 1)
    }
  }

  num.detected = by(df,df$node_id,df.detected.positive)
  return(list(detected=sum(num.detected)/length(num.detected)))
}

df.cor = function(df) {
  return(list(cor=cor(df$true_dnds,df$aln_dnds,method='spearman',use='complete.obs')))
}

df.alignment.accuracy = function(df) {
  df.sps = function(df) {
    return(df[1,]$sum_of_pairs_score)
  }
  df.tcs = function(df) {
    return(df[1,]$total_column_score)
  }
  sps.by.gene = by(df,df$node_id,df.sps)
  tcs.by.gene = by(df,df$node_id,df.tcs)

  df.unfiltered = function(df){return(df[1,]$unfiltered_site_fraction)}
  unfiltered.by.gene = by(df,df$node_id,df.unfiltered)

  return(list(
    sum_of_pairs=mean(sps.by.gene),
    total_column_score=mean(tcs.by.gene),
    aln_entropy=mean(df$aln_entropy),
    true_entropy=mean(df$true_entropy),
    unfiltered_fraction=mean(unfiltered.by.gene)
  ))
}

# Go through each experiment, calculate summary stats, and build a data frame.
summarize.results = function(data,thresh=3.8,paml_thresh=0.95,col.names=c('slrsim_label')) {

  # Paste together some metadata so that we have one ID per experiment.
  #labels <- apply(as.data.frame(data[,col.names]),1,paste,collapse='/')
  #data[,'label'] <- labels
  #attrs = c('label','slrsim_file','alignment_name','filtering_name','alignment_score_threshold','slrsim_ref','sitewise_action','phylosim_insertrate','slrsim_tree_length')

  #ids = rep("",nrow(data))
  #for (attr in attrs) {
  #  ids = paste(ids,data[[attr]],sep=" ")
  #}
  #data$id = ids

  data[,'label'] <- data[,'slrsim_label']

  split.sets = split(data,data$label)
  for (i in 1:length(split.sets)) {
    df = split.sets[[i]]
    print(df[1,]$label)

    # Calculate the summaries.
    stats = df.stats(df,thresh=thresh,paml_thresh=paml_thresh)
    gene.detection.rate = df.gene.detection.rate(df)
    swfwer = df.swfwer(df,thresh=thresh,paml_thresh=paml_thresh)
    cor = df.cor(df)
    acc = df.alignment.accuracy(df)

    # Put them all together into a data frame.
    attr.list = list()
    for (attr in attrs) {
      attr.list[[attr]] = df[[attr]][1]
    }   

    if(!exists('res.df')) {
      res.df = data.frame(attr.list,stats,gene.detection.rate,swfwer,cor,acc)
    } else {
      res.df = rbind(res.df,data.frame(attr.list,stats,gene.detection.rate,swfwer,cor,acc))
    }
  } 
  return(res.df)
}

plot.roc = function(df,color='black',lty='solid',lwd=1,overlay=T) {
  require(RColorBrewer)       

  if (!overlay) {
    plot(c(),xlim=c(0,0.1),ylim=c(0,1))
  }

  a = get.perf(df)
  plot(a,col=color,lty=lty,lwd=lwd,add=overlay,colorize=F)
}

get.perf = function(df,...) {
  if (!is.paml(df)) {
    #print("Signed LRT...")
    # Create a signed LRT if it's SLR-based data.
    df$aln_lrt = sign(df$aln_dnds-1)*df$aln_lrt
  }
  truth = as.integer( df$true_dnds > 1 )
  
  require(ROCR)
  pred = prediction(df$aln_lrt,truth)
  perf = performance(pred,measure="tpr",x.measure="fpr")
  return(perf)
}

get.perf.df = function(label="unlabeled",df) {
  perf = get.perf(df)
  df = data.frame(label=label,fpr=unlist(perf@x.values),tpr=unlist(perf@y.values))
  return(df)
}

empty.plot = function() {
  plot(c(),xlim=c(0,0.1),ylim=c(0,1),xlab='',ylab='')
}

generic.roc.plot = function(data,col.names='label',na.rm=F,plot=T,...) {
  comb.fdrs <- data.frame()
  comb.roc <- data.frame()

  labels <- apply(as.data.frame(data[,col.names]),1,paste,collapse='/')
  data[,'label'] <- labels
  for (lbl in sort(unique(labels))) {
    sub <- subset(data,label==lbl)
    print(paste("subset: ",lbl,nrow(sub),sep="/"))
    if (nrow(sub)==0) {print("Skipping!");next;}
    sub.roc <- slr.roc(sub,na.rm=na.rm)
    comb.roc <- rbind(sub.roc,comb.roc)
    #print(comb.roc[1,])
    fdr.row <- max(which(sub.roc[,'fdr'] <= 0.1))
    if (!is.na(fdr.row)) {
      #print(paste(sub.roc[fdr.row,'tn'],sub.roc[fdr.row,'tp']))
      mark.fdr = sub.roc[fdr.row,]
      mark.fdr$colour = 'black'
      comb.fdrs <- rbind(comb.fdrs,mark.fdr)
    }
    neutral.row <- max(which(sub.roc[,'score'] >= 0))
    if (!is.na(neutral.row)) {
      #print(neutral.row)
      mark.neutral = sub.roc[neutral.row,]
      mark.neutral$colour = 'gray'
      comb.fdrs <- rbind(comb.fdrs,mark.neutral)
    }
  }

  if (plot) {
    plot.roc(comb.roc,plot.unity=T,plot.points=comb.fdrs,...)
  } else {
    return(comb.roc)
  }
}

plot.roc = function(data,plot.x='tn',plot.y='tp',plot.unity=F,fill.below=F,plot.points=NA) {
  data$roc.x = data[,plot.x]
  data$roc.y = data[,plot.y]

  require(ggplot2)
  require(grid)
  print(paste("Plot.roc row count: ",nrow(data)))

  p <- ggplot(data,aes(x=roc.x,y=roc.y,colour=label)) + geom_line()
  no.leg = opts(legend.position="none")

  max.x = max(data$roc.x)

  p1 <- p + no.leg
  p2 <- p + coord_cartesian(xlim=c(0,max.x/20)) + xlim(0,max.x/20) + no.leg
  leg <- p + opts(keep="legend_box")

  if (plot.unity) {
    p1 <- p1 + geom_abline(linetype=2,colour='gray')
    p2 <- p2 + geom_abline(linetype=2,colour='gray')
  }

  if (fill.below) {
    p1 <- p1 + geom_area(data=data,aes(x=roc.x,y=roc.y))
    p2 <- p2 + geom_area(data=data,aes(x=roc.x,y=roc.y))
  }

  if (!is.na(plot.points)) {
    plot.points$roc.x = plot.points[,plot.x]
    plot.points$roc.y = plot.points[,plot.y]
    p1 <- p1 + geom_point(data=subset(plot.points,colour=="black"),aes(x=roc.x,y=roc.y),colour="black")
    p2 <- p2 + geom_point(data=subset(plot.points,colour=="black"),aes(x=roc.x,y=roc.y),colour="black")
    p2 <- p2 + geom_point(data=subset(plot.points,colour=="gray"),aes(x=roc.x,y=roc.y),colour="gray")
  }

  # See http://learnr.wordpress.com/2009/05/26/ggplot2-two-or-more-plots-sharing-the-same-legend/
  lay = grid.layout(nrow=2,ncol=2,
    widths=unit(c(2,0.5),c("null","null")),
    heights=unit(c(1,1),c("null","null")))
  vplayout = function(...) {grid.newpage();pushViewport(viewport(layout=lay))}
  subplot = function(x,y){viewport(layout.pos.row=x,layout.pos.col=y)}    
  vplayout()
  print(p1,vp=subplot(1,1))
  print(p2,vp=subplot(2,1))
  print(leg,vp=subplot(1:2,2))
}

generic.roc = function(data,by='slrsim_label',na.rm=F) {
  comb.roc <- data.frame()
  #labels <- apply(as.data.frame(data[,col.names]),1,paste,collapse='/')
  data[,'label'] <- data[,by]
  labels = data[,'label']
  for (lbl in sort(unique(labels))) {
    print(lbl)
    sub <- subset(data,label==lbl)
    if (nrow(sub)==0) {next}
    sub.roc <- slr.roc(sub,na.rm=na.rm)
    comb.roc <- rbind(sub.roc,comb.roc)
  }
  return(comb.roc)
}

my_summary <- function(df) {
  stats = df.stats(df,thresh=3.8)
  data.frame(
    fpr = stats$fpr,
    fdr = stats$fdr,
    pos = stats$true_pos_count,
    pos_inf = stats$pos_inf,
    pos_true = stats$pos_true
  )
}

summarize.by.labels = function(data,by='slrsim_label') {
  library(plyr)
  a = ddply(data,.(slrsim_label),my_summary)
  print(a)
}

get.fdr.thresholds = function(data,by='slrsim_label',summarize=T) {
  output.df = data.frame()
  data[,'label'] = data[,by]
  for (lbl in unique(data[,'label'])) {
    print(paste("Label subset:",lbl))
    data.subset = subset(data,label==lbl)
    roc = generic.roc(data.subset,by=by)

    for (fdr in c(0.01,0.05,0.1,0.2,0.5)) {
      fdr.ok = roc[which(roc[,'fdr'] <= fdr),]
      if (nrow(fdr.ok) == 0) {next}
      max.row = fdr.ok[nrow(fdr.ok),]
      score.at.fdr = max.row[,'score']
      tpr.at.fdr = max.row[,'tpr']
      max.row[,'fdr.threshold'] = fdr
      max.row[,'pct_below'] = nrow(fdr.ok) / nrow(roc)
      output.df = rbind(output.df,max.row)
    }
  }

  summary = subset(output.df,select=c(slrsim_label,fdr.threshold,fdr,score,pct_below))
  return(summary)
}

facet.roc.plot = function(data,col.names='label',na.rm=F,plot.x='fp',plot.y='tp',facet.x='alignment_name',facet.y='filtering_name',zoom=F) {
  require(ggplot2)
  require(grid)

  data = generic.roc(data=data,col.names=col.names,na.rm=na.rm)

  data$roc.x = data[,plot.x]
  data$roc.y = data[,plot.y]
  data$facet.x = data[,facet.x]
  data$facet.y = data[,facet.y]

  p <- ggplot(data,aes(x=roc.x,y=roc.y,group=label,colour=alignment_score_threshold)) + geom_line()
  p <- p + scale_colour_gradient(low="blue",high="red",limits=c(0,10))
  p <- p + facet_grid(facet.y ~ facet.x)
  if (zoom) {
    max.x = max(data$roc.x)
    p <- p + coord_cartesian(xlim=c(0,max.x/20)) + xlim(0,max.x/20)
  }
  print(p)
}

slr.roc = function(df,na.rm=F) {
  library(doBy)
  library(plyr)

  if (na.rm) {
    n.na = nrow(subset(df,is.na(aln_lrt)))
    df = subset(df,!is.na(aln_lrt))
    print(sprintf("Removed %d NA rows",n.na))
  }

  if (!is.paml(df)) {
    # Create a signed LRT if it's SLR-based data.
    # We'll replace NAs with an extremely low score.
    df[is.na(df$aln_lrt),'aln_lrt'] = 1000
    df[is.na(df$aln_dnds),'aln_dnds'] = 0
    df$score = sign(df$aln_dnds-1)*df$aln_lrt
  } else {
    df$score = df$aln_lrt
  }
  df$truth = as.integer( df$true_dnds > 1 )
  df <- orderBy(~-score,data=df)

  df$tp = cumsum(df$truth)
  df$tn = cumsum(1-df$truth)
  df$count = cumsum(rep(1,nrow(df)))
  df$fp = df$tn

  df$fpr = df$tn / max(df$tn)
  df$tpr = df$tp / max(df$tp)

  df$fdr = df$fp/(df$count)

  return(df)
}

plot.fdr = function(data,plot.y='tpr',plot.x='fdr',facet.x='alignment_name',facet.y='filtering_name',color.by='alignment_score_threshold') {
  data$roc.x = data[,plot.x]
  data$roc.y = data[,plot.y]
  data$facet.x = data[,facet.x]
  data$facet.y = data[,facet.y]
  data$color.by = data[,color.by]

  require(ggplot2)
  require(grid)
  p <- ggplot(data,aes(x=roc.x,y=roc.y)) + geom_line(colour="black",size=1)
  # The 'group.counter' variable is created by the 'fdr.sweep' function to group together FDR chunks.
  p <- p + geom_area(aes(fill=color.by,group=group.counter),position='identity')
  p <- p + scale_fill_continuous(low="white",high="red",limits=c(0,10))
  p <- p + xlab(plot.x) + ylab(plot.y)
  p <- p + facet_grid(facet.y ~ facet.x)
  print(p)
}

fdr.sweep = function(data,step=0.01,min=0,max=1,max.search='tp',facet.names=c('alignment_name','filtering_name'),group.names=c('alignment_score_threshold')) {
  comb.df = data.frame()
  data$max.search = data[,max.search]

  labels <- apply(as.data.frame(data[,facet.names]),1,paste,collapse='/')
  data$label <- labels

  groups <- apply(as.data.frame(data[,group.names]),1,paste,collapse='/')
  data$group <- groups

  group.counter <- 0
  for (lbl in sort(unique(data$label))) {
    group.counter <- group.counter + 1
    sub <- subset(data,label==lbl)
    print(paste("label:",lbl," n:",nrow(sub)))
#     sub = data

    last.row <- NULL
    fdrs <- seq(from=min,to=max,by=step)
    for (i in fdrs) {
      # Find the row which maximizes the max.search value within the given FDR.
      within.fdr <- subset(sub,fdr <= i)
      within.fdr <- orderBy(~-max.search,data=within.fdr)
      row <- within.fdr[1,]
      print(sprintf("fdr:%.3f n:%d best:%s",i,nrow(within.fdr),row$group))
    
      # Some logic to fix up the boundaries between label groups.
      if (!is.null(last.row) && last.row$group != row$group) {
        group.counter <- group.counter + 1
        # Duplicate the row with the next-highest FDR.
        last.row$group.counter <- group.counter
        comb.df <- rbind(comb.df,last.row)
      }
      row$group.counter <- group.counter
  
      comb.df <- rbind(comb.df,row)
      last.row <- row
    }
  }

  data$max.search <- NULL
  comb.df$max.search <- NULL
  return(comb.df)
}

gg.df.cor = function(df) {
  return(cor(df$true_dnds,df$aln_dnds,method='spearman',use='complete.obs'))
}

gg.df.fpr = function(df) {
  df.res = df.stats(df)
  return(df.res[['fpr']])
}

correlation.boxplot = function(data,variable,n=5,title="") {
  png("~/public_html/blank.png")  

  dx = (max(data[[variable]]) - min(data[[variable]])) / (n-1)
  data$var_orig = data[[variable]]
  data$var = round_any(data[[variable]],dx)

  cor = ddply(data,.(var),'gg.df.cor')

  p <- qplot(x=var,y=gg.df.cor,data=cor) + xlab(variable) + ylab("Spearman dN/dS Correlation") + ylim(0,1)
  gr <- geom_rug(data=data,aes(x=var_orig,alpha=0.05,size=0.1))
  p + gr + opts(legend.position="none") + opts(title=title)

  #p + geom_histogram(data=data,aes(x=var, y=..density..))
  #p + stat_summary(fun.data=fun,colour="red",geom="point",size=3)
  #p + geom_boxplot() + xlab(variable) + ylab("Error")
  ggsave("~/public_html/slrsim_cor.png",width=8,height=6)
  dev.off()
}

my.own.roc = function(data, score.field='score', truth.fun=NULL) {
  library(doBy)
  library(plyr)

  #png("~/public_html/test.png")
  #n <- 2000
  #truth <- rbinom(n,1,0.8)
  #score <- rnorm(n)

  truth <- adply(df, 1, truth.fun);
  score <- data[[score.field]]

  data$truth = truth
  data$score = score
  data <- orderBy(~+score,data=data)

  data$tp = cumsum(data$truth)
  data$fp = cumsum(1-data$truth)

  return(data)
}

plot.ncod = function(data,col.names=c('alignment_name')) {
  require(ggplot2)
  labels <- apply(as.data.frame(data[,col.names]),1,paste,collapse='/')
  data[,'label'] <- labels
  print(str(data))
  nodes = sort(unique(data$slrsim_rep))
  for (node in nodes[1]) {
    df = subset(data,slrsim_rep==node)
    p <- ggplot(df,aes(x=seq_position,y=aln_ncod,colour=label)) + geom_line()
    print(p)
  }
}

plot.scatter = function(data,col.names=c('alignment_name','filtering_name')) {
  require(ggplot2)
  labels <- apply(as.data.frame(data[,col.names]),1,paste,collapse='/')
  data[,'label'] <- labels
  p <- ggplot(data, aes(x=true_dnds,y=aln_dnds,colour=aln_lrt)) + xlim(0,2) + ylim(0,2)
  p <- p + opts(aspect.ratio=1) + geom_point()
  p <- p + facet_wrap(~ label)
  print(p)
}

plot.by.columns = function(data,
  base.dir='.',
  script.dir='~/src/greg-ensembl/scripts',
  col.names=c('slrsim_label'),
  skip.plot = F
) {
  source(paste(script.dir,"/aln-tools/aln.tools.R",sep=''))
  source(paste(script.dir,"/aln-tools/phylo.tools.R",sep=''))
  source(paste(script.dir,"/aln-tools/plot.phylo.greg.R",sep=''))
  library(ape)

  labels <- apply(as.data.frame(data[,col.names]),1,paste,collapse='_')
  labels <- sub(" ","",labels)
  data[,'label'] <- labels

  for (lbl in sort(unique(labels))) {
    sub <- subset(data,label==lbl)
    sub <- orderBy(~slrsim_rep,data=sub)
    first.row = sub[1,]
 
    lbl = gsub('[^0-9a-z_]','x',lbl)   
    print(paste(lbl,"..."))

    node_id = first.row$node_id    
    tree_file = paste(base.dir,"/",lbl,'.nh',sep='')
    aln_file = paste(base.dir,"/",lbl,'.fa',sep='')
    sitewise_file = paste(base.dir,"/",lbl,'.csv',sep='')

#    system(paste("perl ~/lib/greg-ensembl/scripts/tree_dump.pl",
    system(paste("perl ~/src/greg-ensembl/scripts/tree_dump.pl",
    " --url=",dbURL,
    " --id=",node_id,
    " --tree=",tree_file,
    " --aln=",aln_file,
    " --sw=",sitewise_file,
    sep=""))

    if (!skip.plot) {
      plot.file = paste(base.dir,"/",first.row$label,".png",sep="")
      png(file=plot.file,width=2000,height=200)
      plot.protein(aln.file=aln_file,tree.file=tree_file,sitewise.file=sitewise_file,remove.files=T)
      dev.off()
    }

  }  
}

plot.protein = function(aln.file,tree.file,sitewise.file,remove.files=F) {
  aln = read.aln(aln.file)
  tree = read.tree(tree.file)
  aln$tree = tree
  
  length = aln$length
  tree.space = length/10
  par(mar=c(0,0,0,0))
  plot.new()
  plot.window(xlim=c(0,length + tree.space),ylim=c(0,aln$num_seqs))
  plot.aln(aln,overlay=T,plot.tree=F,plot.chars=F,
    x.lim=c(tree.space,length+tree.space),
    y.lim=c(0,aln$num_seqs))

#  text(x=tree.space*.9,y=0,labels=label,cex=2,adj=c(0,0))

  rect(xleft=tree.space,xright=length+tree.space,ybottom=0,ytop=aln$num_seqs)
  par(new=T)
  len = tree_length(tree)
  print(len)
  plot.phylo.greg2(tree,x.lim=c(0,(4)/(tree.space/length)),y.lim=c(.5,aln$num_seqs+.5),edge.width=0.5)

  if(remove.files) {
    unlink(c(tree.file,aln.file,sitewise.file))
  }
}