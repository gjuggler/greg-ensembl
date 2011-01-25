library(ggplot2)

get.subsets.for.env = function(env) {
  with(env,{
  # Branch models summary:
  # a: (H, G, others)
  # b: (H#1, G, others)
  # c: (H, G#1, others)
  # d: (H#2, G#1, others)
  # e: (H#1, G#1, others)
  # f: (((H,C),G)#1, others)
  # g: (((H,C),G)$1, others)
  # h: (((H,C),G), others)

  # Which omegas are which for each test?
  # pval.1: fg=b_omega_1, bg=b_omega_0 # human
  # pval.2: fg=c_omega_1, bg=c_omega_0 # gorilla
  # pval.3: fg=d_omega_1, bg=d_omega_0 # gorilla
  # pval.4: fg=d_omega_2, bg=d_omega_0 # human
  # pval.6: fg=d_omega_1, bg=d_omega_2 # sig. hum-gor diff.
  # pval.7: fg=f_omega_1, bg=f_omega_0 # great ape branch
  # pval.8: fg=g_omega_1, bg=g_omega_0 # great ape clade

  # Only look at genes with at least N gorilla subs!!
  #stats.lnl = subset(stats.lnl,(subs_n_Ggor + subs_s_Ggor) >= 4)
  stats.lnl = subset(stats.lnl,((subs_n_Ggor + subs_s_Ggor) >= 3) & ((subs_n_other + subs_s_other) >= 3))

  # Add signed LRTs.
  stats.lnl[,'lrt.1.signed'] = with(stats.lnl,2*pmax(0,(b_lnL-a_lnL)) * sign(b_omega_1-b_omega_0))
  stats.lnl[,'lrt.2.signed'] = with(stats.lnl,2*pmax(0,(c_lnL-a_lnL)) * sign(c_omega_1-c_omega_0))
  stats.lnl[,'lrt.3.signed'] = with(stats.lnl,2*pmax(0,(d_lnL-b_lnL)) * sign(d_omega_1-d_omega_0))
  stats.lnl[,'lrt.4.signed'] = with(stats.lnl,2*pmax(0,(d_lnL-c_lnL)) * sign(d_omega_2-d_omega_0))
  stats.lnl[,'lrt.5.signed'] = with(stats.lnl,2*pmax(0,(e_lnL-a_lnL)) * sign(e_omega_1-e_omega_0))
  stats.lnl[,'lrt.6.signed'] = with(stats.lnl,2*pmax(0,(d_lnL-e_lnL)) * sign(d_omega_1-d_omega_2))
  stats.lnl[,'lrt.7.signed'] = with(stats.lnl,2*pmax(0,(f_lnL-h_lnL)) * sign(f_omega_1-f_omega_0))
  stats.lnl[,'lrt.8.signed'] = with(stats.lnl,2*pmax(0,(g_lnL-h_lnL)) * sign(g_omega_1-g_omega_0))

  # Add the p-values
  stats.lnl[,'pval.1'] = with(stats.lnl,1 - pchisq(2*(b_lnL-a_lnL),df=1))
  stats.lnl[,'pval.2'] = with(stats.lnl,1 - pchisq(2*(c_lnL-a_lnL),df=1))
  stats.lnl[,'pval.3'] = with(stats.lnl,1 - pchisq(2*(d_lnL-b_lnL),df=1))
  stats.lnl[,'pval.5'] = with(stats.lnl,1 - pchisq(2*(e_lnL-a_lnL),df=1))
  stats.lnl[,'pval.4'] = with(stats.lnl,1 - pchisq(2*(d_lnL-c_lnL),df=1))
  stats.lnl[,'pval.6'] = with(stats.lnl,1 - pchisq(2*(d_lnL-e_lnL),df=1))
  stats.lnl[,'pval.7'] = with(stats.lnl,1 - pchisq(2*(f_lnL-h_lnL),df=1))
  stats.lnl[,'pval.8'] = with(stats.lnl,1 - pchisq(2*(g_lnL-h_lnL),df=1))

  method = 'BH'
  stats.lnl[,'pval.1.bh'] = with(stats.lnl,p.adjust(pval.1,method=method))
  stats.lnl[,'pval.2.bh'] = with(stats.lnl,p.adjust(pval.2,method=method))
  stats.lnl[,'pval.3.bh'] = with(stats.lnl,p.adjust(pval.3,method=method))
  stats.lnl[,'pval.4.bh'] = with(stats.lnl,p.adjust(pval.4,method=method))
  stats.lnl[,'pval.5.bh'] = with(stats.lnl,p.adjust(pval.5,method=method))
  stats.lnl[,'pval.6.bh'] = with(stats.lnl,p.adjust(pval.6,method=method))
  stats.lnl[,'pval.7.bh'] = with(stats.lnl,p.adjust(pval.7,method=method))
  stats.lnl[,'pval.8.bh'] = with(stats.lnl,p.adjust(pval.8,method=method))

  stats.lnl[,'fg.1'] = stats.lnl$b_omega_1
  stats.lnl[,'fg.2'] = stats.lnl$c_omega_1
  stats.lnl[,'fg.3'] = stats.lnl$d_omega_1
  stats.lnl[,'fg.4'] = stats.lnl$d_omega_2
  stats.lnl[,'fg.5'] = stats.lnl$e_omega_1
  stats.lnl[,'fg.6'] = stats.lnl$d_omega_1
  stats.lnl[,'fg.7'] = stats.lnl$f_omega_1
  stats.lnl[,'fg.8'] = stats.lnl$g_omega_1

  stats.lnl[,'bg.1'] = stats.lnl$b_omega_0
  stats.lnl[,'bg.2'] = stats.lnl$c_omega_0
  stats.lnl[,'bg.3'] = stats.lnl$d_omega_0
  stats.lnl[,'bg.4'] = stats.lnl$d_omega_0
  stats.lnl[,'bg.5'] = stats.lnl$e_omega_0
  stats.lnl[,'bg.6'] = stats.lnl$d_omega_2
  stats.lnl[,'bg.7'] = stats.lnl$f_omega_0
  stats.lnl[,'bg.8'] = stats.lnl$g_omega_0

  # Raw p-value threshold is 0.05
  for (test in c(1,2,3,4,5,6,7,8)) {
    pval.lbl <- paste('pval.',test,sep='')
    pval.bh.lbl <- paste('pval.',test,'.bh',sep='')
    fg.lbl <- paste('fg.',test,sep='')
    bg.lbl <- paste('bg.',test,sep='')

    up.name <- paste(test,'.up',sep='')
    down.name <- paste(test,'.down',sep='')
    up.bh.name <- paste(test,'.up.bh',sep='')
    down.bh.name <- paste(test,'.down.bh',sep='')

    tmp <- stats.lnl
    
    tmp$sub_pval <- tmp[,pval.lbl]
    tmp$sub_pval_bh <- tmp[,pval.bh.lbl]
    tmp$sub_fg <- tmp[,fg.lbl]
    tmp$sub_bg <- tmp[,bg.lbl]

    t <- 0.05
    sub.up <- subset(tmp,sub_pval < t & sub_fg > sub_bg)
    sub.down <- subset(tmp,sub_pval < t & sub_fg < sub_bg)

    t <- 0.1
    sub.up.bh <- subset(tmp,sub_pval_bh < t & sub_fg > sub_bg)
    sub.down.bh <- subset(tmp,sub_pval_bh < t & sub_fg < sub_bg)
    
    assign(up.name,sub.up)
    assign(down.name,sub.down)
    assign(up.bh.name,sub.up.bh)
    assign(down.bh.name,sub.down.bh)
  }

  ## Parallel h-g shift.
  p.up.genes = intersect(`3.up`$name,`4.up`$name)
  p.up = subset(stats.lnl,name %in% p.up.genes)
  p.up$sub_fg = p.up$d_omega_1
  p.up$sub_pval = p.up$pval.3
  p.up.bh.genes = intersect(`3.up.bh`$name,`4.up.bh`$name)
  p.up.bh = subset(stats.lnl,name %in% p.up.bh.genes)
  p.up.bh$sub_fg = p.up.bh$d_omega_1
  p.up.bh$sub_pval = p.up.bh$pval.3
  p.down.genes = intersect(`3.down`$name,`4.down`$name)
  p.down = subset(stats.lnl,name %in% p.down.genes)
  p.down$sub_fg = p.down$d_omega_1
  p.down$sub_pval = p.down$pval.3
  p.down.bh.genes = intersect(`3.down.bh`$name,`4.down.bh`$name)
  p.down.bh = subset(stats.lnl,name %in% p.down.bh.genes)
  p.down.bh$sub_fg = p.down.bh$d_omega_1
  p.down.bh$sub_pval = p.down.bh$pval.3
  
  rm(t)
  })
}

parallel.genes <- function(human.env,chimp.env) {
  GH <- human.env$`3.up.bh`$name
  GC <- chimp.env$`3.up.bh`$name

  select.cols <- c('name','pval.3.bh','pval.2.bh','d_omega_1','d_omega_0','a_omega_0','go')

  # GH only -- possible GC parallel acceleration.
  GH.only <- setdiff(GH,GC)
  GH.rows <- subset(human.env$`3.up.bh`, name %in% GH.only, select=select.cols)
  GH.rows$label <- '3_up_GH_only'
  # GC only -- possible GH parallel acceleration.
  GC.only <- setdiff(GC,GH)
  GC.rows <- subset(chimp.env$`3.up.bh`, name %in% GC.only, select=select.cols)
  GC.rows$label <- '3_up_GC_only'

  GH <- human.env$`2.up.bh`$name
  GC <- chimp.env$`2.up.bh`$name
  # GH only -- possible GC parallel acceleration.
  GH.only <- setdiff(GH,GC)
  GH2.rows <- subset(human.env$`2.up.bh`, name %in% GH.only, select=select.cols)
  GH2.rows$label <- '2_up_GH_only'
  # GC only -- possible GH parallel acceleration.
  GC.only <- setdiff(GC,GH)
  GC2.rows <- subset(chimp.env$`2.up.bh`, name %in% GC.only, select=select.cols)
  GC2.rows$label <- '2_up_GC_only'

  df <- rbind(GH.rows,GC.rows,GH2.rows,GC2.rows)
  write.csv(df,file="possibly_parallel.csv")
}

env.subsets.genes <- function(env.name) {
  env <- get(env.name)
  env$env.name <- env.name
  outer.df <- with(env, {
     categories = c('up','up.bh','up.bh.above','up.bh.above.15','down','down.bh','down.bh.below','down.bh.below.0.05')
    tests <- c(1,2,3,4,6,7,8,'p')

    df <- data.frame()
    for (test in tests) {
      for (type in categories) {

        cur.row <- c()
        data.string <- paste(test,type,sep='.')
        if (type == 'up.bh.above' || type == 'up.bh.above.15') {
          data.string <- paste(test,'up.bh',sep='.')
        } else if (type == 'down.bh.below' || type == 'down.bh.below.0.05') {
          data.string <- paste(test,'down.bh',sep='.')
        }
        d.sub <- get(data.string)
        if (type == 'up.bh.above') d.sub <- subset(d.sub, sub_fg > 1)
        else if (type == 'up.bh.above.15') d.sub <- subset(d.sub, sub_fg > 1.5)
        else if (type == 'down.bh.below') d.sub <- subset(d.sub, sub_fg < 0.1)
        else if (type == 'down.bh.below.0.05') d.sub <- subset(d.sub, sub_fg < 0.05)
  
        if (nrow(d.sub) == 0) next;
        if (test != 'p' && (type == 'up' || type == 'down')) next;

        # Collect data about the genes in this subset.
        d.sub$label <- paste(test,type,sep='.')
        pval.name <- paste('pval',test,'bh',sep='.')
        if (test != 'p') {
          d.sub$pval.bh <- d.sub[,c(pval.name)]
        } else {
          d.sub$pval.bh <- (d.sub$pval.3 + d.sub$pval.4)
        }
        library(doBy)
        d.sub <- orderBy(~pval.bh,data=d.sub)
        #d.sub <- orderBy(~-grantham_Ggor,data=d.sub)

        #print(str(d.sub))

        cols <- c('label','Hsap_protein','name','pval.bh',
                    'sub_fg','grantham_Ggor','grantham_other','subs_n_Ggor','subs_s_Ggor','subs_n_other','subs_s_other','go')
        cur.row <- subset(d.sub,select=cols)
        #print(paste(test,type))
        df <- rbind(df,cur.row) 
      }
    }
    df$base <- base

    rm(cur.row)
    rm(col.names)
    rm(row.names)
    rm(data.string)
    return(df)
  })
  return(outer.df)
}

env.subsets.matrix <- function(env) {
  outer.df <- with(env, {
    df <- data.frame()

    row.names = c('up','up.bh','up.moderate','up.strong','down','down.bh','down.moderate', 'down.strong')
    col.names = get.subset.prefixes()
  
    for (type in row.names) {
      cur.row <- c()
      for (test in col.names) {
        data.string <- paste(test,type,sep='.')
        if (type == 'up.moderate' || type == 'up.strong') {
          data.string <- paste(test,'up',sep='.')
        } else if (type == 'down.moderate' || type == 'down.strong') {
          data.string <- paste(test,'down',sep='.')
        }
        
        d.sub <- get(data.string)
        
        if (type == 'up.moderate') d.sub <- subset(d.sub, sub_fg > 1)
        else if (type == 'up.strong') d.sub <- subset(d.sub, sub_fg > 1.5)
        else if (type == 'down.moderate') d.sub <- subset(d.sub, sub_fg < 0.2)
        else if (type == 'down.strong') d.sub <- subset(d.sub, sub_fg < 0.1)
  
        cur.row <- c(cur.row,nrow(d.sub))
      }
      df <- rbind(df,cur.row) 
    }
    colnames(df) <- col.names
    rownames(df) <- row.names
    df$total <- nrow(stats.lnl)
    df$subset <- row.names
    df$base <- base
    rm(cur.row)
    rm(col.names)
    rm(row.names)
    return(df)
  })
  return(outer.df)
}

# GO data from http://www.ensembl.org/biomart/martview/202b3d53b5406eb4f9b5f13f4c5dba4c/202b3d53b5406eb4f9b5f13f4c5dba4c/202b3d53b5406eb4f9b5f13f4c5dba4c

compare.tests.by.name <- function(envir,a,b) {
  test.a <- get(a,envir=envir)
  test.b <- get(b,envir=envir)
  res <- compare.tests(test.a,test.b)
  print(paste(a,b,res))
}

compare.tests <- function(test.a,test.b) {
  n.int <- length(intersect(test.a$name,test.b$name))
  n.union <- length(union(test.a$name,test.b$name))
  n.mean <- mean(nrow(test.a),nrow(test.b))
  return(n.int/n.union)
}

correspondence.matrix <- function(h.env,c.env,type='counts') {

  if (type=='counts') {
    # If we're looking at subset counts, make an 'up' and 'down' for each test.
    numbers <- rep(c(1,2,3,4,5,6,7,8,'p'),each=2)
    updown <- c('up','down')
    set.names <- paste(numbers,updown,sep='.')
    h.sets <- lapply(set.names,function(x){
      get(x,envir=h.env)
    })
    c.sets <- lapply(set.names,function(x){
      get(x,envir=c.env)
    })
    all.sets <- c(h.sets,c.sets)
    counts <- matrix(nrow=length(all.sets),ncol=length(all.sets)) 
    rownames(counts) <- c(paste('h_',set.names,sep=''),paste('c_',set.names,sep=''))
    colnames(counts) <- c(paste('h_',set.names,sep=''),paste('c_',set.names,sep=''))
    ret.mat <- counts
  } else {
    # Otherwise, just one row / col per test.
    numbers <- rep(c(1,2,3,4,5,6,7,8),each=1)
    updown <- c('up')    
    set.names <- paste(numbers,updown,sep='.')
    h.stats <- lapply(numbers,function(x){
      lrt.lbl <- paste('lrt',x,'signed',sep='.')
      stats.lnl <- get('stats.lnl',envir=h.env)
      stats.lnl[,'corr.lrt'] <- stats.lnl[,lrt.lbl]
      return(stats.lnl)
    })
    c.stats <- lapply(numbers,function(x){
      lrt.lbl <- paste('lrt',x,'signed',sep='.')
      stats.lnl <- get('stats.lnl',envir=c.env)
      stats.lnl[,'corr.lrt'] <- stats.lnl[,lrt.lbl]
      return(stats.lnl)
    })
    all.sets <- c(h.stats,c.stats)
    correlations <- matrix(nrow=length(all.sets),ncol=length(all.sets))
    rownames(correlations) <- c(paste('h_',set.names,sep=''),paste('c_',set.names,sep=''))
    colnames(correlations) <- c(paste('h_',set.names,sep=''),paste('c_',set.names,sep=''))
    ret.mat <- correlations
  }

  for (i in 1:length(all.sets)) {
    for (j in 1:length(all.sets)) {
      i.set <- all.sets[[i]]
      j.set <- all.sets[[j]]
      if (type == 'counts') {
        ret.mat[i,j] <- compare.tests(i.set,j.set)
      } else {
        ret.mat[i,j] <- compare.test.correlation(i.set,j.set)
      }
    }
  }
  
  return(ret.mat)
}

compare.test.correlation <- function(stats_a,stats_b) {
  stats_a <- stats_a[,c('corr.lrt','name')]
  stats_b <- stats_b[,c('corr.lrt','name')]
  stats_a[,'corr.a'] <- stats_a[,'corr.lrt']
  stats_b[,'corr.b'] <- stats_b[,'corr.lrt']
  merged <- merge(stats_a,stats_b,by='name')
  correlation <- cor(merged$corr.a,merged$corr.b,method='spearman',use='complete.obs')
  return(correlation)
}

# Take a union of all genes to get a combined list.
combined.gene.list <- function(env.list) {
  first <- env.list[[1]]
  first.env <- get(first)
  genes.list <- with(first.env,stats.lnl[,c('Hsap_protein','name')])

  for (i in 2:length(env.list)) {
    other.env <- get(env.list[[i]])
    other.genes <- with(other.env,stats.lnl[,c('Hsap_protein','name')])
    genes.list <- rbind(genes.list,other.genes)
  }

  genes.list <- unique(genes.list)
  return(genes.list)
}

load.env <- function(base) {
  rdata <- paste(base,".Rdata",sep="")
  env <- new.env()
  env$base <- base

  with(env,{load(rdata)})
  get.subsets.for.env(env)
  return(env)
}

load.data.from.directory <- function(dir.name,read.only=F) {
  env.suffixes <- get.env.suffixes()
  var.envs <- c()
  for (env.suffix in env.suffixes) {
    print(env.suffix)
    file <- paste(dir.name,'/','lnl_',env.suffix,sep='')
    cur.env <- load.env(file)
    cur.env$dir.name <- dir.name

    # Save the Rdata with filled p-values etc.
    if (!read.only) {
      save(cur.env,file=paste(file,"_filled.Rdata",sep=''))
    }

    # Plot a correlogram for this data frame.
    stats.lnl <- cur.env$stats.lnl
    if (!read.only) {
      plot.corrgram(stats.lnl,env.suffix)
      # Dump all the subsets.
      dump.subsets.for.env(cur.env)
    }

    # Save the environment in the local scope.
    var.name <- paste(dir.name,'_',env.suffix,sep='')
    assign(var.name,cur.env,envir=.GlobalEnv)
    var.envs <- c(var.envs,var.name)
    print(var.name)
  }  

  # Collect results into test-centric data frames and write to files.
  collate.results(var.envs[1],var.envs[2],read.only=read.only)

  # Dump overlap & correlation matrices.
  d <- correspondence.matrix(get(var.envs[1]),get(var.envs[2]))
  if (!read.only) {
    write.csv(d,file=paste(dir.name,'matrix.counts.csv',sep=''),quote=F)
  }
  var.name <- paste(dir.name,'_overlap',sep='')
  assign(var.name,d,envir=.GlobalEnv)

  d <- correspondence.matrix(get(var.envs[1]),get(var.envs[2]),type='correlation')
  if (!read.only) {
    write.csv(d,file=paste(dir.name,'matrix.correlations.csv',sep=''),quote=F)
  }
  var.name <- paste(dir.name,'_corr',sep='')
  assign(var.name,d,envir=.GlobalEnv)

  # Dump cross-dataset subsets: H,G,C-only genes, and H,G,C-exclusive ones.

}

dump.parallel.subsets <- function(all.data,base.dir) {

  test.fields <- c('3.h','4.h','4.c')
  renamed.fields <- c('g','h','c')

  if (!exists("sig.df")) {
    sig.df <- data.frame()
    for (i in 1:nrow(all.data)) {
      row <- all.data[i,]
      cur.df <- data.frame()
      for (i in 1:length(test.fields)) {
        field <- test.fields[i]
        renamed.field <- renamed.fields[i]
        renamed.lrt <- paste('lrt.',renamed.field,sep='')
        lrt.field <- paste('lrt.',field,sep='')
        p.field <- paste('pval.',field,sep='')
        if (is.na(row[,lrt.field]) || is.na(row[,p.field])) {
          cur.df[1,renamed.field] = 1
          cur.df[1,renamed.lrt] = 0
        } else {
          cur.df[1,renamed.field] = row[1,p.field]
          cur.df[1,renamed.lrt] = row[1,lrt.field]
        }
      }
      cur.df[1,'name'] = row[,'name']
      cur.df[1,'Hsap_protein'] = row[,'Hsap_protein']
     sig.df <- rbind(sig.df,cur.df)
    }
    assign('sig.df',sig.df,envir=.GlobalEnv)
  }

  t <- 0.05
  no.g <- subset(sig.df, `g` > t & `h` < t & `c` < t)
  no.h <- subset(sig.df, `g` < t & `h` > t & `c` < t)
  no.c <- subset(sig.df, `g` < t & `h` < t & `c` > t)
  c.only <- subset(sig.df, `g` > t & `h` > t & `c` < t)
  h.only <- subset(sig.df, `g` > t & `h` < t & `c` > t)
  g.only <- subset(sig.df, `g` < t & `h` > t & `c` > t)
  h.no.c <- subset(sig.df, `h` < t & `c` > t)
  c.no.h <- subset(sig.df, `h` > t & `c` < t)

  out.with.it <- function(genes,file) {
    print(file)

    all.genes <- data.frame(Hsap_protein=all.tests$Hsap_protein,name=all.tests$name)
    merged <- merge(all.genes,genes,by=c('name'),all.x=TRUE)
    merged[is.na(merged$lrt.signed),]$lrt.signed <- 0

    merged <- merged[,c('Hsap_protein.x','name','lrt.signed')]
    names(merged) <- c('Hsap_protein','name','lrt.signed')
    merged <- merged[order(-merged$lrt.signed),]
    write.csv(merged,file=file,quote=F)
  }

  # Human unique
  tmp <- subset(sig.df,lrt.h > 0)
  tmp$lrt.signed <- with(tmp, pmax(0,lrt.h - pmax(lrt.g,lrt.c,0)))
  tmp.b <- subset(sig.df,lrt.h < 0)
  tmp.b$lrt.signed <- with(tmp.b, pmin(0,lrt.h - pmin(lrt.g,lrt.c,0)))
  tmp <- rbind(tmp,tmp.b)
  out.with.it(tmp,file=paste(base.dir,'/human_sorted.csv',sep=''))

  # Chimp unique
  tmp <- subset(sig.df,lrt.c > 0)
  tmp$lrt.signed <- with(tmp, pmax(0,lrt.c - pmax(lrt.g,lrt.h,0)))
  tmp.b <- subset(sig.df,lrt.c < 0)
  tmp.b$lrt.signed <- with(tmp.b, pmin(0,lrt.c - pmin(lrt.g,lrt.h,0)))
  tmp <- rbind(tmp,tmp.b)
  out.with.it(tmp,file=paste(base.dir,'/chimp_sorted.csv',sep=''))

  # Gorilla unique
  tmp <- subset(sig.df,lrt.g > 0)
  tmp$lrt.signed <- with(tmp, pmax(0,lrt.g - pmax(lrt.h,lrt.c,0)))
  tmp.b <- subset(sig.df,lrt.g < 0)
  tmp.b$lrt.signed <- with(tmp.b, pmin(0,lrt.g - pmin(lrt.h,lrt.c,0)))
  tmp <- rbind(tmp,tmp.b)
  out.with.it(tmp,file=paste(base.dir,'/gorilla_sorted.csv',sep=''))

  # Human & Gorilla parallelism
  tmp <- subset(sig.df,lrt.g > 0 & lrt.h > 0)
  tmp$lrt.signed <- with(tmp, pmax(0,pmin(lrt.g,lrt.h) - pmax(lrt.c,0)))
  tmp.b <- subset(sig.df,lrt.g < 0 & lrt.h < 0)
  tmp.b$lrt.signed <- with(tmp.b, pmin(0,pmax(lrt.g,lrt.h) - pmin(lrt.c,0)))
  tmp.c <- subset(sig.df,(lrt.g < 0 & lrt.h > 0) | (lrt.g > 0 & lrt.h < 0))
  tmp.c$lrt.signed <- 0
  tmp <- rbind(tmp, tmp.b, tmp.c)
  out.with.it(tmp,file=paste(base.dir,'/gh_parallel_sorted.csv',sep=''))  

  # Human & Chimp parallelism
  tmp <- subset(sig.df,lrt.c > 0 & lrt.h > 0)
  tmp$lrt.signed <- with(tmp, pmax(0,pmin(lrt.c,lrt.h) - pmax(lrt.g,0)))
  tmp.b <- subset(sig.df,lrt.c < 0 & lrt.h < 0)
  tmp.b$lrt.signed <- with(tmp.b, pmin(0,pmax(lrt.c,lrt.h) - pmin(lrt.g,0)))
  tmp.c <- subset(sig.df,(lrt.c < 0 & lrt.h > 0) | (lrt.c > 0 & lrt.h < 0))
  tmp.c$lrt.signed <- 0
  tmp <- rbind(tmp, tmp.b, tmp.c)
  out.with.it(tmp,file=paste(base.dir,'/hc_parallel_sorted.csv',sep=''))  

  # Gorilla & Chimp parallelism
  tmp <- subset(sig.df,lrt.g > 0 & lrt.c > 0)
  tmp$lrt.signed <- with(tmp, pmax(0,pmin(lrt.g,lrt.c) - pmax(lrt.h,0)))
  tmp.b <- subset(sig.df,lrt.g < 0 & lrt.c < 0)
  tmp.b$lrt.signed <- with(tmp.b, pmin(0,pmax(lrt.g,lrt.c) - pmin(lrt.h,0)))
  tmp.c <- subset(sig.df,(lrt.g < 0 & lrt.c > 0) | (lrt.g > 0 & lrt.c < 0))
  tmp.c$lrt.signed <- 0
  tmp <- rbind(tmp, tmp.b, tmp.c)
  out.with.it(tmp,file=paste(base.dir,'/gc_parallel_sorted.csv',sep=''))  


  tmp <- no.g
  tmp[,'lrt.signed'] <- tmp[,'lrt.h']
  out.with.it(tmp,file=paste(base.dir,'/no_g_sorted.csv',sep=''))

  tmp <- no.h
  tmp[,'lrt.signed'] <- tmp[,'lrt.c']
  out.with.it(tmp,file=paste(base.dir,'/no_h_sorted.csv',sep=''))

  tmp <- no.c
  tmp[,'lrt.signed'] <- tmp[,'lrt.h']
  out.with.it(tmp,file=paste(base.dir,'/no_c_sorted.csv',sep=''))

  tmp <- c.only
  tmp[,'lrt.signed'] <- tmp[,'lrt.h']
  out.with.it(tmp,file=paste(base.dir,'/only_c_sorted.csv',sep=''))

  tmp <- h.only
  tmp[,'lrt.signed'] <- tmp[,'lrt.h']
  out.with.it(tmp,file=paste(base.dir,'/only_h_sorted.csv',sep=''))

  tmp <- g.only
  tmp[,'lrt.signed'] <- tmp[,'lrt.g']
  out.with.it(tmp,file=paste(base.dir,'/only_g_sorted.csv',sep=''))

  tmp <- h.no.c
  tmp[,'lrt.signed'] <- tmp[,'lrt.h']
  out.with.it(tmp,file=paste(base.dir,'/h_no_c_sorted.csv',sep=''))

  tmp <- c.no.h
  tmp[,'lrt.signed'] <- tmp[,'lrt.c']
  out.with.it(tmp,file=paste(base.dir,'/c_no_h_sorted.csv',sep=''))

}

# Run this command to load the datasets and their sig. subsets
#load.data.from.directory("fifth_runs")

library(plyr)
library(doBy)

get.subset.prefixes <- function() {
  return(c(1,2,3,4,5,6,7,8,'p'))
}

get.env.suffixes <- function() {
  env.suffixes <- list(
   'm_Hsap','m_Ptro'
  )
  return(env.suffixes)
}

get.env.names <- function(dir.name) {
  env.suffixes <- get.env.suffixes()
  env.names <- paste(dir.name,env.suffixes,sep='_')
  return(env.names)
}

# This collects all the p-values from the different datasets into a gene-centric list.
collate.results <- function(h.env.name,c.env.name,read.only=F) {
  env.names <- list(h.env.name,c.env.name)
  all.genes <- combined.gene.list(env.names)

  all.df <- all.genes
  print(head(all.df))
  for (test in c(1,2,3,4,5,6,7,8)) {
    test.df <- all.genes
    pval.label <- paste('pval',test,sep='.')
    lrt.label <- paste('lrt',test,'signed',sep='.')
    fg.label <- paste('fg',test,sep='.')
    bg.label <- paste('bg',test,sep='.')
    print(pval.label)
    for (env.name in env.names) {
      this.env <- get(env.name)
      env.short <- 'h'
      if (env.name == c.env.name) {
        env.short <- 'c'
      }
      dir.name <- this.env$dir.name
      if (test == 3 || test == 4) {
        #print(head(this.env$stats.lnl))
        print(fg.label)
      }
      genes.df <- this.env$stats.lnl[,c(
        pval.label,
        lrt.label,
        fg.label,
        bg.label,
        'name'
      )]

      colnames(genes.df) <- c(
        paste('pval',test,env.short,sep='.'),
        paste('lrt',test,env.short,sep='.'),
        paste('fg',test,env.short,sep='.'),
        paste('bg',test,env.short,sep='.'),
        'name'
      )
      test.df <- merge(test.df,genes.df,by='name',all.x=TRUE,all.y=FALSE)
      all.df <- merge(all.df,genes.df,by='name',all.x=TRUE,all.y=FALSE)
    }

    assign(paste('test',test,sep='.'),test.df,envir=.GlobalEnv)
    if (!read.only) {
      write.csv(test.df,file=paste(dir.name,"/collated",'.',test,".csv",sep=''),row.names=F)
    }
  }

  assign('all.tests',all.df,envir=.GlobalEnv)
  if (!read.only) {
    write.csv(all.df,file=paste(dir.name,"/collated.all.csv",sep=''),row.names=F)
  }
}

collect.counts <- function(dir.name) {
  env.names = get.env.names(dir.name)

  all.counts <- data.frame()
  for (env.name in env.names) {
    cur.df <- env.subsets.matrix(env.name)
    rownames(cur.df) <- NULL
    all.counts <- rbind(all.counts,cur.df)
  }
  return(all.counts)
}

plot.corrgram <- function(data,name) {
  library(corrgram)

  data <- data[,c('grantham_Ggor','subs_n_Ggor','subs_s_Ggor','grantham_other','subs_n_other','subs_s_other',
  'pval.1','pval.2','pval.3','pval.4','pval.5','pval.6','pval.7','pval.8','gc_3','gc_cds','gc_genomic')]

  png(file=paste('correlation',name,"png",sep='.'),width=1600,height=1600)

  nums <- unlist(lapply(data,is.numeric))
  print(nums)
  data <- data[,nums]

  corrgram(data,order=NULL,upper.panel=NULL)
  dev.off()
}

dump.subsets.for.env <- function(env) {
  subset.matrix <- env.subsets.matrix(env)
  write.csv(subset.matrix,file=paste(env$base,'_subsets.csv',sep=''),row.names=F,quote=F)

  with(env,{
    write.csv(stats.lnl,file=paste(base,'_all','.csv',sep=''),row.names=F,quote=F)
    write.table(stats.lnl[,c('name')],file=paste(base,'_all_names','.txt',sep=''),row.names=F,quote=F)
    write.table(stats.lnl[,c('Hsap_gene')],file=paste(base,'_all_genes','.txt',sep=''),row.names=F,quote=F)
    write.table(stats.lnl[,c('Hsap_protein')],file=paste(base,'_all_proteins','.txt',sep=''),row.names=F,quote=F)
    for (test in c(1,2,3,4,5,6,7,8,'p')) {

      # Output a list of all genes ordered by the p-value.
      if (test != 'p') {
        all.genes <- stats.lnl
        signed.stat.label <- paste('lrt.',test,'.signed',sep='')
        fg.omega.label <- paste('fg.',test,sep='')
        bg.omega.label <- paste('bg.',test,sep='')
        pval.label <- paste('pval.',test,sep='')
        all.genes[,'lrt.signed'] = all.genes[,signed.stat.label]
        all.genes[,'pval'] = all.genes[,pval.label]
        all.genes[,'fg.omega'] = all.genes[,fg.omega.label]
        all.genes[,'bg.omega'] = all.genes[,bg.omega.label]
        all.genes <- orderBy(~ -lrt.signed,data=all.genes)
        all.genes <- all.genes[,c('name','Hsap_protein','lrt.signed','pval','fg.omega','bg.omega')]
        print(head(all.genes,n=5))
        write.csv(all.genes,file=paste(base,'_',test,'_sorted.csv',sep=''),row.names=F,quote=F)
      }

      # Output the subset list
      for (dir in c('up','down','up.bh','down.bh')) {
        df.name <- paste(test,dir,sep='.')
        print(df.name)
        cur.df <- get(df.name)
        cur.df <- orderBy(~ sub_pval,data=cur.df)
        write.csv(cur.df,file=paste(base,'_',df.name,'_table.csv',sep=''),row.names=F,quote=F)
        write.table(cur.df[,c('Hsap_protein')],file=paste(base,'_',df.name,'_proteins.txt',sep=''),row.names=F,col.names=F,quote=F)
        write.table(cur.df[,c('Hsap_tx')],file=paste(base,'_',df.name,'_txs.txt',sep=''),row.names=F,col.names=F,quote=F)
        write.table(cur.df[,c('Hsap_gene')],file=paste(base,'_',df.name,'_genes.txt',sep=''),row.names=F,col.names=F,quote=F)
        write.table(cur.df[,c('name')],file=paste(base,'_',df.name,'_names.txt',sep=''),row.names=F,col.names=F,quote=F)
      }
    }
  })
}

plot.omegas.new <- function(all.tests,which.tests) {

  rows <- data.frame()

  if (which.tests == 'lineages') {
    tests <- c('3.h','4.h','4.c')
    test.labels <- c('Gorilla','Human','Chimpanzee')
  } else {
    tests <- c('3.h','7.h','8.h')
    test.labels <- c('Gorilla','AGA branch','AGA clade')
    
  }

   for (test in tests) {
    print(test)
    fg <- all.tests[,paste('fg.',test,sep='')]
    lrt <- all.tests[,paste('lrt.',test,sep='')]
    bg <- all.tests[,paste('bg.',test,sep='')]
    dr <- sign(lrt+0.001)

    cur.df <- data.frame(
      fg=fg,
      lrt=lrt,
      bg=bg,
      dr=dr,
      test=test
    )
    cur.df <- subset(cur.df, !is.na(fg) & !is.na(dr))
    rows <- rbind(rows,cur.df)    
  }

  rows$test <- factor(rows$test,labels=test.labels)

  t <- 10
  rows[rows$lrt < -t,]$lrt = -t
  rows[rows$lrt > t,]$lrt = t

  pdf(file="all_lrt.pdf")
  p <- ggplot(data=rows,aes(lrt,fill=test))
  p <- p + geom_bar(position="dodge",aes(group=test),binwidth=0.5)
  print(p)
  dev.off()

  pdf(file="all_fg.pdf")
  p <- ggplot(data=rows,aes(x=fg,fill=test))
  p <- p + geom_bar(position="dodge",aes(group=test),binwidth=0.5)
  p <- p + scale_x_log10()
  print(p)
  dev.off()
}


plot.test.5 <- function(data,dir='up') {
  if (dir == 'up') {
    data <- subset(data,lrt.5.h > 0)
  }  

  sig <- subset(data,pval.5.h < 0.05)
  sig.lm <- lm(lrt.h ~ lrt.g,data=sig)
  nonsig <- subset(data,pval.5.h > 0.05)
  nonsig.lm <- lm(lrt.h ~ lrt.g,data=nonsig)
  print(sig.lm)
  print(nonsig.lm)
  #data$sig.in.five <- as.numeric( data[,'pval.5.h'] < 0.05 )

  p <- ggplot(data,aes(x=pval.g+.01,y=pval.h+.01,colour=pval.c+.01))
  p <- p + geom_point(alpha=0.6,size=2)
  p <- p + scale_fill_brewer()
  p <- p + scale_y_log10() + scale_x_log10()
  pdf(file="test.pdf")
  print(p)
  dev.off()

}

aln.and.subs <- function(env,gene.name) {
  env$gene.name <- gene.name

  with(env, {
    sub.rows <- subset(stats.lnl, name == gene.name)
    if (nrow(sub.rows) == 0) {
      print(paste("Gene ",gene.name," not found!"))
      return()
    }
    row <- sub.rows

    f <- row$filtered_aln_file    
    print(row$filtered_aln_file)    
    alns.file <- paste(dir.name,"/alns/alns.tar",sep="")
    cmd <- paste("tar -x",f,"-f",alns.file)
    print(cmd)
    system(cmd)

    library(phylosim)
    sim <- PhyloSim()
    readAlignment(sim,f)
    pdf(paste(dir.name,'/alns/',f,'.pdf',sep=''))
    plot(sim,plot.chars=F)
    dev.off()

  })
  env$gene.name <- NULL
}

manuscript.table <- function(all.genes,human.env,chimp.env,go.ens) {
  h.stats <- human.env$stats.lnl
  c.stats <- chimp.env$stats.lnl

  all.names <- all.genes[,c('Hsap_protein','name')]
  int.names <- c('ACOX1','ADAMTS2','ADCY10','AC005841.1')
  intr.row <- subset(all.names, name %in% int.names)
#  all.names <- intr.row

  final.df <- ddply(all.names,.(name),function(x) {
    cur.name <- x$name
    print(cur.name)
    all.row <- subset(all.genes,name == cur.name)
    h.row <- subset(h.stats,name == cur.name)
    c.row <- subset(c.stats,name == cur.name)

    if (nrow(h.row) > 0) {
      terms <- as.character(unlist(go.ens[h.row$Hsap_protein]))
    } else if (nrow(c.row) > 0) {
      terms <- as.character(unlist(go.ens[c.row$Hsap_protein]))
    } else {
      terms <- c()
    }
    term.presence.absence <- group.presence.absence(cur.name,terms)

    lrts <- data.frame(
      lrt.score.gorilla = all.row$lrt.3.h,
      lrt.score.human = all.row$lrt.4.h,
      lrt.score.chimpanzee = all.row$lrt.4.c,
      lrt.pval.gorilla = all.row$pval.3.h,
      lrt.pval.human = all.row$pval.4.h,
      lrt.pval.chimpanzee = all.row$pval.4.c,
      dnds.gorilla = all.row$fg.3.h,
      dnds.human = all.row$fg.4.h,
      dnds.chimpanzee = all.row$fg.4.c
    )

    #trms <- paste(as.character(unlist(terms)),collapse=',')
    #go <- data.frame(go.terms = trms)

    human.stats <- data.frame(xyz=1)
    if (nrow(h.row) > 0) {
      human.stats <- data.frame(
        gc.3=h.row$gc_3,
        gc.cds=h.row$gc_cds,
        subs.grantham.gorilla=h.row$grantham_Ggor,
        subs.grantham.human=h.row$grantham_other,
        subs.n.human=h.row$subs_n_other,
        subs.n.gorilla=h.row$subs_n_Ggor,
        subs.s.human=h.row$subs_s_other,
        subs.s.gorilla=h.row$subs_s_Ggor,
        codons=h.row$orig_aln_length,
        ens.protein.id=h.row$Hsap_protein,
        ens.gene.id=h.row$Hsap_gene,
        dnds.all = h.row$a_omega_0
      )
    }

    chimp.stats <- data.frame(xyz=1)
    if (nrow(c.row) > 0) {
      chimp.stats <- data.frame(
        subs.grantham.chimp=c.row$grantham_other,
        subs.n.chimp=c.row$subs_n_other,
        subs.s.chimp=c.row$subs_s_other
      )
    }
    return(data.frame(
      term.presence.absence
      ,lrts
      ,human.stats
      ,chimp.stats
#     ,go
    ))
  })

  sel <- sort(names(final.df))
  final.df <- final.df[,sel]
  final.df <- format.numeric.df(final.df)
  final.df[,'xyz'] <- NULL

  write.csv(final.df,file="manuscript_table.csv",row.names=F,quote=F)

  print(head(final.df))

  assign('final.df',final.df,envir=.GlobalEnv)
}

format.numeric.df <- function(x,digits=3) {  
  nums <- unlist(lapply(x,is.numeric))
  print(nums)
  for (col in names(nums)) {
    x[,col] <- format(x[,col],digits)
  }

  return(x)
}