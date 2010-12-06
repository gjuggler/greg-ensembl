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
  stats.lnl = subset(stats.lnl,(subs_n_Ggor + subs_s_Ggor) >= 4)

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

    row.names = c('up','up.bh','up.bh.above','down','down.bh','down.bh.below')
    col.names = get.subset.prefixes()
  
    for (type in row.names) {
      cur.row <- c()
      for (test in col.names) {
        data.string <- paste(test,type,sep='.')
        if (type == 'up.bh.above' || type == 'up.bh.above.15') {
          data.string <- paste(test,'up.bh',sep='.')
        } else if (type == 'down.bh.below' || type == 'down.bh.below.05') {
          data.string <- paste(test,'down.bh',sep='.')
        }
        
        d.sub <- get(data.string)
        
        if (type == 'up.bh.above') d.sub <- subset(d.sub, sub_fg > 1)
        else if (type == 'up.bh.above.15') d.sub <- subset(d.sub, sub_fg > 1.5)
        else if (type == 'down.bh.below') d.sub <- subset(d.sub, sub_fg < 0.1)
        else if (type == 'down.bh.below.05') d.sub <- subset(d.sub, sub_fg < 0.05)
  
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

test.correspondence <- function(h.env,c.env) {
    compare.tests.by.name(h.env,'2.up','3.up')
    compare.tests.by.name(h.env,'2.down','3.down')
    compare.tests.by.name(h.env,'1.up','4.up')
    compare.tests.by.name(h.env,'1.down','4.down')
    compare.tests.by.name(h.env,'7.up','8.up')
    compare.tests.by.name(h.env,'7.down','8.down')

  h.gor.up <- get('3.up',envir=h.env)
  c.gor.up <- get('3.up',envir=c.env)
  h.gor.down <- get('3.down',envir=h.env)
  c.gor.down <- get('3.down',envir=c.env)
  compare.tests.by.name(environment(),'h.gor.up','c.gor.up')
  compare.tests.by.name(environment(),'h.gor.down','c.gor.down')
}

correspondence.matrix <- function(h.env,c.env) {
  numbers <- rep(c(1,2,3,4,6,7,8),each=2)
  updown <- c('up','down')
  set.names <- paste(numbers,updown,sep='.')
  h.sets <- lapply(set.names,function(x){get(x,envir=h.env)})
  c.sets <- lapply(set.names,function(x){get(x,envir=c.env)})

  all.sets <- c(h.sets,c.sets)

  data <- matrix(nrow=length(all.sets),ncol=length(all.sets))
  rownames(data) <- c(paste('h_',set.names,sep=''),paste('c_',set.names,sep=''))
  colnames(data) <- c(paste('h_',set.names,sep=''),paste('c_',set.names,sep=''))

  for (i in 1:length(all.sets)) {
    for (j in 1:length(all.sets)) {
      i.set <- all.sets[[i]]
      j.set <- all.sets[[j]]
      data[i,j] <- compare.tests(i.set,j.set)
    }
  }
  return(data)
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

load.data.from.directory <- function(dir.name) {
  env.suffixes <- get.env.suffixes()
  var.envs <- c()
  for (env.suffix in env.suffixes) {
    print(env.suffix)
    file <- paste(dir.name,'/','lnl_',env.suffix,sep='')
    cur.env <- load.env(file)
    # Save the Rdata with filled p-values etc.
    save(cur.env,file=paste(file,"_filled.Rdata",sep=''))

    # Dump all the subsets.
    dump.subsets.for.env(cur.env)

    # Save the environment in the local scope.
    var.name <- paste(dir.name,'_',env.suffix,sep='')
    assign(var.name,cur.env,envir=.GlobalEnv)
    var.envs <- c(var.envs,var.name)
    print(var.name)
  }  

  # Dump correspondence matrix.
  d <- correspondence.matrix(get(var.envs[1]),get(var.envs[2]))
  write.csv(d,file=paste(dir.name,'correspondence.matrix.csv',sep=''),quote=F)

}

# Run this command to load the datasets and their sig. subsets
#load.data.from.directory("fifth_runs")

library(plyr)
library(doBy)

get.subset.prefixes <- function() {
  return(c(1,2,3,4,6,7,8,'p'))
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
collate.results <- function(dir.name) {
  env.names <- get.env.names(dir.name)
  # Get collated results for tests 3, 4, 7, and 8
  all.genes <- combined.gene.list(env.names)

  for (test in c(1,2,3,4,6,7,8)) {
    all.df <- all.genes
    pval.label <- paste('pval',test,sep='.')
    lrt.label <- paste('pval',test,sep='.')
    print(pval.label)
    for (env.name in env.names) {
      this.env <- get(env.name)
      genes.df <- this.env$stats.lnl[,c('Hsap_protein',pval.label,lrt.label,'name')]

      colnames(genes.df) <- c('Hsap_protein',paste('pval',env.name,sep='.'),paste('lrt',env.name,sep='.'),'name')

      all.df <- merge(all.df,genes.df,all.x=TRUE,all.y=FALSE)
    }

    all.df <- orderBy(~pval.sixth_runs_gt4_m_Hsap,data=all.df)

    assign(paste('test',test,sep='.'),all.df,envir=.GlobalEnv)
    write.csv(all.df,file=paste("test",test,"csv",sep='.'),row.names=F)
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
  'pval.1','pval.2','pval.3','pval.4','pval.6','pval.7','pval.8','gc_3','gc_cds','gc_genomic')]

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
    for (test in c(1,2,3,4,6,7,8,'p')) {

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

plot.omegas <- function(Hsap.env,Ptro.env) {
  library(ggplot2)

  Hsap <- Hsap.env$stats.lnl
  Ptro <- Ptro.env$stats.lnl

  Hsap[,'lrt.pt.signed'] <- 0
  Hsap[,'fg.pt'] <- 0
  Hsap[,'bg.pt'] <- 0

  Ptro[,'lrt.pt.signed'] <- Ptro[,'lrt.4.signed']
  Ptro[,'fg.pt'] <- Ptro[,'fg.4']
  Ptro[,'bg.pt'] <- Ptro[,'bg.4']

  df <- data.frame()
  for (test in c(1, 2, 3, 'pt', 4, 7, 8)) {
    lrt.lbl <- paste('lrt.',test,'.signed',sep='')
    fg.lbl <- paste('fg.',test,sep='')
    bg.lbl <- paste('bg.',test,sep='')

    rows <- Hsap
    if (test == 'pt') {
      rows <- Ptro
    }

    rows$plot.test <- test

    rows$sig <- abs(rows[,lrt.lbl])
    rows <- subset(rows,!is.na(sig))
    rows$lrt <- rows[,lrt.lbl]
    rows$dr <- sign(rows[,lrt.lbl]-0.0001)
#    rows[rows[,lrt.lbl] > 6,'dr'] <- 2
#    rows[rows[,lrt.lbl] < -6,'dr'] <- -2

    rows$fg <- rows[,fg.lbl]
    rows$bg <- rows[,bg.lbl]
    #rows <- subset(rows, fg > 0.01)

    df <- rbind(df,rows)
  }

  df$plot.test <- factor(df$plot.test,labels=c('hsap','ggor','Ggor','Ptro','Hsap','GA-b','GA-c'))
  df$dr <- factor(df$dr,labels=c('Decelerated','Accelerated'))

  pdf(file="test.pdf")
  p <- ggplot(data=df,aes(x=fg,colour=dr,fill=dr))
  p <- p + geom_bar(position="dodge",aes(group=dr),binwidth=0.2)
  p <- p + scale_x_log10()
  p <- p + facet_grid(plot.test ~ .)

  print(p)
  dev.off()

  df2 <- df
  df2[df2$lrt > 6,]$lrt = 5.9
  df2[df2$lrt < -6,]$lrt = -5.9

  pdf(file="test2.pdf")
  p <- ggplot(data=df2,aes(x=lrt,fill=dr,colour=dr))
  p <- p + geom_bar(position="dodge",aes(group=dr),binwidth=0.5)
  p <- p + scale_x_continuous(limits=c(-6,6))
  p <- p + facet_grid(plot.test ~ .)
  print(p)
  dev.off()

  pdf(file="test3.pdf")

  df2 <- df2[order(df2$sig),]

  p <- ggplot(data=df2,aes(x=bg,y=fg,colour=lrt))
  p <- p + geom_point(size=0.5,alpha=1)
  p <- p + scale_y_log10() + scale_x_log10()
  p <- p + facet_grid(plot.test ~ .)
  print(p)
  dev.off()
}

print.alns <- function(rows) {
  #old.wd <- getwd()
  #setwd('alns')
  for (i in 1:nrow(rows)) {
    row <- rows[i,]
    f <- row$filtered_aln_file
    
    print(row$filtered_aln_file)
    
    cmd <- paste("tar -x",f,"-f ./alns.tar")
    print(cmd)
    system(cmd)

    library(phylosim)
    sim <- PhyloSim()
    readAlignment(sim,f)
    pdf(paste(f,'.pdf',sep=''))
    plot(sim,plot.chars=F)
    dev.off()
  }
  #setwd(old.wd)
}