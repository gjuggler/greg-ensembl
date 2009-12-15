if (!exists('dup_data')) {load('~/scratch/dup_data.Rdata')}
if (!exists('dup_tags')) {load('~/scratch/dup_tags.Rdata')}

neg = c('negative4','negative3','negative2','negative1')
pos = c('positive4','positive3','positive2','positive1')

# Massage the duplication age values.
if (TRUE) {
  tags = as.character(unique(dup_tags$tag))
  dup.df = data.frame(node_id=as.character(unique(dup_tags$node_id)))
  for (t in tags) {
    subs = subset(dup_tags,tag == t)
    tag.df = data.frame(node_id=as.character(subs$node_id))
    tag.df[[t]] = subs$value
    #print(names(tag.df))
    dup.df = merge(dup.df,tag.df)
    #print(dup.df[1:3,])
  }
}

# Collate the site-wise estimates into a few summary values.
if (TRUE) {
  dup_data = as.data.frame(dup_data)
  #a = dup_data[sample(nrow(dup_data),size=100000),]
  a = dup_data
  a = subset(a, c1_omega < 50 & c2_omega < 50 & c1_ncod > 20 & c2_ncod > 20)
  a$node_id = as.character(a$node_id)

  # Overall correlation of c1 and c2 omegas.
  df.cor <- function(x) cor(x$c1_omega,x$c2_omega)
  c = by(a,a$node_id,df.cor)
  d = data.frame(node_id=as.character(names(c)),cor=as.numeric(c))

  # Number of 'shifts' in omega values (normalized by aln length)
  df.shift <- function(x) {
    subs = subset(x,c1_type != '' & c2_type != '')
    shift = subset(subs,(c1_type %in% pos & c2_type %in% neg) | (c1_type %in% neg & c2_type %in% pos))
    return(nrow(shift) / max(1,nrow(subs)))
  }
  c = by(a,a$node_id,df.shift)
  d$shift = as.numeric(c)

  df.count <- function(x) {nrow(x)}
  c = by(a,a$node_id,df.count)
  d$count = as.numeric(c)

  # Create a data frame with rows as unique node_ids.
  d = merge(d,unique(a[,c('node_id','c1_id','c2_id')]),by='node_id')
  # Merge with the per-node tag info.
  data = merge(dup.df,d,by='node_id')
}

# Collect the chromosome location for the 'left' and 'right' duplicates.
if (TRUE) {
  c1_chr = data.frame(node_id=dup.df$node_id,c1_chr=dup.df$dupldiv_human_chr)
  c2_chr = data.frame(node_id=dup.df$node_id,c2_chr=dup.df$dupldiv_human_chr)
  data = merge(data,c1_chr,by.x='c1_id',by.y='node_id')
  data = merge(data,c2_chr,by.x='c2_id',by.y='node_id')
}


# Plots!
if (TRUE) {

  png(file="dupldiv_xyplot.png")
  plot(x=dup_data$c1_omega,y=dup_data$c2_omega,pch='.',col=rgb(0,0,0,.3),log='xy')
  dev.off()

  pdf(file="dupldiv_chrplot.pdf")
  good_chrs = subset(data,c1_chr != '' & c2_chr != '')
  plot.chrs = function(df,threshold) {
    df = subset(df, dupldiv_bl_max <= threshold)
    df$c1_chr = as.numeric(factor(df$c1_chr))
    df$c2_chr = as.numeric(factor(df$c2_chr))
    plot(jitter(df$c1_chr),jitter(df$c2_chr),pch='.',col=rgb(0,0,0,.5))
  }
  avg_dist = median(as.numeric(good_chrs$dupldiv_bl_max))
  plot.chrs(good_chrs,avg_dist*2)
  plot.chrs(good_chrs,avg_dist)
  plot.chrs(good_chrs,avg_dist/2)
  dev.off()

  # Compare mean correlations of duplications from different taxonomic nodes.
  pdf(file="dupldiv_taxon_age.pdf")
  df.meancor = function(df) {mean(df$cor,na.rm=T)}
  c = by(data,data$dupldiv_taxon_name,df.meancor)

  big_dup_data = merge(a,dup.df,by='node_id')
  c = by(big_dup_data,big_dup_data$dupldiv_taxon_name,df.cor)
  d = by(big_dup_data,big_dup_data$dupldiv_taxon_name,function(x){nrow(x)})

  taxon_correlations = data.frame(dupldiv_taxon_name=as.character(names(c)),cor=as.numeric(c),count=as.numeric(d))
  taxon_ages = unique(data[,c('dupldiv_taxon_name','dupldiv_taxon_mya')])
  taxon_data = merge(taxon_correlations,taxon_ages,by='dupldiv_taxon_name')
  taxon_data = subset(taxon_data,count > 1000)
  plot(taxon_data$cor ~ taxon_data$dupldiv_taxon_mya,ylim=c(0,0.5))
  dev.off()

  # Look at correlation of mammalian sub-trees as a reference.
  if (!exists('sites')) {load('~/scratch/ensembl_sites.tsv.Rdata')}
  #if (!exists('sites')) {load('~/scratch/sites_genes.Rdata')}
  good_sites = subset(sites,!is.na(g_omega) & !is.na(l_omega) & !is.na(p_omega))
  good_sites = subset(sites,g_omega < 20 & l_omega < 20)
#  good_sites = sample(nrow(good_sites),n=100000)

  good_dup = subset(big_dup_data,c1_omega < 20 & c2_omega < 20)
  good_dup = subset(good_dup,c1_ncod > 25 & c2_ncod > 25)
#  good_dup = good_dup[1:100000,]
  cols = rainbow(n=4,alpha=0.05)
  col.indices = as.numeric(good_dup$dupldiv_bl_max) * 2
  col.indices = pmin(col.indices,4)
  col.indices = pmax(col.indices,0)
  col.indices = ceiling(col.indices)
  dup.col = cols[col.indices]
  dup.col = rgb(1,0,0,0.05)

  col.indices = (as.numeric(sites$v_ncod) - 20) * 2
  col.indices = pmin(col.indices,4)
  col.indices = pmax(col.indices,0)
  col.indices = ceiling(col.indices)
  site.col = cols[col.indices]
  site.col = rgb(0,0,1,0.05)

  png(file='dupldiv_paralog_ortholog.png',width=960,height=480)
  par(mfrow=c(1,2))
  plot(good_sites$p_omega,good_sites$l_omega,log='xy',col=site.col,pch='.',xlim=c(0.01,5),ylim=c(0.01,5))
  plot(good_dup$c1_omega,good_dup$c2_omega,log='xy',col=dup.col,pch='.',xlim=c(0.01,5),ylim=c(0.01,5))
  dev.off()

}