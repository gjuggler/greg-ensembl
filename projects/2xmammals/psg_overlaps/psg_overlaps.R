get.all.scores = function(merged.genes,score.column='m.dnds.m0') {

  merged.genes$score = merged.genes[,score.column]

  assign('roc.kosiol',compare.kosiol.eslr(merged.genes),pos=.GlobalEnv)
  assign('roc.nielsen',compare.nielsen.eslr(merged.genes),pos=.GlobalEnv)
  assign('roc.clark',compare.clark.eslr(merged.genes),pos=.GlobalEnv)
  assign('roc.sabeti',compare.sabeti.eslr(merged.genes),pos=.GlobalEnv)
}


output.kosiol.locations = function(genes) {
  a = genes[,c('chr.hg19','start.hg19','end.hg19','strand.x')]
  #print(a[1:10,])
  write.table(a,file="kosiol_locs.txt",sep="\t",col.names=F,row.names=F,quote=F)
}

merge.kosiol.lifted = function(kosiol.genes) {
  # Coming from Galaxy liftOver output.
  data = read.delim("kosiol_hg19.bed",sep="\t",header=F)
  lifted = data[,1:4]
  colnames(lifted) = c("chr.hg19","start.hg19","end.hg19","name")

  merged = merge(kosiol.genes,lifted,by.x="transcriptid",by.y="name")

  return(merged)
}

merge.kosiol.mapped = function(genes) {
  mapped <- read.delim(file="kosiol_ensg.txt",header=F)
  genes$ensg.hg19 <- as.vector(mapped[,1])
  return(genes)
}

compare.kosiol.eslr = function(eslr.merged) {

  kosiol.genes.data = read.delim(file="kosiol-gene-definitions.txt",sep="\t",header=T,stringsAsFactors=F)
  kosiol.genes.table = read.delim(file="kosiol-gene-table.tsv",sep="\t",header=T,stringsAsFactors=F)

  kosiol.genes = merge(kosiol.genes.data,kosiol.genes.table,by.x='transcriptid',by.y='name')
  
  kosiol.lifted = merge.kosiol.lifted(kosiol.genes)

  if (!exists('kosiol.mapped')) {
    output.kosiol.locations(kosiol.lifted)  
    cmd = 'perl ~/src/greg-ensembl/scripts/get_genes_by_location.pl < kosiol_locs.txt > kosiol_ensg.txt'
    system(cmd)
    
    kosiol.mapped <- merge.kosiol.mapped(kosiol.lifted)
    assign('kosiol.mapped',kosiol.mapped,pos=.GlobalEnv)
  }
  
  truth <- kosiol.mapped$score > 400
  kosiol.truth <- data.frame(
    id=kosiol.mapped$transcriptid,
    chr=kosiol.mapped$chr.hg19,
    start=kosiol.mapped$start.hg19,
    end=kosiol.mapped$end.hg19,
    ensg=as.character(kosiol.mapped$ensg.hg19),
    truth=as.integer(truth)
    )

  eslr.scores <- get.eslr.score.df(eslr.merged)
  print(paste(nrow(kosiol.truth),nrow(eslr.scores)))

  eslr.kosiol <- merge(kosiol.truth,eslr.scores,by='ensg')
  print(nrow(eslr.kosiol))

  a <- simple.roc(eslr.kosiol)
  return(a)
}

compare.nielsen.eslr = function(eslr.merged) {
  nielsen.genes = read.csv("nielsen-genes.csv")
  
  # Get truth values for Nielsen tests.
  nielsen.genes$truth <- as.integer(nielsen.genes$p.value < 0.2)
  nielsen.genes[is.na(nielsen.genes$truth),]$truth = 0
  nielsen.genes$ensg <- nielsen.genes$ENS.ref
  nielsen.genes <- subset(nielsen.genes,!is.na(ensg))
  print(paste("Number of Nielsen true positives: ",(sum(nielsen.genes$truth))))

  # Get E-SLR scores and merge.
  eslr.scores <- get.eslr.score.df(eslr.merged)
  eslr.nielsen <- merge(nielsen.genes,eslr.scores,by='ensg')

  print(paste("After merging:",sum(eslr.nielsen$truth)))
  a <- simple.roc(eslr.nielsen)
  return(a)
}

compare.clark.eslr = function(eslr.merged) {
  ncbi.gene.info <- get.ncbi.gene.info()
  clark.genes <- read.delim("clark-genes.tsv",header=T,sep="\t")
  clark.genes <- subset(clark.genes,!is.na(LocusLink.ID) & LocusLink.ID != "")

  clark.ncbi <- merge(clark.genes,ncbi.gene.info,by.x="LocusLink.ID",by.y="GeneID")
  #print(clark.ncbi[1:10,c('LocusLink.ID','gene.symbol','ensg')])

  # Get truth values for Clark tests.
  #clark.ncbi$truth = as.integer(clark.ncbi$`Model.2..M2..p.value.chimp` < 0.01 | clark.ncbi$`Model.2..M2..p.value.human` < 0.01)
  clark.ncbi$truth <- as.integer(clark.ncbi$`P.value.model.1..dN.dS..human..one.sided..test.for.pos.selection.` < 0.2
    | clark.ncbi$`P.value.model.1..dN.dS..chimp..one.sided..test.for.pos.selection.` < 0.2)
  
  # Get E-SLR scores and merge.
  eslr.scores <- get.eslr.score.df(eslr.merged)
  eslr.clark <- merge(clark.ncbi,eslr.scores,by='ensg')

  print(paste("True positives after merging:",sum(eslr.clark$truth)))
  a <- simple.roc(eslr.clark)
  return(a)
}

compare.sabeti.eslr = function(eslr.merged) {
  sabeti.genes <- read.delim("sabeti-genes.txt",header=F,stringsAsFactors=F,col.names='id')
  ncbi.gene.info <- get.ncbi.gene.info()
  sabeti.ncbi <- merge(sabeti.genes,ncbi.gene.info,by.x='id',by.y='Symbol')

  eslr.scores <- get.eslr.score.df(eslr.merged)
  eslr.sabeti <- merge(sabeti.ncbi,eslr.scores,by='ensg',all.y=T)

  eslr.sabeti$truth <- as.integer(!is.na(eslr.sabeti$id))

  a <- simple.roc(eslr.sabeti)
  return(a)
}

get.ncbi.gene.info = function() {
  ncbi.gene.info <- read.delim("Homo_sapiens.gene_info",sep="\t",header=F,stringsAsFactors=F)
  names(ncbi.gene.info) <- c("tax_id","GeneID","Symbol","LocusTag","Synonyms","dbXrefs","chromosome","map.location","description","type","symbol.from.nomenclature","full.name.from.nomenclature","nomenclature.status","other.designations","modification.date")

  ncbi.ensgs <- sapply(strsplit(ncbi.gene.info$dbXrefs,"\\|"),function(x){a=grep("Ensembl",x,value=T);a=sub("Ensembl:","",a);a[1]})   
  ncbi.gene.info$ensg <- ncbi.ensgs
  ncbi.gene.info <- subset(ncbi.gene.info,!is.na(ensg) & ensg != "")
  return(ncbi.gene.info)
}

get.eslr.score.df = function(eslr.merged) {
  eslr.merged = subset(eslr.merged,!is.na(human_gene))
  eslr.scores = data.frame(
    ensg = eslr.merged$human_gene,
    score = eslr.merged$score
    )
  return(eslr.scores)
}

simple.roc = function(df,truth='truth',score='score') {
  library(doBy)
  library(plyr)

  df$score = df[[score]]
  df[is.na(df$score),]$score = min(df$score)-1
  df$truth = df[[truth]]
  df <- orderBy(~-score,data=df)

  #print(df[1:10,])
  df$tp = cumsum(df$truth) / sum(df$truth)
  df$tn = cumsum(1-df$truth) / sum(1-df$truth)

  df$p = cumsum(df$truth)
  df$n = cumsum(1-df$truth)
  
  return(df)
}

roc.summary = function(roc) {

  out.df = data.frame()
  for (i in c(0.1,0.25,0.5,0.75,1)) {
    result = roc[max(which(roc$tp <= i)),c('score','tp','tn','p','n')]
    result$thresh = i
    out.df = rbind(out.df,result)
  }
  print(out.df)
}
