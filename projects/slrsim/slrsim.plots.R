source("aln-tools/aln.tools.R")
source("aln-tools/phylo.tools.R")
source("aln-tools/plot.phylo.greg.R")
library(ape)

plot.aln.detail = function(base) {
  t_f = paste(base,".nh",sep="")
  a_f = paste(base,".fa",sep="")
  csv_f = paste(base,".csv",sep="")

  tree = read.tree(t_f)
  aln = read.aln(a_f,seqType='protein')
  csv = read.csv(csv_f)

  sort.aln.by.tree(aln,tree)

  num_rows = 3

  par(mfrow=c(num_rows+1,1),
  mar=c(5,4,0,0))
  xl=c(0,3000)
  xs = csv$aln_position

  plot(y=csv$ncod,x=xs,xlim=xl,type="s")
  plot(y=csv$true_entropy,x=xs,xlim=xl,type="s")

  plot.new()
  plot.window(xlim=xl,ylim=c(0,aln$num_seqs))
  plot.aln(aln,overlay=T)

  return(csv)
}

plot.aln.overview = function(base,len=2000) {
  t_f = paste(base,".nh",sep="")
  a_f = paste(base,".fa",sep="")
  csv_f = paste(base,".csv",sep="")

  tree = read.tree(t_f)
  aln = read.aln(a_f,seqType='protein')
  csv = read.csv(csv_f)
  aln = sort.aln.by.tree(aln,tree)

  row = csv[1,]

  xl = c(-500,len)
  plot.new()
  plot.window(xlim=xl,ylim=c(0,aln$num_seqs))
  plot.aln(aln,overlay=T)

  tx = xl[2]
  ty = aln$num_seqs/2

  txt = paste(
#   "model: ",row$phylosim_insertmodel,"\t",
#   "rate: ",row$phylosim_insertrate,"\t",
#   "root_len: ",row$phylosim_seq_length,"\t",
  "file: ",row$slrsim_file,"\t",
  "aln: ",row$alignment_name,"\t",
  "node: ",row$node_id,

"\n",

#  "size: ",row$leaf_count,"\t",
#  "len: ",row$tree_length,"\t",
  "filter: ",row$filtering_name,"\t",
  "thr: ",row$alignment_score_threshold,"\t",
  "%:",sprintf("%.2f",row$unfiltered_site_fraction),"\t",
  "tcs:",sprintf("%.2f",row$total_column_score),"\t",
  "sps:",sprintf("%.2f",row$sum_of_pairs_score),"\t",
   sep="")

  text(tx,ty,txt,adj=c(1,0.5))

  par(new=T)
  len = tree_length(tree)
  plot.phylo.greg2(tree,x.lim=c(-1,30),y.lim=c(.5,aln$num_seqs+.5),edge.width=0.5)
}

plot.all.tiled = function(dir,len=2000,width=NA,height=NA) {
  files = list.files(path=dir,pattern="fa")
  files = as.numeric(gsub("\\.fa","",files))
  files = sort(files)

  if (is.na(width)) {
    width = 1600
  }
  if (is.na(height)) {
    height = 64 * length(files)
  }

  png(filename=paste(dir,"overview.png",sep=""),width=width,height=height,pointsize=14)
  print(c(dir,width,height))

  n = length(files)

  par(mar=c(0,0,0,0))
  par(mfrow=c(n+1,1))
  for (i in 1:length(files)) {
    f = files[i]
    print(f)
    plot.aln.overview(paste(dir,f,sep=""),len=len)
  }

  dev.off()    
}

plot.all.details = function(dir='.') {
  files = list.files(path=dir,pattern="fa")
  for (i in 1:length(files)) {
    f = files[i]
    f = gsub("\\.fa","",f)
    print(f)

    png(filename=paste(dir,f,".png",sep=""),width=1200,height=2400)
    a = plot.file(paste(dir,f,sep=""))
    dev.off()
  }
}

# Example usage
#base = "NO.backup/2010-03-17/indel_models/"
#plot.all.tiled(base)
#plot.all.details(base)
