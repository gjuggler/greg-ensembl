if (exists('drv')) {
  lapply(dbListConnections(drv),dbDisconnect)
} else {
  library(RMySQL)
  drv <- dbDriver('MySQL')
}
con <- dbConnect(drv, host='mysql-greg.ebi.ac.uk', port=4134, user='slrsim', password='slrsim', dbname='slrsim_anisimova')
url = "mysql://slrsim:slrsim@mysql-greg.ebi.ac.uk:4134/slrsim_anisimova"

source("aln-tools/aln.tools.R")
source("aln-tools/phylo.tools.R")
source("aln-tools/plot.phylo.greg.R")
library(ape)

plot.protein = function(node_id,label='') {
  tree_file = "tree.nh"
  aln_file = "aln.fasta"

  system(paste("perl ~/lib/greg-ensembl/scripts/tree_dump.pl",
    " --url=",url,
    " --id=",node_id,
    " --tree=",tree_file,
    " --aln=",aln_file,
    sep="")
  )

  aln = read.aln(aln_file)
  tree = read.tree(tree_file)
  aln$tree = tree
  
  length = aln$length
  tree.space = length/10
  par(mar=c(0,0,0,0))
  plot.new()
  plot.window(xlim=c(0,length + tree.space),ylim=c(0,aln$num_seqs))
  plot.aln(aln,overlay=T,plot.tree=F,plot.chars=F,
    x.lim=c(tree.space,length+tree.space),
    y.lim=c(0,aln$num_seqs))

  text(x=tree.space*.9,y=0,labels=label,cex=2,adj=c(0,0))

  rect(xleft=tree.space,xright=length+tree.space,ybottom=0,ytop=aln$num_seqs)
  par(new=T)
  len = tree_length(tree)
  print(len)
  plot.phylo.greg2(tree,x.lim=c(0,(4)/(tree.space/length)),y.lim=c(.5,aln$num_seqs+.5),edge.width=0.5)

  unlink(c(tree_file,aln_file))                                 
}

plot.by.column = function(data, column.name="label") {
  source("collect_slrsim.R")

  
}

ran.main = FALSE
main = function() {
  source("collect_slrsim.R")

  data = get.all.data()
  
  # Set up the plotting params.
  pdf("~/public_html/slrsim.pdf")
  par(mfrow=c(10,1))

#  attrs = c('slrsim_file','slrsim_tree_length','phylosim_ins_rate','alignment_name','filtering_name')
   attrs = c('slrsim_file','phylosim_insertmodel','phylosim_insertrate','filtering_name')
  ids = rep("",nrow(data))
  for (attr in attrs) {ids = paste(ids,data[[attr]],sep=" ")} 
  unique_ids = unique(ids)
  print(length(unique_ids))
  for (my_id in sort(unique_ids)) {
    print(my_id)
    df = data[ids==my_id,]
    first.row = df[1,]
    print(first.row)
    plot.protein(first.row$node_id,label=my_id)
  }

  dev.off()

  assign("ran.main",TRUE,pos=.GlobalEnv)
}