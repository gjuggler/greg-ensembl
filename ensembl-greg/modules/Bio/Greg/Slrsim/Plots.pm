package Bio::Greg::Slrsim::Plots;

use strict;
use Cwd;
use POSIX qw(strftime mktime);
use File::Path;
use Bio::Greg::EslrUtils;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::Process');

sub fetch_input {
  my $self = shift;

  $self->load_all_params();

  $self->hash_print($self->params);
  $self->_output_folder();
}

sub run {
  my $self = shift;

  my $experiment_name = $self->param('experiment_name');

  $self->output_params_file;
  $self->dump_sql;

  # Call the method corresponding to the current experiment.
  $self->$experiment_name();
}

sub _output_folder {
  my $self = shift;
  return $self->get_output_folder;
}

sub dump_sql {
  my $self     = shift;
  my $filename = $self->param('output_folder') . '/slrsim.sqldata';
  my $gzip     = $self->param('output_folder') . '/slrsim.sqldata.gz';

  my $dbc  = $self->compara_dba->dbc;
  my $u    = $dbc->username;
  my $p    = $dbc->password;
  my $h    = $dbc->host;
  my $port = $dbc->port;
  my $db   = $dbc->dbname;

  print "Mysqldumping...\n";
  if (!-e $filename && !-e $gzip) {
    #my $cmd = qq^mysqldump -P$port -h$h -u$u -p$p $db > $filename;^;
    #system($cmd);
  }

  print "Gzipping...\n";
  if (!-e $gzip) {
    #my $cmd = qq^gzip $filename;^;
    #system($cmd);
    #unlink($filename);
  }
}

sub script {
  my $self = shift;
  return $self->base . "/projects/slrsim/collect_slrsim.R";
}

sub dump_data {
  my $self = shift;

  my $filename = $self->param('output_folder') . '/slrsim_sites.Rdata';
  my $script   = $self->script;

  my $dbname = $self->dbc->dbname;

  if ( !-e $filename ) {
    my $rcmd = qq^
dbname="$dbname"
source("$script")
data = get.all.data()
save(data,file="${filename}");
^;
    print "$rcmd\n";
    my $params = {};
    Bio::Greg::EslrUtils->run_r( $rcmd, $params );
  }

}

sub slrsim_all {
  my $self = shift;

  my $folder = $self->get_output_folder;
  my $all_file = $self->param('output_folder') . '/df.list.Rdata';
  my $script = $self->script;
  my $dbname = $self->dbc->dbname;
  my $rcmd = qq^
# Try to plot all the alignments in the output folder.
library(R.oo)
library(ape)
library(ggplot2)

dbname="$dbname"
source("$script")

if(!file.exists("${all_file}")) {
  df.list = get.all.by.experiment()
  save(df.list,file="${all_file}")

  for (df in df.list) {
    name = df[1,'experiment_name']
    print(paste("Experiment name:",name))
    data = df
    this.file=paste("$folder",'/',name,'.Rdata',sep='')
    if(!file.exists(this.file)) {
      save(data,file=this.file)
    }
  }
} else {
#  load("${all_file}")
}
^;
  print "$rcmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r( $rcmd, $params );

}

sub slrsim_table {
  my $self = shift;

  my $folder = $self->get_output_folder;
  my $file = "${folder}/df.list.Rdata";
  my $functions = $self->base . "/projects/slrsim/slrsim.functions.R";
  my $table_file = "${folder}/table.csv";
  my $table_one_one = "${folder}/table_1_1.csv";

my $rcmd = qq^
source("${functions}")
load("${file}")

df.f <- function(x) {
  x[,'slrsim_label'] <- paste(x[,'slrsim_tree_file'],x[,'alignment_name'],x[,'filtering_name'])
  table.df <- ddply(x,.(slrsim_label,phylosim_insertrate,tree_mean_path),paper.table)
  return(table.df)
}
paper.df <- ldply(df.list,df.f)
write.csv(paper.df,file="${table_file}",row.names=F)

paper.df <- subset(paper.df, indel == 0.1 & length == 1)
write.csv(paper.df,file="${table_one_one}",row.names=F)

^;
  my $params = {};
  Bio::Greg::EslrUtils->run_r( $rcmd );

}

sub fig_one {
  my $self = shift;

  my $folder = $self->get_output_folder;

  # Call slrsim_all to dump the data.
  $self->slrsim_all;

  # Call slrsim_table to dump the table.
  $self->slrsim_table;

  my $mafft_nofilter_file = "${folder}/1_bglobin_mafft_None.Rdata";
  my $prank_filter_file = "${folder}/1_bglobin_prank_filter_tcoffee.Rdata";
  my $prank_nofilter_file = "${folder}/1_bglobin_prank_None.Rdata";
  my $true_aln_file = "${folder}/1_bglobin_True Alignment_None.Rdata";

  my $file = "${folder}/df.list.Rdata";
  my $f_auc = "${folder}/auc.pdf";
  my $f_fdr = "${folder}/fdr.pdf";
  my $f_sens = "${folder}/sens.pdf";
  my $f_tpr = "${folder}/tpr_at_threshold.pdf";
  my $functions = $self->base . "/projects/slrsim/slrsim.functions.R";
  my $plots = $self->base . "/projects/slrsim/slrsim.plots.R";

  my $table_file = "${folder}/table.csv";

  my $phylosim_dir = $self->base . "/projects/phylosim";

  my $rcmd = qq^
library(ggplot2)
library(fields)
library(RColorBrewer)

plot.f = function(data,palette="Spectral",field='fdr',rev.colors=F,limits=c(0,1),x.lim=c(0,4),y.lim=c(0,0.1),label.species=T,do.plot=T) {
  data[,'z_vals'] = data[,field]
  p <- ggplot(data,aes(x=tree_mean_path,y=phylosim_insertrate*2,z=z_vals))
  p <- p + scale_x_continuous('Mean Path Length')
  p <- p + scale_y_continuous('Indel Rate')
  if (!is.null(x.lim)) {
    p <- p + coord_cartesian(xlim = x.lim,ylim=y.lim)
  }
  p <- p + geom_tile(aes_string(fill=field))
  if (is.null(limits)) {
    limits <- c(min(data[,'z_vals']),max(data[,'z_vals']))
  }
  n <- 10
  breaks <- seq(from=limits[1],to=limits[2],length.out=n+1)
  if (rev.colors) {
    p <- p + scale_fill_gradientn(colours=rev(brewer.pal(n=n,palette)),limits=limits)
  } else {
    p <- p + scale_fill_gradientn(colours=brewer.pal(n=n,palette),limits=limits)
  }
  p <- p + coord_equal(ratio=10)
  if(do.plot) {
    print(p)
  } else {
    return(p)
  }
}

source("${functions}")
source("${plots}")

load("${file}")

na.rm <- TRUE
plot.x <- 'fpr'
plot.y <- 'tpr'
zoom.fpr <- 0.1

df.f <- function(x) {
  slr.threshold <- 0
  a <- summarize.by.labels(x,fn=fig.1.summary,thresh=slr.threshold)
  return(a)
}
summary.df <- ldply(df.list,df.f)
sub.df <- summary.df

art = 'artificial.nh'
bg = 'bglobin.nh'
enc = 'encode.nh'
tree.files <- unique(sub.df[,'slrsim_tree_file'])
tree.labels <- c('6 artificial','17 bglobin', '44 vertebrate')
sub.df[sub.df[,'slrsim_tree_file'] == art,'slrsim_tree_file'] = 0
sub.df[sub.df[,'slrsim_tree_file'] == bg,'slrsim_tree_file'] = 1
sub.df[sub.df[,'slrsim_tree_file'] == enc,'slrsim_tree_file'] = 2
#sub.df[,'slrsim_tree_file'] <- factor(sub.df[,'slrsim_tree_file'],labels=tree.labels)

str(sub.df)
pdf(file="${f_auc}",width=10,height=5)
p <- plot.f(sub.df,field='auc',do.plot=F,label.species=F,rev.colors=F,limits=c(0.4,1))
p <- p + facet_grid(slrsim_tree_file ~ alignment_name)
print(p)
dev.off()

pdf(file="${f_fdr}",width=10,height=5)
p <- plot.f(sub.df,field='fdr',do.plot=F,label.species=F,rev.colors=T)
p <- p + facet_grid(slrsim_tree_file ~ alignment_name)
print(p)
dev.off()

pdf(file="${f_sens}",width=10,height=5)
p <- plot.f(sub.df,field='sens',do.plot=F,label.species=F,rev.colors=F)
p <- p + facet_grid(slrsim_tree_file ~ alignment_name)
print(p)
dev.off()

pdf(file="${f_tpr}",width=10,height=5)
p <- plot.f(sub.df,field='tpr_at_threshold',do.plot=F,label.species=F,rev.colors=F)
p <- p + facet_grid(slrsim_tree_file ~ alignment_name)
print(p)
dev.off()

^;

my $rcmd2 = qq^
source("${functions}")
source("${plots}")

na.rm <- TRUE
plot.x <- 'fpr'
plot.y <- 'tpr'
zoom.fpr <- 0.1


load("${true_aln_file}");
df <- data
load("${prank_nofilter_file}");
df <- rbind(df,data);
load("${prank_filter_file}");
df <- rbind(df,data);

d <- subset(df)
#d <- subset(df, phylosim_insertrate==0.05 & tree_mean_path==1.0)

d[,'slrsim_label'] <- paste(d[,'slrsim_tree_file'],d[,'alignment_name'],d[,'filtering_name'],d[,'alignment_score_threshold'])

# ROC plot of alignments.
f = function(df,thresh) {return(slr.roc(df,na.rm=na.rm))}
comb.roc <- summarize.by.labels(d,f)
max.x <- max(comb.roc[,plot.y]) * zoom.fpr

d[,'slrsim_label'] <- paste(d[,'slrsim_tree_file'],d[,'alignment_name'],d[,'filtering_name'],d[,'alignment_score_threshold'])

sub.roc <- comb.roc
# Remove the tree file name from the labels so the plots only group into align-filter-threshold
sub.roc[,'slrsim_label'] <- paste(sub.roc[,'alignment_name'],sub.roc[,'filtering_name'],sub.roc[,'alignment_score_threshold'])

p <- plot.roc(sub.roc,plot=F,plot.x=plot.x,plot.y=plot.y)
#  p <- p + scale_colour_brewer("Filtering Scheme",palette="Set1")

fdr.line <- 0.2
p <- p + geom_abline(slope=2/fdr.line,colour='gray')
fdr.line <- 0.1
p <- p + geom_abline(slope=2/fdr.line,colour='black')

p <- p + facet_grid(. ~ slrsim_tree_file )

n.trees <- length(unique(d[,'slrsim_tree_file']))
pdf(file=paste("${folder}/filter_roc.pdf",sep=""),width=10*n.trees,height=10)
print(p)
dev.off()

pdf(file=paste("${folder}/filter_roc_zoom.pdf",sep=""),width=10*n.trees,height=10)
p <- p + scale_x_continuous(limits=c(0,max.x))
print(p)
dev.off()

^;

#  Bio::Greg::EslrUtils->run_r($rcmd2,{});
  Bio::Greg::EslrUtils->run_r($rcmd,{});
  
}

sub fig_two {
  my $self = shift;

  # Call slrsim_all to dump the data.
  $self->slrsim_all;

  # Call slrsim_table to dump the table.
  $self->slrsim_table;

  my $folder = $self->get_output_folder;
  my $file = "${folder}/fig_two.Rdata";
  my $roc_zoom = "${folder}/roc_zoom.pdf";
  my $functions = $self->base . "/projects/slrsim/slrsim.functions.R";
  my $plots = $self->base . "/projects/slrsim/slrsim.plots.R";

  my $phylosim_dir = $self->base . "/projects/phylosim";

  my $rcmd = qq^
library(ggplot2)
source("${functions}")
source("${plots}")

#library(R.oo)
#source("${phylosim_dir}/PhyloSimSource.R")
#source("${phylosim_dir}/examples/greg_scripts.R")
#plot.all.in.dir("${folder}/alns");

load("${file}")
print(str(data))
na.rm <- TRUE
plot.x <- 'fpr'
plot.y <- 'tpr'
zoom.fpr <- 0.1

d <- data
labels <- paste(
  data[,'slrsim_tree_file'],
  data[,'alignment_name'],
  data[,'phylosim_insertrate'],
  data[,'tree_mean_path'],
  sep='_'
)
data[,'slrsim_label'] <- factor(labels)

for (lbl in unique(data[,'slrsim_label'])) {
  print(lbl)
  d <- subset(data,slrsim_label == lbl)

  # ROC plot of alignments.
  f = function(df,thresh) {return(slr.roc(df,na.rm=na.rm))}
  d[,'slrsim_label'] <- d[,'filtering_name']
  comb.roc <- summarize.by.labels(d,f)
  max.x <- max(comb.roc[,plot.x]) * zoom.fpr
  sub.roc <- comb.roc

  p <- plot.roc(sub.roc,plot=F,plot.x=plot.x,plot.y=plot.y)

  fdr.line <- 0.2
  p <- p + geom_abline(slope=2/fdr.line,colour='gray')
  fdr.line <- 0.1
  p <- p + geom_abline(slope=2/fdr.line,colour='black')

  pdf(file=paste("${folder}/",lbl,"_roc.pdf",sep=""),width=10,height=10)
  print(p)
  dev.off()

  pdf(file=paste("${folder}/",lbl,"_roc_zoom.pdf",sep=""),width=10,height=10)
  p <- p + scale_x_continuous(limits=c(0,max.x))
  print(p)
  dev.off()
}
^;
  Bio::Greg::EslrUtils->run_r($rcmd,{});


}


sub output_params_file {
  my $self = shift;

  my $filename = $self->param('output_folder') . '/params.txt';

  if ( !-e $filename ) {
    print "$filename\n";
    my $out;
    open( $out, ">$filename" );
    $self->hash_print( $self->params, $out );
    close($out);
  }
}

1;
