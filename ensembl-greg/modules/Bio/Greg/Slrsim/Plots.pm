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

  my $defaults = {};

  $self->load_all_params($defaults);
  $self->_output_folder();
}

sub run {
  my $self = shift;

  my $experiment_name = $self->param('experiment_name');

  $self->output_params_file;
  $self->dump_sql;
  # Tricky Perl: Call the method corresponding to the current experiment.
  eval {
    $self->$experiment_name();
  };
}

sub _output_folder {
  my $self = shift;
  return $self->get_output_folder('/homes/greg/lib/greg-ensembl/projects/slrsim/NO.backup/output');
}

sub dump_sql {
  my $self = shift;
  my $filename = $self->param('output_folder') . '/slrsim.sqldata';
  my $gzip = $self->param('output_folder') . '/slrsim.sqldata.gz';

  my $dbc = $self->compara_dba->dbc;
  my $u = $dbc->username;
  my $p = $dbc->password;
  my $h = $dbc->host;
  my $port = $dbc->port;
  my $db = $dbc->dbname;

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
  my $script = $self->script;

  my $dbname = $self->dbc->dbname;
 

  if (!-e $filename) {
    my $rcmd = qq^
dbname="$dbname"
source("$script")
data = get.all.data()
save(data,file="${filename}");
^;
    print "$rcmd\n";
    my $params = {};
    Bio::Greg::EslrUtils->run_r($rcmd,$params);
  }
  
}

sub slrsim_all {
  my $self = shift;

  my $folder = $self->_output_folder;
  my $all_file = $self->param('output_folder') . '/df.list.Rdata';
  my $script = $self->script;
  my $dbname = $self->dbc->dbname;
  my $rcmd = qq^
# Try to plot all the alignments in the output folder.
library(R.oo)
library(ape)
library(ggplot2)
source("~/lib/phylosim/PhyloSimSource.R")
source("~/lib/phylosim/examples/greg_scripts.R")
#plot.all.in.dir("${folder}/alns");

dbname="$dbname"
source("$script")

if(!file.exists("${all_file}")) {
  df.list = get.all.by.experiment()
  save(df.list,file="${all_file}")
} else {
  load("${all_file}")
}

for (df in df.list) {
  name = df[1,'experiment_name']
  print(paste("Experiment name:",name))
  data = df
  this.file=paste("$folder",'/',name,'.Rdata',sep='')
  if(!file.exists(this.file)) {
    save(data,file=this.file)
  }
}

^;
  print "$rcmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r($rcmd,$params);
  
}

sub fig_two {
  my $self = shift;

  $self->slrsim_all;
  my $folder = $self->_output_folder;
  my $file = "${folder}/fig_two.Rdata";
  my $roc_zoom = "${folder}/roc_zoom.pdf";
  my $functions = $self->base . "/projects/slrsim/slrsim.functions.R";
  my $plots = $self->base . "/projects/slrsim/slrsim.plots.R";

  my $rcmd = qq^
library(ggplot2)
source("${functions}")
source("${plots}")

library(R.oo)
source("~/lib/phylosim/PhyloSimSource.R")
source("~/lib/phylosim/examples/greg_scripts.R")
plot.all.in.dir("${folder}/alns");

load("${file}")
print(str(data))
na.rm <- TRUE
plot.x <- 'fp'
plot.y <- 'tp'
zoom.fpr <- 0.1

n.trees <- length(unique(data[,'slrsim_tree_file']))

for (aln in unique(data[,'alignment_name'])) {
  if (aln == 'True Alignment') {next;}
  d <- subset(data,alignment_name==aln | alignment_name=='True Alignment')

  d[,'slrsim_label'] <- paste(d[,'slrsim_tree_file'],d[,'alignment_name'],d[,'filtering_name'],d[,'alignment_score_threshold'])

  # ROC plot of alignments.
  f = function(df,thresh) {return(slr.roc(df,na.rm=na.rm))}
  comb.roc <- summarize.by.labels(d,f)
  max.x <- max(comb.roc[,plot.y]) * zoom.fpr

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

  pdf(file=paste("${folder}/",aln,"_roc.pdf",sep=""),width=10*n.trees,height=10)
  print(p)
  dev.off()

  pdf(file=paste("${folder}/",aln,"_roc_zoom.pdf",sep=""),width=10*n.trees,height=10)
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

  if (!-e $filename) {
    print "$filename\n";
    my $out;
    open($out,">$filename");
    $self->hash_print($self->params,$out);
    close($out);
  }
}

1;
