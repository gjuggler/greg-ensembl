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

  #$self->hash_print($self->params);
  $self->_output_folder();
}

sub run {
  my $self = shift;

  $self->output_params_file;
  $self->dump_sql;
  $self->dump_data;
#  $self->dump_alns;

  # Call the method corresponding to the current experiment.
  my $experiment_name = $self->param('experiment_name');
  $self->$experiment_name();
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

  if (!-e $filename && !-e $gzip) {
    print "  mysqldumping...\n";
    my $cmd = qq^mysqldump -P$port -h$h -u$u -p$p $db > $filename;^;
    system($cmd);
  } else {
    print "  mysqldump already exists\n";
  }

  if (!-e $gzip) {
    print "  gzipping...\n";
    my $cmd = qq^gzip $filename;^;
    system($cmd);
    unlink($filename);
  } else {
    print "  gzip already exists\n";
  }
}

sub dump_alns {
  my $self = shift;

  my $folder = $self->get_output_folder;

  my $alns_gzip = "$folder/alns.tar.gz";
  if (!-e $alns_gzip) {
    my $cwd = cwd();
    chdir($folder);
    system('tar -cvvf alns.tar data; gzip alns.tar');
    chdir($cwd);
  }
}

sub dump_data {
  my $self = shift;

  my $sites_f = $self->param('output_folder') . '/sites.Rdata';
  my $genes_f = $self->param('output_folder') . '/genes.Rdata';
  my $collect_script   = $self->collect_script;
  my $functions_script   = $self->functions_script;

  my $dbname = $self->dbc->dbname;

  if ( !-e $sites_f || !-e $genes_f) {
    my $rcmd = qq^
dbname="$dbname"
source("${collect_script}")
source("${functions_script}")
sites <- get.all.data()
genes <- get.all.genes()
save(sites, genes, file="${sites_f}");
save(genes, file="${genes_f}");
^;
    #print "$rcmd\n";
    my $params = {};
    Bio::Greg::EslrUtils->run_r( $rcmd, $params );
  }

}

# Slrsim table: generic summary table output.
sub slrsim_table {
  my $self = shift;

  my $folder = $self->get_output_folder;
  my $file = "${folder}/sites.Rdata";
  my $functions = $self->base . "/projects/slrsim/slrsim.functions.R";
  my $table_file = "${folder}/table.csv";

my $rcmd = qq^
### Plots::slrsim_table

source("${functions}")
load("${file}")

sites[, 'slrsim_label'] <- paste( sites[, 'slrsim_label'], sites[, 'analysis'], sep='_')

paper.df <- ddply(sites, .(slrsim_label), paper.table)
paper.df <- format.numeric.df(paper.df, digits=3)
write.csv(paper.df, file="${table_file}", row.names=F, quote=F)
^;
  Bio::Greg::EslrUtils->run_r( $rcmd );

}

sub fig_zero {
  my $self = shift;

  my $folder = $self->get_output_folder;

  # Call slrsim_all to dump the data.
  my $sites_file = "${folder}/sites.Rdata";
  if (!-e $sites_file) {
    $self->slrsim_all;
  }

  # Call slrsim_table to dump the table.
  my $table_file = "${folder}/table.csv";
  if (!-e $table_file) {
    $self->slrsim_table;
  }  

}


sub fig_one_a {
  my $self = shift;

  $self->fig_one;
}

sub fig_one_b {
  my $self = shift;

  $self->fig_one;
}

sub fig_one_c  {
  my $self = shift;

  $self->fig_one;
}

sub fig_one {
  my $self = shift;

  my $folder = $self->get_output_folder;

  # Call slrsim_all to dump the data.
  my $sites_file = "${folder}/sites.Rdata";
  if (!-e $sites_file) {
    $self->slrsim_all;
  }

  # Call slrsim_table to dump the table.
  my $table_file = "${folder}/table.csv";
  if (!-e $table_file) {
    $self->slrsim_table;
  }

  my $f_auc = "${folder}/auc.pdf";
  my $f_auc_full = "${folder}/auc_full.pdf";
  my $f_tpr_at_fdr = "${folder}/tpr_at_fdr.pdf";
  my $f_tpr_at_fpr = "${folder}/tpr_at_fpr.pdf";
  my $f_cor = "${folder}/cor.pdf";

  my $sites_at_default = "${folder}/sites_at_default.Rdata";
  my $rocs_at_default = "${folder}/rocs_at_default.pdf";

  my $indel_sweep_cor = "${folder}/indel_sweep_cor.pdf";
  my $indel_sweep_tpr_fpr = "${folder}/indel_sweep_tpr_fpr.pdf";
  my $indel_sweep_auc = "${folder}/indel_sweep_auc.pdf";

  my $functions = $self->base . "/projects/slrsim/slrsim.functions.R";
  my $plots = $self->base . "/projects/slrsim/slrsim.plots.R";

  my $phylosim_dir = $self->base . "/projects/phylosim";

  my $rcmd = qq^
library(ggplot2)
library(fields)
library(RColorBrewer)

plot.f = function(data,palette="Spectral",field='fdr',rev.colors=F,limits=c(0,1),x.lim=c(0,4),y.lim=c(0,0.1),do.plot=T) {
  data[,'z_vals'] = data[,field]
  p <- ggplot(data,aes(x=length,y=ins_rate,z=z_vals))
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

  length.range = diff(range(data[, 'length']))
  indel.range = diff(range(data[, 'ins_rate']))
  if (indel.range > 0) {
    p <- p + coord_equal(ratio = length.range / indel.range)
  }
  if(do.plot) {
    print(p)
  } else {
    return(p)
  }
}

source("${functions}")
source("${plots}")

summary.df <- read.csv("${table_file}")
tbl <- summary.df

na.rm <- TRUE
plot.x <- 'fpr'
plot.y <- 'tpr'
zoom.fpr <- 0.1

art = 'artificial.nh'
bg = 'bglobin.nh'
enc = 'encode.nh'
tree.files <- unique(tbl[,'tree'])
tree.labels <- c('6 artificial','17 bglobin', '44 vertebrate')
tbl[tbl[,'tree'] == art,'slrsim_tree_file'] = 0
tbl[tbl[,'tree'] == bg,'slrsim_tree_file'] = 1
tbl[tbl[,'tree'] == enc,'slrsim_tree_file'] = 2
#tbl[,'tree'] <- factor(tbl[,'tree'],labels=tree.labels)

w <- 10
h <- 5

# Default range ROC plot
do.default.roc <- FALSE
if (do.default.roc) {
  if (!file.exists("${sites_at_default}")) {
    print("Loading sites subset at default values...")
    load("${sites_file}")
    sub.sites <- subset(sites, tree_length == 1 & ins_rate == 0.05)
    save(sub.sites, file="${sites_at_default}")
  }
  load("${sites_at_default}")
  sub.sites[, 'slrsim_label'] <- paste(
    sub.sites[, 'analysis'], 
    sub.sites[, 'aligner'],
    sep='|'
  )
  f = function(df,thresh) {return(slr.roc(df, na.rm=T))}
  roc <- summarize.by.labels(sub.sites, f)
  p <- plot.roc(roc, plot=F, plot.x='fpr', plot.y='tpr', color.by='aligner')
  p <- p + facet_grid(. ~ analysis)
  full.roc <- p
  max.x <- max(roc[, 'fpr']) * 0.05
  sub.roc <- p + scale_x_continuous(limits=c(0, max.x))
  pdf(file="${rocs_at_default}", width=15, height=10)
  subplot <- function(x, y) viewport(layout.pos.col=x, layout.pos.row=y)
  vplayout <- function(x, y) {
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(y,x)))
  }
  vplayout(1, 2)
  print(sub.roc, vp=subplot(1, 1))
  print(full.roc, vp=subplot(1, 2))
  dev.off()
}

# Single-variable param sweep.
p.sweep <- function(file, field, label) {
  pdf(file=file, width=10, height=8)
  default.rate <- subset(tbl, ins_rate %in% c(0, 0.02, 0.05, 0.08, 0.1))
  default.rate[, 'ins_rate'] <- paste("Ins. Rate = ",default.rate[, 'ins_rate'])
  default.rate[, 'z_value'] <- default.rate[, field]
  p <- ggplot(default.rate, aes(x=length, y=z_value))
  p  <- p + geom_line(aes(colour=analysis), alpha=0.8)
  p <- p + scale_colour_brewer(name="Analysis Method", palette="Set1")
  p <- p + scale_x_continuous(name="Tree Length", limits=c(0.2, 2), breaks=c(0.2, 1.0, 2.0))
  y.lim <- c(0, 1)
  if (field == 'auc') {
    y.lim <- c(0.5, 1)
  }
  p <- p + scale_y_continuous(name=label, limits=y.lim, breaks=c(0, 0.5, 1.0))
  p <- p + facet_grid(ins_rate ~ aligner)
  p <- p + opts(
    panel.margin=unit(rep(0.5, times=4), "lines"),
    title=paste(label," vs Tree Length", sep='')
  )
  p <- p + theme_bw()
  print(p)
  dev.off()
}

p.sweep(file="${indel_sweep_cor}", "cor", "Correlation")
p.sweep(file="${indel_sweep_tpr_fpr}", "tpr_at_fpr", "TPR at FPR")
p.sweep(file="${indel_sweep_auc}", "auc", "AUC_0.2")

# Parameter sweep plots for summary statistics.
pdf(file="${f_auc}",width=w,height=h)
p <- plot.f(tbl,field='auc',do.plot=F,rev.colors=F,limits=c(0.4,1))
p <- p + facet_grid(analysis ~ aligner)
print(p)
dev.off()

pdf(file="${f_auc_full}",width=w,height=h)
p <- plot.f(tbl,field='auc_full',do.plot=F,rev.colors=F,limits=c(0.4,1))
p <- p + facet_grid(analysis ~ aligner)
print(p)
dev.off()

pdf(file="${f_tpr_at_fpr}",width=w,height=h)
p <- plot.f(tbl,field='tpr_at_fpr',do.plot=F,rev.colors=F)
p <- p + facet_grid(analysis ~ aligner)
print(p)
dev.off()

pdf(file="${f_tpr_at_fdr}",width=w,height=h)
p <- plot.f(tbl,field='tpr_at_fdr',do.plot=F,rev.colors=F)
p <- p + facet_grid(analysis ~ aligner)
print(p)
dev.off()

pdf(file="${f_cor}",width=w,height=h)
p <- plot.f(tbl,field='cor',do.plot=F,rev.colors=F)
p <- p + facet_grid(analysis ~ aligner)
print(p)
dev.off()

^;
  Bio::Greg::EslrUtils->run_r($rcmd,{});
  
}

sub fig_two_a {
  my $self = shift;
  $self->fig_two;
}

sub fig_two_b {
  my $self = shift;
  $self->fig_two;
}

sub fig_two_c {
  my $self = shift;
  $self->fig_two;
}

sub fig_two_d {
  my $self = shift;
  $self->fig_two;
}

sub fig_two_e {
  my $self = shift;
  $self->fig_two;
}

sub fig_two_f {
  my $self = shift;
  $self->fig_two;
}

sub fig_two_h {
  my $self = shift;
  $self->fig_two;
}

sub fig_two_clustalw {
  my $self = shift;
  $self->fig_two;
}

sub fig_two_fsa {
  my $self = shift;
  $self->fig_two;
}

sub fig_two {
  my $self = shift;

  my $folder = $self->get_output_folder;

  # Call slrsim_all to dump the data.
  my $sites_file = "${folder}/sites.Rdata";
  if (!-e $sites_file) {
    $self->slrsim_all;
  }

  # Call slrsim_table to dump the table.
  my $table_file = "${folder}/table.csv";
  if (!-e $table_file) {
    $self->slrsim_table;
  }
}

sub ari_indels {
  my $self = shift;

  my $folder = $self->get_output_folder;

  # Call slrsim_all to dump the data.
  my $sites_file = "${folder}/sites.Rdata";
  if (!-e $sites_file) {
    $self->slrsim_all;
  }

  # Call slrsim_table to dump the table.
  my $table_file = "${folder}/table.csv";
  if (!-e $table_file) {
    $self->slrsim_table;
  }
  
}

# Filter table: return a table of results with multiple filtering results on each line.
sub filter_table {
  my $self = shift;

  my $folder = $self->get_output_folder;
  my $file = "${folder}/sites.Rdata";
  my $functions = $self->base . "/projects/slrsim/slrsim.functions.R";
  my $fdr_tpr = "${folder}/fdr_tpr_collated.csv";
  my $fpr_tpr = "${folder}/tpr_tpr_collated.csv";
  my $cor = "${folder}/cor_collated.csv";
  my $auc = "${folder}/auc_collated.csv";

my $rcmd = qq^
source("${functions}")

### Plots::filter_table

load("${file}")

x <- sites
x[, 'slrsim_label'] <- paste(x[, 'slrsim_tree_file'], x[,'alignment_name'], x[,'filtering_name'])
table.df <- ddply(x, .(slrsim_label, phylosim_insertrate, tree_mean_path), paper.table)
paper.df <- table.df

for (col in c('auc', 'fpr_tpr', 'fdr_tpr', 'cor')) {
  filter.df <- ddply(paper.df,
    .(phylosim_insertrate,tree_mean_path,tree,aligner),
    function(x) {
      filter.list <- dlply(x,.(filter),function(y) {y[1,col]})
      return(data.frame(x[1,],filter.list))
    })
  write.csv(filter.df,file=paste("${folder}/",col,"_collated.csv",sep=""),row.names=F,quote=F)
}

^;
  Bio::Greg::EslrUtils->run_r( $rcmd );  
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

sub _output_folder {
  my $self = shift;
  return $self->get_output_folder;
}

sub collect_script {
  my $self = shift;
  return $self->base . "/projects/slrsim/collect_slrsim.R";
}

sub functions_script {
  my $self = shift;
  return $self->base . "/projects/slrsim/slrsim.functions.R";
}

1;
