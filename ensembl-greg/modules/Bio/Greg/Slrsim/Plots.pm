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
  $self->get_output_folder();
}

sub run {
  my $self = shift;

  my $experiment_name = $self->param('experiment_name');

  $self->output_params_file;
  $self->dump_sql;
  $self->dump_data;
  # Tricky Perl: Call the method corresponding to the current experiment.
  $self->$experiment_name();
}

sub get_output_folder {
  my $self = shift;

  if (defined $self->param('output_folder')) {
    return $self->param('output_folder');
  }

  my $date_string = strftime("%Y-%m-%d",localtime);
  my $i=0;
  
  my $filename;
  do {
    $i++;
    $filename = sprintf("NO.backup/output/%s/%s_%.2d",$date_string,$date_string,$i);
  } while (-e $filename);

  print "Output folder: $filename\n";
  $self->param('output_folder',$filename);

  # We'll store this output folder in the meta table, so it will be re-used if this module is run again w/ the same database.
  $self->store_meta({output_folder => $filename});
  mkpath([$filename]);
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

  if (!-e $filename && !-e $gzip) {
    my $cmd = qq^mysqldump -P$port -h$h -u$u -p$p $db > $filename;^;
    system($cmd);
  }

  if (!-e $gzip) {
    my $cmd = qq^gzip $filename;^;
    system($cmd);
    unlink($filename);
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

sub mammals_indel_simulations {
  my $self = shift;

  $self->mammals_simulations;
}

sub mammals_simulations {
  my $self = shift;

  my $filename = $self->param('output_folder') . '/reference_comparison_roc.png';
  my $folder = $self->param('output_folder');
  my $script = $self->base . "/collect_slrsim.R";
  my $rcmd = qq^
dbname="gj1_slrsim"
host = 'ens-research'
port = 3306
user = 'ensadmin'
password = 'ensembl'

source("collect_slrsim.R")
data = get.all.data()

thresh = get.fdr.thresholds(data,fdr=0.05,col.names=c('parameter_set_name'))

#full.roc = generic.roc.plot(data,col.names=c('parameter_set_name'),plot=F)
#
#for (lbl in sort(unique(full.roc[,'label']))) {
#  roc = subset(full.roc,label==lbl)
#  for (fdr in c(0.01,0.05,0.1,0.2)) {
#    fdr.ok = roc[which(roc[,'fdr'] <= fdr),]
#    max.row = fdr.ok[nrow(fdr.ok),]
#    score.at.fdr = max.row[,'score']
#    tpr.at.fdr = max.row[,'tpr']
#    print(paste(lbl,'[',fdr,'=>',score.at.fdr,tpr.at.fdr,']'))
#  }
#}

q()
^;
  print "$rcmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r($rcmd,$params);
}

sub reference_comparison {
  my $self = shift;
  
  my $filename = $self->param('output_folder') . '/reference_comparison_roc.png';
  my $script = $self->script;
  my $rcmd = qq^
source("$script")
data = get.all.data()
png(file="${filename}",width=600,height=600)
generic.roc.plot(data,col.names=c('slrsim_ref','alignment_name'))
dev.off()
q()
^;
  print "$rcmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r($rcmd,$params);
  
}

sub alignment_comparison {
  my $self = shift;

  my $filename = $self->param('output_folder') . '/alignment_roc.png';
  my $script = $self->script;

  my $rcmd = qq^
source("$script")
data = get.all.data()
png(file="${filename}",width=600,height=600,pointsize=24)
generic.roc.plot(data,col.names=c('alignment_name'))
dev.off()
q()
^;
  print "$rcmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r($rcmd,$params);

}

sub filter_order_test {
  my $self = shift;
  
  $self->_plot_proteins;
}

sub filter_test {
  my $self = shift;
  $self->filter_order_test;
  $self->filter_sweep;
}

sub filter_sweeps {
  my $self = shift;
  $self->filter_sweep;
}

sub filter_sweep {
  my $self = shift;

  my $folder = $self->param('output_folder');
  my $filename = $self->param('output_folder') . '/filter_sweep_roc.png';

  my $dbname = $self->compara_dba->dbc->dbname;
  
  $self->_plot_scatter;
  $self->_plot_scores;
  $self->_fdr_sweep;
  $self->_plot_proteins;
  my $script = $self->script;

  my $rcmd = qq^
dbname = "$dbname"
source("$script")

data = get.all.data(dir="$folder")
data[is.na(data[,'alignment_score_threshold']),'alignment_score_threshold'] = 0;

label.names = c('alignment_name','filtering_name','alignment_score_threshold')

png(file="${folder}/roc_facets.png",width=600,height=600)
x = subset(data,alignment_name != "True Alignment")
facet.roc.plot(x,zoom=F,na.rm=F,plot.x='fp',plot.y='tp',facet.x='alignment_name',facet.y='filtering_name',col.names=label.names)
dev.off()

png(file="${folder}/roc_facets_zoom.png",width=600,height=600)
x = subset(data,alignment_name != "True Alignment")
facet.roc.plot(x,zoom=T,na.rm=F,plot.x='fp',plot.y='tp',facet.x='alignment_name',facet.y='filtering_name',col.names=label.names)
dev.off()

png(file="${folder}/roc_abs.png",width=600,height=600)
generic.roc.plot(data,na.rm=F,plot.x='fp',plot.y='tp',col.names=label.names)
dev.off()

png(file="${folder}/roc_fpr.png",width=600,height=600)
generic.roc.plot(data,na.rm=F,plot.x='fpr',plot.y='tpr',col.names=label.names)
dev.off()

png(file="${folder}/roc_fpr_noNA.png",width=600,height=600)
generic.roc.plot(data,na.rm=T,plot.x='fpr',plot.y='tpr',col.names=label.names)
dev.off()

q()
^;
  print "$rcmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r($rcmd,$params);
}

sub _plot_scores {
  my $self = shift;
  my $folder = $self->get_output_folder;
  my $script = $self->script;
  my $rcmd = qq^
source("$script")
data = get.all.data(dir="$folder")
col.names=c('filtering_name','alignment_score_threshold')

png(file="${folder}/ncod_plots.png",width=1800,height=600)
plot.ncod(data,col.names=col.names)
dev.off()
^;
  print "$rcmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r($rcmd,$params);
  
}

sub _plot_scatter {
  my $self = shift;
  my $folder = $self->get_output_folder;
  my $script = $self->script;
  my $rcmd = qq^
source("$script")
data = get.all.data(dir="$folder")
col.names=c('alignment_name','filtering_name','alignment_score_threshold')

png(file="${folder}/scatter_plots.png",width=1800,height=1800)
plot.scatter(data,col.names=col.names)
dev.off()
^;
  print "$rcmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r($rcmd,$params);

}

sub _plot_proteins {
  my $self = shift;
  my $folder = $self->get_output_folder;
  my $script = $self->script;
  my $rcmd = qq^
source("$script")
data = get.all.data(dir="$folder")
col.names=c('alignment_name','filtering_name','alignment_score_threshold','sitewise_filter_order','alignment_score_mask_character_cdna')
plot.by.columns(data,base.dir="${folder}",col.names=col.names)
^;
  print "$rcmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r($rcmd,$params);
}

sub _fdr_sweep {
  my $self = shift;
  my $folder = $self->get_output_folder;
  my $script = $self->script;

my $rcmd = qq^
source("$script")
png(file="${folder}/fdr_sweep.png",width=1200,height=1200)
data = get.all.data(dir="$folder")
label.names = c('alignment_name','filtering_name','alignment_score_threshold')
group.names = c('alignment_name','filtering_name')
roc = generic.roc(data,na.rm=F,col.names=label.names)
roc = subset(roc,alignment_name != "True Alignment")
#roc = subset(roc,filtering_name != "none")
#roc = subset(roc,filtering_name == "tcoffee")

roc.fdr = fdr.sweep(roc)
plot.fdr(roc.fdr,facet.x='alignment_name',facet.y='filtering_name',color.by='alignment_score_threshold');
#write.csv(roc.fdr,file="${folder}/fdr_sweep.csv",row.names=F)
dev.off()
^;
  print "$rcmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r($rcmd,$params);

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
