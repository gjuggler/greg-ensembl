package Bio::Greg::Mammals::OutputMammalsData;

use strict;
use Time::HiRes qw(sleep);

use POSIX qw(strftime mktime);
use Cwd;
use File::Path;

use Bio::EnsEMBL::Hive::Process;

use Bio::Greg::EslrUtils;

use base ('Bio::Greg::Hive::Process');

sub run {
  my $self = shift;

  my $dir = $self->get_output_folder;

  my $c_genes_file = "${dir}/genes_m_c.Rdata";
  my $g_genes_file = "${dir}/genes_m_g.Rdata";

  my $c_sites_file = "${dir}/sites_m_c.Rdata";
  my $g_sites_file = "${dir}/sites_m_g.Rdata";

  if (!-e $c_genes_file) {
    $self->export_genes('genes_c', $c_genes_file);
    $self->export_genes('genes_g', $g_genes_file);
  }

  if (!-e $c_sites_file) {
    $self->export_sites('sites_c', $c_sites_file);
    $self->export_sites('sites_g', $g_sites_file);
  }

  #$self->summarize_sites;

  #  $self->model_sites_distribution;
  #   $self->get_psg_overlaps;
  #  $self->plot_global_distribution;
}

sub export_genes {
  my $self = shift;
  my $table = shift;
  my $output_file = shift;

  my $folder = $self->get_output_folder;
  my $base   = Bio::Greg::EslrUtils->baseDirectory;
  my $csr    = "$base/scripts/collect_sitewise.R";
  my $dbname = $self->hive_dba->dbc->dbname;

  my $cmd = qq^
dbname <- "${dbname}"
source("$csr");

genes <- get.vector(con,"SELECT * from ${table};")
save(genes,file="${output_file}")
^;
  print "$cmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r( $cmd, $params );
}

sub export_sites {
  my $self = shift;
  my $table = shift;
  my $output_file = shift;

  my $folder = $self->get_output_folder;
  my $base   = Bio::Greg::EslrUtils->baseDirectory;
  my $csr    = "$base/scripts/collect_sitewise.R";
  my $dbname = $self->hive_dba->dbc->dbname;

  my $cmd = qq^
dbname <- "${dbname}"
source("$csr");

sites <- get.vector(con,"SELECT * from ${table};")
sites <- factorize(sites, columns=c('chr_name', 'exon_position', 'type', 'pfam_domain'))

save(sites,file="${output_file}")
^;
  print "$cmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r( $cmd, $params );
}


sub summarize_sites {
  my $self = shift;

  my $sites_file = $self->get_output_folder . "/sites.Rdata";
  my $fdr_file   = $self->get_output_folder . "/fdr_thresholds.Rdata";

  my $folder = $self->get_output_folder;

  my $slrsim_r_connection = $self->get_r_dbc_string('gj1_slrsim');
  my $indel_r_connection  = $self->get_r_dbc_string('gj1_indelsim');

  my $cslr = $self->base . "/projects/slrsim/collect_slrsim.R";
  my $csw  = $self->base . "/scripts/collect_sitewise.R";

  my $dbname = $self->compara_dba->dbc->dbname;
  my $cmd    = qq^

if (!file.exists("$fdr_file")) {
  $slrsim_r_connection
  source("${cslr}")
  data = get.all.data()
print(data[1,])
  fdr.thresholds = get.fdr.thresholds(data,col.names=c('parameter_set_name','experiment_name'))

  print(fdr.thresholds[1:2,])
q()

  $indel_r_connection
  source("${cslr}")
  data = get.all.data()
  fdr.indel.thresholds = get.fdr.thresholds(data,col.names=c('parameter_set_name'))
  save(fdr.thresholds,fdr.indel.thresholds,file="$fdr_file")
} else {
  load("$fdr_file")
}

dbname="$dbname"
source("${csw}");
#gcinfo(TRUE)

# Summarize sites.
#summary = sites.summary(sites.dir="${folder}",filter.fn=NULL)
#write.csv(summary,file="${folder}/summary_table.csv",row.names=F)

# Summarize the filtered version.
filter.fn = function(df) {
  return(subset(df,filter_value >= 3))
}
summary = sites.summary(sites.dir="${folder}",filter.fn=filter.fn)
write.csv(summary,file="${folder}/summary_table_filt.csv",row.names=F)

# Summarize the filtered version.
filter.fn = function(df) {
  return(subset(df,!is.na(domain)))
}
summary = sites.summary(sites.dir="${folder}",filter.fn=filter.fn)
write.csv(summary,file="${folder}/summary_table_domains.csv",row.names=F)

^;
  print "$cmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r( $cmd, $params );

}

sub model_sites_distribution {
  my $self = shift;

  my $folder = $self->get_output_folder;
  my $csr    = $self->base . "/scripts/collect_sitewise.R";

  my $rcmd = qq^
source("${csr}")
library(fitdistrplus)

fits = sites.fit.intervals(db="gj1_2x_57",sites.dir="${folder}")
write.csv(fits,file="${folder}/distribution_fits.csv",row.names=F)

filter.fn = function(df) {
  return(subset(df,omega > 0))
}
fits = sites.fit.intervals(db="gj1_2x_57",sites.dir="${folder}",filter.fn=filter.fn)
write.csv(fits,file="${folder}/distribution_fits_above_0.csv",row.names=F)

filter.fn = function(df) {
  return(subset(df,omega < 1 & filter_value >= 3))
}
fits = sites.fit.intervals(db="gj1_2x_57",sites.dir="${folder}",filter.fn=filter.fn)
write.csv(fits,file="${folder}/distribution_fits_filt_below_one.csv",row.names=F)

#filter.fn = function(df) {
#  return(subset(df,!is.na(domain) & filter_value >= 3))
#}
#fits = sites.fit.intervals(db="gj1_2x_57",sites.dir="${folder}",filter.fn=filter.fn)
#write.csv(fits,file="${folder}/distribution_fits_filt_domains.csv",row.names=F)

^;
  print "$rcmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r( $rcmd, $params );

}

sub get_psg_overlaps {
  my $self = shift;

  my $cwd = cwd();
  chdir("psg_overlaps");

  my $folder         = $self->get_output_folder;
  my $base           = $self->base;
  my $psg_overlaps_r = "$base/projects/2xmammals/psg_overlaps/psg_overlaps.R";

  my $rcmd = qq^
source("$psg_overlaps_r");
load("${folder}/genes.Rdata")

get.all.scores(genes,score.column='p_positive_1')
save(roc.kosiol,roc.nielsen,roc.clark,roc.sabeti,roc.summary,file="${folder}/psg_rocs_p_pos1.Rdata")
get.all.scores(genes,score.column='m_positive_1')
save(roc.kosiol,roc.nielsen,roc.clark,roc.sabeti,roc.summary,file="${folder}/psg_rocs_m_pos1.Rdata")

get.all.scores(genes,score.column='p_slr_dnds')
save(roc.kosiol,roc.nielsen,roc.clark,roc.sabeti,roc.summary,file="${folder}/psg_rocs_p_dnds.Rdata")
get.all.scores(genes,score.column='m_slr_dnds')
save(roc.kosiol,roc.nielsen,roc.clark,roc.sabeti,roc.summary,file="${folder}/psg_rocs_m_dnds.Rdata")

genes[,'m_nlog_pval'] = -log(genes[,'m_pval_stouffer'])
genes[,'p_nlog_pval'] = -log(genes[,'p_pval_stouffer'])
get.all.scores(genes,score.column='p_nlog_pval')
save(roc.kosiol,roc.nielsen,roc.clark,roc.sabeti,roc.summary,file="${folder}/psg_rocs_p_pval.Rdata")
get.all.scores(genes,score.column='m_nlog_pval')
save(roc.kosiol,roc.nielsen,roc.clark,roc.sabeti,roc.summary,file="${folder}/psg_rocs_m_pval.Rdata")
^;
  print "$rcmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r( $rcmd, $params );

  chdir($cwd);
}

sub plot_global_distribution {
  my $self = shift;

  my $folder = $self->get_output_folder;

  my $rcmd = qq^
source("../../scripts/collect_sitewise.R");

pdf("${folder}/global_distr.pdf",width=24,height=8,pointsize=24)

#filt.fn = function(data) {return(data[sample(nrow(data),size=100000),])}
sites.plot.distributions(db="gj1_2x_57",sites.dir="${folder}")

dev.off()

^;
  print "$rcmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r( $rcmd, $params );

}

1;
