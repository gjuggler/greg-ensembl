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

  my $base = 'nfs/users/nfs_g/gj1/scratch/mammals/output';
  $self->get_output_folder($base);
  
  $self->export_genes;
#  $self->export_sites;
  $self->summarize_sites;
#  $self->model_sites_distribution;
#   $self->get_psg_overlaps;
#  $self->plot_global_distribution;
}

sub export_genes {
  my $self = shift;

  my $genes_file = $self->get_output_folder . "/genes.Rdata";

  my $folder = $self->get_output_folder;
  my $base = Bio::Greg::EslrUtils->baseDirectory;
  my $csr = "$base/scripts/collect_sitewise.R";

  my $cmd = qq^
source("$csr");
genes <- get.genes.merged(db="gj1_2x_57",exclude.cols=c('tree_newick'))
save(genes,file="${genes_file}")
^;
  print "$cmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r($cmd,$params);

}

sub export_sites {
  my $self = shift;

  my $sth = $self->db_handle->prepare("SELECT max(parameter_set_id) from parameter_set;");
  $sth->execute;
  my @row = $sth->fetchrow_array;
  $sth->finish;
  my $max_pset = 1;
  if (@row) {
    $max_pset = $row[0];
  }

  $max_pset = 1;
  foreach my $i (0 .. $max_pset) {
    my $sites_file = $self->get_output_folder . "/sites_$i.tsv";
    my $sites_psc_rdata = $self->get_output_folder . "/sites_${i}_psc.Rdata";
    my $sites_zipped = $self->get_output_folder . "/sites_$i.tsv.gz";
    my $sites_rdata = $self->get_output_folder . "/sites_$i.Rdata";
    my $csr = $self->base."/scripts/collect_sitewise.R";

    my $mysqlArgs = Bio::Greg::EslrUtils->mysqlArgsFromConnection($self->dbc);
#    if (!-e $sites_zipped) {
      my $cmd = qq^
mysql $mysqlArgs -e "SELECT s.node_id,s.aln_position,s.parameter_set_id,domain,filter_value,omega,omega_lower,omega_upper,lrt_stat,note,type,g.chr_name,g.chr_start from stats_sites s LEFT JOIN sitewise_genome g ON s.node_id=g.node_id AND s.aln_position=g.aln_position WHERE s.parameter_set_id=$i;" > ${sites_file}
^;
      $cmd = qq^
mysql $mysqlArgs -e "SELECT s.node_id,s.aln_position,s.parameter_set_id,domain,filter_value,omega,omega_lower,omega_upper,lrt_stat,note,type,g.chr_name,g.chr_start from stats_sites s LEFT JOIN sitewise_genome g ON s.node_id=g.node_id AND s.aln_position=g.aln_position WHERE s.parameter_set_id=1 LIMIT 500000;" > ${sites_file}
^ if ($i == 0);
      print "$cmd\n";
      my $rc = system($cmd);
    die "Bad return value!" if ($rc);
      $cmd = "rm -f ${sites_zipped}; gzip ${sites_file}";
      print "$cmd\n";
      system($cmd);
#    }

#    if (!-e $sites_rdata) {
      my $rcmd = qq^
source("${csr}");
gcinfo(TRUE)
sites = read.table(gzfile("${sites_zipped}"),header=T,stringsAsFactors=T,na.strings="NULL")

# Add p-value for nonneutral evolution.
#print(sites[1:10,])
p.values = 1 - pchisq(sites[,'lrt_stat'],1)
sites[,'pval'] = p.values

#print(nrow(sites))
print(str(sites))
save(sites,file="${sites_rdata}")

psc.sites = subset(sites,pval < 0.05 & omega > 1)
save(psc.sites,file="${sites_psc_rdata}")
^;
      print "$rcmd\n";
      my $params = {};
      Bio::Greg::EslrUtils->run_r($rcmd,$params);  
#    }
  }
}

sub summarize_sites {
  my $self = shift;

  my $sites_file = $self->get_output_folder . "/sites.Rdata";
  my $fdr_file = $self->get_output_folder . "/fdr_thresholds.Rdata";

  my $folder = $self->get_output_folder;

  my $slrsim_r_connection = $self->get_r_dbc_string('gj1_slrsim');
  my $indel_r_connection = $self->get_r_dbc_string('gj1_indelsim');

  my $cslr = $self->base . "/projects/slrsim/collect_slrsim.R";
  my $csw = $self->base . "/scripts/collect_sitewise.R";

  my $dbname = $self->compara_dba->dbc->dbname;
my $cmd = qq^

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
  Bio::Greg::EslrUtils->run_r($cmd,$params);

}

sub model_sites_distribution {
  my $self = shift;

  my $folder = $self->get_output_folder;
  my $csr = $self->base."/scripts/collect_sitewise.R";

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
  Bio::Greg::EslrUtils->run_r($rcmd,$params);

}

sub get_psg_overlaps {
  my $self = shift;

  my $cwd = cwd();
  chdir("psg_overlaps");

  my $folder = $self->get_output_folder;
  my $base = $self->base;
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
  Bio::Greg::EslrUtils->run_r($rcmd,$params);
  
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
  Bio::Greg::EslrUtils->run_r($rcmd,$params);

}


1;
