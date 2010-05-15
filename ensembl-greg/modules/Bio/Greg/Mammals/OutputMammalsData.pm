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
  
#  $self->export_genes;
  $self->export_sites;
  $self->summarize_sites;
#  $self->model_sites_distribution;
#   $self->get_psg_overlaps;
#  $self->plot_global_distribution;
}

sub export_genes {
  my $self = shift;

  my $genes_file = $self->get_output_folder . "/genes.Rdata";

  my $folder = $self->get_output_folder;

  my $cmd = qq^
source("../../scripts/collect_sitewise.R");
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

  foreach my $i (0 .. $max_pset) {
    my $sites_file = $self->get_output_folder . "/sites_$i.tsv";
    my $sites_zipped = $self->get_output_folder . "/sites_$i.tsv.gz";
    my $sites_rdata = $self->get_output_folder . "/sites_$i.Rdata";

    my $mysqlArgs = Bio::Greg::EslrUtils->mysqlArgsFromConnection($self->dbc);
    if (!-e $sites_zipped) {
      my $cmd = qq^
mysql $mysqlArgs -e "SELECT node_id,parameter_set_id,domain,filter_value,omega,omega_lower,omega_upper,lrt_stat,note,type from stats_sites WHERE parameter_set_id=$i;" > ${sites_file}
^;
      $cmd = qq^
mysql $mysqlArgs -e "SELECT node_id,parameter_set_id,domain,filter_value,omega,omega_lower,omega_upper,lrt_stat,note,type from stats_sites WHERE parameter_set_id=1 LIMIT 500000;" > ${sites_file}
^ if ($i == 0);
      print "$cmd\n";
      system($cmd);
      $cmd = "rm -f ${sites_zipped}; gzip ${sites_file}";
      print "$cmd\n";
      system($cmd);
    }

    if (!-e $sites_rdata) {
      my $rcmd = qq^
source("../../scripts/collect_sitewise.R");
gcinfo(TRUE)
sites = read.table(gzfile("${sites_zipped}"),header=T,stringsAsFactors=T,na.strings="NULL")
print(nrow(sites))
save(sites,file="${sites_rdata}")
^;
      print "$rcmd\n";
      my $params = {};
      Bio::Greg::EslrUtils->run_r($rcmd,$params);  
    }
  }
}

sub summarize_sites {
  my $self = shift;

  my $sites_file = $self->get_output_folder . "/sites.Rdata";
  my $folder = $self->get_output_folder;
my $cmd = qq^
source("../../scripts/collect_sitewise.R");
#gcinfo(TRUE)

# Summarize sites.
#summary = sites.summary(db="gj1_2x_57",sites.dir="${folder}")
#write.csv(summary,file="${folder}/summary_table.csv",row.names=F)

# Summarize the filtered version.
filter.fn = function(df) {
  return(subset(df,filter_value >= 3))
}
summary = sites.summary(db="gj1_2x_57",sites.dir="${folder}",filter.fn=filter.fn)
write.csv(summary,file="${folder}/summary_table_filt.csv",row.names=F)

# Summarize the filtered version.
filter.fn = function(df) {
  return(subset(df,!is.na(domain)))
}
summary = sites.summary(db="gj1_2x_57",sites.dir="${folder}",filter.fn=filter.fn)
write.csv(summary,file="${folder}/summary_table_domains.csv",row.names=F)


^;
  print "$cmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r($cmd,$params);

}

sub model_sites_distribution {
  my $self = shift;

  my $folder = $self->get_output_folder;

my $rcmd = qq^
source("../../scripts/collect_sitewise.R");

fits = sites.fit.distributions(db="gj1_2x_57",sites.dir="${folder}")
write.csv(fits,file="${folder}/distribution_fits.csv",row.names=F)

filter.fn = function(df) {
  return(subset(df,filter_value >= 3))
}
fits = sites.fit.distributions(db="gj1_2x_57",sites.dir="${folder}",filter.fn=filter.fn)
write.csv(fits,file="${folder}/distribution_fits_filt.csv",row.names=F)


filter.fn = function(df) {
  return(subset(df,omega < 1 & filter_value >= 3))
}
fits = sites.fit.distributions(db="gj1_2x_57",sites.dir="${folder}",filter.fn=filter.fn)
write.csv(fits,file="${folder}/distribution_fits_omega_lt_one.csv",row.names=F)

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

  my $rcmd = qq^
source("psg_overlaps.R");
load("${folder}/genes.Rdata")

genes[,'m.cor.psc.count'] = genes[,'m.weak.psc.count'] / genes[,'sitewise_value_count']
get.all.scores(genes,score.column='m.cor.psc.count')
save(roc.kosiol,roc.nielsen,roc.clark,roc.sabeti,roc.summary,file="${folder}/psg_rocs_m_cor.Rdata")

genes[,'p.cor.psc.count'] = genes[,'p.weak.psc.count'] / genes[,'sitewise_value_count']
get.all.scores(genes,score.column='p.cor.psc.count')
save(roc.kosiol,roc.nielsen,roc.clark,roc.sabeti,roc.summary,file="${folder}/psg_rocs_p_cor.Rdata")

get.all.scores(genes,score.column='m.slr.dnds')
save(roc.kosiol,roc.nielsen,roc.clark,roc.sabeti,roc.summary,file="${folder}/psg_rocs_m_slr.Rdata")

get.all.scores(genes,score.column='p.slr.dnds')
save(roc.kosiol,roc.nielsen,roc.clark,roc.sabeti,roc.summary,file="${folder}/psg_rocs_p_slr.Rdata")

get.all.scores(genes,score.column='m.psc.count')
save(roc.kosiol,roc.nielsen,roc.clark,roc.sabeti,roc.summary,file="${folder}/psg_rocs_m_pscs.Rdata")

get.all.scores(genes,score.column='p.psc.count')
save(roc.kosiol,roc.nielsen,roc.clark,roc.sabeti,roc.summary,file="${folder}/psg_rocs_p_pscs.Rdata")

get.all.scores(genes,score.column='m.weak.psc.count')
save(roc.kosiol,roc.nielsen,roc.clark,roc.sabeti,roc.summary,file="${folder}/psg_rocs_m_weak_pscs.Rdata")

get.all.scores(genes,score.column='p.weak.psc.count')
save(roc.kosiol,roc.nielsen,roc.clark,roc.sabeti,roc.summary,file="${folder}/psg_rocs_p_weak_pscs.Rdata")
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
