package Bio::Greg::Slrsim::Plots;

use strict;
use Cwd;
use POSIX qw(strftime mktime);
use Bio::Greg::EslrUtils;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::Process');

sub fetch_input {
  my $self = shift;

  $self->load_all_params();
}

sub run {
  my $self = shift;

  my $old_cwd = cwd();
  chdir("/homes/greg/lib/greg-ensembl/projects/slrsim");

#  my $experiment_name = $self->param('experiment_name');

#  if ($experiment_name eq 'Filter Sweeps') {
    $self->filter_sweep_roc;
#  }

  chdir($old_cwd);
}


sub filter_sweep_roc {
  my $self = shift;

  my $date_string = strftime("%Y-%m-%d",localtime);
  my $i=0;
  
  my $filename;
  do {
    $i++;
    $filename = "output/".$date_string."_".$i.".png";
  } while (-e $filename);

  my $rcmd = qq^
source("collect_slrsim.R")
data = get.all.data()
roc = filter.sweep.roc(data)
png(file="${filename}",width=600,height=600)
plot.roc(roc)
dev.off()
^;
  print "$rcmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r($rcmd,$params);

}

1;
