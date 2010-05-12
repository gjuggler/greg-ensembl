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
  
  my $base = 'nfs/users/nfs_g/gj1/scratch/primate_hiv/output';
  $self->get_output_folder($base);
  
  $self->export_windows;
}

sub export_windows {
  my $self = shift;

  my $genes_file = $self->get_output_folder . "/windows.Rdata";

  my $cmd = qq^
source("../../scripts/collect_sitewise.R");
source("./collect_primate_hiv.R");
genes <- get.genes.merged(db="gj1_hiv_57",exclude.cols=c('tree_newick'))
save(genes,file="${genes_file}")
^;
  print "$cmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r($cmd,$params);
}


1;
