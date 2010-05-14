package Bio::Greg::PrimateHIV::OutputPrimateHIVData;

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
  
  my $base = '/nfs/users/nfs_g/gj1/scratch/primate_hiv/output';
  $self->get_output_folder($base);
  
  $self->export_windows;
}

sub export_windows {
  my $self = shift;

  my $genes_file = $self->get_output_folder . "/gene_windows.Rdata";
  my $genes_csv = $self->get_output_folder . "/gene_windows.csv";

  my $cmd = qq^
source("../../scripts/collect_sitewise.R");
source("./collect_primate_hiv.R");
genes <- get.genes.merged(db="gj1_hiv_57",pset.cols=c('alignment_slice_start'))
save(genes,file="${genes_file}")
write.csv(genes,file="${genes_csv}",row.names=F)

^;
  print "$cmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r($cmd,$params);
}


1;
