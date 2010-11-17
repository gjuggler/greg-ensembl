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
  my $genes_csv  = $self->get_output_folder . "/gene_windows.csv";

  my $output_folder = $self->get_output_folder;

  my $sitewise_script = Bio::Greg::EslrUtils->baseDirectory . "/projects/primate_hiv/collect_primate_hiv.R";

  my $dbname = $self->compara_dba->dbc->dbname;

  my $cmd = qq^
output_folder = "${output_folder}";
dbname = "$dbname"
source("$sitewise_script",echo=T);
^;
  print "$cmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r( $cmd, $params );
}

1;
