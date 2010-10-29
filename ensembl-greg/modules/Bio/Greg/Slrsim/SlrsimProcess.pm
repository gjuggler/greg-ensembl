package Bio::Greg::Slrsim::SlrsimProcess;

use strict;
use Bio::EnsEMBL::Hive::Process;
use Bio::Greg::Hive::Process;
use File::Path;

use base ('Bio::Greg::Hive::Process');

sub get_output_folder {
  my $self = shift;
  my $subfolder = shift;

  my $base = '/homes/greg/lib/greg-ensembl/projects/slrsim/NO.backup/output';

  my $dated_base = $self->SUPER::get_output_folder($base);

  my $full_folder = $dated_base;
  $full_folder = "${dated_base}/$subfolder" if ($subfolder ne '');

  mkpath([$full_folder]);
  return $full_folder;
}

1;
