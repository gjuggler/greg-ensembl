package Bio::Greg::Hive::SystemCommand;

use strict;
use Cwd;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::Process');

sub fetch_input {
  my $self = shift;

  ### DEFAULT PARAMETERS ###
  my $defaults = {
    cmd => qw^echo "No command given!\n"^,
    cwd => ''
  };
  ##########################

  $self->load_all_params($defaults);
}

sub run {
  my $self = shift;

  my $orig_cwd = cwd();
  my $cwd = $self->param('cwd');
  if (defined $cwd) {
    chdir($tmpdir);
  }

  my $cmd = $self->param('cmd');
  if (defined $cmd) {
    system($cmd);
  }

  chdir($cwd);
}

1;
