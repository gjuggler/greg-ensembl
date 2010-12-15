package Bio::Greg::Gorilla::MapGorillaSubs;

use strict;

use base (
  'Bio::Greg::Hive::Process', 'Bio::Greg::StatsCollectionUtils',
  'Bio::Greg::Hive::Align',   'Bio::Greg::Hive::CountSubstitutions'
);

my $TU = 'Bio::EnsEMBL::Compara::TreeUtils';
my $CU = 'Bio::EnsEMBL::Compara::ComparaUtils';
my $AU = 'Bio::EnsEMBL::Compara::AlignUtils';

sub param_defaults {
  return {
  };
}

sub fetch_input {
  my ($self) = @_;

}

sub run {

  my $base           = Bio::Greg::EslrUtils->baseDirectory;

}
