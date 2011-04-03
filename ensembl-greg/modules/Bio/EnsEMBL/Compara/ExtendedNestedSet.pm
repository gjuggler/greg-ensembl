package Bio::EnsEMBL::Compara::ExtendedNestedSet;

use strict;
use Bio::EnsEMBL::Compara::NestedSet;

use base ('Bio::EnsEMBL::Compara::ProteinTree');

sub leaves {
  my $self = shift;

  my @leaves = @{$self->get_all_leaves};
  return @leaves;
}

sub nodes {
  my $self = shift;
  my @nodes = @{$self->get_all_nodes};
  return @nodes;
}

1;
