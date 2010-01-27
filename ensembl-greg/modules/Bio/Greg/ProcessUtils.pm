package Bio::Greg::ProcessUtils;

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;

use Bio::EnsEMBL::Hive;

our @ISA = qw();

sub replace_params {
  my $self = shift;  
  my @params = @_;

  my $final_params = shift @params;
  while (@params) {
    my $p = shift @params;
    if (!ref $p) {
      $p = eval($p);
    }
    
    my $new_params = {};
    foreach my $key (keys %{$final_params}) {
      $new_params->{$key} = $final_params->{$key};
    }
    foreach my $key (keys %{$p}) {
      $new_params->{$key} = $p->{$key};
    }
    $final_params = $new_params;
  }

  return $final_params;
}

sub get_params {
  my $self = shift;
  my $param_string = shift;
  
  my $param_object;
  my $tmp_params = eval($param_string);  
  foreach my $key (keys %$tmp_params) {
    $param_object->{$key} = $tmp_params->{$key};
  }
  return $param_object;
}

1;
