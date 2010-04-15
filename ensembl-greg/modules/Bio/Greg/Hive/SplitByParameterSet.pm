package Bio::Greg::Hive::SplitByParameterSet;

use strict;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::Process');

my $TREE = 'Bio::EnsEMBL::Compara::TreeUtils';

sub fetch_input {
  my $self = shift;

  ### DEFAULT PARAMETERS ###
  my $defaults = {
    flow_parameter_sets => 'all'
  };
  ##########################

  $self->load_all_params($defaults);
}

sub run {
  my $self = shift;

  my $flow_parameter_sets = $self->param('flow_parameter_sets');

  my @param_sets;
  if ($flow_parameter_sets eq 'all') {
    my $query = qq^select distinct(parameter_set_id) FROM parameter_set order by parameter_set_id;^;
    @param_sets = @{$self->db_handle->selectcol_arrayref($query)};
  } else {
    @param_sets = split(",",$flow_parameter_sets);
  }

  foreach my $param_set (@param_sets) {
    my $parameter_set_params = $self->get_params_from_param_set($param_set);
    my $input_params = $self->string_to_hash($self->input_id);
    my $output_params = { %$input_params, %$parameter_set_params };
    $output_params->{parameter_set_id} = $param_set;

    print "Flowing output: \n";
    $self->hash_print($output_params);
    my $job_arr = $self->dataflow_output_id($output_params,1);
    my $job_id = @{$job_arr}[0];
    print "  -> Job id: $job_id \n";
  }
}

sub write_output {
  
}

1;
