package Bio::Greg::Hive::HiveLoaderUtils;

use warnings;
use strict;

use Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor;

use Bio::EnsEMBL::Hive::Utils 'stringify';  # import 'stringify()'

our $analysis_counter = 1;

sub new {
  my ($class,@args) = @_;
  my $self = bless {}, $class;
  
  return $self;
}

sub dba {
  my $self = shift;
  my $dba = shift;

  $self->{_dba} = $dba if (defined $dba);

  return $self->{_dba};
}

sub hive_dba {
  my $self = shift;
  my $hive_dba = shift;

  $self->{_hive_dba} = $hive_dba if (defined $hive_dba);

  return $self->{_hive_dba};
}

sub dbc {
  my $self = shift;
  return $self->dba->dbc;
}

sub init {
  my $self = shift;
  my $url = shift;

  my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-url => $url);
  my $hive_dba = Bio::EnsEMBL::Hive::DBSQL::DBAdaptor->new(-url => $url);
  $self->dba($dba);
  $self->hive_dba($hive_dba);
}


sub add_job_to_analysis {
  my $self = shift;
  my $analysis_name = shift;
  my $input_id_hash = shift;

  my $hive_dba = $self->hive_dba;

  my $analysis_adaptor = $hive_dba->get_AnalysisAdaptor;
  my $analysis         = $analysis_adaptor->fetch_by_logic_name($analysis_name);

  my $job_id = Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor->CreateNewJob(
    -input_id     => $input_id_hash,
    -analysis     => $analysis,
    -input_job_id => 0,
    );

  my $input_string = stringify($input_id_hash);

  print " -> Added job #$job_id to $analysis_name  [$input_string]\n";

  return $job_id;
}

sub create_analysis {
  my $self = shift;
  my $logic_name    = shift;
  my $module        = shift;
  my $params        = shift;
  my $hive_capacity = shift || 100;
  my $batch_size    = shift || 1;

  my $analysis_id = $analysis_counter++;

  print " -> Creating analysis $logic_name\n";

  my $param_string = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
  my $cmd          = qq{REPLACE INTO analysis SET
                 created=now(),
                 analysis_id=$analysis_id,
                 logic_name="$logic_name",
                 module="$module",
                 parameters="$param_string"
                 ;};
  my $dbc = $self->dbc;
  $dbc->do($cmd);
  $cmd = qq{REPLACE INTO analysis_stats SET
              analysis_id=$analysis_id,
              hive_capacity=$hive_capacity,
              batch_size=$batch_size,
              failed_job_tolerance=20000
              ;};
  $dbc->do($cmd);
  return $analysis_id;
}

sub connect_analysis {
  my $self = shift;
  my $from_name   = shift;
  my $to_name     = shift;
  my $branch_code = shift;
  $branch_code = 1 unless ( defined $branch_code );

  my $hive_dba = $self->hive_dba;
  my $dataflow_rule_adaptor = $hive_dba->get_DataflowRuleAdaptor;
  my $analysis_adaptor      = $hive_dba->get_AnalysisAdaptor;

  my $from_analysis = $analysis_adaptor->fetch_by_logic_name($from_name);
  my $to_analysis   = $analysis_adaptor->fetch_by_logic_name($to_name);

  if ( $from_analysis and $to_analysis ) {
    $dataflow_rule_adaptor->create_rule( $from_analysis, $to_analysis, $branch_code );
    warn "Created DataFlow rule: [$branch_code] $from_name -> $to_name\n";
  } else {
    die "Could not fetch analyses $from_analysis -> $to_analysis to create a dataflow rule";
  }
}


sub wait_for {
  my $self = shift;
  my $waiting_name  = shift;
  my $wait_for_list = shift;
  
  my $hive_dba = $self->hive_dba;
  my $ctrl_rule_adaptor = $hive_dba->get_AnalysisCtrlRuleAdaptor;
  my $analysis_adaptor  = $hive_dba->get_AnalysisAdaptor;
  
  my $waiting_analysis = $analysis_adaptor->fetch_by_logic_name($waiting_name);
  
  foreach my $wait_for_name (@$wait_for_list) {
    my $wait_for_analysis = $analysis_adaptor->fetch_by_logic_name($wait_for_name);
    
    if ( $waiting_analysis and $wait_for_analysis ) {
      $ctrl_rule_adaptor->create_rule( $wait_for_analysis, $waiting_analysis );
      warn "Created Control rule: $waiting_name will wait for $wait_for_name\n";
    } else {
      die "Could not fetch $waiting_name -> $wait_for_name to create a control rule";
    }
  }
}

sub clean_hive_tables {
  my $self = shift;
  my $dba = $self->dba;
  my @truncate_tables = qw^
      analysis 
      analysis_ctrl_rule
      analysis_data analysis_description
      analysis_job analysis_job_file
      analysis_stats analysis_stats_monitor
      dataflow_rule
      hive
      monitor
      resource_description
      ^;
  map {
    print "$_\n";
    eval { $dba->dbc->do("truncate table $_"); }
  } @truncate_tables;
}



sub combine_hashes {
  my $self = shift;
  my @hashes = @_;

  my $new_hash = {};
  foreach my $hash (@hashes) {
    foreach my $key (keys %{$hash}) {
      $new_hash->{$key} = $hash->{$key};
    }
  }
  return $new_hash;
}

1;
