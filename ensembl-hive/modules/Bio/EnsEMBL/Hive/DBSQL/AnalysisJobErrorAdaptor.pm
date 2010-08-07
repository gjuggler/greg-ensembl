package Bio::EnsEMBL::Hive::DBSQL::AnalysisJobErrorAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::SqlHelper;

use base ('Bio::EnsEMBL::DBSQL::BaseAdaptor');

sub generic_fetch {
  my ($self) = @_;
  my $obj = ref($self);
  throw('generic_fetch() is not supported from '.$obj.' yet.');
}

=head2 store

  Arg [JOB]   : AnalysisJob instance
  Arg [ERROR] : A reference which represents an error; normally caught from a 
                die signal.
  Example     : $aje->store($aj, $@);
  Description : updates the analysis_job_error table in the hive database
  Returntype  : none
  Exceptions  : If JOB is not of the right type and if the error cannot be 
                written to the DB
  Caller      : AnalysisJob via Worker

=cut

sub store {
  my ($self, @args) = @_;
  
  my ($job, $error) = rearrange([qw(job error)], @args);
  assert_ref($job, 'Bio::EnsEMBL::Hive::AnalysisJob');
  $error = '-NO_ERROR_GIVEN-' if ! defined $error;
  
  my $helper = Bio::EnsEMBL::Utils::SqlHelper->new($self->dbc());
  my $sql = q{insert into analysis_job_error (analysis_job_id, status, 
worker_id, retry_count, error) values (?,?,?,?,?)};
  my $params = [$job->dbID(), $job->status(), $job->worker_id(), 
    $job->retry_count(), $error];
  $helper->execute_update(-SQL => $sql, -PARAMS => $params);
  return;
}

1;

