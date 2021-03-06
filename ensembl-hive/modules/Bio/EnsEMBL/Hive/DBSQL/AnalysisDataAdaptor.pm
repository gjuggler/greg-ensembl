# Perl module for Bio::EnsEMBL::Hive::DBSQL::AnalysisDataAdaptor
#
# Date of creation: 22.03.2004
# Original Creator : Jessica Severin <jessica@ebi.ac.uk>
#
# Copyright EMBL-EBI 2004
#
# You may distribute this module under the same terms as perl itself

=pod

=head1 NAME

  Bio::EnsEMBL::Hive::DBSQL::AnalysisDataAdaptor

=head1 SYNOPSIS

  $dataDBA = $db_adaptor->get_AnalysisDataAdaptor;

=head1 DESCRIPTION

   analysis_data table holds LONGTEXT data for use by the analysis system.
   This data is general purpose and it's up to each analysis to
   determine how to use it.
   This Adaptor module is used to access/store this data.

=head1 CONTACT

  Please contact ehive-users@ebi.ac.uk mailing list with questions/suggestions.

=head1 APPENDIX

  The rest of the documentation details each of the object methods.
  Internal methods are preceded with a _

=cut



package Bio::EnsEMBL::Hive::DBSQL::AnalysisDataAdaptor;

use strict;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base ('Bio::EnsEMBL::DBSQL::BaseAdaptor');


sub fetch_by_dbID {
  my ($self, $data_id) = @_;

  my $sql = "SELECT data FROM analysis_data WHERE analysis_data_id = ?";
  my $sth = $self->prepare($sql);
  $sth->execute($data_id);

  my ($data) = $sth->fetchrow_array();
  $sth->finish();
  return $data;
}

#
# STORE METHODS
#
################

sub store {
  my ($self, $data) = @_;
  my $data_id;
  
  return 0 unless($data);
  
  my $sth2 = $self->prepare("INSERT INTO analysis_data (data) VALUES (?)");
  $sth2->execute($data);
  $data_id = $sth2->{'mysql_insertid'};
  $sth2->finish;

  return $data_id;
}


sub store_if_needed {
  my ($self, $data) = @_;
  my $data_id;

  return 0 unless($data);

  my $sth = $self->prepare("SELECT analysis_data_id FROM analysis_data WHERE data = ?");
  $sth->execute($data);
  ($data_id) = $sth->fetchrow_array();
  $sth->finish;

  if($data_id) {
    # print("data already stored as id $data_id\n");
    return $data_id;
  }

  return $self->store($data);
}

1;





