#
# BioPerl module for DBSQL::Obj
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

=pod

=head1 NAME

Bio::EnsEMBL::Hive::DBSQL::DBAdaptor

=head1 SYNOPSIS

    $db = Bio::EnsEMBL::Hive::DBSQL::DBAdaptor->new(
        -user   => 'root',
        -dbname => 'pog',
        -host   => 'caldy',
        -driver => 'mysql',
        );

=head1 DESCRIPTION

  This object represents the handle for a Hive system enabled database

=head1 CONTACT

  Please contact ehive-users@ebi.ac.uk mailing list with questions/suggestions.

=cut


package Bio::EnsEMBL::Hive::DBSQL::DBAdaptor;

use strict;
use Bio::EnsEMBL::DBSQL::DBConnection;

use base ('Bio::EnsEMBL::DBSQL::DBAdaptor');


#sub get_Queen {
#  my $self = shift;
#
#  return $self->get_QueenAdaptor();
#}

sub get_available_adaptors {
 
    my %pairs =  (
            # Core adaptors extended with Hive stuff:
        'MetaContainer'       => 'Bio::EnsEMBL::Hive::DBSQL::MetaContainer',
        'Analysis'            => 'Bio::EnsEMBL::Hive::DBSQL::AnalysisAdaptor',
            # Hive adaptors:
        'Queen'               => 'Bio::EnsEMBL::Hive::Queen',
        'AnalysisJob'         => 'Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor',
        'AnalysisData'        => 'Bio::EnsEMBL::Hive::DBSQL::AnalysisDataAdaptor',
        'AnalysisStats'       => 'Bio::EnsEMBL::Hive::DBSQL::AnalysisStatsAdaptor',
        'AnalysisCtrlRule'    => 'Bio::EnsEMBL::Hive::DBSQL::AnalysisCtrlRuleAdaptor',
        'DataflowRule'        => 'Bio::EnsEMBL::Hive::DBSQL::DataflowRuleAdaptor',
        'ResourceDescription' => 'Bio::EnsEMBL::Hive::DBSQL::ResourceDescriptionAdaptor',
        'NakedTable'          => 'Bio::EnsEMBL::Hive::DBSQL::NakedTableAdaptor',
        'JobMessage'          => 'Bio::EnsEMBL::Hive::DBSQL::JobMessageAdaptor',
    );
    return (\%pairs);
}
 
1;
