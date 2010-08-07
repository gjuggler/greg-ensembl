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

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Utils::Argument;

@ISA = qw( Bio::EnsEMBL::DBSQL::DBAdaptor );


sub get_Queen {
  my $self = shift;

  return $self->get_QueenAdaptor();
}

sub new {
  my ($class, @args) = @_;

  my ($conf_file, $url, $species) = rearrange(['CONF_FILE', 'URL', 'SPECIES'], @args);

  if ($url and $url =~ /mysql\:\/\/([^\@]+\@)?([^\:\/]+)(\:\d+)?\/(.+)/) {
    my $user_pass = $1;
    my $host = $2;
    my $port = $3;
    my $dbname = $4;

    $user_pass =~ s/\@$//;
    my ($user, $pass) = $user_pass =~ m/([^\:]+)(\:.+)?/;
    $pass =~ s/^\:// if ($pass);
    $port =~ s/^\:// if ($port);
    push(@args, "-user" => $user) if ($user);
    push(@args, "-pass" => $pass) if ($pass);
    push(@args, "-port" => $port) if ($port);
    push(@args, "-host" => $host);
    push(@args, "-dbname" => $dbname);
    if (!$species) {
      push(@args, "-species" => $dbname);
    }
  }

  my $self = $class->SUPER::new(@args);
  return $self;
}

sub get_available_adaptors {
 
    my %pairs =  (
        'MetaContainer'       => 'Bio::EnsEMBL::DBSQL::MetaContainer',
        'Analysis'            => 'Bio::EnsEMBL::DBSQL::AnalysisAdaptor',
        'Queen'               => 'Bio::EnsEMBL::Hive::Queen',
        'AnalysisJob'         => 'Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor',
        'AnalysisJobError'    => 'Bio::EnsEMBL::Hive::DBSQL::AnalysisJobErrorAdaptor',
        'AnalysisData'        => 'Bio::EnsEMBL::Hive::DBSQL::AnalysisDataAdaptor',
        'AnalysisStats'       => 'Bio::EnsEMBL::Hive::DBSQL::AnalysisStatsAdaptor',
        'AnalysisCtrlRule'    => 'Bio::EnsEMBL::Hive::DBSQL::AnalysisCtrlRuleAdaptor',
        'DataflowRule'        => 'Bio::EnsEMBL::Hive::DBSQL::DataflowRuleAdaptor',
        'ResourceDescription' => 'Bio::EnsEMBL::Hive::DBSQL::ResourceDescriptionAdaptor',
    );
    return (\%pairs);
}
 
1;

