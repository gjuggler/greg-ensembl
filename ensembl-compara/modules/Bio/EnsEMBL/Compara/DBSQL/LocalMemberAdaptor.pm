package Bio::EnsEMBL::Compara::DBSQL::LocalMemberAdaptor;

use strict;
use Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor;
use Bio::EnsEMBL::Compara::Member;
use Bio::EnsEMBL::Compara::Attribute;
use Bio::EnsEMBL::Compara::DBSQL::SequenceAdaptor;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

our @ISA = qw(Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor);

sub columns {
    my $self = shift;
    my $base_columns_ref = $self->SUPER::columns;
    my @arr = @{$base_columns_ref};

    push @arr, 'm.cdna_sequence_id';
    return \@arr;
}

sub create_instance_from_rowhash {
  my $self = shift;
  my $rowhash = shift;

  my $member = new Bio::EnsEMBL::Compara::LocalMember;
  $self->init_instance_from_rowhash($member, $rowhash);
  return $member;
}


sub init_instance_from_rowhash {
  my $self = shift;
  my $member = shift;
  my $rowhash = shift;

  $self->SUPER::init_instance_from_rowhash($member,$rowhash);

  $member->cdna_sequence_id($rowhash->{'cdna_sequence_id'});

  return $member;
}



1;
