package Bio::EnsEMBL::Compara::BaseRelation;

use strict;
use Bio::EnsEMBL::Utils::Argument;
use Bio::EnsEMBL::Utils::Exception;

sub new {
  my ($class, @args) = @_;

  my $self = bless {}, $class;

  if (scalar @args) {
    #do this explicitly.
    my ($dbid, $stable_id, $method_link_species_set_id, $method_link_type, $description, $adaptor) = rearrange([qw(DBID STABLE_ID METHOD_LINK_SPECIES_SET_ID METHOD_LINK_TYPE DESCRIPTION  ADAPTOR)], @args);
    
    $dbid && $self->dbID($dbid);
    $stable_id && $self->stable_id($stable_id);
    $description && $self->description($description);
    $method_link_species_set_id && $self->method_link_species_set_id($method_link_species_set_id);
    $method_link_type && $self->method_link_type($method_link_type);
    $adaptor && $self->adaptor($adaptor);
  }
  
  return $self;
}   

=head2 new_fast

  Arg [1]    : hash reference $hashref
  Example    : none
  Description: This is an ultra fast constructor which requires knowledge of
               the objects internals to be used.
  Returntype : 
  Exceptions : none
  Caller     : 

=cut

sub new_fast {
  my ($class, $hashref) = @_;

  return bless $hashref, $class;
}

=head2 dbID

  Arg [1]    : int $dbID (optional)
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub dbID {
  my $self = shift;
  $self->{'_dbID'} = shift if(@_);
  return $self->{'_dbID'};
}

=head2 stable_id

  Arg [1]    : string $stable_id (optional)
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub stable_id {
  my $self = shift;
  $self->{'_stable_id'} = shift if(@_);
  return $self->{'_stable_id'};
}

=head2 description

  Arg [1]    : string $description (optional)
  Example    : 
  Description: 
  Returntype : string
  Exceptions : 
  Caller     : 

=cut

sub description {
  my $self = shift;
  $self->{'_description'} = shift if(@_);
  return $self->{'_description'};
}

=head2 method_link_species_set

  Arg [1]    : MethodLinkSpeciesSet object (optional)
  Example    : 
  Description: getter/setter method for the MethodLinkSpeciesSet for this relation
  Returntype : Bio::EnsEMBL::Compara::MethodLinkSpeciesSet
  Exceptions : 
  Caller     : 

=cut

sub method_link_species_set {
  my $self = shift;

  if(@_) {
    my $mlss = shift;
    unless ($mlss->isa('Bio::EnsEMBL::Compara::MethodLinkSpeciesSet')) {
      throw("Need to add a Bio::EnsEMBL::Compara::MethodLinkSpeciesSet, not a $mlss\n");
    }
    $self->{'_method_link_species_set'} = $mlss;
    $self->{'_method_link_species_set_id'} = $mlss->dbID;
  }

  #lazy load from method_link_species_set_id
  if ( ! defined $self->{'_method_link_species_set'} && defined $self->method_link_species_set_id) {
    my $mlssa = $self->adaptor->db->get_MethodLinkSpeciesSetAdaptor;
    my $mlss = $mlssa->fetch_by_dbID($self->method_link_species_set_id);
    $self->{'_method_link_species_set'} = $mlss;
  }

  return $self->{'_method_link_species_set'};
}

=head2 method_link_species_set_id

  Arg [1]    : integer (optional)
  Example    : 
  Description: 
  Returntype : integer
  Exceptions : 
  Caller     : 

=cut

sub method_link_species_set_id {
  my $self = shift;

  $self->{'_method_link_species_set_id'} = shift if (@_);
  return $self->{'_method_link_species_set_id'};
}

=head2 method_link_type

  Arg [1]    : string $method_link_type (optional)
  Example    : 
  Description: 
  Returntype : string
  Exceptions : 
  Caller     : 

=cut

sub method_link_type {
  my $self = shift;

  $self->{'_method_link_type'} = shift if (@_);
  unless (defined $self->{'_method_link_type'}) {
    my $mlss = $self->method_link_species_set;
    throw("method_link_type needs a valid method_link_species_set") unless($mlss);
    $self->{'_method_link_type'} = $mlss->method_link_type;
  }

  return $self->{'_method_link_type'};
}

=head2 method_link_id

  Arg [1]    : integer (optional)
  Example    : 
  Description: 
  Returntype : integer
  Exceptions : 
  Caller     : 

=cut

sub method_link_id {
  my $self = shift;

  $self->{'_method_link_id'} = shift if (@_);
  unless (defined $self->{'_method_link_id'}) {
    my $mlss = $self->method_link_species_set;
    $self->{'_method_link_id'} = $mlss->method_link_id;
  }

  return $self->{'_method_link_id'};
}

=head2 adaptor

  Arg [1]    : string $adaptor (optional)
               corresponding to a perl module
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub adaptor {
  my $self = shift;
  $self->{'_adaptor'} = shift if(@_);
  return $self->{'_adaptor'};
}

sub add_Member_Attribute {
  my ($self, $member_attribute) = @_;

  my ($member, $attribute) = @{$member_attribute};

  throw("member argument not defined\n") unless($member);
  throw("attribute argument not defined\n") unless($attribute);
  
  unless ($member->isa('Bio::EnsEMBL::Compara::Member')) {
    throw("Need to add a Bio::EnsEMBL::Compara::Member, not a $member\n");
  }
  unless ($attribute->isa('Bio::EnsEMBL::Compara::Attribute')) {
    throw("Need to add a Bio::EnsEMBL::Compara::Attribute, not a $attribute\n");
  }

  if (defined $self->{'_this_one_first'} && $self->{'_this_one_first'} eq $member->stable_id) {
    unshift @{$self->{'_member_array'}}, $member_attribute ;
    unshift @{$self->{'_members_by_source'}{$member->source_name}}, $member_attribute;
    unshift @{$self->{'_members_by_source_taxon'}{$member->source_name."_".$member->taxon_id}}, $member_attribute;
  } else {
    push @{$self->{'_member_array'}}, $member_attribute ;
    push @{$self->{'_members_by_source'}{$member->source_name}}, $member_attribute;
    push @{$self->{'_members_by_source_taxon'}{$member->source_name."_".$member->taxon_id}}, $member_attribute;
  }
}


=head2 get_all_Member_Attribute

  Arg [1]    : None
  Example    : 
  Description: 
  Returntype : array reference of [Bio::EnsEMBL::Compara::Member, Bio::EnsEMBL::Compara::Attribute]
  Exceptions : 
  Caller     : 

=cut

sub get_all_Member_Attribute {
  my ($self) = @_;
  
  unless (defined $self->{'_member_array'}) {

    my $MemberAdaptor = $self->adaptor->db->get_MemberAdaptor();
    my $members = $MemberAdaptor->fetch_by_relation($self);

    $self->{'_member_array'} = [];
    $self->{'_members_by_source'} = {};
    $self->{'_members_by_source_taxon'} = {};
    foreach my $member_attribute (@{$members}) {
      $self->add_Member_Attribute($member_attribute);
    }
  }
  return $self->{'_member_array'}; #should return also attributes
}


=head2 get_all_Members

  Arg [1]    : None
  Example    : 
  Description: 
  Returntype : array reference of Bio::EnsEMBL::Compara::Member
  Exceptions : 
  Caller     : 

=cut

sub get_all_Members {
  my ($self) = @_;

  my $members = [];
  foreach my $member_attribute (@{$self->get_all_Member_Attribute}) {
    my ($member, $attribute) = @$member_attribute;
    push (@$members, $member);
  }

  return $members;
}


=head2 get_Member_Attribute_by_source

  Arg [1]    : string $source_name
               e.g. "ENSEMBLPEP"
  Example    : 
  Description: 
  Returntype : array reference of Bio::EnsEMBL::Compara::Member
  Exceptions : 
  Caller     : 

=cut

sub get_Member_Attribute_by_source {
  my ($self, $source_name) = @_;

  throw("Should give defined source_name as arguments\n") unless (defined $source_name);

  $self->get_all_Member_Attribute;

  $self->{'_members_by_source'}->{$source_name} = [] 
    unless(defined($self->{'_members_by_source'}->{$source_name}));
    
  return $self->{'_members_by_source'}->{$source_name};
}

=head2 get_Member_Attribute_by_source_taxon

  Arg [1]    : string $source_name
  Arg [2]    : int $taxon_id
  Example    : $domain->get_Member_by_source_taxon('ENSEMBLPEP',9606)
  Description: 
  Returntype : array reference of Bio::EnsEMBL::Compara::Member
  Exceptions : 
  Caller     :

=cut

sub get_Member_Attribute_by_source_taxon {
  my ($self, $source_name, $taxon_id) = @_;

  throw("Should give defined source_name and taxon_id as arguments\n") unless (defined $source_name && defined $taxon_id);
  $self->get_all_Member_Attribute;  

  $self->{'_members_by_source_taxon'}->{$source_name."_".$taxon_id} = []
    unless(defined($self->{'_members_by_source_taxon'}->{$source_name."_".$taxon_id}));

  return $self->{'_members_by_source_taxon'}->{$source_name."_".$taxon_id};
}

=head2 Member_count_by_source

  Arg [1]    : string $source_name
               e.g. "ENSEMBLPEP"
  Example    : $domain->Member_count_by_source('ENSEMBLPEP');
  Description: 
  Returntype : int
  Exceptions : 
  Caller     : 

=cut

sub Member_count_by_source {
  my ($self, $source_name) = @_; 
  
  throw("Should give a defined source_name as argument\n") unless (defined $source_name);
  
  return scalar @{$self->get_Member_Attribute_by_source($source_name)};
}

=head2 Member_count_by_source_taxon

  Arg [1]    : string $source_name
  Arg [2]    : int $taxon_id
  Example    : Member_count_by_source_taxon('ENSEMBLPEP',9606);
  Description: 
  Returntype : int
  Exceptions : 
  Caller     : 

=cut

sub Member_count_by_source_taxon {
  my ($self, $source_name, $taxon_id) = @_; 
  
  throw("Should give defined source_name and taxon_id as arguments\n") unless (defined $source_name && defined $taxon_id);

  return scalar @{$self->get_Member_Attribute_by_source_taxon($source_name,$taxon_id)};
}

#
# DEPRECATED METHODS
####################

sub source_id {
  my $self = shift;
  deprecate("source method is deprecated. Calling $self->method_link_id instead\n");

  $self->{'_method_link_id'} = shift if (@_);
  return $self->method_link_id;
}

sub source_name {
  my $self = shift;
  deprecate("source_name method is now deprecated. Calling method_link_type instead.\n");

  $self->{'_method_link_type'} = shift if (@_);
  return $self->method_link_type;
}

sub known_sources {
  my ($self) = @_;
  deprecate();
  throw("Get this data from the Bio::EnsEMBL::Compara::MethodLinkSpeciesSet object\n");
}

1;
