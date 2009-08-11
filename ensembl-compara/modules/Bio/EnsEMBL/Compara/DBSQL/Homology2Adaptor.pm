package Bio::EnsEMBL::Compara::DBSQL::Homology2Adaptor;

use strict;
use Bio::EnsEMBL::Compara::Homology;
use Bio::EnsEMBL::Compara::DBSQL::BaseRelationAdaptor;
use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor;

our @ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 fetch_by_Member

 Arg [1]    : Bio::EnsEMBL::Compara::Member $member of type ENSEMBLGENE
 Example    : $homologies = $HomologyAdaptor->fetch_by_Member($member);
 Description: fetch the homology relationships where the given member is implicated
 Returntype : an array reference of Bio::EnsEMBL::Compara::Homology objects
 Exceptions : none
 Caller     : general

=cut

sub fetch_by_Member {
  # to be renamed fetch_all_by_Member
  my ($self, $member) = @_;

  my $join = [ [['homology_member', 'hm2'], 
                'h.homology_id = hm2.homology_id']
             ];

  my $constraint = "hm2.member_id = " .$member->dbID;;

  # This internal variable is used by add_Member_Attribute method 
  # in Bio::EnsEMBL::Compara::BaseRelation to make sure that the first element
  # of the member array is the one that has been used by the user to fetch the
  # homology object
  $self->{'_this_one_first'} = $member->stable_id;

  return $self->generic_fetch($constraint, $join);
}


sub fetch_all_by_MethodLinkSpeciesSet {
  my ($self, $method_link_species_set) = @_;

  throw("method_link_species_set arg is required\n")
    unless ($method_link_species_set);

  my $constraint = "h.method_link_species_set_id = " . $method_link_species_set->dbID;
  
  return $self->generic_fetch($constraint);
}


=head2 final_clause

  Arg [1]    : <string> SQL clause
  Example    : $adaptor->final_clause("ORDER BY h.description LIMIT 10");
               $homologies = $adaptor->fetch_all;
               $adaptor->final_clause("");
  Description: getter/setter method for specifying an extension to the SQL prior to
               a fetch operation.  Useful final clauses are either 'ORDER BY' or 'LIMIT'
  Returntype : <string>
  Caller     : general

=cut

sub final_clause {
  my $self = shift;
  $self->{'_final_clause'} = shift if(@_);
  return $self->{'_final_clause'};
}


#
# internal methods
#
###################

# internal methods used in multiple calls above to build homology objects from table data  

sub _tables {
  my $self = shift;

  return (['homology', 'h'], ['homology_member', 'hm'], ['member', 'm']);
}

sub _default_where_clause {
  my $self = shift;
  return "h.homology_id = hm.homology_id AND hm.member_id = m.member_id";
}

sub _columns {
  my $self = shift;
  
  return (
          'h.homology_id',
          'h.method_link_species_set_id', 
          'h.stable_id as homology_stable_id', 
          'h.description as homology_description',
          'h.subtype',
          'h.dn',
          'h.ds',
          'h.n',
          'h.s',
          'h.lnl',
          'h.threshold_on_ds',

          'hm.peptide_align_feature_id',
          'hm.cigar_line',
          'hm.cigar_start',
          'hm.cigar_end',
          'hm.perc_cov',
          'hm.perc_id',
          'hm.perc_pos',
          
          @{Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor->columns}
          );

}


sub create_Homology_from_rowhash {
  my ($self, $hashref) = @_;
  
  my $homology = Bio::EnsEMBL::Compara::Homology->new_fast
      ({'_adaptor' => $self,
       '_this_one_first' => $self->{'_this_one_first'},
       '_dbID' => $hashref->{'homology_id'},
       '_stable_id' => $hashref->{'homology_stable_id'},
       '_description' => $hashref->{'homology_description'},
       '_method_link_species_set_id' => $hashref->{'method_link_species_set_id'},
       '_subtype' => $hashref->{'subtype'},
       '_dn' => $hashref->{'dn'},
       '_ds' => $hashref->{'ds'},
       '_n' => $hashref->{'n'},
       '_s' => $hashref->{'s'},
       '_lnl' => $hashref->{'lnl'},
       '_threshold_on_ds' => $hashref->{'threshold_on_ds'},
       });
  return $homology;
}


sub create_homology_attribute_from_rowhash {
  my $self = shift;
  my $rowhash = shift;

  my $attribute = new Bio::EnsEMBL::Compara::Attribute;
  $attribute->member_id($rowhash->{'member_id'});
  $attribute->homology_id($rowhash->{'homology_id'});
  $attribute->peptide_member_id($rowhash->{'peptide_member_id'});
  $attribute->peptide_align_feature_id($rowhash->{'peptide_align_feature_id'});
  $attribute->cigar_line($rowhash->{'cigar_line'});
  $attribute->cigar_start($rowhash->{'cigar_start'});
  $attribute->cigar_end($rowhash->{'cigar_end'});
  $attribute->perc_cov($rowhash->{'perc_cov'});
  $attribute->perc_id($rowhash->{'perc_id'});
  $attribute->perc_pos($rowhash->{'perc_pos'});

  return $attribute;
}


sub _objs_from_sth {
  my ($self, $sth) = @_;

  my %homology_hash;

  while(my $rowhash = $sth->fetchrow_hashref) {
    #printf("create HM: %d - %s(%d)\n",  $rowhash->{'homology_id'},  $rowhash->{'m_stable_id'}, $rowhash->{'member_id'});
    my ($member,$attribute);
    $member = Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor->create_instance_from_rowhash($rowhash);    
    $attribute = $self->create_homology_attribute_from_rowhash($rowhash);
    
    my $homology = $homology_hash{$rowhash->{'homology_id'}};
    unless($homology) {
      #printf("create Homology for id=%d\n", $rowhash->{'homology_id'});
      $homology = $self->create_Homology_from_rowhash($rowhash);
      $homology_hash{$rowhash->{'homology_id'}} = $homology;
    }
    $homology->add_Member_Attribute([$member, $attribute]);
  }
  $sth->finish;
  my @homologies = values(%homology_hash);
  return \@homologies;
}

1;
