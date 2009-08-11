#
# Ensembl module for Bio::EnsEMBL::AssemblyExceptionFeature
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::AssemblyExceptionFeature - A feature that represents an assembly exception

=head1 SYNOPSIS

use Bio::EnsEMBL::AssemblyExceptionFeature;

$feature = Bio::EnsEMBL::AssemblyExceptionFeature->new(-start    => 100,
                                                       -end      => 220,
                                                       -type     => 'HAP',
                                                       -slice    => $slice,
                                                       -adaptor  => $adaptor);

=head1 DESCRIPTION

Certain features, e.g. Haplotypes and PARs, are represented as "exceptions" to the 
normal assembly. This class represents such features.

=head1 AUTHOR - Glenn Proctor

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a _

=cut

use strict;

package Bio::EnsEMBL::AssemblyExceptionFeature;

use vars qw(@ISA);

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw(Bio::EnsEMBL::Feature);


=head2 new

  Arg [TYPE] : The type (e.g. HAP for haplotype, PAR for PAR)
  Arg [...]  : Named arguments passed to superclass
  Example    : $feature = Bio::EnsEMBL::AssemblyExceptionFeature->new
                        (-start           => 1,
                         -end             => 100,
                         -slice           => $slice,
                         -alternate_slice => $alt_slice,
                         -adaptor         => $adaptor,
                         -type            => 'HAP')
  Description: Constructs a new Bio::EnsEMBL::Feature.  Generally subclasses
               of this method are instantiated, rather than this class itself.
  Returntype : Bio::EnsEMBL::Feature
  Exceptions : Thrown on invalid -SLICE arguments
  Caller     : general, subclass constructors
  Status     : Stable

=cut

sub new {

  my $caller = shift;

  # allow this to be called as class or object method
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($type, $alternate_slice) = rearrange(['TYPE', 'ALTERNATE_SLICE'],@_);
  $self->{'type'} = $type;
  $self->{'alternate_slice'} = $alternate_slice;

  return $self;
}

=head2 new_fast

  Arg [1]    : hashref to be blessed
  Description: Construct a new Bio::EnsEMBL::Feature using the hashref.
  Exceptions : none
  Returntype : Bio::EnsEMBL::Feature
  Caller     : general, subclass constructors
  Status     : Stable

=cut  

sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}


=head2 type

  Arg [1]    : (optional) string $value
  Example    : $type = $assembly_exception_feature->type();
  Description: Getter/Setter for the type associated with this
               feature.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub type {

  my $self = shift;

  $self->{'type'} = shift if(@_);

  return $self->{'type'};
}


=head2 alternate_slice

  Arg [1]    : (optional) string $value
  Example    : $alt_slice = $assembly_exception_feature->alternate_slice();
  Description: Getter/Setter for the alternate slice associated with this feature.
               The alternate slice represents the "other side" of the AssemblyExceptionFeature.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub alternate_slice {

  my $self = shift;

  $self->{'alternate_slice'} = shift if(@_);

  return $self->{'alternate_slice'};
}



=head2 display_id

  Arg [1]    : none
  Example    : print $aef->display_id();
  Description: This method returns a string that is considered to be
               the 'display' identifier.  For assembly exception features
               this is the name of the alternate seqregion or '' if the 
               alternate slice is not defined.
  Returntype : string
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub display_id {
  my $self = shift;
  my $slice = $self->{'alternate_slice'};
  return '' if(!$slice);
  return $slice->seq_region_name();
}



1;
