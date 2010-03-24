=head1 LICENSE

  Copyright (c) 1999-2010 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

package Bio::EnsEMBL::IdMapping::StableIdGenerator::AedesAegypti;

# Package that implements incrementing and verification of Aedes aegypti
# stable IDs as used by the VectorBase project.
# Based on Aedes_Aegypti.pm
# Differs from Aedes in that Exon stable ids like Ennnnnn not AAEL.ennnnnn
# and gene/transcript/translation start AGAP not AAEL
# also need to exclude old Ensembl-style ids

use strict;
use warnings;

use base qw(Bio::EnsEMBL::IdMapping::StableIdGenerator::EnsemblGeneric);

sub increment_stable_id {

  # This method will increment a stable ID.  For Anopheles, it will
  # pick out the numerical part of the stable ID (no matter what type of
  # stable ID it is) and increment it by one.  It will then replace the
  # numerical part by the incremented value and return the new stable
  # ID.  The parsing of the stable ID is very naive.

  my ( $self, $stable_id ) = @_;

  if ( !$self->is_valid($stable_id) ) {
    throw("Unknown or missing stable ID: $stable_id.");
  }

  $stable_id =~ /^(\D*)(\d+)(\D*)/;

  my $number_as_string = "$2";
  my $number           = $2 + 1;
  $stable_id = sprintf(
    "%s" . sprintf( "%%0%dd", length($number_as_string) ) . "%s",
    $1, $number, $3 );

  return $stable_id;
}

sub is_valid {

  # A stable ID is a valid Anopheles stable ID if it begins with the
  # character string "AGAP" or (for exons) just "E"
  # explicitly make the exon one  E+digits to exclude old-style ENSANG ids
  # otherwise ENSANGnnn found as higher then AGAPnnn
  # when initial_stable_id method checks archive tables

  my ( $self, $stable_id ) = @_;

  if ( !( defined($stable_id) && ( $stable_id =~ /^AGAP/ || $stable_id =~ /^E\d+$/ ) ) ) {
    return 0;
  }

  return 1;
}

1;
