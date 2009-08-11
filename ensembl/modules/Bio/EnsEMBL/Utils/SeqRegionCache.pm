#
# Ensembl module for Bio::EnsEMBL::Utils::SeqRegionCache
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself

=head1 NAME

Bio::EnsEMBL::Utils::SeqRegionCache - A shared LRU cache of information about
seq_regions

=head1 SYNOPSIS

  use Bio::EnsEMBL::DBSQL::DBAdaptor;

  $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);

  $seq_region_cache = $db->get_SeqRegionCache();

  $key = "$seq_region_name:$coord_system_id";

  $array = $seq_region_cache->{$key};

  if($array) {
    $name = $array->[1];
    $length = $array->[3];
  } else {
    # cache miss, get the info from the database
    ...

    # cache the retrieved information
    $seq_region_cache->{$key} = [$seq_region_id, $seq_region_name,
                                 $coord_system_id, $seq_region_length];
  }



=head1 DESCRIPTION

This module is simply a convenient place to put a cache of sequence
region information which is shared by several adaptors for a given database.

=head1 CONTACT

This module is part of the Ensembl project http://www.ensembl.org

For more information email <ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut

use strict;
use Bio::EnsEMBL::Utils::Cache;

package Bio::EnsEMBL::Utils::SeqRegionCache;

our $SEQ_REGION_CACHE_SIZE = 40000;



sub new {
  my $class = shift;

  my %id_cache;
  my %name_cache;

  #
  # the items to cache should be listrefs to
  # [ sr_id, sr_name, cs_id, sr_length ]
  #
  # The name cache key is "sr_name:cs_id"
  # The id cache is keyed on "sr_id"
  #

  tie(%name_cache, 'Bio::EnsEMBL::Utils::Cache', $SEQ_REGION_CACHE_SIZE);
  tie(%id_cache, 'Bio::EnsEMBL::Utils::Cache', $SEQ_REGION_CACHE_SIZE);

  return bless {'name_cache' => \%name_cache, 
                'id_cache' => \%id_cache}, $class;
}


1;




