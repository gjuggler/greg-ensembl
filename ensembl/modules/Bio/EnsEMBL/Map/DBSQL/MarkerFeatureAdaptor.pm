# EnsEMBL module for MarkerFeatureAdaptor
# Copyright EMBL-EBI/Sanger center 2003
#
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Map::DBSQL::MarkerFeatureAdaptor

=head1 SYNOPSIS


=head1 DESCRIPTION

This object is responisble for all database interaction involving marker
features including the fetching and storing of marker features.

The bulk of this objects methods are inherited from 
Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor


=cut

package Bio::EnsEMBL::Map::DBSQL::MarkerFeatureAdaptor;

use strict;

use Bio::EnsEMBL::Map::MarkerFeature;
use Bio::EnsEMBL::Map::Marker;
use Bio::EnsEMBL::Map::MarkerSynonym;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);



=head2 fetch_all_by_Marker

  Arg [1]    : Bio::EnsEMBL::Map::Marker
  Example    : @ms = @{$marker_feature_adaptor->fetch_by_Marker($mrkr)};
  Description: Retrieves a list of MarkerFeatures for a given marker
  Returntype : listref of Bio::EnsEMBL::MarkerFeatures
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub fetch_all_by_Marker {
  my $self = shift;
  my $marker = shift;

  my $constraint = 'm.marker_id = ' . $marker->dbID;

  return $self->generic_fetch($constraint, @_);
}


=head2 fetch_all_by_Slice_and_priority

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Arg [2]    : (optional) int $priority
  Arg [3]    : (optional) int $map_weight
  Arg [3]    : (optional) string $logic_name
  Example    : @feats = @{$mfa->fetch_all_by_Slice_and_priority($slice,80,2)};
  Description: Retrieves all marker features above a specified threshold 
               priority which overlap the provided slice, below a 
               a specified map_weight.
  Returntype : listref of Bio::EnsEMBL::Map::MarkerFeatures in slice coords
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub fetch_all_by_Slice_and_priority {
  my ($self, $slice, $priority, $map_weight, @args) = @_;

  my $constraint = '';
  if(defined $priority) {
    $constraint = "m.priority > $priority";
  }

  if(defined $map_weight) {
    if($constraint) {
      $constraint .= " AND mf.map_weight < $map_weight";
    } else {
      $constraint = "mf.map_weight < $map_weight";
    }
  }

  return $self->fetch_all_by_Slice_constraint($slice, $constraint, @args);
}



=head2 fetch_all_by_RawContig_and_priority

  Description: DEPRECATED use fetch_all_by_Slice_and_priority instead

=cut

sub fetch_all_by_RawContig_and_priority {
  deprecate('Use fetch_all_by_Slice_and_priority() instead.');
  fetch_all_by_Slice_and_priority(@_);
}


sub fetch_all_by_Slice_and_score {
  my $self = shift;
  throw('fetch_all_by_Slice_and_score should not be' .
        'used to fetch marker_features');
}

sub _columns {
  my $self = shift;

  return ('mf.marker_feature_id', 'mf.marker_id',
          'mf.seq_region_id', 'mf.seq_region_start', 'mf.seq_region_end',
          'mf.analysis_id', 'mf.map_weight',
          'm.left_primer', 'm.right_primer', 'm.min_primer_dist',
          'm.max_primer_dist', 'm.priority', 'm.type', 'ms.marker_synonym_id',
          'ms.name', 'ms.source');
}

sub _tables {
  my $self = shift;

  return (['marker_feature', 'mf'], #primary table
	        ['marker', 'm'],
	        ['marker_synonym', 'ms']);
}

sub _left_join {
  my $self = shift;

  return ( [ 'marker_synonym',
             'm.display_marker_synonym_id = ms.marker_synonym_id' ] );
}

sub _default_where_clause {
  my $self = shift;

  return ('mf.marker_id = m.marker_id');
}

sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;

  my ($marker_feature_id, $marker_id, 
      $seq_region_id, $seq_region_start, $seq_region_end,
      $analysis_id, $map_weight,
      $left_primer, $right_primer, $min_primer_dist, $max_primer_dist,
      $priority, $type, $ms_id, $ms_name, $ms_source);

  #warning: ordering depends on _columns function implementation
  $sth->bind_columns(\$marker_feature_id, \$marker_id,
      \$seq_region_id, \$seq_region_start, \$seq_region_end,
      \$analysis_id, \$map_weight,
      \$left_primer, \$right_primer, \$min_primer_dist, \$max_primer_dist,
      \$priority, \$type, \$ms_id, \$ms_name, \$ms_source);

  my @out = ();

  my %marker_cache;
  my %slice_hash;
#  my %sr_name_hash;
  my %sr_cs_hash;
  my %analysis_cache;
  my $marker_adp = $self->db->get_MarkerAdaptor;
  my $sa  = $self->db->get_SliceAdaptor;
  my $analysis_adp = $self->db->get_AnalysisAdaptor;

  my $asm_cs;
  my $cmp_cs;
  my $asm_cs_vers;
  my $asm_cs_name;
  my $cmp_cs_vers;
  my $cmp_cs_name;
  if($mapper) {
    $asm_cs = $mapper->assembled_CoordSystem();
    $cmp_cs = $mapper->component_CoordSystem();
    $asm_cs_name = $asm_cs->name();
    $asm_cs_vers = $asm_cs->version();
    $cmp_cs_name = $cmp_cs->name();
    $cmp_cs_vers = $cmp_cs->version();
  }

  my $dest_slice_start;
  my $dest_slice_end;
  my $dest_slice_strand;
  my $dest_slice_length;
  if($dest_slice) {
    $dest_slice_start  = $dest_slice->start();
    $dest_slice_end    = $dest_slice->end();
    $dest_slice_strand = $dest_slice->strand();
    $dest_slice_length = $dest_slice->length();
  }

  FEATURE: while($sth->fetch) {
    #create a new marker unless this one has been seen already
    my $marker;
    if(!($marker = $marker_cache{$marker_id})) {
      #create a new marker synonym for the display synonym (if defined)
      my $ms;
      if($ms_id) {
        $ms = Bio::EnsEMBL::Map::MarkerSynonym->new
          ($ms_id, $ms_source, $ms_name);
      }

      #create a new marker
      $marker = Bio::EnsEMBL::Map::Marker->new
        ($marker_id, $marker_adp,
         $left_primer, $right_primer, $min_primer_dist, $max_primer_dist,
         $priority, $type, $ms);
      $marker_cache{$marker_id} = $marker;
    }

    #get the slice object
    my $slice = $slice_hash{$seq_region_id};

    if(!$slice) {
      $slice = $sa->fetch_by_seq_region_id($seq_region_id);
      $slice_hash{$seq_region_id} = $slice;
#      $sr_name_hash{$seq_region_id} = $slice->seq_region_name();
      $sr_cs_hash{$seq_region_id} = $slice->coord_system();
    }

    #retrieve analysis
    my $analysis;
    unless($analysis = $analysis_cache{$analysis_id}) {
      $analysis = $analysis_adp->fetch_by_dbID($analysis_id);
      $analysis_cache{$analysis_id} = $analysis;
    }

    #
    # remap the feature coordinates to another coord system
    # if a mapper was provided
    #
    if($mapper) {
#      my $sr_name = $sr_name_hash{$seq_region_id};
      my $sr_cs   = $sr_cs_hash{$seq_region_id};

     ($seq_region_id,$seq_region_start,$seq_region_end) =
        $mapper->fastmap($slice->seq_region_name(), $seq_region_start, $seq_region_end, 0, $sr_cs);

      #skip features that map to gaps or coord system boundaries
      next FEATURE if(!defined($seq_region_id));

      #get a slice in the coord system we just mapped to
      $slice = $slice_hash{"$seq_region_id"} ||=
	$sa->fetch_by_seq_region_id($seq_region_id);
    }

    #
    # If a destination slice was provided convert the coords
    # If the dest_slice starts at 1 and is foward strand, nothing needs doing
    #
    if($dest_slice) {
      if($dest_slice_start != 1 || $dest_slice_strand != 1) {
	if($dest_slice_strand == 1) {
	  $seq_region_start = $seq_region_start - $dest_slice_start + 1;
	  $seq_region_end   = $seq_region_end   - $dest_slice_start + 1;
	} else {
	  my $tmp_seq_region_start = $seq_region_start;
	  $seq_region_start = $dest_slice_end - $seq_region_end + 1;
	  $seq_region_end   = $dest_slice_end - $tmp_seq_region_start + 1;
	}

	#throw away features off the end of the requested slice
	if($seq_region_end < 1 || $seq_region_start > $dest_slice_length) {
	  next FEATURE;
	}
      }
      $slice = $dest_slice;
    }

    #now create a new marker_feature using the marker
    push @out, Bio::EnsEMBL::Map::MarkerFeature->new
      ($marker_feature_id, $self,
       $seq_region_start, $seq_region_end, $slice,
       $analysis, $marker_id, $map_weight, $marker);
  }

  return \@out;
}




=head2 store

  Arg [1]    : Bio::EnsEMBL::Map::MarkerFeature
  Example    : $marker_feature_adaptor->store(@marker_features);
  Description: Stores a list of marker features in this database.
               The dbID and adaptor of each marker will be set on successful 
               storing.
  Returntype : none
  Exceptions : thrown if not all data needed for storing is populated in the
               marker features
  Caller     : general
  Status     : stable

=cut

sub store {
  my ($self, @mfs) = @_;

  foreach my $mf (@mfs) {

    #
    # Sanity checking!
    #
    if(!ref($mf) || !$mf->isa('Bio::EnsEMBL::Map::MarkerFeature')) {
      $self->throw("Incorrect argument [$mf] to store.  Expected " .
                   'Bio::EnsEMBL::Map::MarkerFeature');
    }

    #don't store this feature if it has already been stored
    if($mf->is_stored($self->db())) {
      warning('MarkerFeature ['.$mf->dbID.'] is already stored in this DB.');
      next;
    }

    # Get/test the marker
    my $marker = $mf->marker;
    if(!$marker || !ref($marker) ||
       !$marker->isa('Bio::EnsEMBL::Map::Marker')) {
      throw('Cannot store MarkerFeature without an associated Marker');
    }
    
    #store the marker if it has not been stored yet 
    if(!$marker->is_stored($self->db())) {
      my $marker_adaptor = $self->db->get_adaptor('Marker');
      $marker_adaptor->store($marker);
    }
    my $marker_id = $marker->dbID ||
        throw('Associated Marker must have dbID to store MarkerFeature');

    # Get/test the analysis
    my $analysis = $mf->analysis;
    if(!$analysis || !ref($analysis) ||
       !$analysis->isa('Bio::EnsEMBL::Analysis')) {
      throw('Cannot store MarkerFeature without an associated Analysis');
    }

    #store the analysis if it has not been stored yet
    if(!$analysis->is_stored($self->db())) {
      my $analysis_adaptor = $self->db->get_adaptor('Analysis');
      $analysis_adaptor->store($mf->analysis());
    }
    my $analysis_id = $analysis->dbID ||
        throw('Associated Analysis must have dbID to store MarkerFeature');

    # Store the marker feature itself
    my $original = $mf;
    my $seq_region_id;
    ($mf, $seq_region_id) = $self->_pre_store($mf);

    my $sth =
      $self->prepare("INSERT INTO marker_feature (marker_id,
                           seq_region_id, seq_region_start, seq_region_end,
                           analysis_id, map_weight)
                      VALUES (?, ?, ?, ?, ?, ?)");
    $sth->execute($marker_id,
                  $seq_region_id, $mf->start, $mf->end,
                  $analysis_id, $mf->map_weight || 0);

    my $dbID = $sth->{'mysql_insertid'};

    $original->dbID($dbID);
    $original->adaptor($self);
  }
}


1;
