package Bio::EnsEMBL::Compara::LocalMember;

use strict;
use Bio::Seq;
use Bio::EnsEMBL::Utils::Argument;
use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Compara::GenomeDB;
use Bio::EnsEMBL::Compara::AlignedMember;
use Bio::EnsEMBL::Compara::Graph::CGObject;
use Bio::EnsEMBL::Compara::NestedSet;

our @ISA = qw(Bio::EnsEMBL::Compara::AlignedMember);

sub new {
  my ($class, @args) = @_;

  my $self = bless {}, $class;
  
  if (scalar @args) {
    #do this explicitly.
    my ($dbid, $stable_id, $description, $source_name, $adaptor, $taxon_id, $genome_db_id, $sequence_id, $cdna_sequence_id) = rearrange([qw(DBID STABLE_ID DESCRIPTION SOURCE_NAME ADAPTOR TAXON_ID GENOME_DB_ID SEQUENCE_ID CDNA_SEQUENCE_ID)], @args);

    $dbid && $self->dbID($dbid);
    $stable_id && $self->stable_id($stable_id);
    $description && $self->description($description);
    $source_name && $self->source_name($source_name);
    $adaptor && $self->adaptor($adaptor);
    $taxon_id && $self->taxon_id($taxon_id);
    $genome_db_id && $self->genome_db_id($genome_db_id);
    $sequence_id && $self->sequence_id($sequence_id);

  }

  return $self;
}

sub copy {
    my $self = shift;

    my $mycopy = $self->SUPER::copy;
    bless $mycopy, ref $self;

    $mycopy->cdna_sequence_id($self->cdna_sequence_id);
    return $mycopy;
}

# GJ 2009-01-15
sub cdna_sequence {
    my $self = shift;

    if(@_) {
	$self->{'_cdna_sequence'} = shift;
    }
    
    if(!defined($self->{'_cdna_sequence'}) and
       $self->cdna_sequence_id() != 0 and     
       defined($self->adaptor))
    {
	$self->{'_cdna_sequence'} = $self->adaptor->_fetch_sequence_by_id($self->cdna_sequence_id);
    }
    
    return $self->{'_cdna_sequence'};
}

sub cdna_sequence_id {
    my $self = shift;
    my $k = '_cdna_sequence_id';

    $self->{$k} = shift if (@_);
    $self->{$k} = 0 if (!defined $self->{$k});
    return $self->{$k};
}

1;
