#
# Ensembl module for Bio::EnsEMBL::Compara::DBSQL::GenomicAlignGroupAdaptor
#
# Copyright Javier Herrero
#
# You may distribute this module under the same terms as perl itself
# 
# POD documentation - main docs before the code
# 

=head1 NAME

Bio::EnsEMBL::Compara::DBSQL::GenomicAlignGroupAdaptor - Object to access data in genomic_align_group table

=head1 SYNOPSIS

  use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor; 
  my $db = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor (
      -host => $host,
      -user => $dbuser,
      -pass => $dbpass,
      -port => $port,
      -dbname => $dbname,
      -conf_file => $conf_file);
  
  my $genomic_align_group_adaptor = $db->get_GenomicAlignGroupAdaptor();

  $genomic_align_group_adaptor->store($genomic_align_group);

  $genomic_align_groups = $genomic_align_group_adaptor->fetch_all_by_GenomicAlign($genomic_align);
  $genomic_align_groups = $genomic_align_group_adaptor->fetch_all_by_genomic_align_id(11223);

  $genomic_align_group = $genomic_align_group_adaptor->fetch_by_GenomicAlign($genomic_align);
  $genomic_align_group = $genomic_align_group_adaptor->fetch_by_genomic_align_id(11223);

=head1 DESCRIPTION

This class is intended to access data in genomic_align_group table

=head1 AUTHOR

Javier Herrero (jherrero@ebi.ac.uk)

This modules is part of the Ensembl project http://www.ensembl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Compara::DBSQL::GenomicAlignGroupAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Compara::GenomicAlignGroup;
use Bio::EnsEMBL::Compara::GenomicAlign;
use Bio::EnsEMBL::Compara::DnaFrag;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 new

  Arg [1]    : list of args to super class constructor
  Example    : $gag_a = new Bio::EnsEMBL::Compara::GenomicAlignGroupAdaptor($dbobj);
  Description: Creates a new GenomicAlignGroupAdaptor. This
               class should be instantiated through the get method on the 
               DBAdaptor rather than calling this method directly.
  Returntype : none
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::DBConnection
  Status     : Stable

=cut

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->{_use_autoincrement} = 1;

  return $self;
}


=head2 store

  Arg  1     : Bio::EnsEMBL::Compara::GenomicAlignGroup
               The things you want to store
  Example    : $gen_ali_grp_adaptor->store($genomic_align_group);
  Description: It stores the given GenomicAlginGroup in the database
  Returntype : Bio::EnsEMBL::Compara::GenomicAlignGroup object
  Exceptions : not stored linked Bio::EnsEMBL::Compara::GenomicAlign objects throw.
  Caller     : general
  Status     : Stable

=cut

sub store {
  my ($self, $genomic_align_group) = @_;

  my $genomic_align_group_sql =
        qq{INSERT INTO genomic_align_group (
                node_id,
                genomic_align_id
        ) VALUES (?,?)};
  
  my @values;
  
  ## CHECKING
  my $all_genomic_aligns = $genomic_align_group->get_all_GenomicAligns;
  foreach my $genomic_align (@$all_genomic_aligns) {
    # check if every GenomicAlgin has a dbID
    if (!defined($genomic_align->dbID)) {
      throw("GenomicAlign [$genomic_align] in GenomicAlignGroup is not in DB");
    }
  }
  
  my $node_id = $genomic_align_group->dbID;
  if (!$node_id) {
    my $method_link_species_set_id;
    # Get common method_link_species_set_id
    foreach my $genomic_align (@$all_genomic_aligns) {
      if (!$genomic_align->method_link_species_set_id()) {
        ## undef value and exit loop if method_link_species_set_id does not match
        $method_link_species_set_id = undef;
        last;
      } elsif (!$method_link_species_set_id) {
        $method_link_species_set_id = $genomic_align->method_link_species_set_id();
      } elsif ($method_link_species_set_id != $genomic_align->method_link_species_set_id()) {
        ## undef value and exit loop if method_link_species_set_id does not match
        $method_link_species_set_id = undef;
        last;
      }
    }
    if ($method_link_species_set_id && !$self->use_autoincrement()) {
      ## Only if method_link_species_set_id is the same for all the GenomicAligns
      my $sql = 
              "SELECT MAX(node_id) FROM genomic_align_group WHERE".
              " node_id > ".$method_link_species_set_id.
              "0000000000 AND node_id < ".
              ($method_link_species_set_id + 1)."0000000000";
      my $sth = $self->prepare($sql);
      $sth->execute();
      $node_id = $sth->fetchrow_array();
      if (defined $node_id) {
        $node_id++;
      } else {
        $node_id = $method_link_species_set_id * 10000000000 + 1;
      }
    }
  }

  ## Stores data, all of them with the same id
  my $sth = $self->prepare($genomic_align_group_sql);
  for (my $i = 0; $i < @$all_genomic_aligns; $i++) {
    my $genomic_align  = $all_genomic_aligns->[$i];
    $sth->execute(
                  ($node_id or "NULL"),
                  $genomic_align->dbID
          );
    if (!$node_id) {$node_id = $sth->{'mysql_insertid'};}

    info("Stored Bio::EnsEMBL::Compara::GenomicAlignGroup ".
          "(".($i+1)."/".scalar(@$all_genomic_aligns).") ".
          ($node_id or "NULL").", ".$genomic_align->dbID, );

  }
  $genomic_align_group->dbID($node_id);
  
  return $genomic_align_group;
}


=head2 fetch_by_dbID

  Arg  1     : integer node_id
  Example    : my $genomic_align_group =
                  $genomic_align_group_adaptor->fetch_by_dbID(12413)
  Description: Returns a Bio::EnsEMBL::Compara::GenomicAlignGroup corresponding
               to the given node_id.
  Returntype : Bio::EnsEMBL::Compara::GenomicAlignGroup
  Exceptions : none
  Caller     : object::methodname
  Status     : Stable

=cut

sub fetch_by_dbID {
  my ($self, $node_id) = @_;
  my $genomic_align_group;

  my $genomic_align_adaptor = $self->db->get_GenomicAlignAdaptor;
  
  my $genomic_align_group_sql = qq{
              SELECT
                  node_id,
                  genomic_align_id
              FROM
                genomic_align_group
              WHERE
                node_id = ?
        };
  
  my @values;

  my $sth = $self->prepare($genomic_align_group_sql);
  $sth->execute($node_id);

  # Group results in order to be able to build Bio::EnsEMBL::Compara::GenomicAlignGroupAdaptor objects
  my $group;
  while (my $values = $sth->fetchrow_arrayref) {
    my ($node_id, $genomic_align_id) = @$values;

    $group->{'node_id'} = $node_id;
    my $this_genomic_align = new Bio::EnsEMBL::Compara::GenomicAlign(
            -dbID => $genomic_align_id,
            -adaptor => $genomic_align_adaptor
        );

    push(@{$group->{'genomic_align_array'}}, $this_genomic_align);
  }

  #return undef if no genomic_align_groups have been found
  if (!defined $group->{'node_id'}) {
      return $genomic_align_group;
  }

  $genomic_align_group = new Bio::EnsEMBL::Compara::GenomicAlignGroup(
          -dbID => $group->{'node_id'},
          -adaptor => $self,
          -genomic_align_array => $group->{'genomic_align_array'}
      );
  foreach my $this_genomic_align (@{$genomic_align_group->genomic_align_array}) {
    $this_genomic_align->genomic_align_group($genomic_align_group);
  }

  return $genomic_align_group;
}

=head2 fetch_by_GenomicAlign

  Arg  1     : Bio::EnsEMBL::Compara::GenomicAlign $genomic_align
  Example    : my $genomic_align_group =
                  $genomic_align_group_adaptor->fetch_by_GenomicAlign(
                          $genomic_align)
  Description: Returns a Bio::EnsEMBL::Compara::GenomicAlignGroup corresponding
               to the given Bio::EnsEMBL::Compara::GenomicAlign and group_type.
  Returntype : Bio::EnsEMBL::Compara::GenomicAlignGroup
  Exceptions : none
  Caller     : object::methodname
  Status     : Stable

=cut

sub fetch_by_GenomicAlign {
  my ($self, $genomic_align) = @_;
  my $genomic_align_group;

  my $genomic_align_adaptor = $self->db->get_GenomicAlignAdaptor;
  
  # Check Bio::EnsEMBL::Compara::GenomicAlign object
  unless($genomic_align && ref $genomic_align && 
        $genomic_align->isa('Bio::EnsEMBL::Compara::GenomicAlign')) {
    throw("[$genomic_align] must be a Bio::EnsEMBL::Compara::GenomicAlign object");
  }

  my $genomic_align_block_sql = qq{
              SELECT
                  b.node_id,
                  b.genomic_align_id
              FROM
                genomic_align_group a, genomic_align_group b
              WHERE
                a.node_id = b.node_id
                AND a.genomic_align_id = ?
        };
  
  my @values;
  
  my $sth = $self->prepare($genomic_align_block_sql);
  $sth->execute($genomic_align->dbID);

  # Group results in order to be able to build Bio::EnsEMBL::Compara::GenomicAlignGroupAdaptor objects
  my $group;
  while (my $values = $sth->fetchrow_arrayref) {
    my ($node_id, $genomic_align_id) = @$values;

    $group->{'node_id'} = $node_id;
    my $this_genomic_align;
    if ($genomic_align_id == $genomic_align->dbID) {
      # Use Bio::EnsEMBL::Compara::GenomicAlign object given by argument if possible
      $this_genomic_align = $genomic_align;
    } else {
      # Create a new Bio::EnsEMBL::Compara::GenomicAlign object otherwise
      $this_genomic_align = new Bio::EnsEMBL::Compara::GenomicAlign(
              -dbID => $genomic_align_id,
              -adaptor => $genomic_align_adaptor
          );
    }

    push(@{$group->{'genomic_align_array'}}, $this_genomic_align);
  }

  $genomic_align_group = new Bio::EnsEMBL::Compara::GenomicAlignGroup(
          -dbID => $group->{'node_id'},
          -adaptor => $self,
          -genomic_align_array => $group->{'genomic_align_array'}
      );
  foreach my $this_genomic_align (@{$genomic_align_group->genomic_align_array}) {
    $this_genomic_align->genomic_align_group($genomic_align_group);
  }

  return $genomic_align_group;
}


=head2 fetch_by_genomic_align_id

  Arg  1     : integer $genomic_align_id
  Example    : my $genomic_align_group =
                  $genomic_align_group_adaptor->fetch_by_genomic_align_id(
                          12322)
  Description: Returns a Bio::EnsEMBL::Compara::GenomicAlignGroup corresponding
               to the given Bio::EnsEMBL::Compara::GenomicAlign defined by the
               $genomic_align_id
  Returntype : Bio::EnsEMBL::Compara::GenomicAlignGroup
  Exceptions : none
  Caller     : object::methodname
  Status     : Stable

=cut

sub fetch_by_genomic_align_id {
  my ($self, $genomic_align_id) = @_;
  my $genomic_align_group;

  my $genomic_align_adaptor = $self->db->get_GenomicAlignAdaptor;
  
  # Get Bio::EnsEMBL::Compara::GenomicAlign object
  my $genomic_align = $genomic_align_adaptor->fetch_by_dbID($genomic_align_id);
  unless($genomic_align && ref $genomic_align && 
        $genomic_align->isa('Bio::EnsEMBL::Compara::GenomicAlign')) {
    throw("[$genomic_align] must be a Bio::EnsEMBL::Compara::GenomicAlign object");
  }

  return $self->fetch_by_GenomicAlign($genomic_align);
}



=head2 use_autoincrement

  [Arg  1]   : (optional)int value
  Example    : $genomic_align_adaptor->use_autoincrement(0);
  Description: Getter/setter for the _use_autoincrement flag. This flag
               is used when storing new objects with no dbID in the
               database. If the flag is ON (default), the adaptor will
               let the DB set the dbID using the AUTO_INCREMENT ability.
               If you unset the flag, then the adaptor will look for the
               first available dbID after 10^10 times the
               method_link_species_set_id.
  Returntype : integer
  Exceptions : 
  Caller     : none
  Status     : Stable

=cut

sub use_autoincrement {
  my ($self, $value) = @_;

  if (defined $value) {
    $self->{_use_autoincrement} = $value;
  }

  return $self->{_use_autoincrement};
}

1;
