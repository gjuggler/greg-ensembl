
#
# BioPerl module for Bio::EnsEMBL::DBSQL::BaseAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::BaseAdaptor - Base Adaptor for DBSQL adaptors

=head1 SYNOPSIS

    # base adaptor provides
    
    # SQL prepare function
    $adaptor->prepare("sql statement");

    # get of root DBAdaptor object
    $adaptor->db();

    # constructor, ok for inheritence
    $adaptor = Bio::EnsEMBL::DBSQL::SubClassOfBaseAdaptor->new($dbobj)

=head1 DESCRIPTION

This is a true base class for Adaptors in the Ensembl DBSQL
system. Original idea from Arne


Adaptors are expected to have the following functions

    $obj = $adaptor->fetch_by_dbID($internal_id);

which builds the object from the primary key of the object. This
function is crucial because it allows adaptors to collaborate
relatively independently of each other - in other words, we can change
the schema under one adaptor without too many knock on changes through
the other adaptors.

Most adaptors will also have

    $dbid = $adaptor->store($obj);

which stores the object. Currently the storing of an object also causes
the objects to set

    $obj->dbID

correctly and attach the adaptor.


Other fetch functions go by the convention of

    @object_array = @{$adaptor->fetch_all_by_XXXX($arguments_for_XXXX)};

sometimes it returns an array ref denoted by the 'all' in the name of the
method, sometimes an individual object. For example

    $gene = $gene_adaptor->fetch_by_stable_id($stable_id);

or

    @fp  = @{$simple_feature_adaptor->fetch_all_by_Slice($slice)};


Occassionally adaptors need to provide access to lists of ids. In this case the
convention is to go list_XXXX, such as

    @gene_ids = @{$gene_adaptor->list_geneIds()};

(note: this method is poorly named)

=head1 CONTACT

Post questions to the EnsEMBL developer mailing list: <ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::DBSQL::BaseAdaptor;
require Exporter;
use vars qw(@ISA @EXPORT);
use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use DBI qw(:sql_types);

@ISA = qw(Exporter);
@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});

=head2 new

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBConnection $dbobj
  Example    : $adaptor = new AdaptorInheritedFromBaseAdaptor($dbobj);
  Description: Creates a new BaseAdaptor object.  The intent is that this
               constructor would be called by an inherited superclass either
               automatically or through $self->SUPER::new in an overridden 
               new method.
  Returntype : Bio::EnsEMBL::DBSQL::BaseAdaptor
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::DBConnection
  Status     : Stable

=cut

sub new {
  my ($class,$dbobj) = @_;
  
  my $self = {};
  bless $self,$class;
  
  if( !defined $dbobj || !ref $dbobj ) {
    throw("Don't have a db [$dbobj] for new adaptor");
  }
  if($dbobj->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')){
    $self->db($dbobj);
    $self->dbc($dbobj->dbc);
  }
  elsif( ref($dbobj) =~ /DBAdaptor$/){
    $self->db($dbobj);
    $self->dbc($dbobj->dbc);
  }
  elsif( ref($dbobj) =~ /DBConnection$/){
    $self->dbc($dbobj);    
  }
  else{
    throw("Don't have a DBAdaptor [$dbobj] for new adaptor");
  }
  
  return $self;
}


=head2 prepare

  Arg [1]    : string $string
               a SQL query to be prepared by this adaptors database
  Example    : $sth = $adaptor->prepare("select yadda from blabla")
  Description: provides a DBI statement handle from the adaptor. A convenience
               function so you dont have to write $adaptor->db->prepare all the
               time
  Returntype : DBI::StatementHandle
  Exceptions : none
  Caller     : Adaptors inherited from BaseAdaptor
  Status     : Stable

=cut

sub prepare{
  my ($self,$string) = @_;
 
# uncomment next line to cancel caching on the sql side. Needed for timing comparisons etc 
#  $string =~ s/SELECT/SELECT SQL_NO_CACHE/i;

  return $self->dbc->prepare($string);
}


=head2 db

  Arg [1]    : (optional) Bio::EnsEMBL::DBSQL::DBAdaptor $obj 
               the database this adaptor is using.
  Example    : $db = $adaptor->db();
  Description: Getter/Setter for the DatabaseConnection that this adaptor is 
               using.
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : Adaptors inherited fro BaseAdaptor
  Status     : Stable

=cut

sub db{
  my $self = shift;
  $self->{'db'} = shift if(@_);

  return $self->{'db'};

}

=head2 dbc

  Arg [1]    : (optional) Bio::EnsEMBL::DBSQL::DBConnection $obj 
               the database this adaptor is using.
  Example    : $db = $adaptor->db();
  Description: Getter/Setter for the DatabaseConnection that this adaptor is 
               using.
  Returntype : Bio::EnsEMBL::DBSQL::DBConnection
  Exceptions : none
  Caller     : Adaptors inherited fro BaseAdaptor
  Status     : Stable

=cut

sub dbc{
  my $self = shift;
  $self->{'dbc'} = shift if(@_);

  return $self->{'dbc'};
}


# list primary keys for a particular table
# args are table name and primary key field
# if primary key field is not supplied, tablename_id is assumed
# returns listref of IDs
sub _list_dbIDs {

  my ($self, $table, $pk, $ordered) = @_;
  if (!defined($pk)) {
    $pk = $table . "_id";
  }

  my @out;
  my $sql = "SELECT " . $pk . "  FROM " . $table;
  if(defined($ordered) and $ordered){
    $sql .= " order by seq_region_id, seq_region_start"
  }	
  my $sth = $self->prepare($sql);
  $sth->execute;

  while (my ($id) = $sth->fetchrow) {
    push(@out, $id);
  }

  $sth->finish;

  return \@out;
}


# _straight_join

#   Arg [1]    : (optional) boolean $new_val
#   Example    : $self->_straight_join(1);
#                $self->generic_fetch($constraint);
#                $self->_straight_join(0);
#   Description: PROTECTED Getter/Setter that turns on/off the use of 
#                a straight join in queries.
#   Returntype : boolean
#   Exceptions : none
#   Caller     : general

sub _straight_join {
  my $self = shift;
  if(@_) {
    $self->{'_straight_join'} = shift;
  }

  return $self->{'_straight_join'};
}


=head2 generic_fetch

  Arg [1]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Arg [2]    : (optional) Bio::EnsEMBL::AssemblyMapper $mapper
               A mapper object used to remap features
               as they are retrieved from the database
  Arg [3]    : (optional) Bio::EnsEMBL::Slice $slice
               A slice that features should be remapped to
  Example    : $fts = $a->generic_fetch('contig_id in (1234, 1235)', 'Swall');
  Description: Performs a database fetch and returns feature objects in
               contig coordinates.
  Returntype : listref of Bio::EnsEMBL::SeqFeature in contig coordinates
  Exceptions : none
  Caller     : BaseFeatureAdaptor, ProxyDnaAlignFeatureAdaptor::generic_fetch
  Status     : Stable

=cut

sub generic_fetch {
  my ($self, $constraint, $mapper, $slice) = @_;

  my @tabs = $self->_tables;
  my $columns = join(', ', $self->_columns());

  my $db = $self->db();

  #
  # Construct a left join statement if one was defined, and remove the
  # left-joined table from the table list
  #
  my @left_join_list = $self->_left_join();
  my $left_join_prefix = '';
  my $left_join = '';
  my @tables;
  if(@left_join_list) {
    my %left_join_hash = map { $_->[0] => $_->[1] } @left_join_list;
    while(my $t = shift @tabs) {
      if( exists $left_join_hash{ $t->[0] } ) {
        my $condition = $left_join_hash{ $t->[0] };
        my $syn = $t->[1];
        $left_join .=
          "\n  LEFT JOIN " . $t->[0] . " $syn ON $condition ) ";
        $left_join_prefix .= '(';
      } else {
        push @tables, $t;
      }
    }
  } else {
    @tables = @tabs;
  }

  my $straight_join = '';

  if($self->_straight_join()) {
    $straight_join = "STRAIGHT_JOIN";
  }

  #construct a nice table string like 'table1 t1, table2 t2'
  my $tablenames = join(', ', map({ join(' ', @$_) } @tables));

  my $sql =
      "SELECT $straight_join $columns\n"
    . "FROM $left_join_prefix ($tablenames) $left_join";

  my $default_where = $self->_default_where_clause;
  my $final_clause = $self->_final_clause;

  #append a where clause if it was defined
  if($constraint) {
    $sql .= "\n WHERE $constraint ";
    if($default_where) {
      $sql .= " AND\n       $default_where ";
    }
  } elsif($default_where) {
    $sql .= "\n WHERE $default_where ";
  }

  #append additional clauses which may have been defined
  $sql .= "\n$final_clause";

  # FOR DEBUG:
  # printf(STDERR "SQL:\n%s\n", $sql);

  my $sth = $db->dbc->prepare($sql);
  $sth->execute;
  my $res = $self->_objs_from_sth($sth, $mapper, $slice);
  $sth->finish();
  return $res;
}


=head2 fetch_by_dbID

  Arg [1]    : int $id
               The unique database identifier for the feature to be obtained
  Example    : $feat = $adaptor->fetch_by_dbID(1234));
               $feat = $feat->transform('contig');
  Description: Returns the feature created from the database defined by the
               the id $id.  The feature will be returned in its native
               coordinate system.  That is, the coordinate system in which it
               is stored in the database.  In order to convert it to a
               particular coordinate system use the transfer() or transform()
               method.  If the feature is not found in the database then
               undef is returned instead
  Returntype : Bio::EnsEMBL::Feature or undef
  Exceptions : thrown if $id arg is not provided
               does not exist
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_dbID{
  my ($self,$id) = @_;

  throw("id argument is required") if(!defined $id);

  #construct a constraint like 't1.table1_id = 123'
  my @tabs = $self->_tables;
  my ($name, $syn) = @{$tabs[0]};
  my $constraint = "${syn}.${name}_id = $id";

  #Should only be one
  my ($feat) = @{$self->generic_fetch($constraint)};

  return undef if(!$feat);

  return $feat;
}


=head2 fetch_all_by_dbID_list

  Arg [1]    : listref of ints $id_list
               The unique database identifiers for the features to be obtained
  Example    : @feats = @{$adaptor->fetch_by_dbID_list([1234, 2131, 982]))};
  Description: Returns the features created from the database defined by the
               the ids in contained in the id list $id_list.  The features 
               will be returned in their native coordinate system. That is, 
               the coordinate system in which they are stored in the database.
               In order to convert the features to a particular coordinate 
               system use the transfer() or transform() method.  If none of the
               features are found in the database a reference to an empty 
               list is returned.
  Returntype : listref of Bio::EnsEMBL::Features
  Exceptions : thrown if $id arg is not provided
               does not exist
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_dbID_list {
  my ($self,$id_list_ref) = @_;

  if(!defined($id_list_ref) || ref($id_list_ref) ne 'ARRAY') {
    throw("id_list list reference argument is required");
  }

  return [] if(!@$id_list_ref);

  my @out;
  #construct a constraint like 't1.table1_id = 123'
  my @tabs = $self->_tables;
  my ($name, $syn) = @{$tabs[0]};

  # mysql is faster and we ensure that we do not exceed the max query size by
  # splitting large queries into smaller queries of 200 ids
  my $max_size = 200;
  my @id_list = @$id_list_ref;

  while(@id_list) {
    my @ids;
    if(@id_list > $max_size) {
      @ids = splice(@id_list, 0, $max_size);
    } else {
      @ids = splice(@id_list, 0);
    }

    my $id_str;
    if(@ids > 1)  {
      $id_str = " IN (" . join(',', @ids). ")";
    } else {
      $id_str = " = " . $ids[0];
    }

    my $constraint = "${syn}.${name}_id $id_str";

    push @out, @{$self->generic_fetch($constraint)};
  }

  return \@out;
}

# might not be a good idea, but for convenience
# shouldnt be called on the BIG tables though

sub fetch_all {
  my $self = shift;
  return $self->generic_fetch();
}


#_tables
#
#  Args       : none
#  Example    : $tablename = $self->_table_name()
#  Description: ABSTRACT PROTECTED
#               Subclasses are responsible for implementing this
#               method.  It should list of [tablename, alias] pairs.
#               Additionally the primary table (with the dbID,
#               analysis_id, and score) should be the first table in
#               the list. e.g:
#               ( ['repeat_feature',   'rf'],
#                 ['repeat_consensus', 'rc']);
#               used to obtain features.  
#  Returntype : list of [tablename, alias] pairs
#  Exceptions : thrown if not implemented by subclass
#  Caller     : BaseFeatureAdaptor::generic_fetch
#

sub _tables {
  throw(   "abstract method _tables not defined "
         . "by implementing subclass of BaseAdaptor" );
}


#_columns
#
#  Args       : none
#  Example    : $tablename = $self->_columns()
#  Description: ABSTRACT PROTECTED
#               Subclasses are responsible for implementing this
#               method.  It should return a list of columns to be
#               used for feature creation.
#  Returntype : list of strings
#  Exceptions : thrown if not implemented by subclass
#  Caller     : BaseFeatureAdaptor::generic_fetch
#

sub _columns {
  throw(   "abstract method _columns not defined "
         . "by implementing subclass of BaseAdaptor" );
}


# _default_where_clause
#
#  Arg [1]    : none
#  Example    : none
#  Description: May be overridden to provide an additional where
#               constraint to the SQL query which is generated to
#               fetch feature records.  This constraint is always
#               appended to the end of the generated where clause
#  Returntype : string
#  Exceptions : none
#  Caller     : generic_fetch
#

sub _default_where_clause { return '' }


# _left_join

#  Arg [1]    : none
#  Example    : none
#  Description: Can be overridden by a subclass to specify any left
#               joins which should occur.  The table name specigfied
#               in the join must still be present in the return
#               values of.
#  Returntype : a {'tablename' => 'join condition'} pair
#  Exceptions : none
#  Caller     : general
#

sub _left_join { return () }


#_final_clause

#  Arg [1]    : none
#  Example    : none
#  Description: May be overriden to provide an additional clause
#               to the end of the SQL query used to fetch feature
#               records.  This is useful to add a required ORDER BY
#               clause to the query for example.
#  Returntype : string
#  Exceptions : none
#  Caller     : generic_fetch

sub _final_clause { return '' }


#_objs_from_sth

#  Arg [1]    : DBI::row_hashref $hashref containing key-value pairs 
#               for each of the columns specified by the _columns method
#  Example    : my @feats = $self->_obj_from_hashref
#  Description: ABSTRACT PROTECTED
#               The subclass is responsible for implementing this
#               method.  It should take in a DBI row hash reference
#               and return a list of created features in contig
#               coordinates.
#  Returntype : list of Bio::EnsEMBL::*Features in contig coordinates
#  Exceptions : thrown if not implemented by subclass
#  Caller     : BaseFeatureAdaptor::generic_fetch

sub _objs_from_sth {
  throw(   "abstract method _objs_from_sth not defined "
         . "by implementing subclass of BaseAdaptor" );
}

1;
