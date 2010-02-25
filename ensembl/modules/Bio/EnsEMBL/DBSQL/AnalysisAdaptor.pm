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

=head1 NAME

Bio::EnsEMBL::DBSQL::AnalysisAdaptor

=head1 SYNOPSIS

  use Bio::EnsEMBL::Registry;

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
  );

  $analysis_adaptor =
    Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "analysis" );

  my $analysis = $analysis_adaptor->fetch_by_logic_name('genscan');

=head1 DESCRIPTION

  Module to encapsulate all db access for persistent class Analysis.
  There should be just one per application and database connection.

=head1 METHODS

=cut


package Bio::EnsEMBL::DBSQL::AnalysisAdaptor;

use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception;

use vars qw(@ISA);
use strict;

@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 new

  Args       : Bio::EnsEMBL::DBSQL::DBAdaptor
  Example    : my $aa = new Bio::EnsEMBL::DBSQL::AnalysisAdaptor();
  Description: Creates a new Bio::EnsEMBL::DBSQL::AnalysisAdaptor object and
               internally loads and caches all the Analysis objects from the 
               database.
  Returntype : Bio::EnsEMBL::DBSQL::AnalysisAdaptor
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::DBAdaptor
  Status     : Stable

=cut

sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  #load and cache all of the Analysis objects
  $self->fetch_all;

  return $self;
}


=head2 fetch_all

  Args       : none
  Example    : my @analysis = @{$analysis_adaptor->fetch_all()};
  Description: fetches all of the Analysis objects from the database and caches
               them internally.
  Returntype : listref of Bio::EnsEMBL::Analysis retrieved from the database
  Exceptions : none
  Caller     : AnalysisAdaptor::new
  Status     : Stable

=cut

sub fetch_all {
  my $self = shift;
  my ( $analysis, $dbID );
  my $rowHashRef;

  $self->{_cache} = {};
  $self->{_logic_name_cache} = {};

  my $sth = $self->prepare( q {
    SELECT analysis.analysis_id, logic_name,
           program, program_version, program_file,
           db, db_version, db_file,
           module, module_version,
           gff_source, gff_feature,
           created, parameters, description, display_label, displayable, web_data
    FROM   analysis
    LEFT JOIN analysis_description
    ON analysis.analysis_id = analysis_description.analysis_id } );
  $sth->execute;

  while( $rowHashRef = $sth->fetchrow_hashref ) {
    my $analysis = $self->_objFromHashref( $rowHashRef  );

    $self->{_cache}->{$analysis->dbID}                    = $analysis;
    $self->{_logic_name_cache}->{lc($analysis->logic_name())} = $analysis;
  }

  my @ana = values %{$self->{_cache}};

  return \@ana;
}


=head2 fetch_all_by_feature_class

  Arg [1]    : string $feature_cless - The name of the feature class
  Example    : my @analyses = @{$analysis_adaptor->fetch_all_by_feature_class('Gene');
  Description: Returns all analyses that correspond to a given 
               feature class; see feature_classes method for a list.
  Returntype : Listref of Bio::EnsEMBL::Analysis
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_feature_class {
  my $self = shift;
  my $feat_class = shift || throw( "Need a feature type, e.g. SimpleFeature" );
 
  my @feature_classes = $self->feature_classes; # List of all feature classes
  my %feat_table_map;
  foreach my $class( @feature_classes ){
    # Map e.g. DnaAlignFeature to dna_align_feature
    my $table = join( "_", map lc, ( $class =~ /([A-Z][a-z]+)/g ) );
    $feat_table_map{$class} = $table;
  }
  $feat_table_map{DensityFeature}='density_type'; # analysis_id in diff table
  my $feat_table = $feat_table_map{$feat_class} || 
      ( warning( "No feature type corresponding to $feat_class" ) &&
        return [] );

  my $sql_t = qq|
SELECT DISTINCT analysis_id FROM %s |;
  
  my $sql = sprintf( $sql_t, $feat_table );
  my $sth = $self->prepare( $sql );
  my $rv  = $sth->execute();
  my $res = $sth->fetchall_arrayref;
  my @analyses;
  foreach my $r( @{$res} ){
    my $analysis = $self->fetch_by_dbID($r->[0]) 
        || throw( "analysis_id $r->[0] from $feat_table table "
                  . "is not in the analysis table!" );
    push @analyses, $analysis;
  }
  return [@analyses];
}


=head2 feature_classes

  Arg [1]    : NONE
  Example    : my @fclasses = $analysis_adaptor->feature_classes;
  Description: Returns a list of the different classes of Ensembl feature 
               object that have an analysis
  Returntype : List of feature classes
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub feature_classes{
  # Can't think of a way to do this programatically, so hard-coded
  return qw(
            DensityFeature
            DnaAlignFeature
            Gene
            MarkerFeature
            PredictionTranscript
            ProteinAlignFeature
            ProteinFeature
            QtlFeature
            RepeatFeature
            SimpleFeature
            );
}

=head2 fetch_by_dbID

  Arg [1]    : int $internal_analysis_id - the database id of the analysis 
               record to retrieve
  Example    : my $analysis = $analysis_adaptor->fetch_by_dbID(1);
  Description: Retrieves an Analysis object from the database via its internal
               id.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $id = shift;

  if( defined $self->{_cache}->{$id} ) {
    return $self->{_cache}->{$id};
  }

  my $query = q{
    SELECT analysis.analysis_id, logic_name,
           program, program_version, program_file,
           db, db_version, db_file,
           module, module_version,
           gff_source, gff_feature,
           created, parameters, description, display_label, displayable, web_data
    FROM   analysis
    LEFT JOIN analysis_description
    ON analysis.analysis_id = analysis_description.analysis_id
    WHERE  analysis.analysis_id = ? };

  my $sth = $self->prepare($query);
  $sth->bind_param(1,$id,SQL_INTEGER);
  $sth->execute();
  my $rowHashRef = $sth->fetchrow_hashref;
  if( ! defined $rowHashRef ) {
    return undef;
  }

  my $anal = $self->_objFromHashref( $rowHashRef );
  $self->{_cache}->{$anal->dbID} = $anal;
  $self->{_logic_name_cache}->{lc($anal->logic_name())} = $anal;
  return $anal;
}


=head2 fetch_by_logic_name

  Arg [1]    : string $logic_name the logic name of the analysis to retrieve
  Example    : my $analysis = $a_adaptor->fetch_by_logic_name('Eponine');
  Description: Retrieves an analysis object from the database using its unique
               logic name.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_logic_name {
  my ( $self, $logic_name ) = @_;

  my $analysis;
  my $rowHash;

  # Check the cache for the logic name
  if ( defined( $self->{_logic_name_cache}{ lc($logic_name) } ) ) {
    return $self->{_logic_name_cache}{ lc($logic_name) };
  }

  my $sth = $self->prepare(
    qq(
SELECT  analysis.analysis_id,
        logic_name,
        program,
        program_version,
        program_file,
        db,
        db_version,
        db_file,
        module,
        module_version,
        gff_source,
        gff_feature,
        created,
        parameters,
        description,
        display_label,
        displayable,
        web_data
FROM    analysis
  LEFT JOIN analysis_description
    ON  ( analysis.analysis_id = analysis_description.analysis_id )
WHERE  LOWER(logic_name) = ?)
  );

  $sth->bind_param( 1, lc($logic_name), SQL_VARCHAR );
  $sth->execute();
  my $rowHashRef = $sth->fetchrow_hashref();

  if ( !defined($rowHashRef) ) { return undef }

  $analysis = $self->_objFromHashref($rowHashRef);

  # place the analysis in the caches, cross referenced by dbID and
  # logic_name
  $self->{_cache}->{ $analysis->dbID() } = $analysis;
  $self->{_logic_name_cache}->{ lc($logic_name) } = $analysis;

  return $analysis;
} ## end sub fetch_by_logic_name


=head2 store

  Arg [1]    : Bio:EnsEMBL::Analysis $analysis
  Example    : $analysis_adaptor->store($analysis);
  Description: Stores $analysis in db.  If the analysis is already stored in
               the database its dbID and adaptor are updated, but the analysis
               is not stored a second time.
               Sets created date if not already set. Sets dbID and adaptor
               inside $analysis. Returns dbID.
  Returntype : int - dbID of stored analysis
  Exceptions : throw on incorrect argument
               throw if analysis argument does not have a logic name
  Caller     : general
  Status     : Stable

=cut

sub store {
  my $self = shift;
  my $analysis = shift;

  if(!ref($analysis) || !$analysis->isa('Bio::EnsEMBL::Analysis')) {
    throw("Bio::EnsEMBL::Analysis argument expected.");
  }

  if($analysis->is_stored($self->db())) {
    return $analysis->dbID();
  }

  if(!$analysis->logic_name()) {
    throw("Analysis cannot be stored without a valid logic_name");
  }
    

  my $rows_inserted = 0;
  my $sth;

  if ( $analysis->created() ) {

    # We use insert IGNORE so that this method can be used in a
    # multi-process environment.  If another process has already written
    # this record then there will not be a problem.

    $sth = $self->prepare(
      q{
      INSERT IGNORE INTO analysis
      SET created = ?,
          logic_name = ?,
          db = ?,
          db_version = ?,
          db_file = ?,
          program = ?,
          program_version = ?,
          program_file = ?,
          parameters = ?,
          module = ?,
          module_version = ?,
          gff_source = ?,
          gff_feature = ?
      }
    );
    $sth->bind_param( 1,  $analysis->created,         SQL_DATETIME );
    $sth->bind_param( 2,  $analysis->logic_name,      SQL_VARCHAR );
    $sth->bind_param( 3,  $analysis->db,              SQL_VARCHAR );
    $sth->bind_param( 4,  $analysis->db_version,      SQL_VARCHAR );
    $sth->bind_param( 5,  $analysis->db_file,         SQL_VARCHAR );
    $sth->bind_param( 6,  $analysis->program,         SQL_VARCHAR );
    $sth->bind_param( 7,  $analysis->program_version, SQL_VARCHAR );
    $sth->bind_param( 8,  $analysis->program_file,    SQL_VARCHAR );
    $sth->bind_param( 9,  $analysis->parameters,      SQL_VARCHAR );
    $sth->bind_param( 10, $analysis->module,          SQL_VARCHAR );
    $sth->bind_param( 11, $analysis->module_version,  SQL_VARCHAR );
    $sth->bind_param( 12, $analysis->gff_source,      SQL_VARCHAR );
    $sth->bind_param( 13, $analysis->gff_feature,     SQL_VARCHAR );

    $rows_inserted = $sth->execute();

  } else {
    $sth = $self->prepare(
      q{
      INSERT IGNORE INTO analysis
      SET created = now(),
          logic_name = ?,
          db = ?,
          db_version = ?,
          db_file = ?,
          program = ?,
          program_version = ?,
          program_file = ?,
          parameters = ?,
          module = ?,
          module_version = ?,
          gff_source = ?,
          gff_feature = ?
       }
    );

    $sth->bind_param( 1,  $analysis->logic_name,      SQL_VARCHAR );
    $sth->bind_param( 2,  $analysis->db,              SQL_VARCHAR );
    $sth->bind_param( 3,  $analysis->db_version,      SQL_VARCHAR );
    $sth->bind_param( 4,  $analysis->db_file,         SQL_VARCHAR );
    $sth->bind_param( 5,  $analysis->program,         SQL_VARCHAR );
    $sth->bind_param( 6,  $analysis->program_version, SQL_VARCHAR );
    $sth->bind_param( 7,  $analysis->program_file,    SQL_VARCHAR );
    $sth->bind_param( 8,  $analysis->parameters,      SQL_VARCHAR );
    $sth->bind_param( 9,  $analysis->module,          SQL_VARCHAR );
    $sth->bind_param( 10, $analysis->module_version,  SQL_VARCHAR );
    $sth->bind_param( 11, $analysis->gff_source,      SQL_VARCHAR );
    $sth->bind_param( 12, $analysis->gff_feature,     SQL_VARCHAR );

    $rows_inserted = $sth->execute();

  } ## end else [ if ( $analysis->created...)]

  my $dbID;
  # If we need to fetch the timestamp, or the insert failed due to
  # existance of an existing entry, we need to retrieve the entry from
  # the database.  Note: $sth->execute() may return 0E0 on error which
  # is zero, but true which is why the $rows_inserted clause was added.
  if ( !$analysis->created() || !$rows_inserted || $rows_inserted == 0 )
  {
    my $new_analysis =
      $self->fetch_by_logic_name( $analysis->logic_name );

    if ( !$new_analysis ) {
      throw("Could not retrieve just stored analysis from database.\n"
          . "Possibly incorrect db permissions or missing analysis table\n"
      );
    }

    $dbID = $new_analysis->dbID();
    $analysis->created( $new_analysis->created() );
  }
  
  $dbID ||= $sth->{'mysql_insertid'};
  $sth->finish();

  # store description and display_label
  if( defined( $analysis->description() ) || defined( $analysis->display_label() )|| defined( $analysis->web_data() )) {
      $sth = $self->prepare( "INSERT IGNORE INTO analysis_description (analysis_id, display_label, description, displayable, web_data) VALUES (?,?,?,?, ?)");

      $sth->bind_param(1,$dbID,SQL_INTEGER);
      $sth->bind_param(2,$analysis->display_label(),SQL_VARCHAR);
      $sth->bind_param(3,$analysis->description,SQL_LONGVARCHAR);
      $sth->bind_param(4,$analysis->displayable,SQL_TINYINT);
      #$sth->bind_param(5,$analysis->web_data(),SQL_LONGVARCHAR);
      my $web_data;
      $web_data = $self->dump_data($analysis->web_data()) if ($analysis->web_data());
      $sth->bind_param(5,$web_data,SQL_LONGVARCHAR);
      $sth->execute();

      $sth->finish();
  }
  


  $self->{_cache}->{$dbID} = $analysis;
  $self->{_logic_name_cache}{lc($analysis->logic_name)} = $analysis;

  $analysis->adaptor( $self );
  $analysis->dbID( $dbID );

  return $dbID;
}



=head2 update

  Arg [1]    : Bio::EnsEMBL::Analysis $anal
  Example    : $adaptor->update($anal)
  Description: Updates this analysis in the database
  Returntype : int 1 if update is performed, undef if it is not
  Exceptions : throw if arg is not an analysis object
  Caller     : ?
  Status     : Stable

=cut

sub update {
  my $self = shift;
  my $a    = shift;

  if (!ref($a) || !$a->isa('Bio::EnsEMBL::Analysis')) {
    throw("Expected Bio::EnsEMBL::Analysis argument.");
  }

  if(!$a->is_stored($self->db())) {
    return undef;
  }

  my $sth = $self->prepare
    ("UPDATE analysis " .
     "SET created = ?, logic_name = ?, db = ?, db_version = ?, db_file = ?, ".
     "    program = ?, program_version = ?, program_file = ?,  ".
     "    parameters = ?, module = ?, module_version = ?, ".
     "    gff_source = ?, gff_feature = ? " .
     "WHERE analysis_id = ?");



  $sth->bind_param(1,$a->created,SQL_DATETIME);
  $sth->bind_param(2,$a->logic_name,SQL_VARCHAR);
  $sth->bind_param(3,$a->db,SQL_VARCHAR);
  $sth->bind_param(4,$a->db_version,SQL_VARCHAR);
  $sth->bind_param(5,$a->db_file,SQL_VARCHAR);
  $sth->bind_param(6,$a->program,SQL_VARCHAR);
  $sth->bind_param(7,$a->program_version,SQL_VARCHAR);
  $sth->bind_param(8,$a->program_file,SQL_VARCHAR);
  $sth->bind_param(9,$a->parameters,SQL_VARCHAR);
  $sth->bind_param(10,$a->module,SQL_VARCHAR);
  $sth->bind_param(11,$a->module_version,SQL_VARCHAR);
  $sth->bind_param(12,$a->gff_source,SQL_VARCHAR);
  $sth->bind_param(13,$a->gff_feature,SQL_VARCHAR);
  $sth->bind_param(14,$a->dbID,SQL_INTEGER);

  $sth->execute();

  $sth->finish();

  # also update description & display label - may need to create these if
  # not already there
  $sth = $self->prepare("SELECT description FROM analysis_description WHERE analysis_id= ?");
  $sth->execute($a->dbID);
  my $web_data; #this is an anonymous reference to a hash, will have to be dumped into string before writing to db
  if ($sth->fetchrow_hashref) { # update if exists
      $web_data = $self->dump_data($a->web_data()) if ($a->web_data());
      $sth = $self->prepare
      ("UPDATE analysis_description SET description = ?, display_label = ?, displayable = ?, web_data = ? WHERE analysis_id = ?");
      $sth->bind_param(1,$a->description,SQL_LONGVARCHAR);     
      $sth->bind_param(2,$a->display_label(),SQL_VARCHAR);
      $sth->bind_param(3,$a->displayable,SQL_TINYINT);
      #      print "after $web_data\n";
      $sth->bind_param(4,$web_data,SQL_LONGVARCHAR);
      $sth->bind_param(5,$a->dbID,SQL_INTEGER);
      $sth->execute();

  } else { # create new entry

    if( $a->description() || $a->display_label() || $a->web_data) {
	$web_data = $self->dump_data($a->web_data()) if ($a->web_data());
      #my $web_data = $self->dump_data($a->web_data());
      $sth = $self->prepare( "INSERT IGNORE INTO analysis_description (analysis_id, display_label, description, displayable, web_data) VALUES (?,?,?,?,?)");
	$sth->bind_param(1,$a->dbID,SQL_INTEGER);	
	$sth->bind_param(2,$a->display_label(),SQL_VARCHAR);
	$sth->bind_param(3,$a->description,SQL_LONGVARCHAR);     
	$sth->bind_param(4,$a->displayable,SQL_TINYINT);
	#my $web_data = $self->dump_data($a->web_data());
	$sth->bind_param(5,$web_data,SQL_LONGVARCHAR);
	$sth->execute();

    }

  }


    $sth->finish();

  # the logic_name cache needs to be re-updated now, since we may have just
  # changed the logic_name
  $self->fetch_all();

  return 1;
}



=head2 remove

  Arg [1]    : Bio::EnsEMBL::Analysis $anal
  Example    : $adaptor->remove($anal)
  Description: Removes this analysis from the database.  This is not really
               safe to execute in a multi process environment, so programs
               should not remove analysis while out on the farm.
  Returntype : none
  Exceptions : thrown if $anal arg is not an analysis object
  Caller     : ?
  Status     : Stable

=cut

sub remove {
  my ($self, $analysis) = @_;

  if (!defined $analysis || !ref $analysis) {
    throw("called remove on AnalysisAdaptor with a [$analysis]");
  }

  if(!$analysis->is_stored($self->db())) {
    return undef;
  }

  my $sth = $self->prepare("DELETE FROM analysis WHERE analysis_id = ?");
  $sth->bind_param(1,$analysis->dbID,SQL_INTEGER);
  $sth->execute();

  $sth = $self->prepare("DELETE FROM analysis_description WHERE analysis_id = ?");
  $sth->execute($analysis->dbID());

  # remove this analysis from the cache
  delete $self->{'_cache'}->{$analysis->dbID()};
  delete $self->{'_logic_name_cache'}->{lc($analysis->logic_name)};


  # unset the adaptor and dbID
  $analysis->dbID(undef);
  $analysis->adaptor(undef);

  return;
}



=head2 exists

  Arg [1]    : Bio::EnsEMBL::Analysis $anal
  Example    : if($adaptor->exists($anal)) #do something
  Description: Tests whether this Analysis already exists in the database
               by checking first if the adaptor and dbID are set and
               secondly by whether it is in this adaptors internal cache.
               Note that this will not actually check the database and will
               not find and analysis which were recently added by other
               processes.  You are better off simply trying to store an
               analysis which will reliably ensure that it is not stored twice
               in the database.
  Returntype : int dbID if analysis is found, otherwise returns undef
  Exceptions : thrown if $anal arg is not an analysis object
  Caller     : store
  Status     : Stable

=cut

sub exists {
  my ($self,$anal) = @_;

  if(!ref($anal) || !$anal->isa("Bio::EnsEMBL::Analysis")) {
    throw("Object is not a Bio::EnsEMBL::Analysis");
  }

  #if this analysis is stored in this db already return its dbID
  if($anal->is_stored($self->db())) {
    return $anal->dbID();
  }

  #this analysis object is not stored but one exactly like it may have been
  foreach my $cacheId (keys %{$self->{_cache}}) {
    if ($self->{_cache}->{$cacheId}->compare($anal) >= 0) {
      # $anal->dbID( $cacheId );
      # $anal->adaptor( $self );
      return $cacheId;
    }
  }

  #no analysis like this one exists in the database
  return undef;
}


=head2 _objFromHashref

  Arg [1]    : hashref $rowHash
  Description: Private helper function generates an Analysis object from a 
               mysql row hash reference.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::AnalsisAdaptor::fetch_* methods
  Status     : Stable

=cut

sub _objFromHashref {
  my $self = shift;
  my $rowHash = shift;

  my $web_data = $rowHash->{web_data} ? $self->get_dumped_data($rowHash->{web_data}) : '';
  my $analysis = Bio::EnsEMBL::Analysis->new(
      -id              => $rowHash->{analysis_id},
      -adaptor         => $self,
      -db              => $rowHash->{db},
      -db_file         => $rowHash->{db_file},
      -db_version      => $rowHash->{db_version},
      -program         => $rowHash->{program},
      -program_version => $rowHash->{program_version},
      -program_file    => $rowHash->{program_file},
      -gff_source      => $rowHash->{gff_source},
      -gff_feature     => $rowHash->{gff_feature},
      -module          => $rowHash->{module},
      -module_version  => $rowHash->{module_version},
      -parameters      => $rowHash->{parameters},
      -created         => $rowHash->{created},
      -logic_name      => $rowHash->{logic_name},
      -description     => $rowHash->{description},
      -display_label   => $rowHash->{display_label},
      -displayable     => $rowHash->{displayable},
					     -web_data        => $web_data,
    );

  return $analysis;
}



1;
