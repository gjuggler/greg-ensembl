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

Bio::EnsEMBL::IdMapping::BaseObject - base object for IdMapping objects

=head1 SYNOPSIS

  # this object isn't instantiated directly but rather extended
  use Bio::EnsEMBL::IdMapping::BaseObject;
  our @ISA = qw(Bio::EnsEMBL::IdMapping::BaseObject);

=head1 DESCRIPTION

This is the base object for some of the objects used in the IdMapping
application. An object that extends BaseObject will have a ConfParser,
Logger and Cache object. BaseObject also implements some useful utility
functions related to file and db access.

This isn't very clean OO design but it's efficient and easy to use...

=head1 METHODS

  new
  get_filehandle
  file_exists
  fetch_value_from_db
  dump_table_to_file
  upload_file_into_table
  logger
  conf
  cache

=cut


package Bio::EnsEMBL::IdMapping::BaseObject;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);


=head2 new

  Arg [LOGGER]: Bio::EnsEMBL::Utils::Logger $logger - a logger object
  Arg [CONF]  : Bio::EnsEMBL::Utils::ConfParser $conf - a configuration object
  Arg [CACHE] : Bio::EnsEMBL::IdMapping::Cache $cache - a cache object
  Example     : my $object = Bio::EnsEMBL::IdMapping::BaseObjectSubclass->new(
                  -LOGGER => $logger,
                  -CONF   => $conf,
                  -CACHE  => $cache
                );
  Description : Constructor
  Return type : implementing subclass type
  Exceptions  : thrown on wrong or missing arguments
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($logger, $conf, $cache) = rearrange(['LOGGER', 'CONF', 'CACHE'], @_);

  unless ($logger and ref($logger) and
          $logger->isa('Bio::EnsEMBL::Utils::Logger')) {
    throw("You must provide a Bio::EnsEMBL::Utils::Logger for logging.");
  }
  
  unless ($conf and ref($conf) and
          $conf->isa('Bio::EnsEMBL::Utils::ConfParser')) {
    throw("You must provide configuration as a Bio::EnsEMBL::Utils::ConfParser object.");
  }
  
  unless ($cache and ref($cache) and
          $cache->isa('Bio::EnsEMBL::IdMapping::Cache')) {
    throw("You must provide configuration as a Bio::EnsEMBL::IdMapping::Cache object.");
  }
  
  my $self = {};
  bless ($self, $class);

  # initialise
  $self->logger($logger);
  $self->conf($conf);
  $self->cache($cache);
  
  return $self;
}


=head2 get_filehandle 

  Arg[1]      : String $filename - filename for filehandle
  Arg[2]      : String $path_append - append subdirectory name to basedir
  Arg[3]      : String $mode - filehandle mode (<|>|>>)
  Example     : my $fh = $object->get_filehandle('mapping_stats.txt', 'stats',
                  '>');
                print $fh "Stats:\n";
  Description : Returns a filehandle to a file for reading or writing. The file
                is qualified with the basedir defined in the configuration and
                an optional subdirectory name.
  Return type : filehandle
  Exceptions  : thrown on missing filename
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_filehandle {
  my $self = shift;
  my $filename = shift;
  my $path_append = shift;
  my $mode = shift;

  throw("Need a filename for this filehandle.") unless (defined($filename));
  
  my $path = $self->conf->param('basedir');
  $path = path_append($path, $path_append) if (defined($path_append));

  $mode ||= '>';
  
  open(my $fh, $mode, "$path/$filename") or
    throw("Unable to open $path/$filename: $!");

  return $fh;
}


=head2 file_exists

  Arg[1]      : String $filename - filename to test
  Arg[2]      : Boolean $path_append - turn on pre-pending of basedir
  Example     : unless ($object->file_exists('gene_mappings.ser', 1)) {
                  $object->do_gene_mapping;
                }
  Description : Tests if a file exists and has non-zero size.
  Return type : Boolean
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub file_exists {
  my $self = shift;
  my $filename = shift;
  my $path_append = shift;

  my $path = $self->conf->param('basedir');
  $path = path_append($path, $path_append) if (defined($path_append));

  return (-s "$path/$filename");
}


=head2 fetch_value_from_db 

  Arg[1]      : DBI::db $dbh - a DBI database handle
  Arg[2]      : String $sql - SQL statement to execute
  Example     : my $num_genes = $object->fetch_value_from_db($dbh,
                  'SELECT count(*) FROM gene');
  Description : Executes an SQL statement on a db handle and returns the first
                column of the first row returned. Useful for queries returning a
                single value, like table counts.
  Return type : Return type of SQL statement
  Exceptions  : thrown on wrong or missing arguments
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub fetch_value_from_db {
  my $self = shift;
  my $dbh = shift;
  my $sql = shift;

  throw("Need a db handle.") unless ($dbh and $dbh->isa('DBI::db'));
  throw("Need an SQL query to execute.") unless ($sql);

  my $sth = $dbh->prepare($sql);
  $sth->execute;
  my ($retval) = $sth->fetchrow_array;

  return $retval;
}


=head2 dump_table_to_file 

  Arg[1]      : String $dbtype - db type (source|target)
  Arg[2]      : String $table - name of table to dump
  Arg[3]      : String $filename - name of dump file
  Arg[4]      : Boolean $check_existing - turn on test for existing dump
  Example     : my $rows_dumped = $object->dump_table_to_file('source',
                  'stable_id_event', 'stable_id_event_existing.txt');
  Description : Dumps the contents of a db table to a tab-delimited file. The
                dump file will be written to a subdirectory called 'tables'
                under the basedir from your configuration.
  Return type : Int - the number of rows dumped
  Exceptions  : thrown on wrong or missing arguments
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub dump_table_to_file {
  my $self = shift;
  my $dbtype = shift;
  my $table = shift;
  my $filename = shift;
  my $check_existing = shift;

  # argument check
  unless (($dbtype eq 'source') or ($dbtype eq 'target')) {
    throw("Missing or unknown db type: $dbtype.");
  }
  throw("Need a table name.") unless ($table);
  throw("Need a filename.") unless ($filename);

  # conditionally check if table was already dumped
  if ($check_existing and $self->file_exists($filename, 'tables')) {
    $self->logger->info("$filename exists, won't dump again.\n");
    return 0;
  }
  
  my $fh = $self->get_filehandle($filename, 'tables');

  my $dba = $self->cache->get_DBAdaptor($dbtype);
  my $dbh = $dba->dbc->db_handle;
  my $sth = $dbh->prepare("SELECT * FROM $table");
  $sth->execute;

  my $i = 0;

  while (my @row = $sth->fetchrow_array) {
    $i++;

    # use '\N' for NULL values
    for (my $j = 0; $j < scalar(@row); $j++) {
      $row[$j] = '\N' unless (defined($row[$j]));
    }
    
    print $fh join("\t", @row);
    print $fh "\n";
  }

  $sth->finish;
  
  return $i;
}


=head2 upload_file_into_table

  Arg[1]      : String $dbtype - db type (source|target)
  Arg[2]      : String $table - name of table to upload the data to
  Arg[3]      : String $filename - name of dump file
  Arg[4]      : Boolean $no_check_empty - don't check if table is empty
  Example     : my $rows_uploaded = $object->upload_file_into_table('target',
                  'stable_id_event', 'stable_id_event_new.txt');
  Description : Uploads a tab-delimited data file into a db table. The data file
                will be taken from a subdirectory 'tables' under your configured
                basedir. If the db table isn't empty and $no_check_empty isn't
                set, no data is uploaded (and a warning is issued).
  Return type : Int - the number of rows uploaded
  Exceptions  : thrown on wrong or missing arguments
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub upload_file_into_table {
  my $self = shift;
  my $dbtype = shift;
  my $table = shift;
  my $filename = shift;
  my $no_check_empty = shift;

  # argument check
  unless (($dbtype eq 'source') or ($dbtype eq 'target')) {
    throw("Missing or unknown db type: $dbtype.");
  }
  throw("Need a table name.") unless ($table);
  throw("Need a filename.") unless ($filename);

  # sanity check for dry run
  if ($self->conf->param('dry_run')) {
    $self->logger->warning("dry_run - skipping db upload for $filename.\n");
    return;
  }
  
  my $file = join('/', $self->conf->param('basedir'), 'tables', $filename);
  my $r = 0;
  
  if (-s $file) {

    $self->logger->debug("$file -> $table\n", 1);
    
    my $dba = $self->cache->get_DBAdaptor($dbtype);
    my $dbh = $dba->dbc->db_handle;

    # check table is empty
    my ($sql, $sth);
    unless ($no_check_empty) {
      $sql = qq(SELECT count(*) FROM $table);
      $sth = $dbh->prepare($sql);
      $sth->execute;
      my ($c) = $sth->fetchrow_array;
      $sth->finish;

      if ($c) {
        $self->logger->warning("Table $table not empty: found $c entries.\n", 1);
        $self->logger->info("Data not uploaded!\n", 1);
        return $r;
      }
    }
    
    # now upload the data
    $sql = qq(LOAD DATA LOCAL INFILE '$file' INTO TABLE $table);
    $sth = $dbh->prepare($sql);
    $r = $sth->execute;
    $sth->finish;

  } else {
    $self->logger->warning("No data found in file $filename.\n", 1);
  }

  return $r;
}


=head2 logger

  Arg[1]      : (optional) Bio::EnsEMBL::Utils::Logger - the logger to set
  Example     : $object->logger->info("Starting ID mapping.\n");
  Description : Getter/setter for logger object
  Return type : Bio::EnsEMBL::Utils::Logger
  Exceptions  : none
  Caller      : constructor
  Status      : At Risk
              : under development

=cut

sub logger {
  my $self = shift;
  $self->{'_logger'} = shift if (@_);
  return $self->{'_logger'};
}


=head2 conf

  Arg[1]      : (optional) Bio::EnsEMBL::Utils::ConfParser - the configuration
                to set
  Example     : my $basedir = $object->conf->param('basedir');
  Description : Getter/setter for configuration object
  Return type : Bio::EnsEMBL::Utils::ConfParser
  Exceptions  : none
  Caller      : constructor
  Status      : At Risk
              : under development

=cut

sub conf {
  my $self = shift;
  $self->{'_conf'} = shift if (@_);
  return $self->{'_conf'};
}


=head2 cache

  Arg[1]      : (optional) Bio::EnsEMBL::IdMapping::Cache - the cache to set
  Example     : $object->cache->read_from_file('source');
  Description : Getter/setter for cache object
  Return type : Bio::EnsEMBL::IdMapping::Cache
  Exceptions  : none
  Caller      : constructor
  Status      : At Risk
              : under development

=cut

sub cache {
  my $self = shift;
  $self->{'_cache'} = shift if (@_);
  return $self->{'_cache'};
}


1;

