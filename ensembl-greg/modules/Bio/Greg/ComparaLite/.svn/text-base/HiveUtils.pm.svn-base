package Bio::Greg::ComparaLite::HiveUtils;

use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

# Define some constants.
my $DBA = "Bio::EnsEMBL::Compara::DBSQL::DBAdaptor";
my $HIVE_SQL;


sub hive_url_to_mysql_args {
  my $class = shift;
  my $url = shift;

  $url =~ m!mysql://(.*?):(.*?)@(.*?)/(.*)!;
  
  # Something like mysql -uensadmin -pensembl -hensdb-2-12 -P5106 -A gj1_compara_53
  # $1 = user, $2 = pw, $3 = server+port, $4 = db
  my ($u,$p,$h,$P,$db) = undef;
  $u = $1;
  $p = $2;
  $h = $3;
  $P = 3306;
  $db = $4;
  if ($h =~ m/(.*?):(.*)/) {
    $h = $1;
    $P = $2;
  }

  my $mysqlargs = sprintf("-u%s -p%s -h%s -P%s -A %s",
			  $u,
			  $p,
			  $h,
			  $P,
			  $db);
  return $mysqlargs;
}

sub hive_url_to_hashref {
  my $class = shift;
  my $url = shift;

  my $args = $class->hive_url_to_mysql_args($url);

  $args =~ m/-u(.*) -p(.*) -h(.*) -P(.*) -A (.*)/;
  my $data = {user => $1,
	      password => $2,
	      host => $3,
	      port => $4,
	      database => $5};
  return $data;
}

# Creates a database $database in the given mysql server at $url.
# Returns the new full URL for accessing the database.
sub hive_create_database {
  my $class = shift;
  my $url = shift;
  my $format_database = shift;

  $format_database = 0 unless (defined $format_database);

  my $obj = $class->hive_url_to_hashref($url);
  $url = sprintf("mysql://%s:%s@%s/",
		 $obj->{'user'},
		 $obj->{'password'},
		 $obj->{'host'}.":".$obj->{'port'}
		 );
  my $database = $obj->{'database'};
  #print "BASE URL: $url\n";
  my $temp_url = $url . "mysql";

  my $dba = $DBA->new(-url => $temp_url);
  eval {
    $dba->dbc->do("drop database if exists $database;") if ($format_database);
  };
  eval {
    $dba->dbc->do("create database if not exists $database;");
  };

}

# Creates the SQL tables necessary for running a Hive.
sub hive_create_hive_tables {
  my $class = shift;
  my $url = shift; 		# URL of the desired MySQL location INCLUDING the database name.
  my $format_database = shift;	# Set to 1 if we should wipe the database clean.

  $format_database = 1 unless (defined $format_database);

  my $dba = $DBA->new(-url => $url);
  # Create the tables.

  open(IN,"/nfs/users/nfs_g/gj1/src/ensembl_svn/ensembl-greg/sql/hive-basic-tables.sql");
  $HIVE_SQL = join("",<IN>);
  close(IN);

  $class->run_sql_commands($dba,$HIVE_SQL);
  # Format if necessary.
  if ($format_database) {
    print "Creating hive database: $url \n";
    my @hive_tables = qw(hive analysis analysis_stats analysis_stats_monitor analysis_ctrl_rule dataflow_rule
			 analysis_job analysis_job_file analysis_data analysis_description monitor
			 );
    map {eval{ $dba->dbc->do("truncate $_;")}} @hive_tables; # Truncate each table in turn.
  }
}

# Creates the SQL tables necessary for housing a ComparaLite database.
sub hive_create_compara_tables {
  my $class = shift;
  my $url = shift;
  my $format_database = shift;

  $format_database = 1 unless (defined $format_database);

  open(IN,"/nfs/users/nfs_g/gj1/src/ensembl_svn/ensembl-greg/sql/compara-basic-tables.sql");
  $COMPARA_SQL = join("",<IN>);
  close(IN);

  my $dba = $DBA->new(-url => $url);
  # Create the tables.
  $class->run_sql_commands($dba,$COMPARA_SQL);
  # Format if necessary.
  if ($format_database) {
    print "Creating Compara database: $url \n";
    my @hive_tables = qw(
			 member sequence protein_tree_member protein_tree_node sitewise_aln
			 );
    map {eval{ $dba->dbc->do("truncate $_;")}} @hive_tables; # Truncate each table in turn.
  }  
}


BEGIN {
  # Define the Hive and Compara SQL tables.

  $HIVE_SQL = qq^
  CREATE TABLE IF NOT EXISTS hive (
  hive_id          int(10) NOT NULL auto_increment,
  analysis_id      int(10) NOT NULL,
  beekeeper        varchar(80) DEFAULT '' NOT NULL,
  host	           varchar(40) DEFAULT '' NOT NULL,
  process_id       varchar(40) DEFAULT '' NOT NULL,
  work_done        int(11) DEFAULT '0' NOT NULL,
  status           enum('READY','GET_INPUT','RUN','WRITE_OUTPUT','DEAD') DEFAULT 'READY' NOT NULL,
  born	           datetime NOT NULL,
  last_check_in    datetime NOT NULL,
  died             datetime DEFAULT NULL,
  cause_of_death   enum('', 'NO_WORK', 'JOB_LIMIT', 'HIVE_OVERLOAD', 'LIFESPAN', 'FATALITY') DEFAULT '' NOT NULL,
  PRIMARY KEY (hive_id),
  INDEX analysis_status (analysis_id, status)
);

CREATE TABLE IF NOT EXISTS dataflow_rule (
  dataflow_rule_id    int(10) unsigned not null auto_increment,
  from_analysis_id    int(10) unsigned NOT NULL,
  to_analysis_url     varchar(255) default '' NOT NULL,
  branch_code         int(10) default 1 NOT NULL,

  PRIMARY KEY (dataflow_rule_id),
  UNIQUE (from_analysis_id, to_analysis_url)
);

CREATE TABLE IF NOT EXISTS analysis_ctrl_rule (
  condition_analysis_url     varchar(255) default '' NOT NULL,
  ctrled_analysis_id         int(10) unsigned NOT NULL,

  UNIQUE (condition_analysis_url, ctrled_analysis_id)
);

CREATE TABLE IF NOT EXISTS analysis_job (
  analysis_job_id           int(10) NOT NULL auto_increment,
  prev_analysis_job_id      int(10) NOT NULL,  #analysis_job which created this from rules
  analysis_id               int(10) NOT NULL,
  input_id                  char(255) not null,
  job_claim                 char(40) NOT NULL default '', #UUID
  hive_id                   int(10) NOT NULL,
  status                    enum('READY','BLOCKED','CLAIMED','GET_INPUT','RUN','WRITE_OUTPUT','DONE','FAILED') DEFAULT 'READY' NOT NULL,
  retry_count               int(10) default 0 not NULL,
  completed                 datetime NOT NULL,
  branch_code               int(10) default 1 NOT NULL,
  runtime_msec              int(10) default 0 NOT NULL, 
  query_count               int(10) default 0 NOT NULL, 

  PRIMARY KEY                  (analysis_job_id),
  UNIQUE KEY input_id_analysis (input_id, analysis_id),
  INDEX claim_analysis_status  (job_claim, analysis_id, status),
  INDEX analysis_status        (analysis_id, status),
  INDEX hive_id                (hive_id)
);

CREATE TABLE IF NOT EXISTS analysis_job_file (
  analysis_job_id         int(10) NOT NULL,
  hive_id                 int(10) NOT NULL,
  retry                   int(10) NOT NULL,
  type                    varchar(16) NOT NULL default '',
  path                    varchar(255) NOT NULL,
  
  UNIQUE KEY job_hive_type  (analysis_job_id, hive_id, type),
  INDEX hive_id             (hive_id)
);

CREATE TABLE IF NOT EXISTS analysis_data (
  analysis_data_id  int(10) NOT NULL auto_increment,
  data              longtext,

  PRIMARY KEY (analysis_data_id),
  KEY data (data(100))
);

CREATE TABLE IF NOT EXISTS analysis_stats (
  analysis_id           int(10) NOT NULL,
  status                enum('BLOCKED', 'LOADING', 'SYNCHING', 'READY', 'WORKING', 'ALL_CLAIMED', 'DONE', 'FAILED')
                          DEFAULT 'READY' NOT NULL,
  batch_size            int(10) default 1 NOT NULL,
  avg_msec_per_job      int(10) default 0 NOT NULL,
  avg_input_msec_per_job      int(10) default 0 NOT NULL,
  avg_run_msec_per_job      int(10) default 0 NOT NULL,
  avg_output_msec_per_job      int(10) default 0 NOT NULL,
  hive_capacity         int(10) default 1 NOT NULL,
  behaviour             enum('STATIC', 'DYNAMIC') DEFAULT 'STATIC' NOT NULL,
  input_capacity        int(10) default 4 NOT NULL,
  output_capacity       int(10) default 4 NOT NULL,
  total_job_count       int(10) NOT NULL,
  unclaimed_job_count   int(10) NOT NULL,
  done_job_count        int(10) NOT NULL,
  max_retry_count       int(10) default 3 NOT NULL,
  failed_job_count      int(10) NOT NULL,
  failed_job_tolerance  int(10) default 0 NOT NULL,
  num_running_workers   int(10) default 0 NOT NULL,
  num_required_workers  int(10) NOT NULL,
  last_update           datetime NOT NULL,
  sync_lock             int(10) default 0 NOT NULL,
  
  UNIQUE KEY   (analysis_id)
);

CREATE TABLE IF NOT EXISTS analysis_stats_monitor (
  time                  datetime NOT NULL default '0000-00-00 00:00:00',
  analysis_id           int(10) NOT NULL,
  status                enum('BLOCKED', 'LOADING', 'SYNCHING', 'READY', 'WORKING', 'ALL_CLAIMED', 'DONE', 'FAILED')
                          DEFAULT 'READY' NOT NULL,
  batch_size            int(10) default 1 NOT NULL,
  avg_msec_per_job      int(10) default 0 NOT NULL,
  avg_input_msec_per_job      int(10) default 0 NOT NULL,
  avg_run_msec_per_job      int(10) default 0 NOT NULL,
  avg_output_msec_per_job      int(10) default 0 NOT NULL,
  hive_capacity         int(10) default 1 NOT NULL,
  behaviour             enum('STATIC', 'DYNAMIC') DEFAULT 'STATIC' NOT NULL,
  input_capacity        int(10) default 4 NOT NULL,
  output_capacity       int(10) default 4 NOT NULL,
  total_job_count       int(10) NOT NULL,
  unclaimed_job_count   int(10) NOT NULL,
  done_job_count        int(10) NOT NULL,
  max_retry_count       int(10) default 3 NOT NULL,
  failed_job_count      int(10) NOT NULL,
  failed_job_tolerance  int(10) default 0 NOT NULL,
  num_running_workers   int(10) default 0 NOT NULL,
  num_required_workers  int(10) NOT NULL,
  last_update           datetime NOT NULL,
  sync_lock             int(10) default 0 NOT NULL
);

CREATE TABLE IF NOT EXISTS monitor (
  time                  datetime NOT NULL default '0000-00-00 00:00:00',
  workers               int(10) NOT NULL default '0',
  throughput            float default NULL,
  per_worker            float default NULL,
  analysis              varchar(255) default NULL
);

CREATE TABLE IF NOT EXISTS meta (

  meta_id                     INT NOT NULL AUTO_INCREMENT,
  species_id                  INT UNSIGNED DEFAULT 1,
  meta_key                    VARCHAR(40) NOT NULL,
  meta_value                  VARCHAR(255) BINARY NOT NULL,

  PRIMARY   KEY (meta_id),
  UNIQUE    KEY species_key_value_idx (species_id, meta_key, meta_value),
            KEY species_value_idx (species_id, meta_value)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

CREATE TABLE IF NOT EXISTS analysis (

  analysis_id                 SMALLINT UNSIGNED NOT NULL AUTO_INCREMENT,
  created                     datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  logic_name                  VARCHAR(40) NOT NULL,
  db                          VARCHAR(120),
  db_version                  VARCHAR(40),
  db_file                     VARCHAR(120),
  program                     VARCHAR(80),
  program_version             VARCHAR(40),
  program_file                VARCHAR(80),
  parameters                  VARCHAR(255),
  module                      VARCHAR(80),
  module_version              VARCHAR(40),
  gff_source                  VARCHAR(40),
  gff_feature                 VARCHAR(40),

  PRIMARY KEY (analysis_id),
  KEY logic_name_idx (logic_name),
  UNIQUE (logic_name)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

CREATE TABLE IF NOT EXISTS analysis_description (

  analysis_id                  SMALLINT UNSIGNED NOT NULL,
  description                  TEXT,
  display_label                VARCHAR(255),
  displayable                  BOOLEAN NOT NULL DEFAULT 1,
  web_data                     TEXT,

  UNIQUE KEY analysis_idx (analysis_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;
^;


  $COMPARA_SQL = qq^

DROP FUNCTION IF EXISTS dubious_dup;
CREATE FUNCTION dubious_dup (n int(20))
RETURNS int(20)
RETURN (SELECT count(*) FROM protein_tree_tag where node_id=n AND tag="dubious_duplication" AND value="1");

DROP FUNCTION IF EXISTS leaf_count;
CREATE FUNCTION leaf_count (n int(20))
RETURNS int(20)
RETURN (SELECT COUNT(*) 
   FROM protein_tree_node n1,protein_tree_node n2,protein_tree_member m 
   WHERE n2.left_index BETWEEN n1.left_index and n1.right_index
# GJ 2009-03-03 - until we redo the leftright-indexing, keep this in:
   AND n2.root_id=n1.root_id
   AND n1.node_id=n
   AND m.node_id=n2.node_id
);

CREATE TABLE IF NOT EXISTS meta (

  meta_id                     INT NOT NULL AUTO_INCREMENT,
  species_id                  INT UNSIGNED DEFAULT 1,
  meta_key                    VARCHAR(40) NOT NULL,
  meta_value                  VARCHAR(255) BINARY NOT NULL,
	
  PRIMARY   KEY (meta_id),
  UNIQUE    KEY species_key_value_idx (species_id, meta_key, meta_value),
  KEY species_value_idx (species_id, meta_value)

) ENGINE=InnoDB;


CREATE TABLE IF NOT EXISTS ncbi_taxa_node (
  taxon_id                        int(10) unsigned NOT NULL,
  parent_id                       int(10) unsigned NOT NULL,

  rank                            char(32) default '' NOT NULL,
  genbank_hidden_flag             boolean default 0 NOT NULL,

  left_index                      int(10) NOT NULL,
  right_index                     int(10) NOT NULL,
  root_id                         int(10) default 1 NOT NULL,
  
  PRIMARY KEY (taxon_id),
  KEY (parent_id),
  KEY (rank)
) ENGINE=InnoDB COLLATE=latin1_swedish_ci;

CREATE TABLE IF NOT EXISTS ncbi_taxa_name (
  taxon_id                    int(10) unsigned NOT NULL,

  name                        varchar(255),
  name_class                  varchar(50),

  KEY (taxon_id),
  KEY (name),
  KEY (name_class)
) ENGINE=InnoDB COLLATE=latin1_swedish_ci;

CREATE TABLE IF NOT EXISTS genome_db (
  genome_db_id                int(10) unsigned NOT NULL auto_increment, # unique internal id
  taxon_id                    int(10) unsigned DEFAULT '0' NOT NULL, # KF taxon.taxon_id
  name                        varchar(40) DEFAULT '' NOT NULL,
  assembly                    varchar(100) DEFAULT '' NOT NULL,
  assembly_default            tinyint(1) DEFAULT 1,
  genebuild                   varchar(100) DEFAULT '' NOT NULL,
  locator                     varchar(255),

  FOREIGN KEY (taxon_id) REFERENCES ncbi_taxa_node(taxon_id),

  PRIMARY KEY (genome_db_id),
  UNIQUE name (name,assembly,genebuild)
) ENGINE=InnoDB COLLATE=latin1_swedish_ci;

CREATE TABLE IF NOT EXISTS sequence (
  sequence_id                 int(10) unsigned NOT NULL auto_increment, # unique internal id
  length                      int(10) NOT NULL,
  sequence                    longtext NOT NULL,

  PRIMARY KEY (sequence_id),
  KEY sequence (sequence(18))
) ENGINE=InnoDB MAX_ROWS = 1000000 AVG_ROW_LENGTH = 19000 COLLATE=latin1_swedish_ci;


CREATE TABLE IF NOT EXISTS `sequence_quality` (
  `sequence_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `length` int(10) NOT NULL,
  `sequence` longtext NOT NULL,
  PRIMARY KEY (`sequence_id`),
  KEY `sequence` (`sequence`(18))
) ENGINE=InnoDB AUTO_INCREMENT=1114901 DEFAULT CHARSET=latin1 AVG_ROW_LENGTH=19000;


CREATE TABLE IF NOT EXISTS member (
  member_id                   int(10) unsigned NOT NULL auto_increment, # unique internal id
  stable_id                   varchar(40) NOT NULL, # e.g. ENSP000001234 or P31946
  version                     int(10) DEFAULT '0', 
  source_name                 varchar(40) NOT NULL,
#  source_name                 ENUM('ENSEMBLGENE','ENSEMBLPEP','Uniprot/SPTREMBL','Uniprot/SWISSPROT') NOT NULL,
  taxon_id                    int(10) unsigned NOT NULL, # FK taxon.taxon_id
  genome_db_id                int(10) unsigned, # FK genome_db.genome_db_id
  sequence_id                 int(10) unsigned, # FK sequence.sequence_id
  cdna_sequence_id            int(10) unsigned, # FK sequence.sequence_id
  gene_member_id              int(10) unsigned, # FK member.member_id
  description                 text DEFAULT NULL,
  chr_name                    char(40),
  chr_start                   int(10),
  chr_end                     int(10),
  chr_strand                  tinyint(1) NOT NULL,
  display_label               varchar(128) default NULL,

#  FOREIGN KEY (taxon_id) REFERENCES ncbi_taxa_node(taxon_id),
  FOREIGN KEY (genome_db_id) REFERENCES genome_db(genome_db_id),
  FOREIGN KEY (sequence_id) REFERENCES sequence(sequence_id),
  FOREIGN KEY (cdna_sequence_id) REFERENCES sequence(sequence_id),
  FOREIGN KEY (gene_member_id) REFERENCES member(member_id),

  PRIMARY KEY (member_id),
  UNIQUE source_stable_id (stable_id, source_name),
  KEY (stable_id),
  KEY (source_name),
  KEY (sequence_id),
  KEY (gene_member_id)
) ENGINE=InnoDB COLLATE=latin1_swedish_ci;


CREATE TABLE IF NOT EXISTS protein_tree_node (
  node_id                         int(10) unsigned NOT NULL auto_increment, # unique internal id
  parent_id                       int(10) unsigned NOT NULL,
  root_id                         int(10) unsigned NOT NULL,
  left_index                      int(10) NOT NULL,
  right_index                     int(10) NOT NULL,
  distance_to_parent              double default 1.0 NOT NULL,

  PRIMARY KEY (node_id),
  KEY (parent_id),
  KEY (root_id),
  KEY (left_index),
  KEY (right_index)
) ENGINE=InnoDB COLLATE=latin1_swedish_ci;


CREATE TABLE IF NOT EXISTS protein_tree_member (
  node_id                     int(10) unsigned NOT NULL,
  member_id                   int(10) unsigned NOT NULL, 
  method_link_species_set_id  int(10) unsigned NOT NULL,
  cigar_line                  mediumtext,
  cigar_start                 int(10),
  cigar_end                   int(10),

  FOREIGN KEY (node_id) REFERENCES protein_tree_node(node_id),

  UNIQUE (node_id),
  KEY (member_id)
) ENGINE=InnoDB COLLATE=latin1_swedish_ci;

CREATE TABLE IF NOT EXISTS protein_tree_member_score LIKE protein_tree_member;

CREATE TABLE IF NOT EXISTS protein_tree_tag (
  node_id                int(10) unsigned NOT NULL,
  tag                    varchar(50),
  value                  mediumtext,

  FOREIGN KEY (node_id) REFERENCES protein_tree_node(node_id),

  UNIQUE tag_node_id (node_id, tag),
  KEY (node_id),
  KEY (tag)
) COLLATE=latin1_swedish_ci;

CREATE TABLE IF NOT EXISTS `sitewise_aln` (
  `sitewise_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `parameter_set_id` int(10) unsigned NOT NULL,
  `aln_position` int(10) unsigned NOT NULL,
  `node_id` int(10) unsigned NOT NULL,
  `tree_node_id` int(10) unsigned NOT NULL,
  `omega` float(10,5) DEFAULT NULL,
  `omega_lower` float(10,5) DEFAULT NULL,
  `omega_upper` float(10,5) DEFAULT NULL,
  `lrt_stat` float(10,5) DEFAULT NULL,
  `ncod` int(10) NOT NULL DEFAULT '0',
  `max_ds_branch` float(10,5) DEFAULT NULL,
  `type` enum('negative1','negative2','negative3','negative4','positive1','positive2','positive3','positive4') DEFAULT NULL,
  `note` enum('all_gaps','constant','synonymous','single_char','random') DEFAULT NULL,
  PRIMARY KEY (`sitewise_id`,`parameter_set_id`,`ncod`),
  UNIQUE KEY `aln_position_node_id_param` (`aln_position`,`node_id`,`ncod`,`parameter_set_id`),
  KEY `tree_node_id` (`tree_node_id`),
  KEY `node_id` (`node_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 
  /*!50100 
	PARTITION BY RANGE (ncod) 
	SUBPARTITION BY KEY (parameter_set_id) 
	SUBPARTITIONS 3 (
		PARTITION p0 VALUES LESS THAN (3) ENGINE = InnoDB,
		PARTITION p1 VALUES LESS THAN (6) ENGINE = InnoDB, 
		PARTITION p2 VALUES LESS THAN (9) ENGINE = InnoDB, 
		PARTITION p3 VALUES LESS THAN (10000) ENGINE = InnoDB
	) 
  */
;

CREATE TABLE IF NOT EXISTS `parameter_set` (
  `parameter_set_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `parameter_name` varchar(40) DEFAULT NULL,
  `parameter_value` mediumtext,
  UNIQUE KEY `parameter_set_parameter` (`parameter_set_id`,`parameter_name`),
  KEY `parameter_set_id` (`parameter_set_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE IF NOT EXISTS node_set (
  node_set_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  name varchar(40) NOT NULL DEFAULT '',
  allows_overlap int(10) unsigned NOT NULL DEFAULT 0,
  PRIMARY KEY (node_set_id),
  KEY node_set_name (name)
) ENGINE=InnoDB;

CREATE TABLE IF NOT EXISTS node_set_member (
  node_set_member_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  node_set_id int(10) unsigned NOT NULL,
  node_id int(10),
  FOREIGN KEY (node_set_id) REFERENCES node_set(node_set_id),
  PRIMARY KEY (node_set_member_id),
  UNIQUE KEY node_set_node_id (node_set_id,node_id)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


# Auto add schema version to database
REPLACE INTO meta (meta_key, meta_value) VALUES ("schema_version", "50");
  ^;

}

sub run_sql_commands {
  my $class = shift;
  my $dba = shift;
  my $string = shift;

  # Remove comment lines.
  $string =~ s/^#.*$//gm;
  $string =~ s/^-.*$//gm;

  print $string."\n";
  # Split the string into separate commands.
  my @sql_array = split(";",$string);
  # Remove whitespace lines.
  @sql_array = grep {$_ !~ m/^\s+$/} @sql_array;
  # Use the map function to run each command in turn.
  map {eval {$dba->dbc->do($_.";")};} @sql_array; # Run the commands in turn.
}




1;
