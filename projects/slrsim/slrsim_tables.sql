# conventions taken from the new clean schema of EnsEMBL
# use lower case and underscores
# internal ids are integers named tablename_id
# same name is given in foreign key relations

#
# Table structure for table 'meta'
#
# This table stores meta information about the compara database
#

CREATE TABLE IF NOT EXISTS meta (

  meta_id                     INT NOT NULL AUTO_INCREMENT,
  species_id                  INT UNSIGNED DEFAULT 1,
  meta_key                    VARCHAR(40) NOT NULL,
  meta_value                  VARCHAR(255) BINARY NOT NULL,
	
  PRIMARY   KEY (meta_id),
  UNIQUE    KEY species_key_value_idx (species_id, meta_key, meta_value),
  KEY species_value_idx (species_id, meta_value)

) ENGINE=InnoDB;


#
# Table structure for tables 'ncbi_taxa_node' and 'ncbi_taxa_name'
#
# Contains all taxa used in this database, which mirror the data and tree structure
# from NCBI Taxonomy database (for more details see ensembl-compara/script/taxonomy/README-taxonomy
# which explain our import process)
#

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

#
# Table structure for table 'genome_db'
#
# Contains information about the version of the genome assemblies used in this database
#

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

#
# Table structure for table 'sequence'
#

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

#
# Table structure for table 'member'
#

CREATE TABLE IF NOT EXISTS member (
  member_id                   int(10) unsigned NOT NULL auto_increment, # unique internal id
  stable_id                   varchar(40) NOT NULL, # e.g. ENSP000001234 or P31946
  version                     int(10) DEFAULT '0', 
  source_name                 varchar(40) NOT NULL,
#  source_name                 ENUM('ENSEMBLGENE','ENSEMBLPEP','Uniprot/SPTREMBL','Uniprot/SWISSPROT') NOT NULL,
  taxon_id                    int(10) unsigned NOT NULL, # FK taxon.taxon_id
  genome_db_id                int(10) unsigned, # FK genome_db.genome_db_id
  sequence_id                 int(10) unsigned, # FK sequence.sequence_id
  cdna_sequence_id            int(10) unsigned,
  gene_member_id              int(10) unsigned, # FK member.member_id
  description                 text DEFAULT NULL,
  chr_name                    char(40),
  chr_start                   int(10),
  chr_end                     int(10),
  chr_strand                  tinyint(1) NOT NULL,
  display_label               varchar(128) default NULL,

  #FOREIGN KEY (taxon_id) REFERENCES ncbi_taxa_node(taxon_id),
  #FOREIGN KEY (genome_db_id) REFERENCES genome_db(genome_db_id),
  FOREIGN KEY (sequence_id) REFERENCES sequence(sequence_id),
  FOREIGN KEY (gene_member_id) REFERENCES member(member_id),

  PRIMARY KEY (member_id),
  UNIQUE source_stable_id (stable_id, source_name),
  KEY (stable_id),
  KEY (source_name),
  KEY (sequence_id),
  KEY (gene_member_id)
) ENGINE=InnoDB COLLATE=latin1_swedish_ci;

------------------------------------------------------------------------------------
--
-- Table structure for table 'protein_tree_node'
--
-- overview:
--   This table holds the protein tree data structure, such as root, relation between 
--   parent and child, leaves
--
-- semantics:
--      node_id               -- PRIMARY node id 
--      parent_id             -- parent node id
--      root_id               -- to quickly isolated nodes of the different rooted tree sets
--      left_index            -- for fast nested set searching
--      right_index           -- for fast nested set searching
--      distance_to_parent    -- distance between node_id and its parent_id

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


------------------------------------------------------------------------------------
--
-- Table structure for table 'protein_tree_member'
--
-- overview:
--   to allow certain nodes (leaves) to have aligned protein members attached to them   
-- semantics:
--    node_id                  -- the id of node associated with this name
--    member_id                -- link to member.member_id in many-1 relation (single member per node)
--    method_link_species_set_id -- foreign key from method_link_species_set table
--    cigar_line               -- compressed alignment information 
--    cigar_start              -- protein start (0 if the whole protein is in the alignment)
--    cigar_end                -- protein end (0 if the whole protein is in the alignment)

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

------------------------------------------------------------------------------------
--
-- Table structure for table 'protein_tree_tag'
--
-- overview: 
--    to allow for tagging nodes.
--    
-- semantics:
--    node_id             -- node_id foreign key from protein_tree_node table
--    tag                 -- tag used to fecth/store a value associated to it
--    value               -- value associated with a particular tag

CREATE TABLE IF NOT EXISTS protein_tree_tag (
  node_id                int(10) unsigned NOT NULL,
  tag                    varchar(50),
  value                  mediumtext,

  FOREIGN KEY (node_id) REFERENCES protein_tree_node(node_id),

  UNIQUE tag_node_id (node_id, tag),
  KEY (node_id),
  KEY (tag)
) COLLATE=latin1_swedish_ci;

# Table sitewise_aln
# This table stores the values of calculating the sitewise dN/dS ratio
#  on node_ids (subtrees) for the GeneTrees. A subtree can also be the
#  root of the tree
# sitewise_id - identifies the sitewise entry
# aln_position - is the position in the whole GeneTree alignment, even
# if it is all_gaps in the subtree
# node_id - is the root of the subtree for which the sitewise is
# calculated
# tree_node_id - is the root of the tree. it will be equal to node_id
# if we are calculating sitewise for the whole tree
# omega is the estimated omega value at the position
# omega_lower is the lower bound of the confidence interval
# omega_upper is the upper bound of the confidence interval
# threshold_on_branch_ds is the used threshold to break a tree into
# subtrees when the dS value of a given branch is too big. This is
# defined in the configuration file for the genetree pipeline
# type is the predicted type for the codon/aminoacid
# (positive4,positive3,positive2,positive1,
#  negative4,negative3,negative2,negative1,
#  constant,all_gaps,single_character,synonymous,default)

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
	PARTITION BY KEY (parameter_set_id) 
	PARTITIONS 10
  */
;

CREATE TABLE IF NOT EXISTS gblocks_aln (
  gblocks_aln_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  parameter_set_id int(10) unsigned NOT NULL,
  node_id int(10) unsigned NOT NULL,
  aln_start int(10) unsigned NOT NULL,
  aln_end int(10) unsigned NOT NULL,
  PRIMARY KEY (gblocks_aln_id,parameter_set_id,node_id),
  UNIQUE KEY `start_end_node_id_param` (aln_start,aln_end,node_id,parameter_set_id),
  KEY node_id (node_id)
) ENGINE=InnoDB;


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

CREATE TABLE IF NOT EXISTS sitewise_genome (
  sitewise_genome_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  node_id                     int(10) unsigned NOT NULL,
  aln_position                int(10) unsigned NOT NULL,
  member_id                   int(10) unsigned NOT NULL,
  chr_name                    char(10) NOT NULL,
  chr_start                   int(10) unsigned NOT NULL,
  chr_end                     int(10) unsigned NOT NULL,
  residue                     char(1) NOT NULL,
  PRIMARY KEY (sitewise_genome_id,node_id,aln_position,chr_name),
  key chr_name (chr_name),
  UNIQUE KEY node_aln_member (node_id,aln_position,member_id,chr_name)
) ENGINE=InnoDB
  /*!50100 
	PARTITION BY KEY (chr_name)
        PARTITIONS 8
  */
;

CREATE table if not exists aln_tag (
  node_id int(10) unsigned NOT NULL,
  aln_position int(10) unsigned NOT NULL,
  member_id int(10) unsigned NOT NULL,
  stable_id varchar(20) default NULL,
  source varchar(20) default NULL,
  tag varchar(16) NOT NULL default '',
  value varchar(16) default NULL,
  feature_position int(10) unsigned default NULL,
  KEY node_id (node_id),
  KEY tag (tag),
#  FOREIGN KEY (node_id, aln_position) REFERENCES sitewise_aln(node_id, aln_position),
  UNIQUE tag_aln_id (node_id,aln_position,tag)
) ENGINE=InnoDB
  /*!50100
     PARTITION BY KEY (tag)
     PARTITIONS 5
  */
;


CREATE TABLE IF NOT EXISTS parameter_set_member (
  parameter_set_member_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  parameter_set_id int(10) unsigned NOT NULL,
  node_id int(10) unsigned NOT NULL,
  parameter_string MEDIUMTEXT NOT NULL,
  PRIMARY KEY (parameter_set_member_id),
  UNIQUE KEY parameter_set_node_id (parameter_set_id,node_id)
) ENGINE=InnoDB;


CREATE TABLE IF NOT EXISTS go_terms (
  go_term_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  `node_id` int(10) unsigned default 0,
  `member_id` int(10) unsigned NOT NULL default 0,
  `stable_id` varchar(40) default NULL,
  `go_term` varchar(20) default NULL,
  PRIMARY KEY (go_term_id),
  UNIQUE KEY `member_go_term` (`go_term`,`member_id`),
  KEY `node_id` (`node_id`),
  KEY `stable_id` (`stable_id`),
  KEY `member_id` (`member_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

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
) ENGINE=InnoDB ;


CREATE TABLE IF NOT EXISTS dataflow_rule (
  dataflow_rule_id    int(10) unsigned not null auto_increment,
  from_analysis_id    int(10) unsigned NOT NULL,
  to_analysis_url     varchar(255) default '' NOT NULL,
  branch_code         int(10) default 1 NOT NULL,

  PRIMARY KEY (dataflow_rule_id),
  UNIQUE (from_analysis_id, to_analysis_url)
)ENGINE=InnoDB ;

CREATE TABLE IF NOT EXISTS analysis (

  analysis_id                 int(10) unsigned NOT NULL auto_increment, # unique internal id
  created                     datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  logic_name                  varchar(40) not null,
  db                          varchar(120),
  db_version                  varchar(40),
  db_file                     varchar(120),
  program                     varchar(80),
  program_version             varchar(40),
  program_file                varchar(80),
  parameters                  MEDIUMTEXT,
  module                      varchar(80),
  module_version              varchar(40),
  gff_source                  varchar(40),
  gff_feature                 varchar(40),

  PRIMARY KEY (analysis_id),
  KEY logic_name_idx( logic_name ),
  UNIQUE (logic_name)

) ENGINE=InnoDB COLLATE=latin1_swedish_ci;


CREATE TABLE IF NOT EXISTS analysis_description (
  analysis_id                int(10) unsigned NOT NULL,
  description                text,
  display_label              varchar(255),
  displayable                boolean not null default 1,
  web_data                   text,

#  FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id),

  UNIQUE KEY analysis_idx( analysis_id )

) ENGINE=InnoDB COLLATE=latin1_swedish_ci;


CREATE TABLE IF NOT EXISTS analysis_ctrl_rule (
  condition_analysis_url     varchar(255) default '' NOT NULL,
  ctrled_analysis_id         int(10) unsigned NOT NULL,

  UNIQUE (condition_analysis_url, ctrled_analysis_id)
) ENGINE=InnoDB;

CREATE TABLE IF NOT EXISTS analysis_job (
  analysis_job_id           int(10) NOT NULL auto_increment,
  prev_analysis_job_id      int(10) NOT NULL,  #analysis_job which created this from rules
  analysis_id               int(10) NOT NULL,
  input_id                  MEDIUMTEXT not null,
  job_claim                 char(40) NOT NULL default '', #UUID
  hive_id                   int(10) NOT NULL,
  status                    enum('READY','BLOCKED','CLAIMED','GET_INPUT','RUN','WRITE_OUTPUT','DONE','FAILED') DEFAULT 'READY' NOT NULL,
  retry_count               int(10) default 0 not NULL,
  completed                 datetime NOT NULL,
  branch_code               int(10) default 1 NOT NULL,
  runtime_msec              int(10) default 0 NOT NULL, 
  query_count               int(10) default 0 NOT NULL, 

  PRIMARY KEY                  (analysis_job_id),
  UNIQUE KEY input_id_analysis (input_id(255), analysis_id),
  INDEX claim_analysis_status  (job_claim, analysis_id, status),
  INDEX analysis_status        (analysis_id, status),
  INDEX hive_id                (hive_id)
) ENGINE=InnoDB MAX_ROWS = 100000000 AVG_ROW_LENGTH = 390 ; 

CREATE TABLE IF NOT EXISTS analysis_job_file (
  analysis_job_id         int(10) NOT NULL,
  hive_id                 int(10) NOT NULL,
  retry                   int(10) NOT NULL,
  type                    varchar(16) NOT NULL default '',
  path                    varchar(255) NOT NULL,
  
  UNIQUE KEY job_hive_type  (analysis_job_id, hive_id, type),
  INDEX hive_id             (hive_id)
) ENGINE=InnoDB;

CREATE TABLE IF NOT EXISTS analysis_data (
  analysis_data_id  int(10) NOT NULL auto_increment,
  data              longtext,

  PRIMARY KEY (analysis_data_id),
  KEY data (data(100))
) ENGINE=InnoDB;

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
) ENGINE=InnoDB;

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
) ENGINE=InnoDB;

CREATE TABLE IF NOT EXISTS monitor (
  time                  datetime NOT NULL default '0000-00-00 00:00:00',
  workers               int(10) NOT NULL default '0',
  throughput            float default NULL,
  per_worker            float default NULL,
  analysis              varchar(255) default NULL
) ENGINE=InnoDB;

# Auto add schema version to database
REPLACE INTO meta (meta_key, meta_value) VALUES ("schema_version", "50");
