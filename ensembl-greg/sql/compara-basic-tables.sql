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

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


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
  KEY (rank),
  KEY (left_index),
  KEY (right_index)
) COLLATE=latin1_swedish_ci;

CREATE TABLE IF NOT EXISTS ncbi_taxa_name (
  taxon_id                    int(10) unsigned NOT NULL,

  name                        varchar(255),
  name_class                  varchar(50),

  KEY (taxon_id),
  KEY (name),
  KEY (name_class)
) COLLATE=latin1_swedish_ci;

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
  locator                     varchar(400),

  FOREIGN KEY (taxon_id) REFERENCES ncbi_taxa_node(taxon_id),

  PRIMARY KEY (genome_db_id),
  UNIQUE name (name,assembly,genebuild)
) COLLATE=latin1_swedish_ci;


#
# Table structure for table 'method_link'
#
# Specifies which kind of link can exist between species
# (dna/dna alignment, synteny regions, homologous gene pairs,...)
#

CREATE TABLE IF NOT EXISTS method_link (
  method_link_id              int(10) unsigned NOT NULL auto_increment, # unique internal id
  type                        varchar(50) DEFAULT '' NOT NULL,
  class                       varchar(50) DEFAULT '' NOT NULL,

  PRIMARY KEY (method_link_id),
  KEY type (type)
) COLLATE=latin1_swedish_ci;




#
# Table structure for table 'species_set'
#
-- KEY species_set_id is a multiple key. It defines a set of species
--   (genome_db_ids).
CREATE TABLE IF NOT EXISTS species_set (
  species_set_id              int(10) unsigned NOT NULL auto_increment,
  genome_db_id                int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (species_set_id,genome_db_id)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


#
# Table structure for table 'species_set_tag'
#
# This table is used to store options on clades and group of species. It
# has been initially developed for the gene tree view.
#

CREATE TABLE IF NOT EXISTS species_set_tag (
  species_set_id              int(10) unsigned NOT NULL, # FK species_set.species_set_id
  tag                         varchar(50) NOT NULL,
  value                       mediumtext,

  # FOREIGN KEY (species_set_id) REFERENCES species_set(species_set_id),

  UNIQUE KEY tag_species_set_id (species_set_id,tag)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


#
# Table structure for table 'method_link_species_set'
#

CREATE TABLE IF NOT EXISTS method_link_species_set (
  method_link_species_set_id  int(10) unsigned NOT NULL AUTO_INCREMENT, # unique internal id
  method_link_id              int(10) unsigned, # FK method_link.method_link_id
  species_set_id              int(10) unsigned NOT NULL default '0',
  name                        varchar(255) NOT NULL default '',
  source                      varchar(255) NOT NULL default 'ensembl',
  url                         varchar(255) NOT NULL default '',

  FOREIGN KEY (method_link_id) REFERENCES method_link(method_link_id),
  # species_set(species_set_id) is not a unique key. Some RDBMS may complain
  # FOREIGN KEY (species_set_id) REFERENCES species_set(species_set_id),

  PRIMARY KEY (method_link_species_set_id),
  UNIQUE KEY method_link_id (method_link_id,species_set_id)
) COLLATE=latin1_swedish_ci;

#
# Table structure for table 'dnafrag'
#
-- Index <name> has genome_db_id in the first place because unless fetching all danfrags
--   or fetching by dnafrag_id, genome_db_id appears always in the WHERE clause
-- Unique key <name> is used to ensure that
--   Bio::EnsEMBL::Compara::DBSQL::DnaFragAdaptor->fetch_by_GenomeDB_and_name
--   will always fetch a single row. This can be used in the EnsEMBL Compara DB
--   because we store top-level dnafrags only.

CREATE TABLE IF NOT EXISTS dnafrag (
  dnafrag_id                  bigint unsigned NOT NULL auto_increment, # unique internal id
  length                      int(11) DEFAULT '0' NOT NULL,
  name                        varchar(40) DEFAULT '' NOT NULL,
  genome_db_id                int(10) unsigned DEFAULT '0' NOT NULL, # FK genome_db.genome_db_id
  coord_system_name           varchar(40) DEFAULT NULL,
  is_reference                tinyint(1) DEFAULT 1,

  FOREIGN KEY (genome_db_id) REFERENCES genome_db(genome_db_id),

  PRIMARY KEY (dnafrag_id),
  UNIQUE name (genome_db_id, name)
) COLLATE=latin1_swedish_ci;


#
# Table structure for table 'genomic_align_block'
#
#    This table indexes the genomic alignments
#
-- All queries in the API uses the primary key as rows are always fetched using
--   the genomic_align_block_id. The key "method_link_species_set_id" is used by
--   MART when fetching all the genomic_align_blocks corresponding to a given
--   method_link_species_set_id.z

CREATE TABLE IF NOT EXISTS genomic_align_block (
  genomic_align_block_id      bigint unsigned NOT NULL AUTO_INCREMENT, # unique internal id
  method_link_species_set_id  int(10) unsigned DEFAULT '0' NOT NULL, # FK method_link_species_set_id.method_link_species_set_id
  score                       double,
  perc_id                     tinyint(3) unsigned DEFAULT NULL,
  length                      int(10),
  group_id                    bigint unsigned DEFAULT NULL,

  # method_link_species_set(method_link_species_set_id) is not a unique key. Some RDBMS may complain
  # FOREIGN KEY (method_link_species_set_id) REFERENCES method_link_species_set(method_link_species_set_id),

  PRIMARY KEY genomic_align_block_id (genomic_align_block_id),
  KEY method_link_species_set_id (method_link_species_set_id)
) COLLATE=latin1_swedish_ci;


#
# Table structure for table 'genomic_align'
#
#   This table stores the sequences belonging to the same genomic_align_block entry
#
-- primary key is used when fecthing by dbID
-- key genomic_align_block_id is used when fetching by genomic_align_block_id
-- key dnafrag is used in all other queries

CREATE TABLE IF NOT EXISTS genomic_align (
  genomic_align_id            bigint unsigned NOT NULL AUTO_INCREMENT, # unique internal id
  genomic_align_block_id      bigint unsigned NOT NULL, # FK genomic_align_block.genomic_align_block_id
  method_link_species_set_id  int(10) unsigned DEFAULT '0' NOT NULL, # FK method_link_species_set_id.method_link_species_set_id
  dnafrag_id                  bigint unsigned DEFAULT '0' NOT NULL, # FK dnafrag.dnafrag_id
  dnafrag_start               int(10) DEFAULT '0' NOT NULL,
  dnafrag_end                 int(10) DEFAULT '0' NOT NULL,
  dnafrag_strand              tinyint(4) DEFAULT '0' NOT NULL,
  cigar_line                  mediumtext,
  level_id                    tinyint(2) unsigned DEFAULT '0' NOT NULL,

  FOREIGN KEY (genomic_align_block_id) REFERENCES genomic_align_block(genomic_align_block_id),
  # method_link_species_set(method_link_species_set_id) is not a unique key. Some RDBMS may complain
  # FOREIGN KEY (method_link_species_set_id) REFERENCES method_link_species_set(method_link_species_set_id),
  FOREIGN KEY (dnafrag_id) REFERENCES dnafrag(dnafrag_id),

  PRIMARY KEY genomic_align_id (genomic_align_id),
  KEY genomic_align_block_id (genomic_align_block_id),
  KEY method_link_species_set_id (method_link_species_set_id),
  KEY dnafrag (dnafrag_id, method_link_species_set_id, dnafrag_start, dnafrag_end)
) MAX_ROWS = 1000000000 AVG_ROW_LENGTH = 60 COLLATE=latin1_swedish_ci;


#
# Table structure for table 'genomic_align_group'
#
#   This table can store several groupings of the genomic aligned sequences
#

CREATE TABLE IF NOT EXISTS genomic_align_group (
  node_id                    bigint unsigned NOT NULL AUTO_INCREMENT, # internal id, groups genomic_align_ids
  genomic_align_id            bigint unsigned NOT NULL, # FK genomic_align.genomic_align_id

  FOREIGN KEY (genomic_align_id) REFERENCES genomic_align(genomic_align_id),

  KEY node_id (node_id),
  UNIQUE KEY genomic_align_id (genomic_align_id)
) COLLATE=latin1_swedish_ci;


#
# Table structure for table 'genomic_align_tree'
#
#   This table stores the phylogenetic relationship between different genomic_align entries.
#   This is required to store the tree alignments, i.e. multiple sequence alignments with
#   ancestral sequence reconstruction. This table stores the tree underlying each tree
#   alignments
#
-- primary key is a foreign key to genomic_align_group.group_id

CREATE TABLE IF NOT EXISTS genomic_align_tree (
  node_id                     bigint(20) unsigned NOT NULL AUTO_INCREMENT, # internal id, FK genomic_align_group.genomic_align_id
  parent_id                   bigint(20) unsigned NOT NULL default '0',
  root_id                     bigint(20) unsigned NOT NULL default '0',
  left_index                  int(10) NOT NULL default '0',
  right_index                 int(10) NOT NULL default '0',
  left_node_id                bigint(10) NOT NULL default '0',
  right_node_id               bigint(10) NOT NULL default '0',
  distance_to_parent          double NOT NULL default '1',

  # genomic_align_group(node_id) is not a unique key. Some RDBMS may complain
  # FOREIGN KEY (node_id) REFERENCES genomic_align_group(node_id),

  PRIMARY KEY node_id (node_id),
  KEY parent_id (parent_id),
  KEY root_id (root_id),
  KEY left_index (left_index),
  KEY right_index (right_index)
) COLLATE=latin1_swedish_ci;


#
# Table structure for table 'synteny_region'
#
# We have now decided that Synteny is inherently pairwise
# these tables hold the pairwise information for the synteny
# regions. We reuse the dnafrag table as a link out for identifiers
# (eg, '2' on mouse).
#

CREATE TABLE IF NOT EXISTS synteny_region (
  synteny_region_id           int(10) unsigned NOT NULL auto_increment, # unique internal id
  method_link_species_set_id  int(10) unsigned NOT NULL, # FK method_link_species_set.method_link_species_set_id

  # method_link_species_set(method_link_species_set_id) is not a unique key. Some RDBMS may complain
  # FOREIGN KEY (method_link_species_set_id) REFERENCES method_link_species_set(method_link_species_set_id),

  PRIMARY KEY (synteny_region_id),
  KEY (method_link_species_set_id)
) COLLATE=latin1_swedish_ci;


#
# Table structure for table 'dnafrag_region'
#

CREATE TABLE IF NOT EXISTS dnafrag_region (
  synteny_region_id           int(10) unsigned DEFAULT '0' NOT NULL, # unique internal id
  dnafrag_id                  bigint unsigned DEFAULT '0' NOT NULL, # FK dnafrag.dnafrag_id
  dnafrag_start               int(10) unsigned DEFAULT '0' NOT NULL,
  dnafrag_end                 int(10) unsigned DEFAULT '0' NOT NULL,
  dnafrag_strand              tinyint(4) DEFAULT '0' NOT NULL,

  FOREIGN KEY (synteny_region_id) REFERENCES synteny_region(synteny_region_id),
  FOREIGN KEY (dnafrag_id) REFERENCES dnafrag(dnafrag_id),

  KEY synteny (synteny_region_id,dnafrag_id),
  KEY synteny_reversed (dnafrag_id,synteny_region_id)
) COLLATE=latin1_swedish_ci;


#
# Table structure for table 'sequence'
#

CREATE TABLE IF NOT EXISTS sequence (
  sequence_id                 int(10) unsigned NOT NULL auto_increment, # unique internal id
  length                      int(10) NOT NULL,
  sequence                    longtext NOT NULL,

  PRIMARY KEY (sequence_id),
  KEY sequence (sequence(18))
) MAX_ROWS = 10000000 AVG_ROW_LENGTH = 19000 COLLATE=latin1_swedish_ci;

-- overview:
--   This table holds the sequence exon boundaries information
CREATE TABLE IF NOT EXISTS sequence_exon_bounded (
  sequence_exon_bounded_id    int(10) unsigned NOT NULL auto_increment, # unique internal id
  member_id                   int(10) unsigned NOT NULL, # unique internal id
  length                      int(10) NOT NULL,
  sequence_exon_bounded       longtext NOT NULL,

  FOREIGN KEY (member_id) REFERENCES member(member_id),

  PRIMARY KEY (sequence_exon_bounded_id),
  KEY (member_id),
  KEY sequence_exon_bounded (sequence_exon_bounded(18))
) MAX_ROWS = 10000000 AVG_ROW_LENGTH = 19000 COLLATE=latin1_swedish_ci;

-- overview:
--   This table holds the sequence cds information
CREATE TABLE IF NOT EXISTS sequence_cds (
  sequence_cds_id             int(10) unsigned NOT NULL auto_increment, # unique internal id
  member_id                   int(10) unsigned NOT NULL, # unique internal id
  length                      int(10) NOT NULL,
  sequence_cds                longtext NOT NULL,

  FOREIGN KEY (member_id) REFERENCES member(member_id),

  PRIMARY KEY (sequence_cds_id),
  KEY (member_id),
  KEY sequence_cds (sequence_cds(64))
) MAX_ROWS = 10000000 AVG_ROW_LENGTH = 60000 COLLATE=latin1_swedish_ci;

#
# Table structure for table 'member'
#

CREATE TABLE IF NOT EXISTS member (
  member_id                   int(10) unsigned NOT NULL auto_increment, # unique internal id
  stable_id                   varchar(128) NOT NULL, # e.g. ENSP000001234 or P31946
  version                     int(10) DEFAULT '0',
  source_name                 varchar(40) NOT NULL,
#  source_name                 ENUM('ENSEMBLGENE','ENSEMBLPEP','Uniprot/SPTREMBL','Uniprot/SWISSPROT','ENSEMBLTRANS','EXTERNALCDS') NOT NULL,
  taxon_id                    int(10) unsigned NOT NULL, # FK taxon.taxon_id
  genome_db_id                int(10) unsigned, # FK genome_db.genome_db_id
  sequence_id                 int(10) unsigned, # FK sequence.sequence_id
  cdna_sequence_id                 int(10) unsigned, # FK sequence.sequence_id
  gene_member_id              int(10) unsigned, # FK member.member_id
  description                 text DEFAULT NULL,
  chr_name                    char(40),
  chr_start                   int(10),
  chr_end                     int(10),
  chr_strand                  tinyint(1) NOT NULL,
  display_label               varchar(128) default NULL,

#  FOREIGN KEY (taxon_id) REFERENCES ncbi_taxa_node(taxon_id),
#  FOREIGN KEY (genome_db_id) REFERENCES genome_db(genome_db_id),
  FOREIGN KEY (sequence_id) REFERENCES sequence(sequence_id),
  FOREIGN KEY (cdna_sequence_id) REFERENCES sequence(sequence_id),
  FOREIGN KEY (gene_member_id) REFERENCES member(member_id),

  PRIMARY KEY (member_id),
  UNIQUE source_stable_id (stable_id, source_name),
  KEY (stable_id),
  KEY (source_name),
  KEY (sequence_id),
  KEY (gene_member_id),
  KEY gdb_name_start_end (genome_db_id,chr_name,chr_start,chr_end)
) MAX_ROWS = 100000000 COLLATE=latin1_swedish_ci;


#
# Table structure for table 'subset'
#

CREATE TABLE IF NOT EXISTS subset (
 subset_id      int(10) NOT NULL auto_increment,
 description    varchar(255),
 dump_loc       varchar(255),

 PRIMARY KEY (subset_id),
 UNIQUE (description)
);

#
# Table structure for table 'subset_member'
#

CREATE TABLE IF NOT EXISTS subset_member (
 subset_id   int(10) NOT NULL,
 member_id   int(10) NOT NULL,

 KEY (member_id),
 UNIQUE subset_member_id (subset_id, member_id)
);


-- ----------------------------------------------------------------------------------
--
-- Table structure for table 'analysis'
--
-- semantics:
-- analysis_id - internal id
-- created   - date to distinguish newer and older versions off the
--             same analysis. Not well maintained so far.
-- logic_name  string to identify the analysis. Used mainly inside pipeline.
-- db, db_version, db_file
--  - db should be a database name, db version the version of that db
--    db_file the file system location of that database,
--    probably wiser to generate from just db and configurations
-- program, program_version,program_file
--  - The binary used to create a feature. Similar semantic to above
-- module, module_version
--  - Perl module names (RunnableDBS usually) executing this analysis.
-- parameters a paramter string which is processed by the perl module
-- gff_source, gff_feature
--  - how to make a gff dump from features with this analysis

CREATE TABLE IF NOT EXISTS analysis (

  analysis_id                 int(10) unsigned NOT NULL auto_increment, # unique internal id
  created                     datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  logic_name                  varchar(40) not null,
  db                          varchar(120),
  db_version                  varchar(40),
  db_file                     varchar(255),
  program                     varchar(255),
  program_version             varchar(40),
  program_file                varchar(255),
  parameters                  TEXT,
  module                      varchar(80),
  module_version              varchar(40),
  gff_source                  varchar(40),
  gff_feature                 varchar(40),

  PRIMARY KEY (analysis_id),
  KEY logic_name_idx( logic_name ),
  UNIQUE (logic_name)

) COLLATE=latin1_swedish_ci;

CREATE TABLE IF NOT EXISTS analysis_description (
  analysis_id                int(10) unsigned NOT NULL,
  description                text,
  display_label              varchar(255),
  displayable                boolean not null default 1,
  web_data                   text,

  FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id),

  UNIQUE KEY analysis_idx( analysis_id )

) COLLATE=latin1_swedish_ci;

-- ----------------------------------------------------------------------------------
--
-- Table structure for table 'peptide_align_feature'
--
-- overview: This tables stores the raw HSP local alignment results
--           of peptide to peptide alignments returned by a BLAST run
--           it is translated from a FeaturePair object
-- semantics:
-- peptide_align_feature_id  - internal id
-- qmember_id                - member.member_id of query peptide
-- hmember_id                - member.member_id of hit peptide
-- qgenome_db_id             - genome_db_id of query peptide (for query optimization)
-- hgenome_db_id             - genome_db_id of hit peptide (for query optimization)
-- qstart                    - start pos in query peptide sequence
-- qend                      - end  pos in query peptide sequence
-- hstart                    - start pos in hit peptide sequence
-- hend                      - end  pos in hit peptide sequence
-- score                     - blast score for this HSP
-- evalue                    - blast evalue for this HSP
-- align_length              - alignment length of HSP
-- identical_matches         - blast HSP match score
-- positive_matches          - blast HSP positive score
-- perc_ident                - percent identical matches in the HSP length
-- perc_pos                  - percent positive matches in the HSP length
-- cigar_line                - cigar string coding the actual alignment

CREATE TABLE IF NOT EXISTS peptide_align_feature (

  peptide_align_feature_id    int(10) unsigned NOT NULL auto_increment, # unique internal id
  qmember_id                  int(10) unsigned NOT NULL, # FK member.member_id
  hmember_id                  int(10) unsigned NOT NULL, # FK member.member_id
  qgenome_db_id               int(10) unsigned NOT NULL, # FK genome.genome_id
  hgenome_db_id               int(10) unsigned NOT NULL, # FK genome.genome_id
  analysis_id                 int(10) unsigned NOT NULL, # FK analysis.analysis_id
  qstart                      int(10) DEFAULT '0' NOT NULL,
  qend                        int(10) DEFAULT '0' NOT NULL,
  hstart                      int(11) DEFAULT '0' NOT NULL,
  hend                        int(11) DEFAULT '0' NOT NULL,
  score                       double(16,4) DEFAULT '0.0000' NOT NULL,
  evalue                      double,
  align_length                int(10),
  identical_matches           int(10),
  perc_ident                  int(10),
  positive_matches            int(10),
  perc_pos                    int(10),
  hit_rank                    int(10),
  cigar_line                  mediumtext,

#  FOREIGN KEY (qmember_id) REFERENCES member(member_id),
#  FOREIGN KEY (hmember_id) REFERENCES member(member_id),
#  FOREIGN KEY (qgenome_db_id) REFERENCES genome_db(genome_db_id),
#  FOREIGN KEY (hgenome_db_id) REFERENCES genome_db(genome_db_id),
#  FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id),

  PRIMARY KEY (peptide_align_feature_id),
  KEY hmember_hit (hmember_id, hit_rank)

#  KEY qmember_id  (qmember_id),
#  KEY hmember_id  (hmember_id),
#  KEY hmember_qgenome  (hmember_id, qgenome_db_id),
#  KEY qmember_hgenome  (qmember_id, hgenome_db_id)
) MAX_ROWS = 300000000 AVG_ROW_LENGTH = 133 COLLATE=latin1_swedish_ci;


#
# Table structure for table 'homology'
#

CREATE TABLE IF NOT EXISTS homology (
  homology_id                 int(10) unsigned NOT NULL auto_increment, # unique internal id
  stable_id                   varchar(40),
  method_link_species_set_id  int(10) unsigned NOT NULL, # FK method_link_species_set.method_link_species_set_id
  description                 ENUM('ortholog_one2one','apparent_ortholog_one2one','ortholog_one2many','ortholog_many2many','within_species_paralog','other_paralog','putative_gene_split','contiguous_gene_split','between_species_paralog','possible_ortholog','UBRH','BRH','MBRH','RHS'),
  subtype                     varchar(40) NOT NULL DEFAULT '',
  dn                          float(10,5),
  ds                          float(10,5),
  n                           float(10,1),
  s                           float(10,1),
  lnl                         float(10,3),
  threshold_on_ds             float(10,5),
  ancestor_node_id            int(10) unsigned NOT NULL,
  tree_node_id                int(10) unsigned NOT NULL,

  # method_link_species_set(method_link_species_set_id) is not a unique key. Some RDBMS may complain
  # FOREIGN KEY (method_link_species_set_id) REFERENCES method_link_species_set(method_link_species_set_id),

  PRIMARY KEY (homology_id),
  KEY (method_link_species_set_id),
  KEY (tree_node_id)
) COLLATE=latin1_swedish_ci;


#
# Table structure for table 'homology_member'
#

CREATE TABLE IF NOT EXISTS homology_member (
  homology_id                 int(10) unsigned NOT NULL, # FK homology.homology_id
  member_id                   int(10) unsigned NOT NULL, # FK member.member_id
  peptide_member_id           int(10) unsigned, # FK member.member_id
  peptide_align_feature_id    int(10) unsigned, # FK peptide_align_feature.peptide_align_feature_id
  cigar_line                  mediumtext,
  cigar_start                 int(10),
  cigar_end                   int(10),
  perc_cov                    int(10),
  perc_id                     int(10),
  perc_pos                    int(10),

  FOREIGN KEY (homology_id) REFERENCES homology(homology_id),
  FOREIGN KEY (member_id) REFERENCES member(member_id),
  FOREIGN KEY (peptide_member_id) REFERENCES member(member_id),
#  FOREIGN KEY (peptide_align_feature_id) REFERENCES peptide_align_feature(peptide_align_feature_id),

  UNIQUE homology_member_id (homology_id,member_id),
  KEY (homology_id),
  KEY (member_id),
  KEY (peptide_member_id),
  KEY (peptide_align_feature_id)
) MAX_ROWS = 300000000 COLLATE=latin1_swedish_ci;


#
# Table structure for table 'family'
#

CREATE TABLE IF NOT EXISTS family (
  family_id                   int(10) unsigned NOT NULL auto_increment, # unique internal id
  stable_id                   varchar(40) NOT NULL, # unique stable id, e.g. 'ENSFM'.'0053'.'1234567890'
  version                     INT UNSIGNED NOT NULL,# version of the stable_id (changes only when members move to/from existing families)
  method_link_species_set_id  int(10) unsigned NOT NULL, # FK method_link_species_set.method_link_species_set_id
  description                 varchar(255),
  description_score           double,

  # method_link_species_set(method_link_species_set_id) is not a unique key. Some RDBMS may complain
  # FOREIGN KEY (method_link_species_set_id) REFERENCES method_link_species_set(method_link_species_set_id),

  PRIMARY KEY (family_id),
  UNIQUE (stable_id),
  KEY (method_link_species_set_id),
  KEY (description)
) COLLATE=latin1_swedish_ci;


#
# Table structure for table 'family_member'
#

CREATE TABLE IF NOT EXISTS family_member (
  family_id                   int(10) unsigned NOT NULL, # FK family.family_id
  member_id                   int(10) unsigned NOT NULL, # FK member.memeber_id
  cigar_line                  mediumtext,

  FOREIGN KEY (family_id) REFERENCES family(family_id),
  FOREIGN KEY (member_id) REFERENCES member(member_id),

  UNIQUE family_member_id (family_id,member_id),
  KEY (family_id),
  KEY (member_id)
) COLLATE=latin1_swedish_ci;


#
# Table structure for table 'domain'
#

CREATE TABLE IF NOT EXISTS domain (
  domain_id                   int(10) unsigned NOT NULL auto_increment, # unique internal id
  stable_id                   varchar(40) NOT NULL,
#  source_id                   int(10) NOT NULL,
  method_link_species_set_id  int(10) unsigned NOT NULL, # FK method_link_species_set.method_link_species_set_id
  description                 varchar(255),

  # method_link_species_set(method_link_species_set_id) is not a unique key. Some RDBMS may complain
  # FOREIGN KEY (method_link_species_set_id) REFERENCES method_link_species_set(method_link_species_set_id),

  PRIMARY KEY (domain_id),
  UNIQUE (stable_id, method_link_species_set_id)
) COLLATE=latin1_swedish_ci;


#
# Table structure for table 'domain_member'
#

CREATE TABLE IF NOT EXISTS domain_member (
  domain_id                   int(10) unsigned NOT NULL, # FK domain.domain_id
  member_id                   int(10) unsigned NOT NULL, # FK member.member_id
  member_start                int(10),
  member_end                  int(10),

  FOREIGN KEY (domain_id) REFERENCES domain(domain_id),
  FOREIGN KEY (member_id) REFERENCES member(member_id),

  UNIQUE (domain_id,member_id,member_start,member_end),
  UNIQUE (member_id,domain_id,member_start,member_end)
) COLLATE=latin1_swedish_ci;

#
# Table structure for table 'conservation_score'
#

CREATE TABLE IF NOT EXISTS conservation_score (
  genomic_align_block_id bigint unsigned not null,
  window_size            smallint unsigned not null,
  position	         int unsigned not null,
  expected_score         blob,
  diff_score             blob,

  FOREIGN KEY (genomic_align_block_id) REFERENCES genomic_align_block(genomic_align_block_id),

  KEY (genomic_align_block_id, window_size)
) MAX_ROWS = 15000000 AVG_ROW_LENGTH = 841 COLLATE=latin1_swedish_ci;

#
# Table structure for table 'constrained_element'
CREATE TABLE IF NOT EXISTS constrained_element (
  constrained_element_id bigint(20) unsigned NOT NULL,
  dnafrag_id bigint unsigned NOT NULL,
  dnafrag_start int(12) unsigned NOT NULL,
  dnafrag_end int(12) unsigned NOT NULL,
  dnafrag_strand int(2),
  method_link_species_set_id int(10) unsigned NOT NULL,
  p_value mediumtext,
  taxonomic_level mediumtext,
  score double NOT NULL default '0',

  FOREIGN KEY (dnafrag_id) REFERENCES dnafrag(dnafrag_id),
  FOREIGN KEY (method_link_species_set_id) REFERENCES method_link_species_set(method_link_species_set_id),

  KEY constrained_element_id_idx (constrained_element_id),
  KEY mlssid_idx (method_link_species_set_id),
  KEY mlssid_dfId_dfStart_dfEnd_idx (method_link_species_set_id,dnafrag_id,dnafrag_start,dnafrag_end),
  KEY mlssid_dfId_idx (method_link_species_set_id,dnafrag_id)
) COLLATE=latin1_swedish_ci;

-- ----------------------------------------------------------------------------------
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
--      clusterset_id         -- node id of the set of clusters
--      left_index            -- for fast nested set searching
--      right_index           -- for fast nested set searching
--      distance_to_parent    -- distance between node_id and its parent_id

CREATE TABLE IF NOT EXISTS protein_tree_node (
  node_id                         int(10) unsigned NOT NULL auto_increment, # unique internal id
  parent_id                       int(10) unsigned NOT NULL,
  root_id                         int(10) unsigned NOT NULL,
  clusterset_id                   int(10) unsigned NOT NULL,
  left_index                      int(10) NOT NULL,
  right_index                     int(10) NOT NULL,
  distance_to_parent              double default 1.0 NOT NULL,

  PRIMARY KEY (node_id),
  KEY (parent_id),
  KEY (root_id),
  KEY (left_index),
  KEY (right_index)
) COLLATE=latin1_swedish_ci;

-- overview:
--   This table holds the member information for the nc trees
CREATE TABLE IF NOT EXISTS nc_tree_node LIKE protein_tree_node;

-- overview:
--   This table holds the few big super protein tree alignments that are then broken down.
CREATE TABLE IF NOT EXISTS super_protein_tree_node LIKE protein_tree_node;

-- ----------------------------------------------------------------------------------
--
-- Table structure for table 'protein_tree_member'
--
-- overview:
--   to allow certain nodes (leaves) to have aligned protein members attached to them
-- semantics:
--    node_id                  -- the id of node associated with this name
--    root_id                  -- the id of tree root node
--    member_id                -- link to member.member_id in many-1 relation (single member per node)
--    method_link_species_set_id -- foreign key from method_link_species_set table
--    cigar_line               -- compressed alignment information
--    cigar_start              -- protein start (0 if the whole protein is in the alignment)
--    cigar_end                -- protein end (0 if the whole protein is in the alignment)

CREATE TABLE IF NOT EXISTS protein_tree_member (
  node_id                     int(10) unsigned NOT NULL,
  root_id                     int(10) unsigned NOT NULL,
  member_id                   int(10) unsigned NOT NULL,
  method_link_species_set_id  int(10) unsigned NOT NULL,
  cigar_line                  mediumtext,
  cigar_start                 int(10),
  cigar_end                   int(10),

  FOREIGN KEY (node_id) REFERENCES protein_tree_node(node_id),
  FOREIGN KEY (root_id) REFERENCES protein_tree_node(root_id),

  UNIQUE (node_id),
  KEY (member_id)
) COLLATE=latin1_swedish_ci;

-- overview:
--   This table holds the member information for the nc trees
CREATE TABLE IF NOT EXISTS nc_tree_member LIKE protein_tree_member;

-- overview:
--   This table holds the few big super protein tree alignments that are then broken down.
CREATE TABLE IF NOT EXISTS super_protein_tree_member LIKE protein_tree_member;

-- ----------------------------------------------------------------------------------
--
-- Table structure for table 'protein_tree_member_score'
--
-- overview:
--   to allow certain nodes (leaves) to have aligned protein member_scores attached to them
-- semantics:
--    node_id                  -- the id of node associated with this name
--    root_id                  -- the id of the tree root
--    member_id                -- link to member.member_id in many-1 relation (single member per node)
--    method_link_species_set_id -- foreign key from method_link_species_set table
--    cigar_line               -- string with the alignment score values 
--    cigar_start              -- protein start (0 if the whole protein is in the alignment)
--    cigar_end                -- protein end (0 if the whole protein is in the alignment)

CREATE TABLE IF NOT EXISTS protein_tree_member_score (
  node_id                     int(10) unsigned NOT NULL,
  root_id                     int(10) unsigned NOT NULL,
  member_id                   int(10) unsigned NOT NULL,
  method_link_species_set_id  int(10) unsigned NOT NULL,
  cigar_line                  mediumtext,
  cigar_start                 int(10),
  cigar_end                   int(10),

  FOREIGN KEY (node_id) REFERENCES protein_tree_node(node_id),
  FOREIGN KEY (root_id) REFERENCES protein_tree_node(root_id),

  UNIQUE (node_id),
  KEY (member_id)
) COLLATE=latin1_swedish_ci;


-- ----------------------------------------------------------------------------------
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

-- overview:
--   This table holds the tag information for the nc trees
CREATE TABLE IF NOT EXISTS nc_tree_tag LIKE protein_tree_tag;

CREATE TABLE IF NOT EXISTS nc_profile (
  model_id                    varchar(10) NOT NULL,
  name                        varchar(40),
  type                        varchar(40) DEFAULT 'ncrna' NOT NULL,
  hc_profile                  mediumtext,

  PRIMARY KEY (model_id)
) COLLATE=latin1_swedish_ci;


-- overview:
--   This table holds the few big super protein tree alignments that are then broken down.
CREATE TABLE IF NOT EXISTS super_protein_tree_tag LIKE protein_tree_tag;

-- ------------------------------------------------------------------------------------
--
-- Table structure for table 'protein_tree_stable_id'
--
-- overview:
--     to allow protein trees have trackable stable_ids.
--
-- semantics:
--    node_id           - node_id of the root of the tree
--    stable_id         - the main part of the stable_id ( follows the pattern: label(5).release_introduced(4).unique_id(10) )
--    version           - numeric version of the stable_id (changes only when members move to/from existing trees)

CREATE TABLE IF NOT EXISTS protein_tree_stable_id (
    node_id   INT(10) UNSIGNED NOT NULL,
    stable_id VARCHAR(40)  NOT NULL, # unique stable id, e.g. 'ENSGT'.'0053'.'1234567890'
    version   INT UNSIGNED NOT NULL, # version of the stable_id (changes only when members move to/from existing trees)
    PRIMARY KEY ( node_id ),
    UNIQUE KEY ( stable_id )
);


-- ------------------------------------------------------------------------------------
--
-- Table structure for table 'mapping_session'
--
-- overview:
--      A single mapping_session is the event when mapping between two given releases
--      for a particular class type ('family' or 'tree') is loaded.
--      The whole event is thought to happen momentarily at 'when_mapped' (used for sorting in historical order).

CREATE TABLE IF NOT EXISTS mapping_session (
    mapping_session_id INT UNSIGNED NOT NULL AUTO_INCREMENT,
    type               ENUM('family', 'tree'),
    when_mapped        TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    rel_from           INT UNSIGNED,
    rel_to             INT UNSIGNED,
    PRIMARY KEY ( mapping_session_id ),
    UNIQUE KEY  ( type, rel_from, rel_to )
);


-- ------------------------------------------------------------------------------------
--
-- Table structure for table 'stable_id_history'
--
-- overview:
--      'stable_id_history' table keeps the history of stable_id changes from one release to another.
--
--      The primary key 'object' describes a set of members migrating from stable_id_from to stable_id_to.
--      Their volume (related to the 'shared_size' of the new class) is reflected by the fractional 'contribution' field.
--
--      Since both stable_ids are listed in the primary key,
--      they are not allowed to be NULLs. We shall treat empty strings as NULLs.
--
--      If stable_id_from is empty, it means these members are newcomers into the new release.
--      If stable_id_to is empty, it means these previously known members are disappearing in the new release.
--      If both neither stable_id_from nor stable_id_to is empty, these members are truly migrating.

CREATE TABLE IF NOT EXISTS stable_id_history (
    mapping_session_id INT UNSIGNED NOT NULL,
    stable_id_from     VARCHAR(40) NOT NULL DEFAULT '',
    version_from       INT UNSIGNED NULL DEFAULT NULL,
    stable_id_to       VARCHAR(40) NOT NULL DEFAULT '',
    version_to         INT UNSIGNED NULL DEFAULT NULL,
    contribution       FLOAT,
    PRIMARY KEY ( mapping_session_id, stable_id_from, stable_id_to )
);


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

CREATE TABLE IF NOT EXISTS sitewise_aln (
  sitewise_id                 int(10) unsigned NOT NULL auto_increment, # unique internal id
  aln_position                int(10) unsigned NOT NULL,
  node_id                     int(10) unsigned NOT NULL,
  tree_node_id                int(10) unsigned NOT NULL,
  omega                       float(10,5),
  omega_lower                 float(10,5),
  omega_upper                 float(10,5),
  optimal                     float(10,5),
  ncod                        int(10),
  threshold_on_branch_ds      float(10,5),
  type                        ENUM('single_character','random','all_gaps','constant','default','negative1','negative2','negative3','negative4','positive1','positive2','positive3','positive4','synonymous') NOT NULL,

  FOREIGN KEY (node_id) REFERENCES protein_tree_node(node_id),

  UNIQUE aln_position_node_id_ds (aln_position,node_id,threshold_on_branch_ds),
  PRIMARY KEY (sitewise_id),
  KEY (tree_node_id),
  KEY (node_id)
) COLLATE=latin1_swedish_ci;


CREATE TABLE IF NOT EXISTS `sitewise_omega` (
  `data_id` int unsigned NOT NULL,
  `node_id` int unsigned NOT NULL,
  `parameter_set_id` tinyint unsigned NOT NULL,
  `aln_position` mediumint unsigned NOT NULL,

  `aln_position_fraction` float(10,5) DEFAULT NULL,
  `ncod` smallint unsigned NOT NULL,

  `omega` float(10,5) DEFAULT NULL,
  `omega_lower` float(10,5) DEFAULT NULL,
  `omega_upper` float(10,5) DEFAULT NULL,
  `lrt_stat` float(10,5) DEFAULT NULL,
  `type` enum('negative1','negative2','negative3','negative4','positive1','positive2','positive3','positive4') DEFAULT NULL,
  `note` enum('all_gaps','constant','synonymous','single_char','random') DEFAULT NULL,

#  FOREIGN KEY (node_id) REFERENCES protein_tree_node(node_id) ON DELETE CASCADE ON UPDATE CASCADE,
#  FOREIGN KEY (parameter_set_id) REFERENCES parameter_set(parameter_set_id) ON DELETE CASCADE ON UPDATE CASCADE,

  UNIQUE (node_id,data_id,parameter_set_id,aln_position),
  KEY(ncod),
  KEY(node_id),
  KEY(data_id)
);


-- ----------------------------------------------------------------------------------
--
-- Table structure for table 'lr_index_offset'
--
-- overview:
--   Used to store the current maximum left right index for a given table. Table
--   name is unique and lr_index should be equal to the SQL statement
--   select max(right_index) from table_name
-- semantics:
--   table_name      -- name of the table this lr_index corresponds to
--   lr_index        -- max right index for the given table

CREATE TABLE IF NOT EXISTS lr_index_offset (
	table_name  varchar(64) NOT NULL,
	lr_index    int(10) unsigned NOT NULL,

	PRIMARY KEY (table_name)
) COLLATE=latin1_swedish_ci;

-- ----------------------------------------------------------------------------------
--
-- Table structure for table 'protein_tree_hmmprofile'
--
-- overview:
--   to allow nodes to have hmm profiles attached to them
-- semantics:
--    node_id                  -- the id of the root node associated with this hmm profile
--    type                     -- type of hmm profile (eg: 'hmmls','hmmfs','hmms','hmmsw')
--    hmmprofile               -- foreign key from method_link_species_set table

CREATE TABLE IF NOT EXISTS protein_tree_hmmprofile (
  node_id                     int(10) unsigned NOT NULL,
  type                        varchar(40) DEFAULT '' NOT NULL,
  hmmprofile                  mediumtext,

  KEY(node_id),
  UNIQUE KEY type_node_id (type, node_id)
) COLLATE=latin1_swedish_ci;


# Auto add schema version to database (this will override whatever hive puts there)
REPLACE INTO meta (species_id, meta_key, meta_value) VALUES (NULL, "schema_version", "60");

#Add schema type
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, "schema_type", "compara");
