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
  meta_value                  mediumtext NOT NULL,
	
  PRIMARY   KEY (meta_id),
  UNIQUE    KEY species_key_value_idx (species_id, meta_key, meta_value(40)),
  KEY species_value_idx (species_id, meta_value(40))

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
  taxon_id                    int(10) unsigned DEFAULT '0', # KF taxon.taxon_id
  name                        varchar(40) DEFAULT '' NOT NULL,
  assembly                    varchar(100) DEFAULT '' NOT NULL,
  assembly_default            tinyint(1) DEFAULT 1,
  genebuild                   varchar(100) DEFAULT '' NOT NULL,
  locator                     varchar(255),

  FOREIGN KEY (taxon_id) REFERENCES ncbi_taxa_node(taxon_id) ON DELETE SET NULL ON UPDATE CASCADE,

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
  chr_name                    varchar(40),
  chr_start                   int(10),
  chr_end                     int(10),
  chr_strand                  tinyint(1) NOT NULL,
  display_label               varchar(128) default NULL,

  #FOREIGN KEY (taxon_id) REFERENCES ncbi_taxa_node(taxon_id),
  #FOREIGN KEY (genome_db_id) REFERENCES genome_db(genome_db_id),
  FOREIGN KEY (sequence_id) REFERENCES sequence(sequence_id) ON DELETE SET NULL ON UPDATE CASCADE,
  FOREIGN KEY (gene_member_id) REFERENCES member(member_id) ON DELETE SET NULL ON UPDATE CASCADE,

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

  FOREIGN KEY (node_id) REFERENCES protein_tree_node(node_id) ON DELETE CASCADE ON UPDATE CASCADE,

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

  FOREIGN KEY (node_id) REFERENCES protein_tree_node(node_id) ON DELETE SET NULL ON UPDATE CASCADE,

  UNIQUE tag_node_id (node_id, tag),
  KEY (node_id),
  KEY (tag)
) COLLATE=latin1_swedish_ci;

CREATE TABLE IF NOT EXISTS `parameter_set` (
  `parameter_set_id` tinyint unsigned NOT NULL AUTO_INCREMENT,
  `parameter_name` varchar(40) DEFAULT NULL,
  `parameter_value` mediumtext,
#  PRIMARY KEY (parameter_set_id),
  UNIQUE (parameter_set_id,parameter_name)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


#CREATE TABLE IF NOT EXISTS parameter_set_member (
#  parameter_set_member_id int(10) unsigned NOT NULL AUTO_INCREMENT,
#  parameter_set_id int(10) unsigned NOT NULL,
#  node_id int(10) unsigned NOT NULL,
#  parameter_string MEDIUMTEXT NOT NULL,
#  PRIMARY KEY (parameter_set_member_id),
#  UNIQUE KEY parameter_set_node_id (parameter_set_id,node_id)
#) ENGINE=InnoDB;


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

CREATE TABLE IF NOT EXISTS `sitewise_omega` (
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

  UNIQUE (node_id,parameter_set_id,aln_position),
  KEY(ncod)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 
  /*!50100 
	PARTITION BY KEY (parameter_set_id) 
	PARTITIONS 10
  */
;


CREATE TABLE IF NOT EXISTS `sitewise_pfam` (
  `node_id` int unsigned NOT NULL,
  `parameter_set_id` tinyint unsigned NOT NULL,
  `aln_position` mediumint unsigned NOT NULL,
  `pfam_id` char(12) NOT NULL DEFAULT '',

  `pf_position` mediumint unsigned NOT NULL,
  `score` float(10,5) DEFAULT NULL,

#  FOREIGN KEY (node_id) REFERENCES protein_tree_node(node_id) ON DELETE CASCADE ON UPDATE CASCADE,
#  FOREIGN KEY (parameter_set_id) REFERENCES parameter_set(parameter_set_id) ON DELETE CASCADE ON UPDATE CASCADE,

  UNIQUE (node_id,parameter_set_id,aln_position,pfam_id),
  KEY `pfam_id` (`pfam_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 
;

CREATE table if not exists sitewise_tag (
  node_id int unsigned NOT NULL,
  parameter_set_id tinyint unsigned NOT NULL,
  aln_position mediumint unsigned NOT NULL,

  tag varchar(16) NOT NULL default '',  # Tag
  value varchar(16) default NULL,       # Value, i.e. type of modification, secondary structure, accessibility value
  source varchar(8) default NULL,      # Source, i.e. PDB, Uniprot, etc.

#  FOREIGN KEY (node_id) REFERENCES protein_tree_node(node_id) ON DELETE CASCADE ON UPDATE CASCADE,
#  FOREIGN KEY (parameter_set_id) REFERENCES parameter_set(parameter_set_id) ON DELETE CASCADE ON UPDATE CASCADE,

  UNIQUE (node_id,parameter_set_id,aln_position,tag),
  KEY (source)
) ENGINE=InnoDB
  /*!50100
     PARTITION BY KEY (tag)
     PARTITIONS 5
  */
;

CREATE TABLE IF NOT EXISTS sitewise_genome (
  node_id                     int unsigned NOT NULL,
  parameter_set_id            tinyint unsigned NOT NULL,
  aln_position                mediumint unsigned NOT NULL,
  member_id                   int unsigned NOT NULL,
  chr_name                    char(10) NOT NULL,
  chr_start                   int unsigned NOT NULL,
  chr_end                     int unsigned NOT NULL,

#  FOREIGN KEY (node_id) REFERENCES protein_tree_node(node_id) ON DELETE CASCADE ON UPDATE CASCADE,
#  FOREIGN KEY (member_id) REFERENCES member(member_id) ON DELETE CASCADE ON UPDATE CASCADE,
#  FOREIGN KEY (parameter_set_id) REFERENCES parameter_set(parameter_set_id) ON DELETE CASCADE ON UPDATE CASCADE,

  UNIQUE (node_id,chr_name,aln_position,parameter_set_id,member_id),
  key (chr_name)
) ENGINE=InnoDB
  /*!50100 
	PARTITION BY KEY (chr_name)
        PARTITIONS 8
  */
;

CREATE TABLE IF NOT EXISTS node_set (
  node_set_id tinyint unsigned NOT NULL AUTO_INCREMENT,
  name varchar(40) NOT NULL DEFAULT '',
  allows_overlap tinyint unsigned NOT NULL DEFAULT 0,
  PRIMARY KEY (node_set_id),
  UNIQUE (node_set_id,name)
) ENGINE=InnoDB;

CREATE TABLE IF NOT EXISTS node_set_member (
  node_set_id tinyint unsigned NOT NULL,
  node_id int unsigned NOT NULL,

  FOREIGN KEY (node_set_id) REFERENCES node_set(node_set_id) ON DELETE CASCADE ON UPDATE CASCADE,
  UNIQUE (node_set_id,node_id)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


CREATE TABLE IF NOT EXISTS go_terms (
  go_term_id int unsigned NOT NULL AUTO_INCREMENT,
  `node_id` int unsigned default 0,
  `member_id` int unsigned NOT NULL default 0,
  `source_taxon` int unsigned NOT NULL default 0,
  `stable_id` varchar(40) default NULL,
  `go_term` varchar(20) default NULL,
  `evidence_code` varchar(10) default NULL,
  PRIMARY KEY (go_term_id),
  UNIQUE KEY `member_go_term` (`go_term`,`member_id`,`source_taxon`),
  KEY `node_id` (`node_id`),
  KEY `stable_id` (`stable_id`),
  KEY `member_id` (`member_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Add schema version to database
REPLACE INTO meta (meta_key, meta_value) VALUES ("schema_version", "57");

DROP FUNCTION IF EXISTS dubious_dup;
CREATE FUNCTION dubious_dup (n int(20))
RETURNS int(20)
READS SQL DATA
RETURN (SELECT count(*) FROM protein_tree_tag where node_id=n AND tag="dubious_duplication" AND value="1");

DROP FUNCTION IF EXISTS  leaf_count;
CREATE FUNCTION leaf_count (n int(20))
RETURNS int(20)
READS SQL DATA
RETURN (SELECT COUNT(*) 
   FROM protein_tree_node n1,protein_tree_node n2,protein_tree_member m 
   WHERE n2.left_index BETWEEN n1.left_index and n1.right_index
# GJ 2009-03-03 - until we redo the leftright-indexing, keep this in:
   AND n2.root_id=n1.root_id
   AND n1.node_id=n
   AND m.node_id=n2.node_id
);

DROP FUNCTION IF EXISTS  node_count;
CREATE FUNCTION node_count (n int(20))
RETURNS int(20)
READS SQL DATA
RETURN (SELECT COUNT(*) 
   FROM protein_tree_node n1,protein_tree_node n2
   WHERE n2.left_index BETWEEN n1.left_index and n1.right_index
# GJ 2009-03-03 - until we redo the leftright-indexing, keep this in:
   AND n2.root_id=n1.root_id
   AND n1.node_id=n
);


DROP FUNCTION IF EXISTS  psc_count;
CREATE FUNCTION psc_count (n int(20),pset int(20))
RETURNS int(20)
READS SQL DATA
RETURN (SELECT COUNT(*) 
   FROM sitewise_aln sa
   WHERE sa.parameter_set_id=pset
   AND sa.node_id=n
   AND sa.type='positive4'
   GROUP BY node_id
);

DROP FUNCTION IF EXISTS  num_human_genes;
CREATE FUNCTION num_human_genes (n int(20))
RETURNS int(20)
READS SQL DATA
RETURN (SELECT COUNT(*)
   FROM protein_tree_node n, protein_tree_node n2, protein_tree_member ptm, member m
   WHERE n.node_id=n
   AND n2.left_index BETWEEN n.left_index and n.right_index
   AND ptm.node_id=n2.node_id
   AND m.member_id=ptm.member_id
   AND m.taxon_id = 9606
);

DROP FUNCTION IF EXISTS human_protein_list;
CREATE FUNCTION human_protein_list (n int(20))
RETURNS TEXT
READS SQL DATA
RETURN (SELECT GROUP_CONCAT(m.stable_id SEPARATOR ',')
   FROM protein_tree_node n, protein_tree_node n2, protein_tree_member ptm, member m
   WHERE n.node_id=n
   AND n2.left_index BETWEEN n.left_index and n.right_index
   AND ptm.node_id=n2.node_id
   AND m.member_id=ptm.member_id
   AND m.taxon_id = 9606
);

DROP FUNCTION IF EXISTS human_gene_list;
CREATE FUNCTION human_gene_list (n int(20))
RETURNS TEXT
READS SQL DATA
RETURN (SELECT GROUP_CONCAT(gm.stable_id SEPARATOR ',')
   FROM protein_tree_node n, protein_tree_node n2, protein_tree_member ptm, member m, member gm
   WHERE n.node_id=n
   AND n2.left_index BETWEEN n.left_index and n.right_index
   AND ptm.node_id=n2.node_id
   AND m.member_id=ptm.member_id
   AND m.taxon_id = 9606
   AND m.gene_member_id=gm.member_id
);


DROP FUNCTION IF EXISTS  human_protein_under_node;
CREATE FUNCTION human_protein_under_node (n int(20))
RETURNS varchar(40)
READS SQL DATA
RETURN (SELECT m.stable_id
   FROM protein_tree_node n, protein_tree_node n2, protein_tree_member ptm, member m
   WHERE n.node_id=n
   AND n2.left_index BETWEEN n.left_index and n.right_index
   AND ptm.node_id=n2.node_id
   AND m.member_id=ptm.member_id
   AND m.taxon_id = 9606
   LIMIT 1
);

DROP FUNCTION IF EXISTS  human_gene_under_node;
CREATE FUNCTION human_gene_under_node (n int(20))
RETURNS varchar(40)
READS SQL DATA
RETURN (SELECT gm.stable_id
   FROM protein_tree_node n, protein_tree_node n2, protein_tree_member ptm, member m, member gm
   WHERE n.node_id=n
   AND n2.left_index BETWEEN n.left_index and n.right_index
   AND ptm.node_id=n2.node_id
   AND m.member_id=ptm.member_id
   AND gm.member_id=m.gene_member_id
   AND m.taxon_id = 9606
   LIMIT 1
);

DROP PROCEDURE IF EXISTS sequences_under_node;
DELIMITER //
CREATE PROCEDURE sequences_under_node(IN n INT(20))
READS SQL DATA
   BEGIN 
   	SELECT s.sequence_id,m.stable_id,s.sequence
	   FROM protein_tree_node n, protein_tree_node n2, protein_tree_member ptm, member m, sequence s
	   WHERE n.node_id=n
	   AND n2.left_index BETWEEN n.left_index and n.right_index
	   AND ptm.node_id=n2.node_id
	   AND m.member_id=ptm.member_id
	   AND s.sequence_id = m.sequence_id
	;
   END //
DELIMITER ;


DROP PROCEDURE IF EXISTS genes_under_node;
DELIMITER //
CREATE PROCEDURE genes_under_node(IN n INT(20))
READS SQL DATA
   BEGIN 
   	SELECT distinct gm.stable_id
	   FROM protein_tree_node n, protein_tree_node n2, protein_tree_member ptm, member m, member gm
	   WHERE n.node_id=n
	   AND n2.left_index BETWEEN n.left_index and n.right_index
	   AND ptm.node_id=n2.node_id
	   AND m.member_id=ptm.member_id
	   AND gm.member_id=m.gene_member_id
	;
   END //
DELIMITER ;

DROP PROCEDURE IF EXISTS members_under_node;
DELIMITER //
CREATE PROCEDURE members_under_node(IN n INT(20))
READS SQL DATA
   BEGIN 
   	SELECT distinct ptm.node_id,m.member_id,m.sequence_id,m.stable_id,ptm.cigar_line
	   FROM protein_tree_node n, protein_tree_node n2, protein_tree_member ptm, member m
	   WHERE n.node_id=n
	   AND n2.left_index BETWEEN n.left_index and n.right_index
	   AND ptm.node_id=n2.node_id
	   AND m.member_id=ptm.member_id
	   AND m.source_name="ENSEMBLPEP"
	;
   END //
DELIMITER ;

DROP FUNCTION IF EXISTS  sitewise_node_for_peptide;
CREATE FUNCTION sitewise_node_for_peptide (stable_id VARCHAR(40),parameter_set_id INT(20))
RETURNS INT(20)
READS SQL DATA
RETURN (SELECT IFNULL(sa.node_id,0)
   FROM sitewise_aln sa, protein_tree_member ptm, member m, protein_tree_node ptn, protein_tree_node ptn2
   WHERE m.stable_id=stable_id
   AND m.source_name!="ENSEMBLGENE"
   AND ptm.member_id=m.member_id
   AND ptn.node_id=ptm.node_id
   AND ptn.left_index BETWEEN ptn2.left_index AND ptn2.right_index
   AND sa.node_id=ptn2.node_id
   AND sa.parameter_set_id=parameter_set_id
   LIMIT 1
);

DROP FUNCTION IF EXISTS  bak_sitewise_node_for_peptide;
CREATE FUNCTION bak_sitewise_node_for_peptide (stable_id VARCHAR(40),parameter_set_id INT(20))
RETURNS INT(20)
READS SQL DATA
   RETURN (SELECT sa.node_id
       FROM bak_sitewise_aln sa, protein_tree_member ptm, member m, protein_tree_node ptn, protein_tree_node ptn2
       WHERE m.stable_id=stable_id
       AND m.source_name="ENSEMBLPEP"
       AND ptm.member_id=m.member_id
       AND ptn.node_id=ptm.node_id
       AND ptn.left_index BETWEEN ptn2.left_index AND ptn2.right_index
       AND sa.node_id=ptn2.node_id
       AND sa.parameter_set_id=parameter_set_id
   LIMIT 1
);

DROP FUNCTION IF EXISTS  peptide_for_gene;
CREATE FUNCTION peptide_for_gene (stable_id VARCHAR(40))
RETURNS VARCHAR(40)
READS SQL DATA
RETURN (SELECT pep.stable_id
   FROM protein_tree_member ptm, member pep, member gene
   WHERE ptm.member_id=pep.member_id
   AND pep.gene_member_id = gene.member_id
   AND gene.stable_id=stable_id
   LIMIT 1
);


DROP FUNCTION IF EXISTS  sitewise_node_for_gene;
CREATE FUNCTION sitewise_node_for_gene (stable_id VARCHAR(40),parameter_set_id INT(20))
RETURNS INT(20)
READS SQL DATA
RETURN (sitewise_node_for_peptide(peptide_for_gene(stable_id),parameter_set_id));

DROP FUNCTION IF EXISTS  root_node_for_peptide;
CREATE FUNCTION root_node_for_peptide (stable_id VARCHAR(40))
RETURNS INT(20)
READS SQL DATA
RETURN (SELECT ptn.node_id
   FROM protein_tree_node ptn, protein_tree_node ptn2, protein_tree_member ptm, member m WHERE
   m.stable_id=stable_id AND
   ptm.member_id=m.member_id AND
   ptn2.node_id = ptm.node_id AND
   ptn2.left_index BETWEEN ptn.left_index AND ptn.right_index AND
   ptn.parent_id=1
   LIMIT 1
);

DROP FUNCTION IF EXISTS  num_dups_under_node;
CREATE FUNCTION num_dups_under_node (n INT(20))
RETURNS INT(20)
READS SQL DATA
RETURN (SELECT count(*)
   FROM protein_tree_node ptn1, protein_tree_node ptn2, protein_tree_tag t1
   WHERE ptn1.node_id=n
   AND ptn2.left_index BETWEEN ptn1.left_index AND ptn1.right_index
   AND t1.node_id=ptn2.node_id
   AND t1.tag="Duplication"
   AND t1.value > 0
);

DROP FUNCTION IF EXISTS  avg_omega_for_node;
CREATE FUNCTION avg_omega_for_node (n INT(20), pset INT(20))
RETURNS FLOAT(10,5)
READS SQL DATA
RETURN (SELECT AVG(omega)
   FROM sitewise_aln sa
   WHERE sa.node_id=n
   AND sa.parameter_set_id=pset
   AND sa.omega_upper > sa.omega
   AND sa.ncod >= 4
   AND sa.note != 'random'
   GROUP BY sa.node_id
);

DROP FUNCTION IF EXISTS  avg_lrt_for_node;
CREATE FUNCTION avg_lrt_for_node (n INT(20), pset INT(20))
RETURNS FLOAT(10,5)
READS SQL DATA
RETURN (SELECT AVG(lrt_stat)
   FROM sitewise_aln sa
   WHERE sa.node_id=n
   AND sa.parameter_set_id=pset
   AND sa.omega_upper > sa.omega
   AND sa.ncod >= 4
   AND sa.note != 'random'
   GROUP BY sa.node_id
);


DROP FUNCTION IF EXISTS one_sided_gappiness;
DELIMITER //
CREATE FUNCTION one_sided_gappiness (n_id INT(20), aln_pos INT(20), window_size INT(20), direction INT(20))
RETURNS FLOAT(10,5)
READS SQL DATA
BEGIN
IF direction> 0 THEN
  RETURN (SELECT  1 - AVG(LENGTH(REPLACE(SUBSTR(cigar_line,aln_pos,window_size),"-",""))/window_size)
  FROM protein_tree_member_score ptm, protein_tree_node ptn, protein_tree_node ptn2
  WHERE ptn2.node_id=ptm.node_id and ptn2.left_index between ptn.left_index and ptn.right_index
   AND ptn.node_id=n_id );
ELSE
  RETURN (SELECT 1 - AVG(LENGTH(REPLACE(SUBSTR(cigar_line,aln_pos-window_size,window_size),"-",""))/window_size) 
  FROM protein_tree_member_score ptm, protein_tree_node ptn, protein_tree_node ptn2
  WHERE ptn2.node_id=ptm.node_id and ptn2.left_index between ptn.left_index and ptn.right_index
   AND ptn.node_id=n_id );
END IF;
END //
DELIMITER ;

DROP FUNCTION IF EXISTS gappiness_around_site;
DELIMITER //
CREATE FUNCTION gappiness_around_site (n_id INT(20), aln_pos INT(20), window_size INT(20))
RETURNS FLOAT(10,5)
READS SQL DATA
BEGIN
	DECLARE val_right FLOAT(10,5);
	DECLARE val_left FLOAT(10,5);
	SET val_right = one_sided_gappiness(n_id,aln_pos,window_size/2,1);
	SET val_left = one_sided_gappiness(n_id,aln_pos,window_size/2,-1);
	IF val_right > val_left THEN RETURN(val_right);
	ELSE RETURN (val_left);
	END IF;
END //
DELIMITER ;

DROP FUNCTION IF EXISTS gappiness;
CREATE FUNCTION gappiness (n_id INT(20))
RETURNS FLOAT(10,5)
READS SQL DATA
     RETURN (SELECT  1 - AVG(LENGTH(REPLACE(cigar_line,"-",""))/length(cigar_line))
     FROM protein_tree_member_score ptm, protein_tree_node ptn, protein_tree_node ptn2
     WHERE ptn2.node_id=ptm.node_id and ptn2.left_index between ptn.left_index and ptn.right_index
       AND ptn.node_id=n_id
);


DROP FUNCTION IF EXISTS aln_length;
DELIMITER //
CREATE FUNCTION aln_length (n_id INT(20))
RETURNS INT(20)
READS SQL DATA
BEGIN
        RETURN (SELECT length(cigar_line) FROM protein_tree_member_score ptm, protein_tree_node ptn, protein_tree_node ptn2
                WHERE ptn2.node_id=ptm.node_id AND ptn2.left_index BETWEEN ptn.left_index AND ptn.right_index AND ptn.node_id=n_id
                LIMIT 1);
END //
DELIMITER ;

DROP FUNCTION IF EXISTS within_gblocks;
CREATE FUNCTION within_gblocks (node_id INT(20), aln_pos INT(20))
RETURNS BOOLEAN
READS SQL DATA
RETURN (
       SELECT count(*) FROM
       gblocks_aln ga
       WHERE ga.node_id=node_id AND
       aln_pos BETWEEN ga.aln_start AND ga.aln_end
       LIMIT 1
);

DROP FUNCTION IF EXISTS good_sitewise_pos;
CREATE FUNCTION good_sitewise_pos (node_id INT(20), aln_pos INT(20), parameter_set_id INT(20))
RETURNS BOOLEAN
READS SQL DATA
RETURN (
       SELECT count(*) FROM
       sitewise_aln sa
       WHERE sa.node_id=node_id AND
       sa.aln_position = aln_pos AND
       sa.parameter_set_id=parameter_set_id AND
       sa.ncod >= 4 AND
       sa.omega_upper > sa.omega AND
       sa.note != 'random'
       LIMIT 1
);