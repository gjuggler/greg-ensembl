# A table for storing longer strings for site-wise stuff.
CREATE TABLE IF NOT EXISTS sitewise_stats (
  node_id int unsigned NOT NULL,
  parameter_set_id tinyint unsigned NOT NULL,
  aln_position mediumint unsigned NOT NULL,

  stats MEDIUMTEXT default NULL,       # Value

#  FOREIGN KEY (node_id) REFERENCES protein_tree_node(node_id) ON DELETE CASCADE ON UPDATE CASCADE,
#  FOREIGN KEY (parameter_set_id) REFERENCES parameter_set(parameter_set_id) ON DELETE CASCADE ON UPDATE CASCADE,

  UNIQUE (node_id,parameter_set_id,aln_position)
) ENGINE=InnoDB;
