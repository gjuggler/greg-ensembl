create table if not exists aln_tag (
  node_id int(10) unsigned NOT NULL,
  aln_position int(10) unsigned NOT NULL,
  member_id int(10) unsigned NOT NULL,
  stable_id varchar(20) default NULL,
  source varchar(20) default NULL,
  tag varchar(16) NOT NULL default '',
  value varchar(16) default NULL,
  feature_position int(10) unsigned default NULL,
  FOREIGN KEY (node_id, aln_position) REFERENCES sitewise_aln(node_id, aln_position),
  UNIQUE tag_aln_id (node_id,aln_position,tag)
);