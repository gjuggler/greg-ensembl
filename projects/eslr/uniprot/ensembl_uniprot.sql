create table if not exists ensembl_uniprot (
  `stable_id` varchar(30) NOT NULL,
  `uniprot_id` varchar(30) NOT NULL,
  PRIMARY KEY `stable_uniprot_id` (`stable_id`,`uniprot_id`)
);