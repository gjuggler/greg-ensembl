#!/usr/bin/perl

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::OntologyDBAdaptor;

Bio::EnsEMBL::Registry->load_registry_from_multiple_dbs(
  {
    -host => 'ens-livemirror',
    -user => 'ensro',
  });
my $compara_dba = Bio::EnsEMBL::Registry->get_DBAdaptor('Multi', 'Compara');
my $ontology_dba = Bio::EnsEMBL::Registry->get_DBAdaptor('Multi', 'Ontology');

my $mba = $compara_dba->get_MemberAdaptor;
my $goa = $ontology_dba->get_OntologyTermAdaptor;
my @peps = @{$mba->fetch_all_by_source_taxon('ENSEMBLPEP', 9606)};

foreach my $pep (@peps) {
  my $gene = $pep->get_Gene;
  my @dblinks = @{$gene->get_all_DBLinks('GO')};
  while (my $dblink = shift @dblinks){
    my $term = $goa->fetch_by_accession($dblink->display_id);    
    print join("\t", $pep->stable_id, $term->accession, @{$dblink->get_all_linkage_types}[0])."\n";
  }
}
