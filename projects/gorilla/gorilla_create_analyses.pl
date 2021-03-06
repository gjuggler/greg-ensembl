#!/usr/bin/env perl
use warnings;
use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::Greg::Hive::ComparaHiveLoaderUtils;

my $url = 'mysql://ensadmin:ensembl@ens-research:3306/gj1_gor2';
my $clean = 1;

my $h = new Bio::Greg::Hive::ComparaHiveLoaderUtils;
$h->init($url);

# Clean up our mess.
if ($clean) {
  $h->clean_compara_analysis_tables;
  $h->clean_hive_tables;
}

# Node sets.
node_sets();

### Genomewide omegas track.
tree_stats();
split_by_subtrees();
paml_omegas();
### End genomewide omegas.

### Branch-model tests track.
filter_one_to_one();
#lnl_c();
#lnl_p();
lnl_m();
#lnl_mp();
#lnl_nf();
### End branch-model tests.

### Collect everything at the end. Run this one off-hive with an empty input_id
output_gorilla_data();

map_gorilla_subs();
$h->add_inputs_to_analysis("MapGorillaSubs",[{}]);  
$h->connect_analysis("MapGorillaSubs","MapGorillaSubs",99);

### Genomewide omega track.
# Add a single trigger job to SplitBySubtrees.
$h->add_inputs_to_analysis("SplitBySubtrees",[{}]);  
$h->wait_for("SplitBySubtrees",["NodeSets"]);
$h->connect_analysis("SplitBySubtrees","PamlOmegas",1);

### Branch models track.
$h->connect_analysis("NodeSets","FilterOneToOne",1);
#$h->connect_analysis("FilterOneToOne","lnl_c",1);
#$h->connect_analysis("FilterOneToOne","lnl_p",1);
$h->connect_analysis("FilterOneToOne","lnl_m",1);
#$h->connect_analysis("FilterOneToOne","lnl_mp",1);
#$h->connect_analysis("FilterOneToOne","lnl_mf",1);

#add_all_nodes();
add_genes();

sub add_genes {
  my @genes = ();

  open INPUT, "< gorilla_genes.txt";
  my @lines = <INPUT>;
  close INPUT;
  @genes =  @lines;

  $h->add_genes_to_analysis("NodeSets",\@genes);
}

sub add_all_nodes {
  my $cmd = "SELECT node_id FROM protein_tree_node WHERE parent_id=1;";
  my @nodes = _select_node_ids($cmd);
  my $logic_name = "NodeSets";
  _add_nodes_to_analysis($logic_name,\@nodes);
}

sub node_sets {
  my $logic_name = "NodeSets";
  my $module = "Bio::Greg::Hive::NodeSets";
  my $params = {
    flow_node_set => 'Primates'
  };
  my $id = $h->create_analysis($logic_name,$module,$params,150,1);
}

###
sub tree_stats {
  my $logic_name = "TreeStats";
  my $module = "Bio::Greg::Hive::CollectTreeStats";
  my $params = {
  };
  my $id = $h->create_analysis($logic_name,$module,$params,150,1);
}
sub split_by_subtrees {
  my $logic_name = "SplitBySubtrees";
  my $module = "Bio::Greg::Hive::SplitBySubtrees";
  my $params = {
    seed_species => '9606',
    out_to => 'Mammalia'
  };
  $h->create_analysis($logic_name,$module,$params,100,1);
}
sub paml_omegas {
  my $logic_name = "PamlOmegas";
  my $module = "Bio::Greg::Hive::GenomewideOmegas";
  my $params = {
    method => 'paml',
    foreground_species => '9606,9593'
  };
  $h->create_analysis($logic_name,$module,$params,100,1);
}
###

###
sub collect_duplications {
  my $logic_name = "CollectDuplications";
  my $module = "Bio::Greg::Gorilla::CollectDuplicationStats";
  my $params = {
    collect_duplication_species => '9606,9598,9593,9600'
  };
  $h->create_analysis($logic_name,$module,$params,100,1);
}
###

###
sub filter_one_to_one {
  my $logic_name = "FilterOneToOne";
  my $module = "Bio::Greg::Hive::FilterOneToOneOrthologs";
  my $params = {
    one_to_one_taxon_list => '9606,9598,9593,9600,9544,9483'
  };
  $h->create_analysis($logic_name,$module,$params,150,1);
}

sub lnl_c {
  my $logic_name = "lnl_c";
  my $module = "Bio::Greg::Gorilla::LikelihoodTests";
  my $params = {
    output_table => 'lnl_c',
    aln_type => 'compara'
  };
  $h->create_analysis($logic_name,$module,$params,300,1);
}

sub lnl_p {
  my $logic_name = "lnl_p";
  my $module = "Bio::Greg::Gorilla::LikelihoodTests";
  my $params = {
    output_table => 'lnl_p',
    aln_type => 'genomic_primates',
  };
  $h->create_analysis($logic_name,$module,$params,300,1);
}

sub lnl_m {
  my $logic_name = "lnl_m";
  my $module = "Bio::Greg::Gorilla::LikelihoodTests";
  my $params = {
    output_table => 'lnl_m',
    aln_type => 'genomic_mammals',
  };
  $h->create_analysis($logic_name,$module,$params,300,1);
}

sub lnl_mp {
  my $logic_name = "lnl_mp";
  my $module = "Bio::Greg::Gorilla::LikelihoodTests";
  my $params = {
    output_table => 'lnl_mp',
    aln_type => 'genomic_mammals',
    species_taxon_ids => '9606, 9593, 9598, 9600, 9544, 9483, 9478, 30608, 30611'
  };
  $h->create_analysis($logic_name,$module,$params,300,1);
}

sub lnl_nf {
  my $logic_name = "lnl_nf";
  my $module = "Bio::Greg::Gorilla::LikelihoodTests";
  my $params = {
    output_table => 'lnl_nf',
    aln_type => 'genomic_mammals',
    quality_threshold => 0,
    likelihood_filter_substitution_runs => 0
  };
  $h->create_analysis($logic_name,$module,$params,300,1);
}

sub output_gorilla_data {
  my $logic_name = "OutputGorillaData";
  my $module = "Bio::Greg::Gorilla::OutputGorillaData";
  my $params = {
  };
  $h->create_analysis($logic_name,$module,$params,1,1);
}

sub map_gorilla_subs {
  my $logic_name = "MapGorillaSubs";
  my $module = "Bio::Greg::Gorilla::MapGorillaSubs";
  my $params = {
  };
  $h->create_analysis($logic_name,$module,$params,100,1);
}


########*********########
#-------~~~~~~~~~-------#
########*********########

# Subroutines to return a list of taxon IDs with specific features.
sub clade_taxon_ids {
  my $clade = shift || 1;

  my $dba = $h->dba;
  my @genomes = Bio::EnsEMBL::Compara::ComparaUtils->get_genomes_within_clade($dba,$clade);
  my @taxon_ids = map {$_->taxon_id} @genomes;
  return @taxon_ids;
}
sub coverage_taxon_ids {
  my $coverage = shift;
  
  my $dba = $h->dba;
  my @output;
  my @all_gdb = Bio::EnsEMBL::Compara::ComparaUtils->get_genomes_within_clade($dba,1);
  foreach my $gdb (@all_gdb) {
    # This is finicky: we need to call the "db_adaptor" method to get the Bio::EnsEMBL::DBSQL::DBAdaptor object, and then the meta container.
    my $meta = $gdb->db_adaptor->get_MetaContainer;
    my $str = @{$meta->list_value_by_key('assembly.coverage_depth')}[0];
    push @output, $gdb->taxon_id if ($str eq $coverage);
  }  return @output;
}
sub subtract {
  my $list_a = shift;
  my @remove_us = @_;
  my $hash;
  map {$hash->{$_}=1} @$list_a;
  foreach my $list_b (@remove_us) {
    map {delete $hash->{$_}} @$list_b;
}
  return keys %$hash;
}

sub _select_node_ids {
  my $cmd = shift;
  if ( !defined $cmd ) {
    $cmd = "SELECT node_id FROM protein_tree_node WHERE parent_id=1 AND root_id=1";
  }
  my $dbc = $h->dbc;
  my $sth = $dbc->prepare($cmd);
  $sth->execute();

  my $array_ref = $sth->fetchall_arrayref( [0] );
  my @node_ids = @{$array_ref};
  @node_ids =
    map { @{$_}[0] } @node_ids;    # Some weird mappings to unpack the numbers from the arrayrefs.
  $sth->finish;
  return @node_ids;
}

sub _add_nodes_to_analysis {
  my $logic_name = shift;
  my $node_id_arrayref = shift;

  my @node_ids = @{$node_id_arrayref};
  foreach my $node_id (@node_ids) {
    my $input_id = {
      node_id => $node_id
    };
    $h->add_job_to_analysis($logic_name,$input_id);
  }
}
