#!/usr/bin/env perl
use warnings;
use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::Greg::Hive::ComparaHiveLoaderUtils;

#my ($url) = undef;
#GetOptions('url=s' => \$url);

my $url = 'mysql://ensadmin:ensembl@ens-research:3306/gj1_gor_58';
my $clean = 1;

my $h = new Bio::Greg::Hive::ComparaHiveLoaderUtils;
$h->init($url);

# Clean up our mess.
if ($clean) {
  $h->clean_compara_analysis_tables;
  $h->clean_hive_tables;
}

# Define parameters (species sets, filtering options, etc).
parameter_sets();

# Node sets.
node_sets();

### Genomewide omegas track.
tree_stats();
split_by_subtrees();
paml_omegas();
### End genomewide omegas.

### Duplication counts track.
collect_duplications();
### End duplication counts.

### Branch-model tests track.
filter_one_to_one();
count_sites();
count_sites_outgroup();
likelihood_tests();
likelihood_stats();
collect_likelihoods();
collect_go();
### End branch-model tests.

### Collect everything at the end. Run this one off-hive with an empty input_id
output_all();

# Connect the dots.
### Genomewide omega track.
$h->connect_analysis("NodeSets","TreeStats",1);
$h->connect_analysis("TreeStats","SplitBySubtrees",1);
$h->connect_analysis("SplitBySubtrees","PamlOmegas",1);
###
### Duplication counts track.
$h->connect_analysis("NodeSets","CollectDuplications",1);
###
### Branch models track.
$h->connect_analysis("NodeSets","FilterOneToOne",1);
$h->connect_analysis("FilterOneToOne","CountSites",1);
$h->connect_analysis("FilterOneToOne","CountSitesOutgroup",1);
$h->connect_analysis("FilterOneToOne","LikelihoodTests",1);
$h->connect_analysis("LikelihoodTests","LikelihoodStats",1);
$h->connect_analysis("LikelihoodStats","CollectLikelihoods",1);
$h->connect_analysis("NodeSets","CollectGo",1);

add_all_nodes();

sub add_genes {
  my @genes = ();

  open INPUT, "<candidate_genes.txt";
  my @lines = <INPUT>;
  close INPUT;
  #@genes =  @lines;
  @genes = ('PAX2','CYBB','MB','BRCA2');

  $h->add_genes_to_analysis("NodeSets",\@genes);
}

sub add_all_nodes {
  my $cmd = "SELECT node_id FROM protein_tree_node WHERE parent_id=1;";
  my @nodes = _select_node_ids($cmd);
  my $logic_name = "NodeSets";
  _add_nodes_to_analysis($logic_name,\@nodes);
}

sub parameter_sets {
  my $params;
  my $base_params={};

  my @all_arr = clade_taxon_ids();
  my @mammals_arr = clade_taxon_ids("Eutheria");
  my @primates_arr = clade_taxon_ids("Primates");
  my @homininae_arr = clade_taxon_ids("Homininae");
  my @hominidae_arr = clade_taxon_ids("Hominidae");

  my $everything = join(",",clade_taxon_ids());

  my $primates = join(",",clade_taxon_ids("Primates"));
  my $homininae = join(",",clade_taxon_ids("Homininae"));
  my $hominidae = join(",",clade_taxon_ids("Hominidae"));
  my $non_homininae = join(",",subtract(\@primates_arr,\@homininae_arr));
  my $non_hominidae = join(",",subtract(\@primates_arr,\@hominidae_arr));
  my $non_gorilla = join(",",subtract(\@hominidae_arr,[9593]));
  my $simiiformes = join(",",clade_taxon_ids("Simiiformes"));
  my $haplorrhini = join(",",clade_taxon_ids("Haplorrhini"));

  my $glires = join(",",clade_taxon_ids("Glires"));
  my $laurasiatheria = join(",",clade_taxon_ids("Laurasiatheria"));
  my $afrotheria = join(",",clade_taxon_ids("Afrotheria"));

  my $mammals = join(",",clade_taxon_ids("Eutheria"));

  my $non_primates = join(",",subtract(\@mammals_arr,\@primates_arr));
 
 $params = {
    parameter_set_name => "Homininae",
    parameter_set_shortname => 'hmn',
    keep_species => $homininae
  };
  $h->add_parameter_set($params);

  $params = {
    parameter_set_name => "Hominidae",
    parameter_set_shortname => 'hmd',
    keep_species => $hominidae
  };
  $h->add_parameter_set($params);

  $params = {
    parameter_set_name => "NonGorillaHominidae",
    parameter_set_shortname => 'nongor',
    keep_species => $non_gorilla
  };
  $h->add_parameter_set($params);

  $params = {
    parameter_set_name => "NonHominidPrimates",
    parameter_set_shortname => 'nonhom',
    keep_species => $non_hominidae
  };
  $h->add_parameter_set($params);

  $params = {
    parameter_set_name => "Haplorrhini",
    parameter_set_shortname => 'hpl',
    keep_species => $haplorrhini
  };
  $h->add_parameter_set($params);

  $params = {
    parameter_set_name => "Primates",
    parameter_set_shortname => 'p',
    keep_species => $primates
  };
  $h->add_parameter_set($params);
 
  $params = {
    parameter_set_name => "Glires",
    parameter_set_shortname => 'g',
    keep_species => $glires
  };
  $h->add_parameter_set($params);

  $params = {
    parameter_set_name => "Laurasiatheria",
    parameter_set_shortname => 'l',
    keep_species => $laurasiatheria
  };
  $h->add_parameter_set($params);

  $params = {
    parameter_set_name => "Afrotheria",
    parameter_set_shortname => 'a',
    keep_species => $afrotheria
  };
  $h->add_parameter_set($params);

  $params = {
    parameter_set_name => "Mammals",
    parameter_set_shortname => 'm',
    keep_species => $mammals
  };
  $h->add_parameter_set($params);

  $params = {
    parameter_set_name => "NonPrimateMammals",
    parameter_set_shortname => 'nonprm',
    keep_species => $non_primates
  };
  $h->add_parameter_set($params);

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

sub count_sites {
  my $logic_name = "CountSites";
  my $module = "Bio::Greg::Gorilla::CountSites";
  my $params = {
    gorilla_count_species => '9593,9598,9606',
    counts_sites_table => 'counts_sites',
    counts_genes_table => 'counts_genes'
  };
  $h->create_analysis($logic_name,$module,$params,50,1);
}

sub count_sites_outgroup {
  my $logic_name = "CountSitesOutgroup";
  my $module = "Bio::Greg::Gorilla::CountSites";
  my $params = {
    gorilla_count_species => '9593,9598,9606,9600',
    counts_sites_table => 'outgroup_sites',
    counts_genes_table => 'outgroup_genes'
  };
  $h->create_analysis($logic_name,$module,$params,50,1);
}

sub likelihood_tests {
  my $logic_name = "LikelihoodTests";
  my $module = "Bio::Greg::Gorilla::LikelihoodTests";
  my $params = {
  };
  $h->create_analysis($logic_name,$module,$params,300,1);
}

sub likelihood_stats {
  my $logic_name = "LikelihoodStats";
  my $module = "Bio::Greg::Gorilla::CollectLikelihoodStats";
  my $params = {
  };
  $h->create_analysis($logic_name,$module,$params,50,1);
}

sub collect_go {
  my $logic_name = "CollectGo";
  my $module = "Bio::Greg::Hive::CollectGO";
  my $params = {
    go_taxon_ids => '9606,9593',
  };
  $h->create_analysis($logic_name,$module,$params,200,1);
}

sub collect_likelihoods {
  my $logic_name = "CollectLikelihoods";
  my $module = "Bio::Greg::Gorilla::CollectGorillaStats";
  my $params = {
  };
  $h->create_analysis($logic_name,$module,$params,80,1);
}

sub output_all {
  my $logic_name = "OutputAll";
  my $module = "Bio::Greg::Gorilla::OutputGorillaData";
  my $params = {
  };
  $h->create_analysis($logic_name,$module,$params,1,1);
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
