#!/usr/bin/env perl

use warnings;
use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::ComparaUtils;
use File::Path;
use File::Basename;

use Bio::Greg::Hive::HiveLoaderUtils;
use Bio::Greg::Hive::ComparaHiveLoaderUtils;

use Bio::TreeIO;

my $tree_file = "trees/mammals44tree.nh";
my $treeio = Bio::TreeIO->new(-file => $tree_file);
my $tree = $treeio->next_tree;

map {$_->id(lc($_->id))} $tree->leaves;

my @eutherian = qw(alpaca armadillo bushbaby cat chimp cow dog dolphin elephant guinea_pig hedgehog horse human kangaroo_rat megabat microbat mouse mouse_lemur pika rabbit rat rhesus rock_hyrax shrew sloth squirrel tarsier tenrec treeshrew);  
my @full = qw(chicken chimp cow dog fugu guinea_pig horse human lizard monodelphis mouse platypus rat rhesus stickleback tetraodon zebrafish zebra_finch);
my @full_mammal = qw(chimp cow dog guinea_pig horse human mouse rat rhesus);
my @hmrd = qw(human mouse rat dog);

my $subtree;
$subtree = $tree->slice_by_ids(@eutherian);
print "Eutherian:".$subtree->total_branch_length."\n";
print $subtree->max_distance_to_leaf."\n";
$subtree->to_file("trees/mammals_eutherian.nh");

$subtree = $tree->slice_by_ids(@full_mammal);
print "Full mammal:".$subtree->total_branch_length."\n";
print $subtree->max_distance_to_leaf."\n";
$subtree->to_file("trees/mammals_full_genomes.nh");

$subtree = $tree->slice_by_ids(@hmrd);
print "HMRD:".$subtree->total_branch_length."\n";
print $subtree->max_distance_to_leaf."\n";
$subtree->to_file("trees/mammals_hmrd.nh");

my ( $experiment_name ) = undef;
GetOptions(
  'experiment=s' => \$experiment_name
);
die("No experiment name given!") unless (defined $experiment_name);

# First create a database for the experiment.
my $mysql_base;
if (Bio::Greg::EslrUtils->is_ebi) {
  $mysql_base = 'mysql://slrsim:slrsim@mysql-greg.ebi.ac.uk:4134/';
} else {
  $mysql_base  = 'mysql://ensadmin:ensembl@ens-research/';
}
print "  MySQL base [$mysql_base]\n";

my $url = "${mysql_base}gj1_slrsim";
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new( -url => $url );
$dba->dbc->do("create database if not exists gj1_${experiment_name};");

$url = $mysql_base . "gj1_${experiment_name}";
my $h = new Bio::Greg::Hive::ComparaHiveLoaderUtils;
$h->init($url);

print "  cleaning tables...\n";
$h->init_compara_tables;

$h->clean_hive_tables;
$h->clean_compara_analysis_tables;
$h->clean_compara_tree_tables;
my @truncate_tables = qw^
      aln aln_scores omega
      sites genes merged slrsim_results
  ^;
$h->truncate_tables(\@truncate_tables);

load_trees();
load_tree();
slrsim();
calculate_results();
plots();

$h->connect_analysis( "LoadTrees",   "LoadTree", 1 );
$h->connect_analysis( "LoadTrees",   "CalculateResults", 2 );
$h->connect_analysis( "LoadTree",   "Slrsim" );
$h->wait_for("CalculateResults",["LoadTrees","LoadTree", "Slrsim"]);
$h->wait_for("Plots",["LoadTrees","LoadTree", "Slrsim", "CalculateResults"]);

sub load_trees {
  my $logic_name  = "LoadTrees";
  my $module      = "Bio::Greg::Slrsim::LoadTrees";
  my $params      = { experiment_name => $experiment_name };
  $h->create_analysis( $logic_name, $module, $params, 1, 1 );

  $h->add_job_to_analysis( "LoadTrees", {} );    # Add a dummy job to run and load the trees.
}

sub load_tree {
  my $logic_name  = "LoadTree";
  my $module      = "Bio::Greg::Slrsim::LoadTree";
  my $params      = {};
  $h->create_analysis( $logic_name, $module, $params, 50, 1 );
}

sub slrsim {
  my $logic_name = "Slrsim";
  my $module     = "Bio::Greg::Slrsim::Slrsim";
  my $params     = {};

  if (Bio::Greg::EslrUtils->is_ebi) {  
    $h->create_analysis( $logic_name, $module, $params, 250, 1 );
  } else {
    $h->create_analysis( $logic_name, $module, $params, 800, 1 );
  }
}

sub calculate_results {
  my $logic_name = "CalculateResults";
  my $module     = "Bio::Greg::Slrsim::CalculateResults";
  my $params     = {};
  $h->create_analysis( $logic_name, $module, $params, 30, 1 );

  $params = {
  };
}

sub plots {
  my $logic_name = "Plots";
  my $module     = "Bio::Greg::Slrsim::Plots";
  my $params     = {};
  $h->create_analysis( $logic_name, $module, $params, 1, 1 );

  $params = {
    experiment_name => $experiment_name
  };
  $h->add_job_to_analysis( "Plots", $params );    # Add a dummy job to plot at the end.
}

########*********########
#-------~~~~~~~~~-------#
########*********########
