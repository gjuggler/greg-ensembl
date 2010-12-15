#!/usr/bin/env perl

use warnings;
use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::Greg::EslrUtils;
use File::Path;
use File::Basename;

use Bio::Greg::Hive::HiveLoaderUtils;
use Bio::Greg::Hive::ComparaHiveLoaderUtils;

my ( $clean, $url, $experiment_name ) = undef;
GetOptions(
  'clean' => \$clean,
  'url=s' => \$url,
  'experiment=s' => \$experiment_name
);
Bio::EnsEMBL::Registry->no_version_check(1);

die("No experiment name given!") unless (defined $experiment_name);

$url = 'mysql://slrsim:slrsim@mysql-greg.ebi.ac.uk:4134/slrsim_anisimova' if ( !$url );
my $h = new Bio::Greg::Hive::ComparaHiveLoaderUtils;
$h->init($url);

if ($clean) {
  $h->clean_hive_tables;
  $h->clean_compara_analysis_tables;
  $h->clean_compara_tree_tables;
  my @truncate_tables = qw^
      aln aln_scores omega
      meta
  ^;
  $h->truncate_tables(\@truncate_tables);
}
#output_dir();

load_trees();
load_tree();
simulate();
align();
scores();
dump_alignments();
omegas();
collect_stats();
plots();

$h->connect_analysis( "LoadTrees",   "LoadTree" );
$h->connect_analysis( "LoadTree",   "PhyloSim" );
$h->connect_analysis( "PhyloSim",    "Align" );
$h->connect_analysis( "Align",       "AlignScores" );
$h->connect_analysis( "AlignScores", "DumpAlignments" );
$h->connect_analysis( "AlignScores", "Omegas" );
$h->connect_analysis( "Omegas",      "CollectStats" );

$h->wait_for("Plots",["Omegas","CollectStats"]);
$h->wait_for("PhyloSim",["LoadTrees","LoadTree"]);
#$h->wait_for("Align",["PhyloSim"]);
$h->wait_for("AlignScores",["PhyloSim"]);

sub output_dir {
  $h->dba->dbc->do("insert into meta(meta_key,meta_value) VALUES ('hive_output_dir','/homes/greg/hive_temp');");
}

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
  my $params      = { experiment_name => $experiment_name };
  $h->create_analysis( $logic_name, $module, $params, 100, 1 );
}

sub simulate {
  my $logic_name = "PhyloSim";
  my $module     = "Bio::Greg::Hive::PhyloSim";
  my $params     = {};
  $h->create_analysis( $logic_name, $module, $params, 200, 1 );
}

sub align {
  my $logic_name = "Align";
  my $module     = "Bio::Greg::Hive::Align";
  my $params     = {
    # These params will be filled in by the LoadTree simulation definitions.
  };
  $h->create_analysis( $logic_name, $module, $params, 300, 1 );
}

sub scores {
  my $logic_name = "AlignScores";
  my $module     = "Bio::Greg::Hive::AlignmentScores";
  my $params     = {

    # These params will be filled in by the LoadTree simulation definitions.
  };
  $h->create_analysis( $logic_name, $module, $params, 300, 1 );
}

sub dump_alignments {
  my $logic_name = "DumpAlignments";
  my $module     = "Bio::Greg::Slrsim::DumpSlrsimAlignment";
  my $params     = {};
  $h->create_analysis( $logic_name, $module, $params, 80, 1 );
}

sub omegas {
  my $logic_name = "Omegas";
  my $module     = "Bio::Greg::Hive::PhyloAnalysis";
  my $params     = {

    # These params will be filled in by the LoadTree simulation definitions.
  };
  $h->create_analysis( $logic_name, $module, $params, 500, 1 );
}

sub collect_stats {
  my $logic_name = "CollectStats";
  my $module     = "Bio::Greg::Slrsim::CollectSlrsimStats";
  my $params     = {};
  $h->create_analysis( $logic_name, $module, $params, 100, 1 );
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
