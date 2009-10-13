#!/usr/bin/env perl

use warnings;
use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::Greg::EslrUtils;
use Bio::Greg::EslrPlots;
use Bio::Greg::ComparaLite::HiveUtils;
use Bio::EnsEMBL::Hive::DBSQL::DBAdaptor;
use File::Path;
use File::Basename;
use Cwd;

my ($url,$input,$output,$search) = undef;
GetOptions('url=s' => \$url,
           'input=s' => \$input,
           'output=s' => \$output,
           'search=s' => \$search
	   );

$url = 'mysql://greg:TMOqp3now@mysql-greg.ebi.ac.uk:4134/gj1_slrsim_1' if (!$url);
my ($mysql_args,$database,$project_base) = undef;
$mysql_args = Bio::Greg::ComparaLite::HiveUtils->hive_url_to_mysql_args($url);
$database = Bio::Greg::ComparaLite::HiveUtils->hive_url_to_hashref($url)->{'database'};
$project_base = getcwd();

# Load the adaptors.
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-url => $url);
my $dbc = $dba->dbc;
my $nsa = $dba->get_NestedSetAdaptor();
my $mba = $dba->get_MemberAdaptor();
my $pta = $dba->get_ProteinTreeAdaptor();

my $action = shift @ARGV;
exit(0) if (!$action);
plot() if ($action eq 'plot');
prune() if ($action eq 'prune');
plot_sweep('filter') if ($action eq 'plot_filter');
plot_sweep('params') if ($action eq 'plot_params');

sub plot {
  my @searches = ();
  $searches[0] = $search if ($search);
  while (my $arg = shift @ARGV) {
    push @searches, $arg;
  }

  my @outputs;
  foreach my $search (@searches) {
    my $cmd = qq^SELECT node_id,value FROM protein_tree_tag WHERE value LIKE "%$search%" LIMIT 1;^;
    my $sth = $dbc->prepare($cmd);
    $sth->execute();
    my @arr = $sth->fetchrow_array();
    my $node_id = $arr[0];
    my $value = $arr[1];
    print "$value ($node_id)\n";

    my $cur_out = $output."_${search}.pdf";
    my $params = {
      node_id => $node_id,
      max_tree_length => 4,
      dba => $dba,
      input_table => 'aln_mcoffee',
      alignment_score_filtering => 1,
      alignment_score_table => 'aln_mcoffee_prank',
      alignment_score_threshold => 5,
      quiet => 0,
      cleanup_temp => 0
    };
    Bio::Greg::EslrPlots->plotTree($cur_out,$params);
    push @outputs, $cur_out;
  }

  my $imgs = join(" ",@outputs);
  my $n = scalar(@outputs);
  my $tile_s = "-tile 1x${n}";
  `montage $imgs $tile_s -geometry 4800x400 $output`;

  map {unlink($_)} @outputs;
}

sub plot_sweep {
  my $what_to_sweep = shift;
  my $search = shift @ARGV;
  my @outputs;
  
  my $cmd = qq^SELECT node_id,value FROM protein_tree_tag WHERE value LIKE "%$search%" LIMIT 1;^;
  print $cmd."\n";
  my $sth = $dbc->prepare($cmd);
  $sth->execute();
  my @arr = $sth->fetchrow_array();
  my $node_id = $arr[0];
  my $value = $arr[1];
  print "$value ($node_id)\n";

  if ($what_to_sweep eq 'filter') {
    foreach my $threshold (1,3,5,7,9) {
      my $cur_out = $output."_${search}_${threshold}.pdf";
      my $params = {
        node_id => $node_id,
        max_tree_length => 4,
        dba => $dba,
        input_table => 'aln_mcoffee',
        alignment_score_filtering => 1,
        alignment_score_table => 'aln_mcoffee_prank',
        alignment_score_threshold => $threshold,
        quiet => 1
      };
      Bio::Greg::EslrPlots->plotTree($cur_out,$params);
      push @outputs, $cur_out;
    }
  } elsif ($what_to_sweep eq 'params') {
    foreach my $parameter_set_id (1..6) {
      my $cur_out = $output."_${search}_${parameter_set_id}.pdf";
      my $params = {
        parameter_set_id => $parameter_set_id,
        node_id => $node_id,
        max_tree_length => 4,
        dba => $dba,
        input_table => 'aln_mcoffee',
        quiet => 1
      };
      Bio::Greg::EslrPlots->plotTree($cur_out,$params);
      push @outputs, $cur_out;
    }    
  }

  my $imgs = join(" ",@outputs);
  my $n = scalar(@outputs);
  my $tile_s = "-tile 1x${n}";
  `montage $imgs $tile_s -geometry 4800x400 $output`;
  map {unlink($_)} @outputs;
}

sub prune {
  # Prunes everything outside of the species listed.
  $input = shift @ARGV if (!$input);
  my $tree = Bio::EnsEMBL::Compara::TreeUtils->from_file($input);
  my @keepers = @ARGV;
  my $keep_hash;
  map {$keep_hash->{$_} = 1} @keepers;
  foreach my $leaf ($tree->leaves) {
    if (!exists $keep_hash->{$leaf->name}) {
      Bio::EnsEMBL::Compara::TreeUtils->delete_lineage($tree,$leaf);
      $tree->minimize_tree;
    }
  }
  print Bio::EnsEMBL::Compara::TreeUtils->to_newick($tree)."\n";
}
