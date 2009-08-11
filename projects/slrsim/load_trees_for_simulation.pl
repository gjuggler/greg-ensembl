#!/usr/bin/env perl

use warnings;
use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::Greg::ComparaLite::HiveUtils;
use Bio::EnsEMBL::Hive::DBSQL::DBAdaptor;
use File::Path;
use File::Basename;
use Cwd;

my ($url,$folder,$config,$clean) = undef;
GetOptions('url=s' => \$url,
	   'folder=s' => \$folder,
	   'config=s' => \$config,
	   'clean' => \$clean
	   );
$url = 'mysql://ensadmin:ensembl@ensdb-2-12:5106/gj1_slrsim' if (!$url);
$folder = 'trees' if (!$folder);

my ($mysql_args,$database,$project_base) = undef;
$mysql_args = Bio::Greg::ComparaLite::HiveUtils->hive_url_to_mysql_args($url);
$database = Bio::Greg::ComparaLite::HiveUtils->hive_url_to_hashref($url)->{'database'};
$project_base = getcwd();


#if ($config) {
#  
#}

# Load the adaptors.
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-url => $url);
my $nsa = $dba->get_NestedSetAdaptor();
my $mba = $dba->get_MemberAdaptor();
my $pta = $dba->get_ProteinTreeAdaptor();

my $tree_table = "protein_tree";

if ($clean) {
  # Delete all trees and reset the counter.
  $dba->dbc->do("truncate table  ${tree_table}_member;");
  $dba->dbc->do("truncate table ${tree_table}_node;");
  $dba->dbc->do("truncate table  ${tree_table}_tag;");
  $dba->dbc->do("truncate table member;");
}

# The list of simulation replicates to use. These will be stored as tags in the XYZ_tag table,
# and accessible using the $node->get_tagvalue("tag_key") method.
my @simulation_sets = grep {$_ =~ /simset_/} keys %{$params};

# The number of replicates to simulate.
my $replicates = $params->{'simulation_replicates'};

my $tree_dir = $params->{'tree_dir'};
my @files = <$tree_dir/*.{nh,nhx,txt}>;
foreach my $file (sort @files) {    
  open(IN,"$file");
  my $newick_str = join("",<IN>);
  close(IN);
  
  foreach my $sim_set (@simulation_sets) {
    my $sim_param_str = $params->{$sim_set};
    foreach my $sim_rep (1 .. $replicates) {
      print "$file $sim_set $sim_rep\n";
      #my $node = Bio::EnsEMBL::Compara::TreeUtils->tree_from_newick($newick_str);

      # Go through each leaf and store the member objects.
      #foreach my $leaf ($node->leaves) {
	#$leaf->stable_id($leaf->name);
	#$leaf->source_name("ENSEMBLGENE");
	#$mba->store($leaf);
      #}

      # Store the tree in the XYZ_node table.
      #$pta->store($node);

      # Store the tags.
      #$node->store_tag("input_file",$file);
      #$node->store_tag("sim_set",$sim_set);
      #$node->store_tag("sim_rep",$sim_rep);
      #$node->store_tag("sim_params",$sim_param_str);
    }
  }
}

foreach my $node (@{$pta->fetch_all_roots}) {
  print $node->get_tagvalue("input_file")."\t".$node->node_id."\t" . scalar(@{$node->get_all_leaves}) . "\n";
}
