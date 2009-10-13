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
	   'config=s' => \$config,
	   'clean' => \$clean
	   );

#$url = 'mysql://ensadmin:ensembl@ensdb-2-12:5106/gj1_slrsim' if (!$url);
$url = 'mysql://greg:TMOqp3now@mysql-greg.ebi.ac.uk:4134/gj1_slrsim_1';

my ($mysql_args,$database,$project_base) = undef;
$mysql_args = Bio::Greg::ComparaLite::HiveUtils->hive_url_to_mysql_args($url);
$database = Bio::Greg::ComparaLite::HiveUtils->hive_url_to_hashref($url)->{'database'};
$project_base = getcwd();

# Load the adaptors.
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-url => $url);
my $nsa = $dba->get_NestedSetAdaptor();
my $mba = $dba->get_MemberAdaptor();
my $pta = $dba->get_ProteinTreeAdaptor();

my $params = {
  'tree_dir' => 'trees',
};

sub replace {
  my $p_one = shift;
  my $p_two = shift;
  return Bio::EnsEMBL::Compara::ComparaUtils->replace_params($p_one,$p_two);
}

sub create_simsets {
  my $base_p = {
    simulation_program => 'indelible',
    simulation_replicates => 20,
    tree_file => 'artificial.nh',
    tree_length => 1,
    seq_length => 400,
#    meanlog => -4.079,
#    sdlog => 1.23,
    ins_rate => 0.05,
    del_rate => 0.05,
    omega_distribution => 'M3',
    p0 => 0.386,
    p1 => 0.535,
    p2 => 0.079,
    w0 => 0.018,
    w1 => 0.304,
    w2 => 1.691,
  };
  
  my $i=1;
  my $base_length = 1;
  my $s_params;

  foreach my $j (0.1,1,5,10,20) {
#  foreach my $j (10) {
    $s_params = {
      w2 => 1.691,
      tree_file => 'artificial.nh',
      tree_length => $base_length * $j
    };
    $params->{sprintf 'simset_artificial_%02d',$i++} = replace($base_p,$s_params);
  }

  $i=1;
  foreach my $j (0.1,1,5,10,20) {
#  foreach my $j (10) {
    $s_params = {
      w2 => 4.739,
      tree_file => 'artificial.nh',
      tree_length => $base_length * $j
    };
    $params->{sprintf 'simset_artificial_strong_%02d',$i++} = replace($base_p,$s_params);
  }

  $i=1;
  $base_length = 1;
  foreach my $j (1,5,10,20,50,100,200) {
#  foreach my $j (10) {
    $s_params = {
      w2 => 1.691,
      tree_file => '2xmammals.nh',
      tree_length => $base_length * $j
    };
#    $params->{sprintf 'simset_2xmammal_%02d',$i++} = replace($base_p,$s_params);
  }

  $i=1;
  $base_length = 1;
  foreach my $j (1,5,10,20,50,100,200) {
    $s_params = {
      w2 => 4.739,
      tree_file => '2xmammals.nh',
      tree_length => $base_length * $j
    };
#    $params->{sprintf 'simset_2xmammal_strong_%02d',$i++} = replace($base_p,$s_params);
  }

}

create_simsets();
  
my $tree_table = "protein_tree";
if ($clean) {
  # Delete all trees and reset the counter.
  $dba->dbc->do("truncate table  ${tree_table}_member;");
  $dba->dbc->do("truncate table ${tree_table}_node;");
  $dba->dbc->do("truncate table  ${tree_table}_tag;");
  $dba->dbc->do("truncate table member;");
  $dba->dbc->do("truncate table sequence;");
  $dba->dbc->do("truncate table sitewise_aln;");
  $dba->dbc->do("truncate table omega_tr;");
  $dba->dbc->do("truncate table omega_mc;");
}

# The list of simulation replicates to use. These will be stored as tags in the XYZ_tag table,
# and accessible using the $node->get_tagvalue("tag_key") method.
my @simulation_sets = sort grep {$_ =~ /simset_/} keys %{$params};

my $tree_dir = $params->{'tree_dir'};
foreach my $sim_set (@simulation_sets) {
  my $p = $params->{$sim_set};
  my $replicates = $p->{'simulation_replicates'};
  my $tree_length = $p->{'tree_length'};
  my $file = $tree_dir.'/'.$p->{'tree_file'};
  
  open(IN,"$file");
  my $newick_str = join("",<IN>);
  close(IN);
  
  my $sim_param_str = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($p);
  foreach my $sim_rep (1 .. $replicates) {
    print "$file $sim_set $sim_rep\n";
    my $node = Bio::EnsEMBL::Compara::TreeUtils->from_newick($newick_str);
    my $bl = Bio::EnsEMBL::Compara::TreeUtils->total_distance($node);
    my $n = scalar $node->leaves;
    print "$bl $n \n";
    if ($tree_length) {
      $node = Bio::EnsEMBL::Compara::TreeUtils->scale_to($node,$tree_length);
      print "Tree length: $tree_length\n";
    }
    #print "Rescaled: ". Bio::EnsEMBL::Compara::TreeUtils->to_newick($node)."\n";
    
    # Go through each leaf and store the member objects.
    foreach my $leaf ($node->leaves) {
      $leaf->stable_id($leaf->name);
      #$leaf->source_name("ENSEMBLGENE");
      $leaf->source_name($sim_set.'_'.$sim_rep);
      $mba->store($leaf);
    }
    
    # Store the tree in the XYZ_node table.
    $pta->store($node);
    
    # Store the tags.
    $node->store_tag("input_file",$file);
    $node->store_tag("sim_set",$sim_set);
    $node->store_tag("sim_rep",$sim_rep);
    $node->store_tag("sim_params",$sim_param_str);
  }
}

foreach my $node (@{$pta->fetch_all_roots}) {
  print $node->get_tagvalue("input_file")."\t".$node->node_id."\t" . scalar(@{$node->get_all_leaves}) . "\n";
}
