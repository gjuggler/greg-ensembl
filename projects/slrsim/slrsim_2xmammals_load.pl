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
$url = 'mysql://greg:TMOqp3now@mysql-greg.ebi.ac.uk:4134/gj1_slrsim_2xmammals';

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


create_simsets();

sub create_simsets {
  my $base_p = {
    simulation_program => 'indelible',
    simulation_replicates => 100,
    tree_length => 1,
    seq_length => 400,
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

  $s_params = {
    w2 => 3,
    tree_file => '2x_nox.nh',
    tree_length => 0.787
  };
  $params->{sprintf 'simset_nox'} = replace($base_p,$s_params);

  $s_params = {
    w2 => 3,
    tree_file => '2x_v.nh',
    tree_length => 2.15
  };
  $params->{sprintf 'simset_v'} = replace($base_p,$s_params);

  $s_params = {
    w2 => 3,
    tree_file => '2x_p.nh',
    tree_length => 0.351
  };
  $params->{sprintf 'simset_p'} = replace($base_p,$s_params);

  $s_params = {
    w2 => 3,
    tree_file => '2x_g.nh',
    tree_length => 0.735
  };
  $params->{sprintf 'simset_g'} = replace($base_p,$s_params);

  $s_params = {
    w2 => 3,
    tree_file => '2x_l.nh',
    tree_length => 0.776
  };
  $params->{sprintf 'simset_l'} = replace($base_p,$s_params);

  
}
  
my $tree_table = "protein_tree";
if ($clean) {
  # Delete all trees and reset the counter.
  $dba->dbc->do("truncate table ${tree_table}_member;");
  $dba->dbc->do("truncate table ${tree_table}_node;");
  $dba->dbc->do("truncate table  ${tree_table}_tag;");
  $dba->dbc->do("truncate table member;");
  $dba->dbc->do("truncate table sequence;");
  $dba->dbc->do("truncate table sitewise_aln;");
  $dba->dbc->do("create table if not exists omega_tr like sitewise_aln;");
  $dba->dbc->do("truncate table omega_tr;");
  $dba->dbc->do("create table if not exists omega_mc like sitewise_aln;");
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
    print $node->newick_format()."\n";
    my $bl = Bio::EnsEMBL::Compara::TreeUtils->total_distance($node);
    my $n = scalar $node->leaves;
    if ($tree_length) {
      $node = Bio::EnsEMBL::Compara::TreeUtils->scale_to($node,$tree_length);
      print "  $bl -> $tree_length\n";
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

sub replace {
  my $p_one = shift;
  my $p_two = shift;
  return Bio::EnsEMBL::Compara::ComparaUtils->replace_params($p_one,$p_two);
}
