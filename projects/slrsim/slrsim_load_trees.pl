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
$url = 'mysql://greg:TMOqp3now@mysql-greg.ebi.ac.uk:4134/gj1_slrsim_1' if (!$url);

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
my @simulation_sets = ();

sub replace {
  my $result = {};
  foreach my $p (@_) {
    $result = Bio::EnsEMBL::Compara::ComparaUtils->replace_params($result,$p);
  }
  return $result;
}

create_simsets();

sub create_simsets {
  my $i=1;
  my $base_length = 1;
  my $s_params;

  my $base_p = {
    simulation_program => 'indelible',
    simulation_replicates => 100,
    sim_file => 'artificial.nh',
    sim_length => 1,
    ins_rate => 0,
    del_rate => 0
  };

  my $anisimova_02_M3 = {
    seq_length => 500,
    omega_distribution => 'M3',
    p0 => 0.386,
    p1 => 0.535,
    p2 => 0.079,
    w0 => 0.018,
    w1 => 0.304,
    w2 => 1.691
  };

  my $neutral = {
    seq_length => 500,
    omega_distribution => 'M3',
    p0 => 0.5,
    p1 => 0.5,
    w0 => 1,
    w1 => 1
  };

  my $sim_counter=0;

  foreach my $j (0.11, 1.1, 5.5, 11) {
    foreach my $indel (0,0.01,0.05,0.1,0.2) {
      $s_params = {
        sim_name => 'art',
        sim_file => 'anisimova_01_artificial.nh',
        sim_ref => '1',
        sim_length => $j,
        ins_rate => $indel,
        del_rate => $indel
      };
      push @simulation_sets, replace($base_p,$anisimova_02_M3,$s_params);
    }
  }

  foreach my $j (0.11, 1.1, 5.5, 11) {
    $s_params = {
      sim_name => 'art_neutral',
      sim_file => 'anisimova_01_artificial.nh',
      sim_ref => '1',
      sim_length => $j
    };
    push @simulation_sets, replace($base_p,$neutral,$s_params); # Note: we're using the neutral evolution params here.
  }

#  foreach my $j (0.38,2.11,16.88,33.76) {
#    foreach my $indel (0, 0.05, 0.1) {
#      $s_params = {
#        sim_name => 'bglobin_human',
#        sim_file => 'anisimova_01_bglobin.nh',
#        sim_ref => 'human',
#        sim_length => $j,
#        ins_rate => $indel,
#        del-rate => $indel
#      };
#      push @simulation_sets, replace($base_p,$anisimova_02_M3,$s_params);
#    }
#  }

#  foreach my $j (0.38,2.11,16.88,33.76) {
#    foreach my $indel (0, 0.05, 0.1) {
#      $s_params = {
#        sim_name => 'bglobin_xen',
#        sim_file => 'anisimova_01_bglobin.nh',
#        sim_ref => 'xenlaev',
#        sim_length => $j,
#        ins_rate => $indel,
#        del_rate => $indel
#      };
#      push @simulation_sets, replace($base_p,$anisimova_02_M3,$s_params);
#    }
#  }

}

  
my $tree_table = "protein_tree";
if ($clean) {
  # Delete all trees and reset the counter.
  $dba->dbc->do("truncate table protein_tree_member;");
  $dba->dbc->do("truncate table protein_tree_node;");
  $dba->dbc->do("truncate table protein_tree_tag;");
  $dba->dbc->do("truncate table member;");
  $dba->dbc->do("truncate table sequence;");
}

# The list of simulation replicates to use. These will be stored as tags in the XYZ_tag table,
# and accessible using the $node->get_tagvalue("tag_key") method.
# my @simulation_sets = sort grep {$_ =~ /simset_/} keys %{$params};

my $tree_dir = $params->{'tree_dir'};
foreach my $params (@simulation_sets) {
  my $sim_set = $params->{'sim_name'};
  my $replicates = $params->{'simulation_replicates'};
  my $tree_length = $params->{'sim_length'};
  my $file = $tree_dir.'/'.$params->{'sim_file'};
  
  open(IN,"$file");
  my $newick_str = join("",<IN>);
  close(IN);
  
  my $sim_param_str = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
  foreach my $sim_rep (1 .. $replicates) {
    print "$file $sim_set $sim_rep\n";
    my $node = Bio::EnsEMBL::Compara::TreeUtils->from_newick($newick_str);
    my $bl = Bio::EnsEMBL::Compara::TreeUtils->total_distance($node);
    my $n = scalar $node->leaves;
    if ($tree_length) {
      $node = Bio::EnsEMBL::Compara::TreeUtils->scale_to($node,$tree_length);
    }
    my $length = Bio::EnsEMBL::Compara::TreeUtils->total_distance($node);

    #print "Rescaled: ". Bio::EnsEMBL::Compara::TreeUtils->to_newick($node)."\n";
    
    # Go through each leaf and store the member objects.
    foreach my $leaf ($node->leaves) {
      $leaf->stable_id($leaf->name);
      $leaf->source_name($sim_set.'_'.$sim_rep);
      $mba->store($leaf);
    }
    
    # Store the tree in the XYZ_node table.
    $pta->store($node);
    
    # Store the tags.
    foreach my $tag (grep {$_ =~ m/sim_/g} keys %{$params}) {
      $node->store_tag($tag,$params->{$tag});
    }
    $node->store_tag("sim_rep",$sim_rep);
    $node->store_tag("params_slrsim",$sim_param_str);
  }
}

foreach my $node (@{$pta->fetch_all_roots}) {
  print $node->get_tagvalue("input_file")."\t".$node->node_id."\t" . scalar(@{$node->get_all_leaves}) . "\n";
}

