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
use Digest::MD5 qw(md5_hex);
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

#create_simsets();
x_simsets();

sub create_simsets {
  my $i=1;
  my $base_length = 1;
  my $s_params;

  my $base_p = {
    simulation_program => 'indelible',
    simulation_replicates => 100,
    sim_file => 'artificial.nh',
    sim_length => 1,
    seq_length => 500,
    ins_rate => 0,
    del_rate => 0
  };

  my $tree_anisimova_bglobin = 'anisimova_01_bglobin.nh';
  my $tree_anisimova_artificial = 'anisimova_01_artificial.nh';
  my $tree_primates = '2x_p.nh';
  my $tree_vertebrates = '2x_v.nh';

  my $anisimova_02_M3 = {
    omega_distribution => 'M3',
    p0 => 0.386,
    p1 => 0.535,
    p2 => 0.079,
    w0 => 0.018,
    w1 => 0.304,
    w2 => 1.691
  };
  my $anisimova_02_M3_hi = replace($anisimova_02_M3,{w2 => 4.739});

  my $massingham_05_A = {
    omega_distribution => 'M8',
    p0 => 0.9432,
    p => 0.572,
    q => 2.172,
    w => 2.081
  };
  # the same set of params was used by anisimova et al:
  my $anisimova_02_M8 = $massingham_05_A;

  my $massingham_05_B = {
    omega_distribution => 'M3',
    p0 => 0.75,
    p1 => 0.25,
    w0 => 0.5,
    w1 => 1.5
  };

  my $neutral = {
    omega_distribution => 'constant',
    w => 1
  };

  foreach my $j (0.11, 1.1, 5.5, 11) {
    foreach my $indel (0,0.01,0.05,0.1,0.2) {
      $s_params = {
        sim_name => 'art_anisimova',
        sim_file => $tree_anisimova_artificial,
        sim_ref => '1',
        sim_length => $j,
        ins_rate => $indel,
        del_rate => $indel
      };
      push @simulation_sets, replace($base_p,$anisimova_02_M3,$s_params);

      $s_params = replace($s_params, {sim_name => 'art_anisimova_hi'});
      push @simulation_sets, replace($base_p,$anisimova_02_M3_hi,$s_params);
    }
  }

  foreach my $j (0.11, 1.1, 5.5, 11) {
    $s_params = {
      sim_name => 'art_neutral',
      sim_file => $tree_anisimova_artificial,
      sim_ref => '1',
      sim_length => $j
    };
    push @simulation_sets, replace($base_p,$neutral,$s_params);
  }

  foreach my $j (0.38,2.11,16.88,33.76) {
    $s_params = {
      sim_name => 'bglobin_neutral',
      sim_file => $tree_anisimova_bglobin,
      sim_ref => 'human',
      sim_length => $j
    };
    push @simulation_sets, replace($base_p,$neutral,$s_params);
  }

  $s_params = {
    sim_name => 'mass_05_A',
    sim_file => $tree_anisimova_artificial,
    sim_ref => '1',
    sim_length => 1.1
  };
  push @simulation_sets, replace($base_p,$massingham_05_A,$s_params);

  $s_params = {
    sim_name => 'mass_05_B',
    sim_file => $tree_anisimova_artificial,
    sim_ref => '1',
    sim_length => 1.1
  };
  push @simulation_sets, replace($base_p,$massingham_05_B,$s_params);

  foreach my $j (0.38, 2.11, 16.88) {
    foreach my $indel (0,0.01,0.05,0.1,0.2) {
      $s_params = {
      sim_name => 'bglobin_anisimova',
      sim_file => $tree_anisimova_bglobin,
      sim_ref  => 'human',
      sim_length => $j,
      ins_rate => $indel,
      del_rate => $indel
      };

      push @simulation_sets, replace($base_p,$anisimova_02_M3,$s_params);
    }
  }

  foreach my $indel (0, 0.01, 0.02, 0.05) {
    $s_params = {
      sim_name => 'bglobin_ref_hum',
      sim_file => $tree_anisimova_bglobin,
      sim_ref => 'human',
      sim_length => 16.88,
      ins_rate => $indel,
      del_rate => $indel
    };
    push @simulation_sets, replace($base_p,$anisimova_02_M3_hi,$s_params);

    $s_params = replace($s_params, {sim_name => 'bglobin_ref_xen',
                                    sim_ref => 'xenlaev'});
    push @simulation_sets, replace($base_p,$anisimova_02_M3_hi,$s_params);
  }

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

sub x_simsets {
  my $base_p = {
    simulation_program => 'indelible',
    simulation_replicates => 100,
    sim_file => 'artificial.nh',
    seq_length => 500,
    tree_mult => 1,
    ins_rate => 0.05,
    del_rate => 0.05
  };

  my $full = '2x_v.nh';
  my $nox = '2x_nox.nh';
  my $primates = '2x_p.nh';
  my $glires = '2x_g.nh';
  my $fortyfourmammals = '44mammals.nh';
  my $ensembl = 'ensembl.nh';

  my $distr_lognormal = {
    omega_distribution => 'lognormal',
    meanlog => -4.079,
    sdlog => 1.23
  };

  my $s_params;

  foreach my $mult (0.1, 0.5, 1, 2, 5) {
    $base_p = replace($base_p,{tree_mult => $mult});

    $s_params = {
      sim_name => '2x full',
      sim_file => $full,
      sim_ref => 'human'
    };
    push @simulation_sets, replace($base_p,$distr_lognormal,$s_params);
    
    $s_params = {
      sim_name => '2x no2x',
      sim_file => $nox,
      sim_ref => 'human'
    };
    push @simulation_sets, replace($base_p,$distr_lognormal,$s_params);
    
    $s_params = {
      sim_name => '2x primates',
      sim_file => $primates,
      sim_ref => 'human'
    };
    push @simulation_sets, replace($base_p,$distr_lognormal,$s_params);
    
    $s_params = {
      sim_name => '2x glires',
      sim_file => $glires,
      sim_ref => 'mouse'
    };
    push @simulation_sets, replace($base_p,$distr_lognormal,$s_params);
  }

}

# The list of simulation replicates to use. These will be stored as tags in the XYZ_tag table,
# and accessible using the $node->get_tagvalue("tag_key") method.
# my @simulation_sets = sort grep {$_ =~ /simset_/} keys %{$params};

my $tree_dir = $params->{'tree_dir'};
foreach my $params (@simulation_sets) {
  my $sim_set = $params->{'sim_name'};
  my $replicates = $params->{'simulation_replicates'};
  my $tree_length = $params->{'sim_length'};
  my $tree_mult = $params->{'tree_mult'};
  my $file = $tree_dir.'/'.$params->{'sim_file'};  
  open(IN,"$file");
  my $newick_str = join("",<IN>);
  close(IN);
  
  my $sim_param_str = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
  foreach my $sim_rep (1 .. $replicates) {
    print "$file $sim_set $sim_rep\n";
    
    my $md5 = Digest::MD5->new;
    $md5->add($params);
    $md5->add($sim_rep);
    my $unique_string = substr($md5->hexdigest,20);

    my $node = Bio::EnsEMBL::Compara::TreeUtils->from_newick($newick_str);
    my $bl = Bio::EnsEMBL::Compara::TreeUtils->total_distance($node);
    my $n = scalar $node->leaves;
    if ($tree_length) {
      $node = Bio::EnsEMBL::Compara::TreeUtils->scale_to($node,$tree_length);
    }
    if ($tree_mult) {
      $node = Bio::EnsEMBL::Compara::TreeUtils->scale($node,$tree_mult);
    }
    my $length = Bio::EnsEMBL::Compara::TreeUtils->total_distance($node);
    print "  -> $length\n";

    # Go through each leaf and store the member objects.
    foreach my $leaf ($node->leaves) {
      $leaf->stable_id($leaf->name);
      $leaf->source_name($unique_string);
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

