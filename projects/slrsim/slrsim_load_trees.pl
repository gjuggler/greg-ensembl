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

create_simsets();

sub create_simsets {
  my $i=1;
  my $base_length = 1;

  my $base_params = {
    slrsim_replicates => 50,
    slrsim_file => 'artificial.nh',
    slrsim_tree_mult => 1,

    phylosim_simulation_program => 'indelible',
    phylosim_seq_length => 400,
    phylosim_ins_rate => 0,
    phylosim_del_rate => 0
  };

  my $tree_anisimova_bglobin = 'anisimova_01_bglobin.nh';
  my $tree_anisimova_artificial = 'anisimova_01_artificial.nh';
  my $full = '2x_v.nh';
  my $nox = '2x_nox.nh';
  my $primates = '2x_p.nh';
  my $glires = '2x_g.nh';
  my $fortyfourmammals = '44mammals.nh';
  my $ensembl = 'ensembl.nh';

  my $anisimova_02_M3 = {
    slrsim_scheme_name => "Anisimova 2002 M3",
    phylosim_omega_distribution => 'M3',
    phylosim_p0 => 0.386,
    phylosim_p1 => 0.535,
    phylosim_p2 => 0.079,
    phylosim_w0 => 0.018,
    phylosim_w1 => 0.304,
    phylosim_w2 => 1.691
  };
  my $anisimova_02_M3_hi = replace($anisimova_02_M3,{w2 => 4.739});

  my $massingham_05_A = {
    slrsim_scheme_name => "Massingham 2005-A (M8)",    
    phylosim_omega_distribution => 'M8',
    phylosim_p0 => 0.9432,
    phylosim_p => 0.572,
    phylosim_q => 2.172,
    phylosim_w => 2.081
  };
  # the same set of params was used by anisimova et al:
  my $anisimova_02_M8 = $massingham_05_A;  

  my $massingham_05_B = {
    slrsim_scheme_name => "Massingham 2005-B (M3)",    
    phylosim_omega_distribution => 'M3',
    phylosim_p0 => 0.75,
    phylosim_p1 => 0.25,
    phylosim_w0 => 0.5,
    phylosim_w1 => 1.5
  };

  my $neutral = {
    slrsim_scheme_name => "Neutral evolution",
    phylosim_omega_distribution => 'constant',
    phylosim_w => 1
  };

  my $lognormal = {
    # These lognormal parameters should result in p(w>1) = 0.0835
    slrsim_scheme_name => "2xmammals Lognormal",
    phylosim_omega_distribution => 'lognormal',
    phylosim_meanlog => -1.7,
    phylosim_sdlog => 1.23
  };

  my $slr = {
    sitewise_name => 'SLR',
    sitewise_action => 'slr',
  };
  my $paml = {
    sitewise_name => 'PAML M8A/M8B',
    sitewise_action => 'paml_sitewise paml_lrt',
    paml_model_a => 'M8a',
    paml_model_b => 'M8'
  };
  my $paml_alt = {
    sitewise_name => 'PAML M2/M3',
    sitewise_action => 'paml_sitewise paml_lrt',
    paml_model_a => 'M2',
    paml_model_b => 'M3'
  };
  

  # Looking at bglobin and difference between reference sequences
  my $sim_p = replace($massingham_05_A,{
    slrsim_tree_mult => 1,
  });
  my @sim_params;
  foreach my $indel (0.01,0.02,0.05) {
    push @sim_params, replace($sim_p,{phylosim_ins_rate => $indel,phylosim_del_rate => $indel});
  }
  my @tree_params = map {tree_param($_)} ($tree_anisimova_bglobin);
  my @aln_params = map {aln_param($_)} ('mcoffee');
  my @filter_params = map {filter_param($_)} ('none');
  my @species_params = map {species_param($_)} ('human','mouse','xenlaev');
  my @sitewise_params = ($slr,$paml,$paml_alt);
  foreach my $tr (@tree_params) {
    foreach my $sim (@sim_params) {
      foreach my $aln (@aln_params) {
        foreach my $f (@filter_params) {
          foreach my $sp (@species_params) {
            foreach my $sw (@sitewise_params) {
              my $p = replace($base_params,$tr,$sim,$aln,$f,$sp,$sw);
              verify_params($p);
              push @simulation_sets,$p;
            }
          }
        }
      }
    }  
  }


  $sim_p = replace($lognormal,{
    slrsim_tree_mult => 2,
  });
  foreach my $indel (0.05) {
    push @sim_params, replace($sim_p,{phylosim_ins_rate => $indel,phylosim_del_rate => $indel});
  }
  @tree_params = map {tree_param($_)} ($full);
  @aln_params = map {aln_param($_)} ('mcoffee');
  @filter_params = map {filter_param($_)} ('none','prank','trimal','indelign','gblocks');
  @species_params = map {species_param($_)} ('human');
  @sitewise_params = ($slr,$paml_alt);
  foreach my $tr (@tree_params) {
    foreach my $sim (@sim_params) {
      foreach my $aln (@aln_params) {
        foreach my $f (@filter_params) {
          foreach my $sp (@species_params) {
            foreach my $sw (@sitewise_params) {
              my $p = replace($base_params,$tr,$sim,$aln,$f,$sp,$sw);
              verify_params($p);
              push @simulation_sets,$p;
            }
          }
        }
      }
    }  
  }

}  

sub verify_params {
  my $p = shift;

  # Check that we have proper filter and alignment tables.
  $p->{alignment_score_table} = $p->{alignment_table}.'_'.$p->{alignment_scores_action};
  create_filter_table($p->{alignment_score_table});

  create_aln_table($p->{alignment_table});
  create_aln_table($p->{alignment_table}.'_score');
  create_omega_table($p->{omega_table});
}

sub tree_param {
  my $tree = shift;
  return {slrsim_file => $tree};
}

sub aln_param {
  my $aln = shift;
  my $aln_params = {
    alignment_name => $aln,
    alignment_method => $aln,
    alignment_table => 'aln',
    omega_table => 'omega'
  };
  if ($aln eq 'none') {
    $aln_params = {
      alignment_name => "True Alignment",
      alignment_method => 'none',
      alignment_table => 'protein_tree_member',
      omega_table => 'omega_true'
    };
  }
  return $aln_params;
}

sub filter_param {
  my $filter = shift;
  
  my $f = {
    filtering_name => $filter,
    alignment_scores_action => $filter,
    alignment_score_filtering => 1,
    alignment_score_threshold => 5
  };
  $f->{alignment_score_threshold} = 7 if ($filter eq 'prank');
  if ($filter eq 'mcoffee') {
    $f->{alignment_scores_action} = 'score';
  }
  if ($filter eq 'none') {
    $f = {
      filtering_name => 'None',
      alignment_scores_action => 'score',
      alignment_score_filtering => 0
    };    
  }
  return $f;
}

sub species_param {
  my $species = shift;
  my $p = {
    species_name => $species,
    slrsim_ref => $species
  };
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
  my $sim_set = $params->{'slrsim_scheme_name'};
  my $replicates = $params->{'slrsim_replicates'};
  my $tree_length = $params->{'slrsim_tree_length'};
  my $tree_mult = $params->{'slrsim_tree_mult'};
  my $file = $tree_dir.'/'.$params->{'slrsim_file'};
  open(IN,"$file");
  my $newick_str = join("",<IN>);
  close(IN);
  
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
    my $final_length = Bio::EnsEMBL::Compara::TreeUtils->total_distance($node);
    print "  -> $final_length\n";

    # Go through each leaf and store the member objects.
    foreach my $leaf ($node->leaves) {
      $leaf->stable_id($leaf->name);
      $leaf->source_name($unique_string);
      $mba->store($leaf);
    }
    
    # Store the tree in the XYZ_node table.
    $pta->store($node);
    
    # Store all parameters as tags.
    $params->{'slrsim_rep'} = $sim_rep;
    $params->{'slrsim_tree_length'} = $final_length;
    my $sim_param_str = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
    foreach my $tag (keys %{$params}) {
      $node->store_tag($tag,$params->{$tag});
    }
    # Store the whole parameter set as a string.
    $node->store_tag("params_slrsim",$sim_param_str);
  }
}

foreach my $node (@{$pta->fetch_all_roots}) {
  print $node->get_tagvalue("input_file")."\t".$node->node_id."\t" . scalar(@{$node->get_all_leaves}) . "\n";
}


sub create_aln_table {
  my $table = shift;
  $dba->dbc->do("CREATE TABLE IF NOT EXISTS $table LIKE protein_tree_member");
  $dba->dbc->do("TRUNCATE TABLE $table");
}
sub create_omega_table {
  my $table = shift;
  $dba->dbc->do("CREATE TABLE IF NOT EXISTS $table LIKE sitewise_omega");
  $dba->dbc->do("TRUNCATE TABLE $table");
}
sub create_filter_table {
  my $table = shift;
  create_aln_table($table);
}
