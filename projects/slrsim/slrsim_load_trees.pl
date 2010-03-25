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

$url = 'mysql://slrsim:slrsim@mysql-greg.ebi.ac.uk:4134/slrsim_anisimova' if (!$url);

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
    slrsim_replicates => 1,
    slrsim_file => 'artificial.nh',
    slrsim_tree_mult => 1,
    
    phylosim_simulation_program => 'indelible',
    phylosim_seq_length => 400,
    phylosim_insertrate => 0,
    phylosim_deleterate => 0,
    
    sitewise_name => 'None',
    sitewise_action => ''
  };
  
  my $tree_anisimova_bglobin = 'anisimova_01_bglobin.nh';
  my $tree_anisimova_artificial = 'anisimova_01_artificial.nh';
  my $full = '2x_v.nh';
  my $nox = '2x_nox.nh';
  my $primates = '2x_p.nh';
  my $glires = '2x_g.nh';
  my $fortyfourmammals = '44mammals.nh';
  my $twoxmammals = '2xmammals.nh';
  my $ensembl = 'ensembl.nh';
  my $ensembl_two = 'ensembl_2.nh';

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

  my $massingham_05_A2 = {
    slrsim_scheme_name => "Massingham 2005-A2 (M8, 15% pos-sel)",
    phylosim_omega_distribution => 'M8',
    phylosim_p0 => 0.85,
    phylosim_p => 0.572,
    phylosim_q => 2.172,
    phylosim_w => 2.081
  };

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

  my $sitewise_none = {
    sitewise_name => 'None',
    sitewise_action => ''
  };
  my $slr = {
    sitewise_name => 'SLR',
    sitewise_action => 'slr'
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
  
  my $sim_p = replace($massingham_05_A, {
    slrsim_tree_mult => 1,
    phylosim_insertmodel => 'NB 0.18 2',
    phylosim_deletemodel => 'NB 0.18 2'
                   });

  my $sim_p_b = replace($massingham_05_B, {
    slrsim_tree_mult => 1,
    phylosim_insertmodel => 'POW 2.2 30',
    phylosim_deletemodel => 'POW 2.2 30'
                         });
  
  sub indel_models {
    my @sets = ();
    my @sim_params = ();
    my @tree_params = 
    
    my $cur_params = replace($base_params,$sitewise_none,aln_param('true'),{
      slrsim_ref => 'Human',
      omega_table => 'sitewise_omega'
                             }); 
    my $indel = 0.05;
    # Power-law sweep.
    my $max_length = 100;
    foreach my $a (1.6, 1.8, 2.0, 2.2) {
      my $p = replace($sim_p_b,{
        phylosim_insertrate => $indel,
        phylosim_deleterate => $indel,
        phylosim_insertmodel => "POW $a $max_length",
        phylosim_deletemodel => "POW $a $max_length"
                      });
      push @sim_params, $p;
    }    
    # Negative binomial sweep.
    foreach my $a (0.2, 0.3, 0.4) {
      my $p = replace($sim_p_b, {
        phylosim_insertrate => $indel,
        phylosim_deleterate => $indel,
        phylosim_insertmodel => "NB $a 2",
        phylosim_deletemodel => "NB $a 2",        
                      });
      push @sim_params, $p;
    }

    foreach my $tr (map {tree_param($_)} ($fortyfourmammals,$tree_anisimova_artificial)) {
      foreach my $sp (@sim_params) {
        my $p = replace($cur_params,$tr,$sp);
        Bio::EnsEMBL::Compara::ComparaUtils->hash_print($p);
        verify_params($p);
        push @sets,$p;
      }
    }  
    return @sets;
  }

  sub indel_rates {
    my @sets = ();
    my @sim_params = ();
    
    my $cur_params = replace($base_params,$sitewise_none,aln_param('true'),{
      slrsim_ref => 'Human',
      omega_table => 'sitewise_omega'
                             }); 
    foreach my $indel (0.01, 0.02, 0.05, 0.1, 0.2) {
      my $p = replace($sim_p_b,{
        phylosim_insertrate => $indel,
        phylosim_deleterate => $indel,
        phylosim_insertmodel => "POW 1.8 50",
        phylosim_deletemodel => "POW 1.8 50"
                      });
      push @sim_params, $p;
    }    

    foreach my $mult (0.5, 1, 2) {
      my $p = replace($sim_p_b,{
        phylosim_insertrate => 0.05,
        phylosim_deleterate => 0.05,
        phylosim_insertmodel => "POW 1.8 50",
        phylosim_deletemodel => "POW 1.8 50",
        slrsim_tree_mult => $mult
                      });
      push @sim_params, $p;      
    }
    
    foreach my $tr (map {tree_param($_)} ($fortyfourmammals,$tree_anisimova_bglobin,$tree_anisimova_artificial)) {
      foreach my $sp (@sim_params) {
        my $p = replace($cur_params,$tr,$sp);
        Bio::EnsEMBL::Compara::ComparaUtils->hash_print($p);
        verify_params($p);
        push @sets,$p;
      }
    }  
    return @sets;
  }
  
  sub all_trees {
    my @sets = ();
    my $cur_params = replace($base_params,$sim_p_b,$sitewise_none,aln_param('true'),{
      slrsim_ref => 'Human',
      omega_table => 'sitewise_omega',
      phylosim_insertmodel => 'POW 1.8 50',
      phylosim_deletemodel => 'POW 1.8 50',
      phylosim_insertrate => 0.05,
      phylosim_deleterate => 0.05,
      slrsim_tree_mult => 1
                             }); 
    
    my @tree_params = map {tree_param($_)} (
      $ensembl,
      $ensembl_two,
      $twoxmammals,
      $tree_anisimova_artificial,
      $tree_anisimova_bglobin,
      $full,
      $primates,
      $glires,
      $nox,
      $fortyfourmammals
    );
    foreach my $tr (@tree_params) {
      my $p = replace($cur_params,$tr);
      verify_params($p);
      push @sets, $p;
    }
    return @sets;
  }
  
  sub filtering_sim {
    my @sets = ();
    my $cur_params = replace($base_params,
                             $sim_p_b,
                             $sitewise_none,{
                             slrsim_ref => 'Human',
                             omega_table => 'sitewise_omega',
                             phylosim_insertmodel => 'POW 1.8 50',
                             phylosim_deletemodel => 'POW 1.8 50',
                             phylosim_insertrate => 0.05,
                             phylosim_deleterate => 0.05,
                             slrsim_tree_mult => 1
      });

    my @filter_params = ();
     foreach my $threshold (0,2,5,7) {
      my $p = replace(filter_param('prank_treewise'),{alignment_score_threshold => $threshold});
      push @filter_params,$p;
    }
 
    foreach my $threshold (0,2,5,7) {
      my $p = replace(filter_param('prank_mean'),{alignment_score_threshold => $threshold});
      push @filter_params,$p;
    }

    foreach my $threshold (0,2,5,7) {
      my $p = replace(filter_param('indelign'),{alignment_score_threshold => $threshold});
      push @filter_params,$p;
    }
    
    foreach my $threshold (0,5) { # TrimAl is on-or-off, no need to sweep here.
      my $p = replace(filter_param('trimal'),{alignment_score_threshold => $threshold});
      push @filter_params,$p;
    }

    foreach my $threshold (0,2,5,7) {
      my $p = replace(filter_param('tcoffee'),{alignment_score_threshold => $threshold});
      push @filter_params,$p;
    }

    my @tree_params = map {tree_param($_)} ($fortyfourmammals);
    my @aln_params = map {aln_param($_)} ('true','mcoffee','prank','prank_f','prank_codon_f');
    foreach my $tr (@tree_params) {
      foreach my $aln (@aln_params) {
        foreach my $fp (@filter_params) {
          my $p = replace($cur_params,$aln,$tr,$fp);
          verify_params($p);
          push @sets, $p;
        }
      }
    }
    return @sets;
  }  
    
  sub sitewise_sim {
    my @sets = ();
    my $cur_params = replace($base_params,
                             $massingham_05_A2, # Sitewise sim.
                             $slr, # Sitewise analysis.
                             {
                               slrsim_ref => 'Human',
                               phylosim_insertmodel => 'POW 1.8 50',
                               phylosim_deletemodel => 'POW 1.8 50',
                               phylosim_insertrate => 0.05,
                               phylosim_deleterate => 0.05,
                               alignment_score_threshold => 5,
                               slrsim_tree_mult => 1,
                               slrsim_replicates => 100
                             });

    my @tree_params = map {tree_param($_)} ($fortyfourmammals);
    my @filter_params = map {filter_param($_)} ('none','prank_treewise','prank_mean','indelign','trimal','tcoffee');
    my @aln_params = map {aln_param($_)} ('true','mcoffee','prank_f');#,'prank_f','prank_codon_f');
    foreach my $tr (@tree_params) {
      foreach my $aln (@aln_params) {
        foreach my $fp (@filter_params) {
          my $p = replace($cur_params,$aln,$tr,$fp);
          verify_params($p);
          push @sets, $p;
        }
      }
    }
    return @sets;
  }  

  sub ref_sim {
    my @sets = ();
    my $cur_params = replace($base_params,
                             $massingham_05_A2, # Sitewise sim.
                             $slr, # Sitewise analysis.
                             {
                               phylosim_insertmodel => 'POW 1.8 50',
                               phylosim_deletemodel => 'POW 1.8 50',
                               phylosim_insertrate => 0.05,
                               phylosim_deleterate => 0.05,
                               alignment_score_threshold => 5,
                               slrsim_tree_mult => 1,
                               slrsim_replicates => 1,
                               slrsim_file => $fortyfourmammals
                             });
    my @species_params = map {species_param($_)} ('Human','Mouse','Xtropicalis','Stickleback');
    my @aln_params = map {aln_param($_)} ('true','mcoffee','prank_f');
    my @filter_params = map {filter_param($_)} ('none','tcoffee');
    foreach my $sp (@species_params) {
      foreach my $aln (@aln_params) {
        foreach my $fp (@filter_params) {
          my $p = replace($cur_params,$sp,$aln,$fp);
          verify_params($p);
          push @sets, $p;
        }
      }
    }
    return @sets;
  }
  
  @simulation_sets = ref_sim();
  #@simulation_sets = sitewise_sim();
  #@simulation_sets = filtering_sim();
  #@simulation_sets = all_trees();
  #@simulation_sets = indel_models();
  #@simulation_sets = indel_rates();
}


sub verify_params {
  my $p = shift;

  $p->{alignment_table} = 'aln' unless defined ($p->{alignment_table});
  $p->{alignment_scores_action} = '' unless defined ($p->{alignment_scores_action});
  $p->{omega_table} = 'sitewise_omega' unless defined ($p->{omega_table});

  # Check that we have proper filter and alignment tables.
  #$p->{alignment_score_table} = $p->{alignment_table}.'_'.$p->{alignment_scores_action};
  $p->{alignment_score_table} = $p->{alignment_table}.'_scores';
  create_filter_table($p->{alignment_score_table});

  create_aln_table($p->{alignment_table});
  create_omega_table($p->{omega_table});
}

sub mult_param {
  my $mult = shift;
  return {slrsim_tree_mult => $mult};
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
  if ($aln eq 'none' || $aln eq 'true') {
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
  $newick_str =~ s/\n//g;
  close(IN);
  
  my $base_node = Bio::EnsEMBL::Compara::TreeUtils->from_newick($newick_str);
  foreach my $sim_rep (1 .. $replicates) {
    print "$file $sim_rep\n";
    my $node = $base_node->copy;

    my $md5 = Digest::MD5->new;
    $md5->add($params);
    $md5->add($sim_rep);
    my $unique_string = substr($md5->hexdigest,20);

    my $bl = Bio::EnsEMBL::Compara::TreeUtils->total_distance($node);
    my $n = scalar $node->leaves;
    if ($tree_length) {
      $node = Bio::EnsEMBL::Compara::TreeUtils->scale_to($node,$tree_length);
    }
    if ($tree_mult) {
      $node = Bio::EnsEMBL::Compara::TreeUtils->scale($node,$tree_mult);
    }
    my $final_length = Bio::EnsEMBL::Compara::TreeUtils->total_distance($node);
    my $mean_length = $final_length / scalar($node->nodes);
#    printf "  -> total: %.3f  mean: %.3f\n",$final_length,$mean_length;

    # Go through each leaf and store the member objects.
    foreach my $leaf ($node->leaves) {
      $leaf->stable_id($leaf->name);
      $leaf->source_name($unique_string);
      $mba->store($leaf);
      $leaf->adaptor(undef);
    }

    # Store the tree in the XYZ_node table.
    $pta->store($node);

    printf "  > %.50s\n",$node->newick_format;
    
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
  printf "  > %.50s\n",$node->newick_format;
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
