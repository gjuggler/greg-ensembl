package Bio::Greg::Slrsim::LoadTrees;

use strict;
use Cwd;
use Time::HiRes;
use POSIX qw(strftime mktime);

use Bio::EnsEMBL::Hive::Process;
use Bio::EnsEMBL::Utils::Exception;

use base ('Bio::Greg::Hive::Process');

sub fetch_input {
  my $self = shift;

  ### DEFAULT PARAMETERS ###
  my $defaults = {
    experiment_name => 'filter_sweeps',
    tree_root_dir   => $self->base . '/projects/slrsim/trees'
  };
  ##########################

  $self->load_all_params($defaults);

  $self->load_simulation_params();
  $self->load_simulation_sets();
}

sub run {
  my $self = shift;

  my @simsets;

  my $experiment_name = $self->param('experiment_name');
  if ( defined $experiment_name ) {

    # Tricky Perl: Call the method corresponding to the current experiment.
    @simsets = @{ $self->$experiment_name() };
  }

  printf "Creating %d simulation sets...\n", scalar(@simsets);

  my $creation_tags = {
    timestamp       => strftime( "%Y-%m-%d", localtime ),
    experiment_name => $experiment_name
  };
  $self->store_meta($creation_tags);

  throw("No simulation sets for experiment $experiment_name !") unless ( defined @simsets );

  my $tree_dir = $self->param('tree_root_dir') || '.';

  $self->disconnect_when_inactive(0);

  foreach my $params (@simsets) {
    $self->verify_params($params);

    print "Loading params:\n";
    $self->hash_print($params);

    my $replicates = $params->{slrsim_replicates};
    $replicates = 1 if ( $replicates < 1 );

    my $base_node;
    if ( defined $params->{slrsim_tree} ) {
      $base_node = $params->{slrsim_tree};
    } elsif ( defined $params->{slrsim_tree_file} ) {
      $base_node = $self->_get_tree_from_file( $params->{'slrsim_tree_file'} );
    } elsif ( defined $params->{slrsim_tree_newick} ) {
      print "Loading from newick...\n";
      $base_node = Bio::EnsEMBL::Compara::TreeUtils->from_newick( $params->{slrsim_tree_newick} );
    }

    my @tree_lengths;
    if ( defined $params->{slrsim_tree_lengths} ) {
      @tree_lengths = @{ $params->{slrsim_tree_lengths} };

      #print join(",",@tree_lengths)."\n";
    }

    foreach my $sim_rep ( 1 .. $replicates ) {

      # Take the next tree length if we're provided with a list of lengths.
      if ( defined @tree_lengths ) {
        my $length = $tree_lengths[ $sim_rep - 1 ];
        $params->{slrsim_tree_length} = $length;
      }
      $self->load_tree_into_database( $base_node, $sim_rep, $params );
      sleep(0.2);
    }
  }

  #  my $pta = $self->pta;
  #  foreach my $node ( @{ $pta->fetch_all_roots } ) {
  #    print $node->get_tagvalue("input_file") . "\t"
  #      . $node->node_id . "\t"
  #      . scalar( @{ $node->get_all_leaves } ) . "\n";
  #    printf "  > %.50s\n", $node->newick_format;
  #  }

}

sub _get_tree_from_file {
  my $self     = shift;
  my $filename = shift;

  my $tree_dir = $self->param('tree_root_dir') || '.';
  my $file = $tree_dir . '/' . $filename;
  print " -> Tree file: $file\n";
  open( IN, "$file" );
  my $newick_str = join( "", <IN> );
  $newick_str =~ s/\n//g;
  close(IN);
  print " -> Tree string: $newick_str\n";
  my $tree = Bio::EnsEMBL::Compara::TreeUtils->from_newick($newick_str);
  return $tree;
}

sub load_tree_into_database {
  my $self    = shift;
  my $tree    = shift;
  my $sim_rep = shift;
  my $params  = shift;

  my $mba = $self->compara_dba->get_MemberAdaptor;
  my $pta = $self->compara_dba->get_ProteinTreeAdaptor;

  my $node = $tree->copy;

  my $md5 = Digest::MD5->new;
  $md5->add($params);
  $md5->add($sim_rep);
  my $unique_string = substr( $md5->hexdigest, 20 );

  my $tree_mult = $params->{'slrsim_tree_mult'};
  if ($tree_mult) {
    Bio::EnsEMBL::Compara::TreeUtils->scale( $node, $tree_mult );
  }
  my $tree_length = $params->{'slrsim_tree_length'};
  if ($tree_length) {
    $node = Bio::EnsEMBL::Compara::TreeUtils->scale_to( $node, $tree_length );
  }
  my $tree_max_path = $params->{'slrsim_tree_max_path'};
  if ($tree_max_path) {
    print "Tree max path: $tree_max_path\n";
    $node = Bio::EnsEMBL::Compara::TreeUtils->scale_max_to( $node, $tree_max_path );
  }
  my $tree_mean_path = $params->{'slrsim_tree_mean_path'};
  if ($tree_mean_path) {
    print "Tree mean path: $tree_mean_path\n";
    $node = Bio::EnsEMBL::Compara::TreeUtils->scale_mean_to( $node, $tree_mean_path );
  }

  my $final_length = Bio::EnsEMBL::Compara::TreeUtils->total_distance($node);

  print "Final lenth: $final_length\n";

  # Go through each leaf and store the member objects.
  foreach my $leaf ( $node->leaves ) {

    #print "leaf: ".$leaf->name."\n";
    $leaf->stable_id( $leaf->name );
    $leaf->source_name($unique_string);
    $mba->store($leaf);
    $leaf->adaptor(undef);
  }

  # Store the tree in the XYZ_node table.
  $pta->store($node);

  my $node_id = $node->node_id;
  print " -> Node ID: $node_id\n";
  $self->param( 'node_id', $node_id );

  #  printf "  > %.50s\n", $node->newick_format;

  # Store all parameters as tags.
  $params->{'slrsim_rep'}         = $sim_rep;
  $params->{'slrsim_tree_length'} = $final_length;
  delete $params->{'slrsim_tree_lengths'};
  my $sim_param_str = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
  foreach my $tag ( sort keys %{$params} ) {
    $self->store_tag( $tag, $params->{$tag} );
    sleep(0.1);
  }

  # Store the whole parameter set as a string.
  $self->store_tag( "params_slrsim", $sim_param_str );

  my $output_params = {
    node_id          => $node_id,
    parameter_set_id => $self->parameter_set_id,
    experiment_name  => $params->{'experiment_name'}
  };
  my $data_id = $self->new_data_id($output_params);
  $output_params->{data_id} = $data_id;

  my ($job_id) = @{ $self->dataflow_output_id( $output_params, 1 ) };
  print "  -> Created job: $job_id \n";

  $node->release_tree;
}

sub load_simulation_params {
  my $self = shift;

  # Trees.
  my $trees = {
    'anisimova_bglobin'    => 'anisimova_01_bglobin.nh',
    'anisimova_artificial' => 'anisimova_01_artificial.nh',
    '2x_full'              => '2x_v.nh',
    '2x_nox'               => '2x_nox.nh',
    '2x_primates'          => '2x_p.nh',
    '2x_glires'            => '2x_g.nh',
    '44mammals'            => '44mammals.nh',
    '2xmammals'            => '2xmammals.nh',
    'ensembl_a'            => 'ensembl.nh',
    'ensembl_b'            => 'ensembl_2.nh',
    'mammals_test'         => 'mammals_test.nh'
  };
  $self->param( 'trees', $trees );

  # Omega distributions.
  my $omega_distributions = {
    anisimova_02_M3 => {
      omega_distribution_name     => "Anisimova 2002 M3",
      phylosim_omega_distribution => 'M3',
      phylosim_p0                 => 0.386,
      phylosim_p1                 => 0.535,
      phylosim_p2                 => 0.079,
      phylosim_w0                 => 0.018,
      phylosim_w1                 => 0.304,
      phylosim_w2                 => 1.691,
    },
    massingham_05_A => {
      omega_distribution_name     => "Massingham 2005-A (M8)",
      phylosim_omega_distribution => 'M8',
      phylosim_p0                 => 0.9432,
      phylosim_p                  => 0.572,
      phylosim_q                  => 2.172,
      phylosim_w                  => 2.081
    },
    massingham_05_A2 => {
      omega_distribution_name     => "Massingham 2005-A2 (M8, 15% pos-sel)",
      phylosim_omega_distribution => 'M8',
      phylosim_p0                 => 0.85,
      phylosim_p                  => 0.572,
      phylosim_q                  => 2.172,
      phylosim_w                  => 2.081
    },
    massingham_05_B => {
      omega_distribution_name     => "Massingham 2005-B (M3)",
      phylosim_omega_distribution => 'M3',
      phylosim_p0                 => 0.75,
      phylosim_p1                 => 0.25,
      phylosim_w0                 => 0.5,
      phylosim_w1                 => 1.5
    },
    neutral => {
      phylosim_omega_distribution => 'constant',
      phylosim_w                  => 1
    },
    purifying => {
      phylosim_omega_distribution => 'constant',
      phylosim_w                  => 0.2
    },
    lognormal_wide_2 => {
      phylosim_omega_distribution => 'lognormal',
      phylosim_meanlog            => -1.25,
      phylosim_sdlog              => 0.45
    },
    lognormal_wide => {
      phylosim_omega_distribution => 'lognormal',
      phylosim_meanlog            => -1.5,
      phylosim_sdlog              => 0.84
    },
    lognormal => {
      phylosim_omega_distribution => 'lognormal',
      phylosim_meanlog            => -1.864,
      phylosim_sdlog              => 1.2007
    },
    lognormal_narrow => {
      phylosim_omega_distribution => 'lognormal',
      phylosim_meanlog            => -2.25,
      phylosim_sdlog              => 1.49
    },
    lognormal_narrow_2 => {
      phylosim_omega_distribution => 'lognormal',
      phylosim_meanlog            => -3,
      phylosim_sdlog              => 1.925
    },
    lognormal_low => {
      phylosim_omega_distribution => 'lognormal',
      phylosim_meanlog            => -2.803,
      phylosim_sdlog              => 1.0
    },
    lognormal_low_2 => {
      phylosim_omega_distribution => 'lognormal',
      phylosim_meanlog            => -4,
      phylosim_sdlog              => 1.84
    },
    lognormal_high => {
      phylosim_omega_distribution => 'lognormal',
      phylosim_meanlog            => -1.819,
      phylosim_sdlog              => 1.5
    },
    lognormal_high_2 => {
      phylosim_omega_distribution => 'lognormal',
      phylosim_meanlog            => -1.0,
      phylosim_sdlog              => 0.784
    },
    uniform_neutral => {
      phylosim_omega_distribution => 'constant',
      phylosim_w                  => 1
    },
    uniform_low => {
      phylosim_omega_distribution => 'constant',
      phylosim_w                  => 0.3
    },
    uniform_low2 => {
      phylosim_omega_distribution => 'constant',
      phylosim_w                  => 0.1
    },
    uniform_pos => {
      phylosim_omega_distribution => 'constant',
      phylosim_w                  => 1.5
    },
    uniform_pos2 => {
      phylosim_omega_distribution => 'constant',
      phylosim_w                  => 2
    },

  };
  $self->param( 'omega_distributions', $omega_distributions );

  ## Phylo analysis settings.
  my $phylo_analyses = {
    none => {
      phylo_analysis_name => "None",
      analysis_action     => 'none',
    },
    slr => {
      phylo_analysis_name => "SLR",
      analysis_action     => "slr",
    },
    paml_m8 => {
      phylo_analysis_name => "PAML M8A/M8B",
      analysis_action     => "paml_sitewise paml_lrt",
      paml_model_a        => 'M8a',
      paml_model_b        => 'M8',
    },
    paml_m2 => {
      phylo_analysis_name => 'PAML M2/M3',
      analysis_action     => 'paml_sitewise paml_lrt',
      paml_model_a        => 'M2',
      paml_model_b        => 'M3',
    },
  };
  $self->param( 'phylo_analyses', $phylo_analyses );

  my $indel_models = {
    pow => {
      phylosim_insertmodel => 'POW 1.7 50',
      phylosim_deletemodel => 'POW 1.7 50'
    },
    nb => {
      phylosim_insertmodel => 'NB 0.35 2',
      phylosim_deletemodel => 'NB 0.35 2'
    },
    none => {
      phylosim_insertmodel => 'NB 0.35 2',
      phylosim_deletemodel => 'NB 0.35 2',
      phylosim_insertrate  => 0,
      phylosim_deleterate  => 0,
    }
  };
  $self->param( 'indel_models', $indel_models );
}

sub load_simulation_sets {
  my $self = shift;

}

sub mammals_simulations {
  my $self = shift;
  my $use_indels = shift || 0;

  my $n_replicates = 100;

  my $params = {
    slrsim_replicates => $n_replicates,
    experiment_name   => "mammals_simulations",

    #slrsim_tree_length => 1,
    phylosim_seq_length => 1000,
    slrsim_ref          => 'Homo_sapiens'
  };

  my $indel    = $self->param('indel_models')->{'none'};
  my $distr    = $self->param('omega_distributions')->{'lognormal'};
  my $filter   = $self->filter_param('none');
  my $analysis = $self->param('phylo_analyses')->{'slr'};
  my $aln      = $self->aln_param('true');

  my $base_params = $self->replace_params( $params, $indel, $distr, $filter, $analysis, $aln );

  my @sets = ();

# Go through each parameter set, get the median branch length, and scale our simulations by that parameter.
  my $dbname      = $self->compara_dba->dbc->dbname;
  my $compara_dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
    -url => 'mysql://ensadmin:ensembl@ens-research/gj1_2x_57' );
  my $tree_f       = $self->base . "/projects/slrsim/trees/2xmammals.nh";
  my $species_tree = Bio::EnsEMBL::Compara::TreeUtils->from_file($tree_f);

  my @param_sets = $self->get_parameter_sets('gj1_2x_57');

  #  @param_sets = @param_sets[1..2];
  #  @param_sets = ($param_sets[0]);
  foreach my $pset (@param_sets) {
    my $i      = $pset->{parameter_set_id};
    my $params = eval $pset->{params};
    my $subtree =
      Bio::EnsEMBL::Compara::ComparaUtils->get_species_subtree( $compara_dba, $species_tree,
      $params );

    my $tree_newick = $subtree->newick_format;
    my @tree_lengths = $self->_get_distribution_of_tree_lengths( $i, $n_replicates );

    my $shortname = $params->{parameter_set_shortname};
    print "Pset $shortname\n";

    my $pset_params = {
      parameter_set_name  => $params->{parameter_set_shortname},
      slrsim_tree_lengths => \@tree_lengths,
      slrsim_tree_newick  => $tree_newick
    };

    my $cur_params = $self->replace_params( $base_params, $pset_params );

    # Neutral sim.
    my $neutral_params = $self->clone($cur_params);
    $neutral_params =
      $self->replace_params( $neutral_params, $self->param('omega_distributions')->{'neutral'} );
    $neutral_params->{slrsim_label} = "$shortname noindel neutral";
    push @sets, $neutral_params;

    # No indels.
    my $noindel_params = $self->clone($cur_params);
    $noindel_params->{slrsim_label} = "$shortname noindel lnorm normal";
    push @sets, $noindel_params;

    # No indels, but wider omega distr.
    my $wide_params =
      $self->replace_params( $cur_params, $self->param('omega_distributions')->{'lognormal_wide'} );
    $wide_params->{slrsim_label} = "$shortname noindel lnorm wide";
    push @sets, $wide_params;

    # Very wide omega distr.
    my $v_wide_params =
      $self->replace_params( $cur_params,
      $self->param('omega_distributions')->{'lognormal_extreme'} );
    $v_wide_params->{slrsim_label} = "$shortname noindel lnorm vwide";
    push @sets, $v_wide_params;

    # Narrow omega distr.
    my $narrow_params =
      $self->replace_params( $cur_params,
      $self->param('omega_distributions')->{'lognormal_narrow'} );
    $narrow_params->{slrsim_label} = "$shortname noindel lnorm narrow";
    push @sets, $narrow_params;

    my $indel = $self->param('indel_models')->{'pow'};
    my $aln   = $self->aln_param('prank_f');

    # Low indels.
    my $indel_lo = $self->replace_params( $cur_params, $indel, $aln );
    $indel_lo->{phylosim_insertrate} = 0.01;
    $indel_lo->{phylosim_deleterate} = 0.01;
    $indel_lo->{slrsim_label}        = "$shortname indel 01";
    push @sets, $indel_lo;

    # Med indels.
    my $indel_mid = $self->replace_params( $cur_params, $indel, $aln );
    $indel_mid->{phylosim_insertrate} = 0.025;
    $indel_mid->{phylosim_deleterate} = 0.025;
    $indel_mid->{slrsim_label}        = "$shortname indel 025";
    push @sets, $indel_mid;

    # Med-hi indels.
    #my $indel_midhi = $self->replace_params($cur_params,$indel,$aln);
    #$indel_midhi->{phylosim_insertrate} = 0.03;
    #$indel_midhi->{phylosim_deleterate} = 0.03;
    #$indel_midhi->{slrsim_label} = "indel_03";
    #push @sets, $indel_midhi;

    # Med-hi indels.
    #$indel_midhi = $self->replace_params($cur_params,$indel,$aln);
    #$indel_midhi->{phylosim_insertrate} = 0.04;
    #$indel_midhi->{phylosim_deleterate} = 0.04;
    #$indel_midhi->{slrsim_label} = "indel_04";
    #push @sets, $indel_midhi;

    # High indels.
    #my $indel_hi = $self->replace_params($cur_params,$indel,$aln);
    #$indel_hi->{phylosim_insertrate} = 0.05;
    #$indel_hi->{phylosim_deleterate} = 0.05;
    #$indel_hi->{slrsim_label} = "indel_05";
    #push @sets, $indel_hi;

  }

# Store the base simulation parameters in the meta table. This will be later dumped by the Plots.pm script.
  $self->store_meta($base_params);

  return \@sets;
}

# Gets the median tree length from all the gene trees in a compara database.
sub _get_distribution_of_tree_lengths {
  my $self             = shift;
  my $parameter_set_id = shift;
  my $n_samples        = shift;

  use Bio::Greg::EslrUtils;
  my $script = Bio::Greg::EslrUtils->baseDirectory . "/scripts/collect_sitewise.R";

  my $rcmd = qq^
dbname = "gj1_2x_57"
host = 'ens-research'
port = 3306
user = 'ensro'
password = ''
script = "$script"
source(script)

genes = get.genes(parameter.set.id=$parameter_set_id)
good.genes = subset(genes,duplication_count < 3 & tree_length < 8 & tree_length > 0.01 & tree_length < 10 & !is.na(tree_length))
tree.lengths = good.genes[,'tree_length']
library(MASS)

# Fit a gamma distribution to the empirical set of branch lengths.
gamma.fit = fitdistr(tree.lengths,dgamma,list(shape=1.5,rate=1),lower=0.001)
est = as.vector(gamma.fit[['estimate']])
#print(est[1])
#print(est[2])
#median.length = median(tree.lengths)
#print(median(tree.lengths))

# Randomly sample N values from the gamma distribution.
print(rgamma($n_samples,est[1],est[2]))
#print(median.length)
#print(sd(tree.lengths))
^;
  my @values = Bio::Greg::EslrUtils->get_r_values( $rcmd, $self->worker_temp_directory );
  print "$parameter_set_id: [" . join( ",", @values ) . "]\n";
  return @values;
}

sub alignment_comparison {
  my $self = shift;

  my @sets = ();

  my $params = {
    slrsim_replicates   => 30,
    experiment_name     => "alignment_comparison",
    slrsim_tree_file    => $self->param('trees')->{'anisimova_bglobin'},
    slrsim_tree_mult    => 1.3,
    phylosim_seq_length => 400,
    slrsim_ref          => 'human'
  };

  my $indel       = $self->param('indel_models')->{'power_law'};
  my $distr       = $self->param('omega_distributions')->{'massingham_05_A2'};
  my $filter      = $self->filter_param('none');
  my $analysis    = $self->param('phylo_analyses')->{'slr'};
  my $base_params = $self->replace_params( $params, $indel, $distr, $filter, $analysis );

  # First, add the true alignment.
  my $p = $self->replace( $base_params, $self->aln_param('true') );
  push @sets, $p;

  # Now add the alignment sets.
  my @alns =
    map { $self->aln_param($_) }
    ( 'clustalw', 'muscle', 'cmcoffee', 'prank', 'prank_f', 'prank_codon', 'papaya' );
  foreach my $aln (@alns) {
    my $p = $self->replace( $base_params, $aln );
    push @sets, $p;
  }

# Store the overall simulation parameters in the meta table. This will be later dumped by the Plots.pm script.
  $self->store_meta($base_params);

  return \@sets;
}

sub reference_comparison {
  my $self = shift;

  my @sets = ();

  my $params = {
    slrsim_replicates   => 5,
    experiment_name     => "reference_comparison",
    slrsim_tree_file    => $self->param('trees')->{'2xmammals'},
    slrsim_tree_mult    => 1,
    phylosim_seq_length => 400,
    slrsim_ref          => 'human'
  };

  my $indel       = $self->param('indel_models')->{'pow'};
  my $distr       = $self->param('omega_distributions')->{'massingham_05_A2'};
  my $analysis    = $self->param('phylo_analyses')->{'slr'};
  my $filter      = $self->filter_param('none');
  my $base_params = $self->replace_params( $params, $indel, $distr, $analysis, $filter );

  # First, add the true alignment.
  my $p = $self->replace( $base_params, $self->species_param('human'), $self->aln_param('true') );
  push @sets, $p;

  # Now add the aligned sets, with different reference species.
  my @ref_params =
    map { $self->species_param($_) } ( 'human', 'mm9', 'loxAfr2', 'monDom4', 'xenTro2' );
  my @aln_params = map { $self->aln_param($_) } ( 'cmcoffee', 'prank_codon' );
  foreach my $ref (@ref_params) {
    foreach my $aln (@aln_params) {
      my $p = $self->replace( $base_params, $ref, $aln );
      push @sets, $p;
    }
  }

# Store the overall simulation parameters in the meta table. This will be later dumped by the Plots.pm script.
  $self->store_meta($base_params);

  return \@sets;
}

sub generic_sim_params {
  my $self = shift;

  my $base = {
    slrsim_replicates   => 10,
    experiment_name     => 'generic_sim',
    slrsim_tree_file    => $self->param('trees')->{'anisimova_bglobin'},
    slrsim_tree_mult    => 1,
    phylosim_seq_length => 400,
    slrsim_ref          => 'human'
  };

  my $indel    = $self->param('indel_models')->{'pow'};
  my $distr    = $self->param('omega_distributions')->{'massingham_05_A2'};
  my $analysis = $self->param('phylo_analyses')->{'slr'};

  my $aln  = $self->aln_param('fmcoffee');
  my $filt = $self->filter_param('tcoffee');

  my $params = $self->replace_params( $base, $indel, $distr, $analysis, $aln, $filt );
  return $params;
}

sub filter_test {
  my $self = shift;
  my @sets = ();
  my $p;

  my $base = $self->generic_sim_params;

  $p = $self->replace( $base, $self->aln_param('true'), $self->filter_param('none'), );
  push @sets, $p;

  $p = $self->replace(
    $base,
    $self->filter_thresh(0),
    $self->filter_order('before'),
    { alignment_score_mask_character_cdna => 'N' }
  );
  push @sets, $p;

  $p = $self->replace(
    $base,
    $self->filter_thresh(5),
    $self->filter_order('before'),
    { alignment_score_mask_character_cdna => 'N' }
  );
  push @sets, $p;

  $p = $self->replace(
    $base,
    $self->filter_thresh(9),
    $self->filter_order('before'),
    { alignment_score_mask_character_cdna => 'N' }
  );
  push @sets, $p;

  return \@sets;
}

sub filter_order_test {
  my $self   = shift;
  my @sets   = ();
  my $params = {
    slrsim_replicates   => 20,
    experiment_name     => 'filter_order_test',
    slrsim_tree_file    => $self->param('trees')->{'2x_primates'},
    slrsim_tree_mult    => 4,
    phylosim_seq_length => 400,
    slrsim_ref          => 'human'
  };
  my $indel    = $self->param('indel_models')->{'pow'};
  my $distr    = $self->param('omega_distributions')->{'massingham_05_A2'};
  my $analysis = $self->param('phylo_analyses')->{'slr'};
  $params = $self->replace_params( $params, $indel, $distr, $analysis );

  my $p = $self->replace( $params, $self->aln_param('true'), $self->filter_param('none') );
  push @sets, $p;

  my @aln_params = map { $self->aln_param($_) } ('cmcoffee');
  my @filter_params =
    map { $self->filter_param($_) }
    ('tcoffee');    #, 'tcoffee', 'indelign', 'prank_treewise', 'prank_mean' );
  my @filter_orders = map { { 'sitewise_filter_order' => $_ } } ( 'before', 'after' );
  my @threshold_params = map { { 'alignment_score_threshold' => $_ } } ( 0, 3, 5, 8, 9 );

  foreach my $aln (@aln_params) {
    foreach my $fp (@filter_params) {
      foreach my $fo (@filter_orders) {
        foreach my $tp (@threshold_params) {
          my $p = $self->replace( $params, $aln, $fp, $fo, $tp );
          push @sets, $p;
        }
      }
    }
  }

# Store the overall simulation parameters in the meta table. This will be later dumped by the Plots.pm script.
  $self->store_meta($params);

  return \@sets;
}

sub filter_sweep {
  my $self = shift;

  my @sets = ();

  my $final_params = {
    slrsim_replicates    => 20,
    experiment_name      => "filter_sweep",
    slrsim_tree_file     => $self->param('trees')->{'mammals_test'},
    slrsim_tree_max_path => 0.7,
    phylosim_seq_length  => 500,
    slrsim_ref           => 'human'
  };

  my $indel       = $self->param('indel_models')->{'pow'};
  my $distr       = $self->param('omega_distributions')->{'lognormal'};
  my $analysis    = $self->param('phylo_analyses')->{'slr'};
  my $base_params = $self->replace_params( $indel, $distr, $analysis );

  my $p = $self->replace(
    $base_params,
    $self->aln_param('true'),
    $self->filter_param('none'),
    $final_params
  );
  push @sets, $p;

  my @aln_params =
    map { $self->aln_param($_) } ('clustalw');    #,'prank', 'papaya', 'prank_f', 'prank_codon' );
  my @filter_params =
    map { $self->filter_param($_) }
    ('prank_column')
    ;    #, 'indelign', 'prank_mean' ); #, 'tcoffee', 'indelign', 'prank_treewise', 'prank_mean' );

  foreach my $aln (@aln_params) {
    foreach my $fp (@filter_params) {
      foreach my $threshold ( 0, 2, 4, 6, 8 ) {
        my $p =
          $self->replace( $base_params, $aln, $fp, { alignment_score_threshold => $threshold, },
          $final_params );
        push @sets, $p;
      }
    }
  }

# Store the overall simulation parameters in the meta table. This will be later dumped by the Plots.pm script.
  $self->store_meta( $self->replace( $base_params, $final_params ) );

  return \@sets;
}

####
####
####

### SLRSIM DEFAULT PARAMS
sub slrsim_defaults {
  my $self = shift;

  my $tree = $self->tree_param('anisimova_bglobin');
  $tree->{slrsim_tree_mean_path} = 1;

  my $indel_rate  = $self->indel_rate_param(0.1);
  my $indel_shape = $self->indel_shape_param('pow');

  my $reps   = $self->reps_param(20);
  my $length = $self->length_param(500);

  my $aln      = $self->aln_param('muscle');
  my $omegas   = $self->omega_param('lognormal_high');
  my $analysis = $self->analysis_param('slr');

  my $filter        = $self->filter_param('none');
  my $filter_thresh = $self->filter_thresh_param(0);

  return $self->replace(
    $tree, $indel_rate, $indel_shape, $reps,   $length,
    $aln,  $omegas,     $analysis,    $filter, $filter_thresh
  );
}

sub slrsim_all {
  my $self = shift;

  my @all_sets;
  my @cur_set;

  foreach my $sub_experiment ('slrsim_one_c') {
    @cur_set = @{ $self->$sub_experiment() };
    $self->apply_experiment_name_to_sets( \@cur_set, $sub_experiment );
    push @all_sets, @cur_set;
  }

  return \@all_sets;
}

# SLRsim one: indel rate and tree length sweep.
sub slrsim_one_a {
  my $self = shift;

  return $self->_generic_tree_indel_sweep( '44mammals', 1 );
}

sub slrsim_one_b {
  my $self = shift;

  return $self->_generic_tree_indel_sweep( 'anisimova_artificial', 1 );
}

sub slrsim_one_c {
  my $self = shift;

  return $self->_generic_tree_indel_sweep( 'anisimova_bglobin', 1 );
}

sub slrsim_one_d {
  my $self = shift;

  my $subset = [
    'Human',   'Chimp',      'Gorilla',  'Orangutan', 'Rhesus', 'Marmoset',
    'Tarsier', 'Mouselemur', 'Bushbaby', 'TreeShrew'
  ];
  my $tree_obj = $self->_get_tree_from_file( $self->tree_param('44mammals')->{slrsim_tree_file} );
  print $tree_obj->newick_format() . "\n";
  my $subtree = Bio::EnsEMBL::Compara::TreeUtils->subtree( $tree_obj, $subset );
  print "sub: " . $subtree->newick_format() . "\n";

  return $self->_generic_tree_indel_sweep( $subtree, 1 );
}

sub _generic_tree_indel_sweep {
  my $self             = shift;
  my $tree             = shift;    # Either a tree filename OR a NestedSet object.
  my $base_path_length = shift;

  my $reps        = $self->reps_param(20);
  my $tree_params = $self->tree_param($tree) if ( !ref $tree );
  my $defaults    = $self->slrsim_defaults;

  my $base_params = $self->replace( $defaults, $reps, $tree_params );

  my $tree_obj;
  if ( ref $tree && $tree->isa('Bio::EnsEMBL::Compara::NestedSet') ) {
    $tree_obj = $tree;
  } else {
    $tree_obj = $self->_get_tree_from_file( $base_params->{slrsim_tree_file} );
  }

  my @sets;
  my $params;
  my @lengths = map { $_ * $base_path_length / 5 } 1 .. 100;
  my @indels  = map { $_ * .025 } 0 .. 40;
  foreach my $path_length (@lengths) {
    $params = $self->replace( $base_params, { slrsim_tree_mean_path => $path_length } );
    foreach my $indel (@indels) {
      if ( $indel * $path_length > 0.4 || $path_length == 0 ) {
        print "Too much: $path_length $indel\n";
        next;
      }

      $tree_obj = Bio::EnsEMBL::Compara::TreeUtils->scale_mean_to( $tree_obj, $path_length );
      my $total_length = Bio::EnsEMBL::Compara::TreeUtils->total_distance($tree_obj);

      $params = $self->replace(
        $params, {
          phylosim_insertrate => $indel,
          phylosim_deleterate => $indel,
          slrsim_label        => $path_length . '|' . $indel,
          slrsim_tree         => $tree_obj,
        }
      );
      push @sets, $params;
    }
  }
  my $n = scalar(@sets);
  print "Tree count: $n\n";

  $self->store_meta($base_params);
  return \@sets;
}

# SLRsim two: Adding mammalian species to the tree.
sub slrsim_two {
  my $self = shift;

  my $reps     = $self->reps_param(10);
  my $defaults = $self->slrsim_defaults;

  my $base_params = $self->replace( $defaults, $reps );

  my $tree_param = $self->tree_param('44mammals');
  my $tree_obj   = $self->_get_tree_from_file( $tree_param->{slrsim_tree_file} );
  print $tree_obj->newick_format() . "\n";

  my @subsets;

  # Ingroup addition.
  map { push @subsets, [ 'Human', $_ ] } (

    # Ingroups.
    'Chimp', 'Gorilla', 'Orangutan', 'Rhesus', 'Marmoset', 'Tarsier', 'Mouselemur', 'Bushbaby',
    'TreeShrew',

    # Outgroups.
    'Lamprey',  'Zebrafish', 'Medaka', 'Tetraodon', 'Zebrafinch', 'Platypus', 'Opossum', 'Tenrec',
    'Hedgehog', 'Microbat',  'Dog',    'Mouse',     'Rabbit'
  );

  my @sets;
  foreach my $subset (@subsets) {
    $tree_obj = $tree_obj->copy;
    my $subtree = Bio::EnsEMBL::Compara::TreeUtils->subtree( $tree_obj, $subset );
    print "sub: " . $subtree->newick_format() . "\n";
    foreach my $indel_var ( 0, 0.05, 0.1, 0.2 ) {
      my $indel_p = $self->indel_rate_param($indel_var);
      push @sets,
        $self->replace(
        $base_params,
        $indel_p, {
          slrsim_tree_newick => $subtree->newick_format(),
          slrsim_label       => join( ',', @$subset, $indel_var )
        }
        );
    }
  }
  return \@sets;
}

sub slrsim_two_a {
  my $self = shift;

  return $self->slrsim_two_ab(0);
}

sub slrsim_two_b {
  my $self = shift;

  return $self->slrsim_two_ab(1);
}

sub slrsim_two_ab {
  my $self    = shift;
  my $reverse = shift;

  my $defaults = $self->slrsim_defaults;
  my $reps     = $self->reps(50);

  my $base_params = $self->replace( $defaults, $reps );

  my $tree_param = $self->tree_param('44mammals');
  my $tree_obj   = $self->_get_tree_from_file( $tree_param->{slrsim_tree_file} );
  print $tree_obj->newick_format() . "\n";

  my @subsets;

  # Ingroup addition.
  my @add_me_in_order = (
    'Chimp',      'Gorilla',  'Orangutan',   'Rhesus',    'Marmoset',  'Tarsier',
    'Mouselemur', 'Bushbaby', 'TreeShrew',   'Mouse',     'Squirrel',  'Rabbit',
    'Cow',        'Dog',      'Shrew',       'Tenrec',    'Opossum',   'Platypus',
    'Chicken',    'Lizard',   'Xtropicalis', 'Tetraodon', 'Zebrafish', 'Lamprey'
  );

  @add_me_in_order = reverse @add_me_in_order if ( $reverse == 1 );

  my @sets;
  my @subset_collector = ('Human');
  foreach my $addition (@add_me_in_order) {
    push @subset_collector, $addition;

    $tree_obj = $tree_obj->copy;
    my $subtree = Bio::EnsEMBL::Compara::TreeUtils->subtree( $tree_obj, \@subset_collector );
    print "sub: " . $subtree->newick_format() . "\n";
    foreach my $indel_var ( 0, 0.1 ) {
      my $indel_p = $self->indel_rate_param($indel_var);
      push @sets,
        $self->replace(
        $base_params,
        $indel_p, {
          slrsim_tree_newick => $subtree->newick_format(),
          slrsim_label       => join( ',', $indel_var, scalar( $subtree->leaves ), $addition )
        }
        );
    }
  }
  return \@sets;
}

# SLRsim 3: indel length distributions.
sub slrsim_three {
  my $self = shift;

  my $defaults = $self->slrsim_defaults;
  my $tree     = $self->tree_param('anisimova_bglobin');
  $tree->{slrsim_tree_mean_path} = 2;
  my $reps = {};    #$self->reps_param(3);

  my $base_params = $self->replace( $defaults, $tree, $reps );

  my @indel_params = (
    'NB .5 3',
    'NB .5 2',
    'NB .5 1',
    'NB .4 1',
    'NB .3 1',
    'NB .2 1',
    'POW 3 40',
    'POW 2 40',
    'POW 1.8 40',
    'POW 1.6 40',
    'POW 1.4 40',
    'POW 1.2 40',
    'POW 1.2 10',
    'POW 1.2 5'
  );

  my @sim_sets;
  foreach my $indel_shape (@indel_params) {
    my $indel_param = {
      phylosim_insertmodel => $indel_shape,
      phylosim_deletemodel => $indel_shape,
      slrsim_label         => $indel_shape
    };
    my $final = $self->replace( $base_params, $indel_param );
    push @sim_sets, $final;
  }
  return \@sim_sets;
}

# SLRsim four: different omega distributions.
sub slrsim_four {
  my $self = shift;

  my $defaults = $self->slrsim_defaults;
  my $reps     = {};                                       #$self->reps_param(20);
  my $tree     = $self->tree_param('anisimova_bglobin');

  my $base_params = $self->replace( $defaults, $reps, $tree );

  my @omega_dists = (
    'lognormal_narrow_2', 'lognormal_narrow', 'lognormal',       'lognormal_wide',
    'lognormal_wide_2',   'lognormal_low',    'lognormal_low_2', 'lognormal_high',
    'lognormal_high_2',   'uniform_neutral',  'uniform_low',     'uniform_low2',
    'uniform_pos',        'uniform_pos2'
  );

  my @sets;
  foreach my $omega_dist (@omega_dists) {
    my $omega_params = $self->omega_param($omega_dist);

    my $params = $self->replace( $base_params, $omega_params, { slrsim_label => $omega_dist } );
    push @sets, $params;
  }

  $self->store_meta($base_params);
  return \@sets;
}

# SLRsim five: different reference sequences.
sub slrsim_five {
  my $self = shift;

  my $defaults = $self->slrsim_defaults;

  my $reps = {};                                       #$self->reps(30);
  my $tree = $self->tree_param('anisimova_bglobin');

  my $base_params = $self->replace( $defaults, $reps, $tree );

  my @sets;
  foreach my $ref ( 'xenlaev', 'duck', 'rat', 'bushbaby', 'human' ) {

    #  foreach my $ref ( 'Human', 'Bushbaby', 'Mouse', 'Tenrec', 'Lizard', 'Xtropicalis', 'Fugu' ) {
    my $ref_params = $self->ref_param($ref);

    my $params = $self->replace( $base_params, $ref_params, { slrsim_label => $ref } );
    push @sets, $params;
  }

  $self->store_meta($base_params);
  return \@sets;
}

# SLRsim six: different alignment programs.
sub slrsim_six {
  my $self = shift;

  my $defaults = $self->slrsim_defaults;
  my $reps     = {};                       #$self->reps(20);

  my $base_params = $self->replace( $defaults, $reps );

  # In order from bad to worse.
  my @aligners =
    ( 'true', 'prank_f', 'prank', 'probcons', 'mcoffee', 'mafft', 'muscle', 'clustalw' );

  my @sets;
  foreach my $aligner (@aligners) {
    my $aln_params = $self->aln_param($aligner);
    my $params = $self->replace( $base_params, $aln_params, { slrsim_label => $aligner } );
    push @sets, $params;
  }
  return \@sets;
}

# SLRsim seven: filters.
sub slrsim_seven {
  my $self = shift;

  my $defaults = $self->slrsim_defaults;
  my $reps     = {};                       #$self->reps(20);

  my $base_params = $self->replace( $defaults, $reps );

  my @filters = ( 'none', 'indelign', 'prank_mean', 'tcoffee', 'trimal', 'gblocks', 'oracle' );

  my @sets;
  foreach my $filter (@filters) {
    my $filter_params = $self->filter_param($filter);
    my $filter_thresh = $self->filter_thresh_param(7);
    my $params =
      $self->replace( $base_params, $filter_params, $filter_thresh, { slrsim_label => $filter } );
    push @sets, $params;
  }

  # Add one item for the true alignment / no filtering.
  my $true_aln = $self->aln_param('true');
  my $true_params = $self->replace( $base_params, $true_aln, { slrsim_label => 'True Alignment' } );
  push @sets, $true_params;

  return \@sets;
}

####
####
####

sub apply_experiment_name_to_sets {
  my $selft           = shift;
  my $sets_arrayref   = shift;
  my $experiment_name = shift;

  my @new_sets = map { $_->{experiment_name} = $experiment_name } @$sets_arrayref;
  return \@new_sets;
}

sub tree_param {
  my $self      = shift;
  my $tree_name = shift;

  my $tree = $self->param('trees')->{$tree_name};

  return { slrsim_tree_file => $tree };
}

sub omega_param {
  my $self       = shift;
  my $omega_name = shift;

  my $distr_params = $self->param('omega_distributions')->{$omega_name};
  $self->hash_print($distr_params);
  return $distr_params;
}

sub analysis_param {
  my $self          = shift;
  my $analysis_name = shift;

  my $analysis_params = $self->param('phylo_analyses')->{$analysis_name};

  return $analysis_params;
}

sub aln_param {
  my $self       = shift;
  my $aln        = shift;
  my $aln_params = {
    alignment_name   => $aln,
    alignment_method => $aln,
    alignment_table  => 'aln'
  };
  if ( $aln eq 'none' || $aln eq 'true' ) {
    $aln_params = {
      alignment_name   => "True Alignment",
      alignment_method => 'none'
    };
  }
  return $aln_params;
}

sub reps_param {
  my $self = shift;
  my $reps = shift;

  return { slrsim_replicates => $reps };
}

sub length_param {
  my $self   = shift;
  my $length = shift;

  return { phylosim_seq_length => $length };
}

sub indel_rate_param {
  my $self       = shift;
  my $indel_rate = shift;

  return {
    phylosim_insertrate => $indel_rate / 2,
    phylosim_deleterate => $indel_rate / 2
  };
}

sub indel_shape_param {
  my $self       = shift;
  my $indel_name = shift;

  my $distr_params = $self->param('indel_models')->{$indel_name};

  if ( !defined $distr_params ) {
    $distr_params = {
      phylosim_insertmodel => $indel_name,
      phylosim_deletemodel => $indel_name
    };
  }
  return $distr_params;
}

sub filter_order_param {
  my $self  = shift;
  my $order = shift;

  return { sitewise_filter_order => $order };
}

sub filter_thresh_param {
  my $self   = shift;
  my $thresh = shift;

  return { alignment_score_threshold => $thresh };
}

sub filter_param {
  my $self   = shift;
  my $filter = shift;

  my $f = {
    filtering_name            => $filter,
    alignment_scores_action   => $filter,
    alignment_score_filtering => 1
  };

  if ( $filter eq 'none' ) {
    $f = {
      filtering_name            => 'None',
      alignment_scores_action   => 'none',
      alignment_score_filtering => 0,
      alignment_score_threshold => 0
    };
  }

  if ( $filter =~ m/after/i ) {
    $f->{sitewise_filter_order} = 'after';
  } else {
    $f->{sitewise_filter_order} = 'before';
  }
  return $f;
}

sub ref_param {
  my $self    = shift;
  my $species = shift;
  my $p       = { slrsim_ref => $species };
  return $p;
}

sub verify_params {
  my $self = shift;
  my $p    = shift;

  $p->{alignment_table}         = 'aln'   unless defined( $p->{alignment_table} );
  $p->{alignment_scores_action} = 'none'  unless defined( $p->{alignment_scores_action} );
  $p->{omega_table}             = 'omega' unless defined( $p->{omega_table} );

  # Check that we have proper filter and alignment tables.
  $p->{alignment_score_table} = $p->{alignment_table} . '_scores';
  $self->create_filter_table( $p->{alignment_score_table} );

  $self->create_aln_table( $p->{alignment_table} );
  $self->create_omega_table( $p->{omega_table} );
}

sub create_aln_table {
  my $self  = shift;
  my $table = shift;
  my $dba   = $self->compara_dba;
  $dba->dbc->do("CREATE TABLE IF NOT EXISTS $table LIKE protein_tree_member");
  $dba->dbc->do("TRUNCATE TABLE $table");
  print "CREATED TABLE $table\n";

}

sub create_omega_table {
  my $self  = shift;
  my $table = shift;
  my $dba   = $self->compara_dba;
  $dba->dbc->do("CREATE TABLE IF NOT EXISTS $table LIKE sitewise_omega");
  $dba->dbc->do("TRUNCATE TABLE $table");
}

sub create_filter_table {
  my $self  = shift;
  my $table = shift;
  $self->create_aln_table($table);
}

1;
