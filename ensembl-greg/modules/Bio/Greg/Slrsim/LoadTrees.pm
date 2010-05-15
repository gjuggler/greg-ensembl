package Bio::Greg::Slrsim::LoadTrees;

use strict;
use Cwd;
use Bio::EnsEMBL::Hive::Process;
use POSIX qw(strftime mktime);

use base ('Bio::Greg::Hive::Process');

sub fetch_input {
  my $self = shift;

  ### DEFAULT PARAMETERS ###
  my $defaults = {
    experiment_name => 'filter_sweeps',
    tree_root_dir   => '/homes/greg/lib/greg-ensembl/projects/slrsim/trees'
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

    throw("Experiment $experiment_name not found!") unless ( defined @simsets );

  my $tree_dir = $self->param('tree_root_dir') || '.';
  foreach my $params (@simsets) {
    $self->verify_params($params);

    my $replicates = $params->{'slrsim_replicates'};

    my $file = $tree_dir . '/' . $params->{'slrsim_tree_file'};
    print " -> Tree file: $file\n";
    open( IN, "$file" );
    my $newick_str = join( "", <IN> );
    $newick_str =~ s/\n//g;
    close(IN);
    print " -> Tree string: $newick_str\n";
    my $base_node = Bio::EnsEMBL::Compara::TreeUtils->from_newick($newick_str);

    foreach my $sim_rep ( 1 .. $replicates ) {
      $self->load_tree_into_database( $base_node, $sim_rep, $params );
    }
  }

  my $pta = $self->compara_dba->get_ProteinTreeAdaptor;
  foreach my $node ( @{ $pta->fetch_all_roots } ) {
    print $node->get_tagvalue("input_file") . "\t"
      . $node->node_id . "\t"
      . scalar( @{ $node->get_all_leaves } ) . "\n";
    printf "  > %.50s\n", $node->newick_format;
  }

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
  my $final_length = Bio::EnsEMBL::Compara::TreeUtils->total_distance($node);

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

  printf "  > %.50s\n", $node->newick_format;

  # Store all parameters as tags.
  $params->{'slrsim_rep'}         = $sim_rep;
  $params->{'slrsim_tree_length'} = $final_length;
  my $sim_param_str = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
  foreach my $tag ( keys %{$params} ) {
    $self->store_tag( $tag, $params->{$tag} );
  }

  # Store the whole parameter set as a string.
  $self->store_tag( "params_slrsim", $sim_param_str );

  my $output_params = {node_id => $node_id, parameter_set_id => $self->parameter_set_id, experiment_name => $self->param('experiment_name')};
  my $data_id = $self->new_data_id($output_params);

  my ($job_id) = @{$self->dataflow_output_id($output_params, 1)};
  print "  -> Created job: $job_id \n";

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
    '2xmammals'  => '2xmammals.nh',
    'ensembl_a'            => 'ensembl.nh',
    'ensembl_b'            => 'ensembl_2.nh'
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
      omega_distribution_name     => "Neutral evolution",
      phylosim_omega_distribution => 'constant',
      phylosim_w                  => 1
    },
    purifying => {
      omega_distribution_name     => "Purifying evolution",
      phylosim_omega_distribution => 'constant',
      phylosim_w                  => 0.2
    },
    lognormal => {
      omega_distribution_name     => "2xmammals Lognormal",
      phylosim_omega_distribution => 'lognormal',
      phylosim_meanlog            => -5.39,
      phylosim_sdlog              => 4.01

    }
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
    power_law => {
      phylosim_insertmodel => 'POW 1.8 50',
      phylosim_deletemodel => 'POW 1.8 50',
      phylosim_insertrate  => 0.05,
      phylosim_deleterate  => 0.05,
    },
    negative_binomial => {
      phylosim_insertmodel => 'NB 0.35 2',
      phylosim_deletemodel => 'NB 0.35 2',
      phylosim_insertrate  => 0.05,
      phylosim_deleterate  => 0.05,
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
  
  my $params = {
    slrsim_replicates => 100,
    experiment_name   => "mammals_simulations",
    slrsim_tree_mult => 1,
    phylosim_seq_length => 1000,
    slrsim_ref => 'human'
  };

  my $indel       = $self->param('indel_models')->{'none'};
  my $distr       = $self->param('omega_distributions')->{'lognormal'};
  my $filter      = $self->filter_param('none');
  my $analysis    = $self->param('phylo_analyses')->{'slr'};
  my $aln         = $self->aln_param('true');
  my $base_params = $self->replace_params($params, $indel, $distr, $filter, $analysis, $aln );

  my @sets = ();

  # First, add the true alignment.
  push @sets, $self->replace_params($base_params,$self->tree_param('2x_full'));
  push @sets, $self->replace_params($base_params,$self->tree_param('2x_nox'));
  push @sets, $self->replace_params($base_params,$self->tree_param('2x_primates'));
  push @sets, $self->replace_params($base_params,$self->tree_param('2x_glires'));

  # Store the overall simulation parameters in the meta table. This will be later dumped by the Plots.pm script.
  $self->store_meta($base_params);

  return \@sets;
}

sub alignment_comparison {
  my $self = shift;

  my @sets = ();

  my $params = {
    slrsim_replicates => 30,
    experiment_name   => "alignment_comparison",
    slrsim_tree_file  => $self->param('trees')->{'anisimova_bglobin'},
    slrsim_tree_mult => 1.3,
    phylosim_seq_length => 400,
    slrsim_ref => 'human'
  };
  
  my $indel       = $self->param('indel_models')->{'power_law'};
  my $distr       = $self->param('omega_distributions')->{'massingham_05_A2'};
  my $filter      = $self->filter_param('none');
  my $analysis    = $self->param('phylo_analyses')->{'slr'};
  my $base_params = $self->replace_params($params, $indel, $distr, $filter, $analysis );

  # First, add the true alignment.
  my $p = $self->replace(
    $base_params,
    $self->aln_param('true')
  );
  push @sets, $p;

  # Now add the alignment sets.
  my @alns = map { $self->aln_param($_) } ('clustalw', 'muscle', 'cmcoffee', 'prank', 'prank_f', 'prank_codon', 'papaya');
  foreach my $aln (@alns) {
    my $p = $self->replace($base_params, $aln);
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
    slrsim_replicates => 5,
    experiment_name   => "reference_comparison",
    slrsim_tree_file  => $self->param('trees')->{'2xmammals'},
    slrsim_tree_mult => 1,
    phylosim_seq_length => 400,
    slrsim_ref => 'human'
  };

  my $indel       = $self->param('indel_models')->{'power_law'};
  my $distr       = $self->param('omega_distributions')->{'massingham_05_A2'};
  my $analysis    = $self->param('phylo_analyses')->{'slr'};
  my $filter      = $self->filter_param('none');
  my $base_params = $self->replace_params($params, $indel, $distr, $analysis, $filter );

  # First, add the true alignment.
  my $p = $self->replace(
    $base_params,
    $self->species_param('human'),
    $self->aln_param('true')
  );
  push @sets, $p;

  # Now add the aligned sets, with different reference species.
  my @ref_params = map { $self->species_param($_) } ('human', 'mm9', 'loxAfr2', 'monDom4', 'xenTro2');
  my @aln_params = map { $self->aln_param($_) } ('cmcoffee','prank_codon');
  foreach my $ref (@ref_params) {
    foreach my $aln (@aln_params) {
      my $p = $self->replace($base_params, $ref, $aln);
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
    slrsim_replicates => 10,
    experiment_name => 'generic_sim',
    slrsim_tree_file => $self->param('trees')->{'anisimova_bglobin'},
    slrsim_tree_mult => 1,
    phylosim_seq_length => 400,
    slrsim_ref => 'human'
  };

  my $indel = $self->param('indel_models')->{'power_law'};
  my $distr       = $self->param('omega_distributions')->{'massingham_05_A2'};
  my $analysis    = $self->param('phylo_analyses')->{'slr'};

  my $aln = $self->aln_param('fmcoffee');
  my $filt = $self->filter_param('tcoffee');

  my $params = $self->replace_params( $base, $indel, $distr, $analysis, $aln, $filt);
  return $params;
}

sub filter_test {
  my $self = shift;
  my @sets = ();
  my $p;

  my $base = $self->generic_sim_params;

  $p = $self->replace(
    $base,
    $self->aln_param('true'),
    $self->filter_param('none'),
    );
  push @sets,$p;

  $p = $self->replace(
    $base,
    $self->filter_thresh(0),
    $self->filter_order('before'),
    {
      alignment_score_mask_character_cdna => 'N'
    }
    );
  push @sets,$p;

  $p = $self->replace(
    $base,
    $self->filter_thresh(5),
    $self->filter_order('before'),
    {
      alignment_score_mask_character_cdna => 'N'
    }
    );
  push @sets,$p;

  $p = $self->replace(
    $base,
    $self->filter_thresh(9),
    $self->filter_order('before'),
    {
      alignment_score_mask_character_cdna => 'N'
    }
    );
  push @sets,$p;

  return \@sets;
}

sub filter_order_test {
  my $self = shift;  
  my @sets = ();
  my $params = {
    slrsim_replicates => 20,
    experiment_name   => 'filter_order_test',
    slrsim_tree_file  => $self->param('trees')->{'2x_primates'},
    slrsim_tree_mult => 4,
    phylosim_seq_length => 400,
    slrsim_ref => 'human'
  };
  my $indel       = $self->param('indel_models')->{'power_law'};
  my $distr       = $self->param('omega_distributions')->{'massingham_05_A2'};
  my $analysis    = $self->param('phylo_analyses')->{'slr'};
  $params = $self->replace_params( $params, $indel, $distr, $analysis );

  my $p = $self->replace(
    $params,
    $self->aln_param('true'),
    $self->filter_param('none')
  );
  push @sets, $p;

  my @aln_params = map { $self->aln_param($_) } ( 'cmcoffee' );
  my @filter_params =
    map { $self->filter_param($_) } ( 'tcoffee' ); #, 'tcoffee', 'indelign', 'prank_treewise', 'prank_mean' );
  my @filter_orders = map {{'sitewise_filter_order' => $_}} ('before','after');
  my @threshold_params = map{{'alignment_score_threshold' => $_}} (0, 3, 5, 8, 9);

  foreach my $aln (@aln_params) {
    foreach my $fp (@filter_params) {
      foreach my $fo (@filter_orders) {
        foreach my $tp (@threshold_params) {
          my $p = $self->replace( $params, $aln, $fp, $fo, $tp);
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
    slrsim_replicates => 30,
    experiment_name   => "filter_sweep",
    slrsim_tree_file  => $self->param('trees')->{'anisimova_bglobin'},
    slrsim_tree_mult => 1.8,
    phylosim_seq_length => 400,
    slrsim_ref => 'human'
  };

  my $indel       = $self->param('indel_models')->{'power_law'};
  my $distr       = $self->param('omega_distributions')->{'massingham_05_A2'};
  my $analysis    = $self->param('phylo_analyses')->{'slr'};
  my $base_params = $self->replace_params( $indel, $distr, $analysis );

  my $p = $self->replace(
    $base_params,
    $self->aln_param('true'),
    $self->filter_param('none'),
    $final_params
  );
  push @sets, $p;

  my @aln_params = map { $self->aln_param($_) } ('fmcoffee'); #,'prank', 'papaya', 'prank_f', 'prank_codon' );
  my @filter_params =
    map { $self->filter_param($_) } ( 'indelign','tcoffee','trimal','prank_treewise_after'); #, 'indelign', 'prank_mean' ); #, 'tcoffee', 'indelign', 'prank_treewise', 'prank_mean' );

  foreach my $aln (@aln_params) {
    foreach my $fp (@filter_params) {
      foreach my $threshold (0, 2, 4, 6, 8) {
        my $p = $self->replace( $base_params, $aln, $fp, { alignment_score_threshold => $threshold, },
          $final_params );
        push @sets, $p;
      }
    }
  }

  # Store the overall simulation parameters in the meta table. This will be later dumped by the Plots.pm script.
  $self->store_meta($self->replace($base_params,$final_params));

  my $creation_tags = {
    timestamp => strftime("%Y-%m-%d",localtime),
  };
  $self->store_meta($creation_tags);

  return \@sets;
}

sub tree_param {
  my $self = shift;
  my $tree_name = shift;

  my $tree = $self->param('trees')->{$tree_name};
  
  return {slrsim_tree_file => $tree};
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

sub filter_order {
  my $self = shift;
  my $order = shift;

  return {sitewise_filter_order => $order};
}

sub filter_thresh {
  my $self = shift;
  my $thresh = shift;

  return {alignment_score_threshold => $thresh};
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
      alignment_score_filtering => 0
    };
  }

  if ( $filter =~ m/after/i) {
    $f->{sitewise_filter_order} = 'after';
  } else {
    $f->{sitewise_filter_order} = 'before';
  }
  return $f;
}

sub species_param {
  my $self    = shift;
  my $species = shift;
  my $p       = {
    slrsim_ref   => $species
  };
  return $p;
}

sub verify_params {
  my $self = shift;
  my $p    = shift;

  $p->{alignment_table}         = 'aln'            unless defined( $p->{alignment_table} );
  $p->{alignment_scores_action} = 'none'           unless defined( $p->{alignment_scores_action} );
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
