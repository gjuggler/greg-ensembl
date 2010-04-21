package Bio::Greg::Slrsim::LoadTrees;

use strict;
use Cwd;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::Process');

sub fetch_input {
  my $self = shift;

  ### DEFAULT PARAMETERS ###
  my $defaults = {
    experiment_name => 'Slrsim Experiment',
    simulation_set  => 'slrsim_test',
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

  my $sim_set_name = $self->param('simulation_set');
  if ( defined $sim_set_name ) {
    @simsets = @{ $self->param($sim_set_name) };
    throw("Simulation set $sim_set_name not found!") unless ( defined @simsets );
  }

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

  my $tree_length = $params->{'slrsim_tree_length'};
  if ($tree_length) {
    $node = Bio::EnsEMBL::Compara::TreeUtils->scale_to( $node, $tree_length );
  }
  my $tree_mult = $params->{'slrsim_tree_mult'};
  if ($tree_mult) {
    $node = Bio::EnsEMBL::Compara::TreeUtils->scale( $node, $tree_mult );
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

  my $output_params = {node_id => $node_id};
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
    '2xmammals_autosomal'  => '2xmammals.nh',
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
      phylosim_meanlog            => -1.7,
      phylosim_sdlog              => 1.23

    }
  };
  $self->param( 'omega_distributions', $omega_distributions );

  ## Phylo analysis settings.
  my $phylo_analyses = {
    none => {
      phylo_analysis_name => "None",
      sitewise_action     => 'none',
    },
    slr => {
      phylo_analysis_name => "SLR",
      sitewise_action     => "slr",
    },
    paml_m8 => {
      phylo_analysis_name => "PAML M8A/M8B",
      sitewise_action     => "paml_sitewise paml_lrt",
      paml_model_a        => 'M8a',
      paml_model_b        => 'M8',
    },
    paml_m2 => {
      phylo_analysis_name => 'PAML M2/M3',
      sitewise_action     => 'paml_sitewise paml_lrt',
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
      phylosim_deleterate  => 0.05
    }
  };
  $self->param( 'indel_models', $indel_models );
}

sub load_simulation_sets {
  my $self = shift;

  $self->param( 'filter_sweeps', $self->filter_sweeps() );
}

sub filter_sweeps {
  my $self = shift;
  
  my @sets = ();

  my $final_params = {
    slrsim_replicates => 1,
    experiment_name   => "Filter Sweeps",
    slrsim_tree_file  => $self->param('trees')->{'44mammals'},
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

  my @aln_params = map { $self->aln_param($_) } ( 'muscle' );
  my @filter_params =
    map { $self->filter_param($_) } ( 'tcoffee', 'indelign', 'prank_treewise', 'prank_mean' );

  foreach my $aln (@aln_params) {
    foreach my $fp (@filter_params) {
      foreach my $threshold ( 5 ) {
        my $p = $self->replace( $base_params, $aln, $fp, { alignment_score_threshold => $threshold, },
          $final_params );
        push @sets, $p;
      }
    }
  }
  return \@sets;
}

sub aln_param {
  my $self       = shift;
  my $aln        = shift;
  my $aln_params = {
    alignment_name   => $aln,
    alignment_method => $aln,
    alignment_table  => 'aln',
    omega_table      => 'omega'
  };
  if ( $aln eq 'none' || $aln eq 'true' ) {
    $aln_params = {
      alignment_name   => "True Alignment",
      alignment_method => 'none',
      alignment_table  => 'protein_tree_member',
      omega_table      => 'omega_true'
    };
  }
  return $aln_params;
}

sub filter_param {
  my $self   = shift;
  my $filter = shift;

  my $f = {
    filtering_name            => $filter,
    alignment_scores_action   => $filter,
    alignment_score_filtering => 1,
    alignment_score_threshold => 5
  };

  if ( $filter eq 'none' ) {
    $f = {
      filtering_name            => 'None',
      alignment_scores_action   => 'none',
      alignment_score_filtering => 0
    };
  }
  return $f;
}

sub species_param {
  my $self    = shift;
  my $species = shift;
  my $p       = {
    species_name => $species,
    slrsim_ref   => $species
  };
  return $p;
}

sub verify_params {
  my $self = shift;
  my $p    = shift;

  $p->{alignment_table}         = 'aln'            unless defined( $p->{alignment_table} );
  $p->{alignment_scores_action} = 'none'           unless defined( $p->{alignment_scores_action} );
  $p->{omega_table}             = 'sitewise_omega' unless defined( $p->{omega_table} );

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
