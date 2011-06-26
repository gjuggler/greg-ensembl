package Bio::Greg::Slrsim::LoadTrees;

use strict;
use Cwd;
use Time::HiRes;
use POSIX qw(strftime mktime);

use Bio::EnsEMBL::Hive::Process;
use Bio::EnsEMBL::Utils::Exception;

use File::Path;

use base ('Bio::Greg::Hive::Process');

sub param_defaults {
  my $self = shift;
  return {
    experiment_name => 'filter_sweeps',
    tree_root_dir   => $self->base . '/projects/slrsim/trees'
  };
}

sub fetch_input {
  my $self = shift;

  $self->load_all_params();
  $self->load_simulation_params();

  $self->_output_folder();
  mkpath( [ $self->_output_folder ]);

}

sub _output_folder {
  my $self = shift;
  return $self->get_output_folder('/homes/greg/lib/greg-ensembl/projects/slrsim/NO.backup/output');
}

sub run {
  my $self = shift;

  my @simsets;

  my $experiment_name = $self->param('experiment_name');
  if ( defined $experiment_name ) {

    # Tricky Perl: Call the method corresponding to the current experiment.
    @simsets = @{ $self->$experiment_name() };
  } else {
    throw("No experiment defined!");
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

  my $node_id = 1;
  foreach my $params (@simsets) {
    $self->verify_params($params);

    print "Loading params:\n";
    $self->hash_print($params);

    my $replicates = $params->{slrsim_replicates};
    $replicates = 1 if ( $replicates < 1 );

    my $base_node;
    if ( defined $params->{slrsim_tree} ) {
      $base_node = $params->{slrsim_tree};
      delete $params->{slrsim_tree};
    } elsif ( defined $params->{slrsim_tree_file} ) {
      $base_node = $self->_get_tree_from_file( $params->{'slrsim_tree_file'} );
    } elsif ( defined $params->{slrsim_tree_newick} ) {
      $base_node = Bio::EnsEMBL::Compara::TreeUtils->from_newick( $params->{slrsim_tree_newick} );
      delete $params->{slrsim_tree_newick};
    }

    my @tree_lengths;
    if ( defined $params->{slrsim_tree_lengths} ) {
      @tree_lengths = @{ $params->{slrsim_tree_lengths} };
    }

    my $tree_string = $base_node->newick_format;
    print "$tree_string\n";

    # Take the next tree length if we're provided with a list of lengths.
    if ( defined @tree_lengths ) {
      my $length = $tree_lengths[ 0 ];
      $params->{slrsim_tree_length} = $length;
    }
    
    my $p = {
      node_id => $node_id++,
      slrsim_replicates       => $replicates,
      tree_string      => $tree_string,
      parameter_set_id => $self->parameter_set_id,
    };
    my $data_id = $self->new_data_id($p);
    $p->{data_id} = $data_id;

    my $combined_params = $self->replace( $params, $p );
    
    my ($job_id) = @{ $self->dataflow_output_id( $combined_params, 1 ) };
    print "  -> Created slrsim job: $job_id \n";    
    sleep(0.05);
  }

  # For each unique slrsim_label, flow out a collect_slrsim job.
  my $label_hash;
  map {$label_hash->{$_->{slrsim_label}} = 1} @simsets;
  foreach my $label (keys %$label_hash) {
    my $params = {
      slrsim_label => $label
    };
    my ($job_id) = @{ $self->dataflow_output_id( $params, 2)};
    print "  -> Created collection job: $job_id \n";
  }

  my $len_ins_hash;
  map {$len_ins_hash->{$_->{slrsim_tree_length}.'_'.$_->{phylosim_insertrate}} = 1} @simsets;
  foreach my $len_ins (keys %$len_ins_hash) {
    my @toks = split('_', $len_ins);
    my $len = $toks[0];
    my $ins = $toks[1];
    my $params = {
      slrsim_tree_length => $len,
      phylosim_insertrate => $ins
    };
    my ($job_id) = @{ $self->dataflow_output_id( $params, 2)};
  }

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
    'bglobin'    => 'bglobin.nh',
    'artificial' => 'artificial.nh',
    'encode'     => 'encode.nh'
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
    anisimova_02_M3_hi => {
      omega_distribution_name     => "Anisimova 2002 M3",
      phylosim_omega_distribution => 'M3',
      phylosim_p0                 => 0.386,
      phylosim_p1                 => 0.535,
      phylosim_p2                 => 0.079,
      phylosim_w0                 => 0.018,
      phylosim_w1                 => 0.304,
      phylosim_w2                 => 4.739,
    },
    lognormal => {
      phylosim_omega_distribution => 'lognormal',
      phylosim_meanlog            => -1.864,
      phylosim_sdlog              => 1.2007
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
      slrsim_analysis_name => "None",
      analysis_action     => 'none',
    },
    slr => {
      slrsim_analysis_name => "SLR Sitewise",
      analysis_action     => "slr",
    },
    paml_m8 => {
      slrsim_analysis_name => "PAML Sitewise M8",
      analysis_action     => "paml_m8"
    },
    paml_m2 => {
      slrsim_analysis_name => 'PAML Sitewise M2',
      analysis_action     => 'paml_m2'
    },
  };
  $self->param( 'phylo_analyses', $phylo_analyses );

  my $indel_models = {
    pow => {
      phylosim_insertmodel => 'POW 1.8 40',
      phylosim_deletemodel => 'POW 1.8 40'
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

####
### SLRSIM DEFAULT PARAMS
sub slrsim_defaults {
  my $self = shift;

  my $tree = $self->tree_param('anisimova_bglobin');
  $tree->{slrsim_tree_mean_path} = 1;

  my $indel_rate = $self->indel_rate_param(0.1)
    ;    # The *combined* total indel rate (half goes to insertions, half to deletions)
  my $indel_shape = $self->indel_shape_param('pow');

  my $reps   = $self->reps_param(20);
  my $length = $self->length_param(500);

  my $aln      = $self->aln_param('muscle');
  my $omegas   = $self->omega_param('lognormal_high');
  my $analysis = $self->analysis_param('slr');

  my $filter        = $self->filter_param('none', 0.6);

  return $self->replace(
    $tree, $indel_rate, $indel_shape, $reps,   $length,
    $aln,  $omegas,     $analysis,    $filter
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

sub scheme_a {
  my $self = shift;

  my $tree_param  = $self->tree_param('artificial');
  my $omega_param = $self->omega_param('lognormal');
  my $ref_param = $self->ref_param('Human');
  my $params      = $self->replace( $self->slrsim_defaults, $tree_param, $ref_param, $omega_param );
  return $params;
}

sub scheme_b {
  my $self = shift;

  my $tree_param  = $self->tree_param('bglobin');
  my $omega_param = $self->omega_param('lognormal');
  my $ref_param = $self->ref_param('human');
  my $params      = $self->replace( $self->slrsim_defaults, $tree_param, $ref_param, $omega_param );
  return $params;
}

sub scheme_c {
  my $self = shift;

  my $tree_param  = $self->tree_param('encode');
  my $omega_param = $self->omega_param('lognormal');
  my $ref_param = $self->ref_param('Human');
  my $params      = $self->replace( $self->slrsim_defaults, $tree_param, $ref_param, $omega_param );
  return $params;
}

sub fig_zero {
  my $self = shift;

  my @sets;

  push @sets, @{$self->_fig_zero($self->scheme_a, 'slr')};
  push @sets, @{$self->_fig_zero($self->scheme_a, 'paml_m8')};
  push @sets, @{$self->_fig_zero($self->scheme_a, 'paml_m2')};
  push @sets, @{$self->_fig_zero($self->scheme_b, 'slr')};
  push @sets, @{$self->_fig_zero($self->scheme_b, 'paml_m8')};
  push @sets, @{$self->_fig_zero($self->scheme_b, 'paml_m2')};
  push @sets, @{$self->_fig_zero($self->scheme_c, 'slr')};
  push @sets, @{$self->_fig_zero($self->scheme_c, 'paml_m8')};
  push @sets, @{$self->_fig_zero($self->scheme_c, 'paml_m2')};

  my $n = scalar(@sets);
  print "Tree count: $n\n";

  $self->store_meta($sets[0]);
  return \@sets;
}

sub _fig_zero {
  my $self = shift;
  my $scheme = shift;
  my $analysis_type = shift;

  my $params = $scheme;
  $params = $self->replace($params, $self->analysis_param($analysis_type));
  $params = $self->replace($params, $self->reps_param(50));

  my @sets;
  my @lengths = map { $_ * 1 / 5 } (0.25, 0.5, 1 .. 10);
  foreach my $path_length (@lengths) {
    my $cur_params = $self->replace( $params, $self->mean_path_param($path_length) );    
    my $tree = $cur_params->{tree_name};
    my $tree_obj;
    if ( ref $tree && $tree->isa('Bio::EnsEMBL::Compara::NestedSet') ) {
      $tree_obj = $tree;
    } else {
      $tree_obj = $self->_get_tree_from_file( $cur_params->{slrsim_tree_file} );
    }    
    $tree_obj = Bio::EnsEMBL::Compara::TreeUtils->scale_mean_to( $tree_obj, $path_length );
    my $total_length = Bio::EnsEMBL::Compara::TreeUtils->total_distance($tree_obj);

    $cur_params = $self->replace(
      $cur_params, {
        phylosim_insertrate => 0,
        phylosim_deleterate => 0,
        slrsim_label        => "${tree}_${analysis_type}_${path_length}",
        slrsim_tree         => $tree_obj
      }
      );
    push @sets, $cur_params;
  }
  return \@sets;
}

sub fig_one_a {
  my $self = shift;

  my @sets;

  push @sets, @{$self->_fig_one($self->scheme_a, 'slr')};
  return \@sets;
}

sub fig_one_b {
  my $self = shift;

  my @sets;

  push @sets, @{$self->_fig_one($self->scheme_b, 'slr')};
  return \@sets;
}

sub fig_one_c {
  my $self = shift;

  my @sets;

  push @sets, @{$self->_fig_one($self->scheme_c, 'slr')};
  return \@sets;
}

sub _fig_one {
  my $self = shift;
  my $scheme = shift;
  my $analysis_type = shift;

  my @array;
  foreach my $aln ('true', 'clustalw', 'mafft', 'prank', 'prank_codon') {
    my $p = {};
    if ( $aln eq 'mafft_filter' ) {
      my $aln_p         = $self->aln_param('mafft');
      my $filter_p      = $self->filter_param('tcoffee', 0.6);
      $p = $self->replace( $scheme, $aln_p, $filter_p);
      $p->{alignment_name} = 'mafft_filter';
    } elsif ( $aln eq 'mafft_optimal' ) {
      my $aln_p         = $self->aln_param('mafft');
      my $filter_p      = $self->filter_param('oracle', 0.6);
      $p = $self->replace( $scheme, $aln_p, $filter_p);
      $p->{alignment_name} = 'mafft_optimal';
    } else {
      my $aln_p = $self->aln_param($aln);
      $p = $self->replace( $scheme, $aln_p );
    }

    # Add the analysis_action parameter.
    $p = $self->replace($p, $self->analysis_param($analysis_type));
    
    push @array, @{ $self->_fig_one_indel_sweep($p) };
  }
  return \@array;
}

sub _fig_one_indel_sweep {
  my $self   = shift;
  my $scheme = shift;

  my $reps        = $self->reps_param(100);
  my $base_params = $self->replace( $scheme, $reps);

  my $tree = $base_params->{tree_name};
  my $aln  = $base_params->{alignment_name};
  my $analysis = $base_params->{analysis_action};
  my $filt = $base_params->{filtering_name};

  my $tree_obj;
  if ( ref $tree && $tree->isa('Bio::EnsEMBL::Compara::NestedSet') ) {
    $tree_obj = $tree;
  } else {
    $tree_obj = $self->_get_tree_from_file( $base_params->{slrsim_tree_file} );
  }

  my @sets;
  my $params;
  my @lengths = map { $_ * 1 / 5 } 1 .. 10;
  my @ins_rates  = map { $_ / 100 } 0 .. 10;

  push @lengths, 0.1;
  push @lengths, 0.05;

  foreach my $path_length (@lengths) {
    $params = $self->replace( $base_params, $self->mean_path_param($path_length) );
    foreach my $ins_rate (@ins_rates) {
      $tree_obj = Bio::EnsEMBL::Compara::TreeUtils->scale_mean_to( $tree_obj, $path_length );
      my $total_length = Bio::EnsEMBL::Compara::TreeUtils->total_distance($tree_obj);

      $params = $self->replace(
        $params, {
          phylosim_insertrate => $ins_rate,
          phylosim_deleterate => $ins_rate,
          slrsim_label        => "${tree}_${analysis}_${aln}_${path_length}_${ins_rate}",
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

sub fig_two_short {
  my $self = shift;
  return $self->fig_two(1);

}

sub fig_two_a {
  my $self = shift;  
  my @sets;
  push @sets, @{$self->fig_two('mafft', $self->scheme_b)};
  return \@sets;
}

sub fig_two_b {
  my $self = shift;
  my @sets;
  push @sets, @{$self->fig_two('prank', $self->scheme_b)};
  return \@sets;
}

sub fig_two_c {
  my $self = shift;  
  my @sets;
  push @sets, @{$self->fig_two('prank_codon', $self->scheme_b)};
  return \@sets;
}

sub fig_two_d {
  my $self = shift;  
  my @sets;
  push @sets, @{$self->fig_two('probcons', $self->scheme_b)};
  return \@sets;
}

sub fig_two_e {
  my $self = shift;  
  my @sets;
  push @sets, @{$self->fig_two('pagan_codon', $self->scheme_b)};
  return \@sets;
}

sub fig_two_clustalw {
  my $self = shift;  
  my @sets;
  push @sets, @{$self->fig_two('clustalw', $self->scheme_b)};
  return \@sets;
}

sub fig_two_fsa {
  my $self = shift;  
  my @sets;
  push @sets, @{$self->fig_two('fsa', $self->scheme_b)};
  return \@sets;
}

sub fig_two {
  my $self = shift;
  my $aln = shift;
  my $scheme = shift;

  my @sets;
  my @subsets_array;
  my @pairs = ([1, 0.05],
               [1, 0.1],
               [2, 0.1],
               [2, 0.2]);
  foreach my $pair (@pairs) {
    my ($length, $indel) = @$pair;
    push @subsets_array, {
      aln => $aln,
      scheme => $scheme,
      length => $length,
      indel => $indel
    };
  }

  foreach my $subset ( @subsets_array ) {
    $self->hash_print($subset);

    my $scheme = $subset->{scheme};
    my $aln    = $self->aln_param( $subset->{aln} );
    my $aln_s = $subset->{aln};
    my $length = $self->mean_path_param( $subset->{length} );
    my $indel  = $self->indel_rate_param( $subset->{indel} );

    foreach my $filter ( 'true', 'none', 'gblocks', 'branchlength_lo', 'branchlength_hi', 'optimal_lo', 'optimal_hi', 'tcoffee_lo', 'tcoffee_hi', 'guidance_lo', 'guidance_hi' ) {
      my $tree = $scheme->{tree_name};
      my $len_indel_s = '_'.$subset->{length} . '_'. $subset->{indel};
      
      my $max_mask = 0.6;
      if ($filter =~ m/lo/i) {
        $max_mask = 0.2;
      }
      
      my $params;
      if ( $filter eq 'true' ) {
        $params = $self->replace(
          $scheme,
          $self->aln_param('true'),
          $self->filter_param('true', $max_mask),
          { slrsim_label => "${tree}_True_none_${len_indel_s}" }
          );
        $params->{filtering_name} = 'true';
        $params->{alignment_name} = $aln->{alignment_name};
      } else {
        $params = $self->replace(
          $scheme, $aln,
          $self->filter_param($filter, $max_mask),
          { slrsim_label => "${tree}_${aln_s}_${filter}_${len_indel_s}" }
          );
      }
      
      $params =
        $self->replace( 
          $params, 
          $self->reps_param(50),
          $indel, 
          $length,
          { experiment_name => 'fig_two' } );
      push @sets, $params;
      
      $self->store_meta($params) if ( scalar(@sets) == 1 );
    }
  }
  return \@sets;
}

sub fig_three_a {
  my $self = shift;  
  my @sets;
  push @sets, @{$self->fig_three('prank_codon', $self->scheme_b)};
  return \@sets;
}

sub fig_three_b {
  my $self = shift;  
  my @sets;
  push @sets, @{$self->fig_three('clustalw', $self->scheme_b)};
  return \@sets;
}

sub fig_three {
  my $self = shift;
  my $aln = shift;
  my $scheme = shift;

  my @sets;
  my @subsets_array;

  my $params;
  my @lengths = map { $_ * 1 / 2.5 } 1 .. 5;
  my @ins_rates  = map { $_ / 50 } 0 .. 5;

  my @pairs;
  foreach my $length (@lengths) {
    foreach my $ins_rate (@ins_rates) {
      push @pairs, [$length, $ins_rate];
    }
  }

  foreach my $pair (@pairs) {
    my ($length, $ins_rate) = @$pair;
    push @subsets_array, {
      aln => $aln,
      scheme => $scheme,
      length => $length,
      ins_rate => $ins_rate
    };
  }

  foreach my $subset ( @subsets_array ) {
    $self->hash_print($subset);

    my $scheme = $subset->{scheme};
    my $aln    = $self->aln_param( $subset->{aln} );
    my $aln_s = $subset->{aln};
    my $length = $self->mean_path_param( $subset->{length} );
    my $indel = {
      phylosim_insertrate => $subset->{ins_rate},
      phylosim_deleterate => $subset->{ins_rate}
    };

    foreach my $filter ( 'true', 'none', 'gblocks', 'optimal_a', 'optimal_b', 'optimal_c', 'tcoffee', 'guidance', 'branchlength' ) {
#    foreach my $filter ( 'true', 'none', 'gblocks', 'branchlength_lo', 'optimal_lo', 'tcoffee_lo') {
      my $tree = $scheme->{tree_name};
      my $len_indel_s = $subset->{length} . '_'. $subset->{ins_rate};
      
      my $max_mask = 0.2;
      $max_mask = 0.05 if ($filter eq 'optimal_a');
      $max_mask = 0.1 if ($filter eq 'optimal_b');
      $max_mask = 0.2 if ($filter eq 'optimal_c');
      
      my $params;
      $params = $self->replace(
        $scheme, $aln,
        $self->filter_param($filter, $max_mask),
        { slrsim_label => "${tree}_${aln_s}_${filter}_${len_indel_s}" }
        );
      
      if ( $filter eq 'true' ) {
        $params = $self->replace($params, $self->aln_param('true'));
      }
      
      $params =
        $self->replace(
          $params,
          $self->reps_param(100),
          $indel,
          $length,
          { experiment_name => 'fig_three' } );
      push @sets, $params;
      
      $self->store_meta($params) if ( scalar(@sets) == 1 );
    }
  }
  return \@sets;
}


sub ari_indels {
  my $self = shift;

  my $scheme = $self->scheme_b;
  my $analysis_type = 'slr';

  my @array;
  foreach my $aln ('true', 
                   'mafft', 
                   'probcons',
                   'fsa',
                   'fsa_careful',
                   'prank', 
                   'prank_codon', 
                   'pagan', 
                   'pagan_codon',
                   ) {
    my $aln_p = $self->aln_param($aln);
    my $p = $self->replace( $scheme, $aln_p );

    # Add the analysis_action parameter.
    $p = $self->replace($p, $self->analysis_param($analysis_type));
    push @array, @{ $self->_ari_indel_sweep($p) };
  }
    
  return \@array;
}

sub _ari_indel_sweep {
  my $self   = shift;
  my $scheme = shift;

  my $reps        = $self->reps_param(40);
  my $base_params = $self->replace( $scheme, $reps);

  my $tree = $base_params->{tree_name};
  my $aln  = $base_params->{alignment_name};
  my $analysis = $base_params->{analysis_action};

  my $tree_obj;
  if ( ref $tree && $tree->isa('Bio::EnsEMBL::Compara::NestedSet') ) {
    $tree_obj = $tree;
  } else {
    $tree_obj = $self->_get_tree_from_file( $base_params->{slrsim_tree_file} );
  }

  my @sets;
  my $params;
  my @pairs = ([0.5, 0.2],
               [1, 0.05],
               [1, 0.1],
               [2, 0.1],
               [2, 0.2],
               [4, 0.05]);

  foreach my $pair (@pairs) {
    my ($length, $indel) = @$pair;
    $params = $self->replace( $base_params, $self->mean_path_param($length) );
    $tree_obj = Bio::EnsEMBL::Compara::TreeUtils->scale_mean_to( $tree_obj, $length );
    my $total_length = Bio::EnsEMBL::Compara::TreeUtils->total_distance($tree_obj);

    $params = $self->replace(
      $params, {
        phylosim_insertrate => $indel,
        phylosim_deleterate => $indel,
        slrsim_label        => "${tree}_${analysis}_${aln}_${length}_${indel}",
        slrsim_tree         => $tree_obj,
      }
    );
    push @sets, $params;
  }
  my $n = scalar(@sets);
  print "Tree count: $n\n";

  $self->store_meta($base_params);
  return \@sets;
}

sub mammals_sim {
  my $self = shift;

  my @sets;
  push @sets, $self->_mammals_sets('mammals_eutherian.nh', 'eutherian');
  push @sets, $self->_mammals_sets('mammals_full_genomes.nh', 'full_genomes');
  push @sets, $self->_mammals_sets('mammals_hmrd.nh', 'hmrd');
  return \@sets;
}

sub _mammals_sets {
  my $self = shift;
  my $tree_file = shift;
  my $label = shift;

  my $tree_obj = $self->_get_tree_from_file( $tree_file );

  my $mean_path = Bio::EnsEMBL::Compara::TreeUtils->mean_path($tree_obj);

  my @sets;
  foreach my $indel (0, 0.05) {
    my $scheme = $self->scheme_b;
    my $params = $self->replace(
      $scheme,
      $self->reps_param(100),
      $self->length_param(500),
      $self->aln_param('prank_codon'),
      $self->indel_rate_param($indel), {
        slrsim_tree => $tree_obj,
        slrsim_tree_file => $tree_file,
        slrsim_label => "${label}_${indel}",
        tree_mean_path => $mean_path
      }
      );
    push @sets, $params;
  }
  
  print $tree_obj->ascii;
  return @sets;
}

sub phd_anim {
  my $self = shift;

  my @sets;

  my $scheme = $self->scheme_b;
  my @lengths = map { $_ * 1 / 10 } 1 .. 20;
  my @indels = (0, 0.01, 0.02, 0.1);

  foreach my $indel (@indels) {
    foreach my $length(@lengths) {
      foreach my $aln ('true', 'mafft') {
        my $params = $self->replace($scheme,
                                    $self->analysis_param('none'),
                                    $self->reps_param(1),
                                    $self->length_param(100),
                                    $self->aln_param($aln),
                                    $self->indel_rate_param($indel),
                                    $self->mean_path_param($length)
          );
        $params->{random_seed} = 1;
        $params->{aln_output_id} = "${aln}_${indel}_${length}";
        push @sets, $params;
      }
    }
  }
  return \@sets;
}

sub phd_short {
  my $self = shift;

  my $subsets = [
      {
        aln    => 'clustalw',
        scheme => $self->scheme_c,
        length => 0.6,
        indel  => 0.06
      }, {
        aln    => 'clustalw',
        scheme => $self->scheme_c,
        length => 1.0,
        indel  => 0.1
      }, {
        aln    => 'clustalw',
        scheme => $self->scheme_c,
        length => 1.5,
        indel  => 0.14
      }
      ];

  my @sets;

  foreach my $subset ( @{$subsets} ) {
    $self->hash_print($subset);

    my $scheme = $subset->{scheme};
    my $aln    = $self->aln_param( $subset->{aln} );
    my $aln_s = $subset->{aln};
    my $length = $self->mean_path_param( $subset->{length} );
    my $indel  = $self->indel_rate_param( $subset->{indel} );

    my $params = $self->replace($scheme,$aln,$length,$indel,
                                    $self->analysis_param('none'),
                                    $self->reps_param(1),
                                    $self->length_param(50)
      );
    $params->{random_seed} = 1;
    $params->{aln_output_id} = $subset->{aln}.'_'.$subset->{indel}.'_'.$subset->{length};
    
    push @sets, $params;
  }    
  return \@sets;
}

sub apply_experiment_name_to_sets {
  my $self           = shift;
  my $sets_arrayref   = shift;
  my $experiment_name = shift;

  my @new_sets = map { $_->{experiment_name} = $experiment_name } @$sets_arrayref;
  return \@new_sets;
}

sub tree_param {
  my $self      = shift;
  my $tree_name = shift;

  my $tree = $self->param('trees')->{$tree_name};

  return {
    slrsim_tree_file => $tree,
    tree_name        => $tree_name
  };
}

sub omega_param {
  my $self       = shift;
  my $omega_name = shift;

  my $distr_params = $self->param('omega_distributions')->{$omega_name};
  #$self->hash_print($distr_params);
  return $distr_params;
}

sub mean_path_param {
  my $self = shift;
  my $path = shift;

  return { slrsim_tree_mean_path => $path };
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
    aligner => $aln,
    alignment_name   => $aln,
    alignment_method => $aln,
    alignment_table  => 'aln'
  };
  if ( $aln eq 'none' || $aln eq 'true' ) {
    $aln_params = {
      aligner => 'none',
      alignment_name   => "True_Alignment",
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
  my $max_mask = shift;

  die("Need to specify max mask fraction") unless (defined $max_mask);
  
  my $f = {
    filter => $filter,
    alignment_scores_action   => $filter,
    alignment_score_filtering => 1,
    maximum_mask_fraction => $max_mask
  };

  if ( $filter eq 'none' || $filter eq 'true' ) {
    $f = {
      filter => $filter,
      alignment_scores_action   => 'none',
      alignment_score_filtering => 0,
      alignment_score_threshold => 0,
      maximum_mask_fraction => 0
    };
  }

  if ( $filter eq 'gblocks' ) {
    $f->{maximum_mask_fraction} = 1;
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
