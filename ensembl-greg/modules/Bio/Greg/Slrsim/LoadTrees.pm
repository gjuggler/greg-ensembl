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
  mkpath([$self->_output_folder."/alns"]);

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

    foreach my $sim_rep ( 1 .. $replicates ) {
      # Take the next tree length if we're provided with a list of lengths.
      if ( defined @tree_lengths ) {
        my $length = $tree_lengths[ $sim_rep - 1 ];
        $params->{slrsim_tree_length} = $length;
      }

      my $p = {
        slrsim_rep => $sim_rep,
        tree_string => $tree_string,
        parameter_set_id => $self->parameter_set_id,
      };
      my $data_id = $self->new_data_id($p);
      $p->{data_id} = $data_id;

      my $combined_params = $self->replace($params,$p);

      my ($job_id) = @{ $self->dataflow_output_id( $combined_params, 1 ) };
      print "  -> Created job: $job_id \n";
      
      #$self->load_tree_into_database( $base_node, $sim_rep, $params );
      sleep(0.2);
    }
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
    'bglobin' => 'bglobin.nh',
    'artificial' => 'artificial.nh',
    'encode' => 'encode.nh'
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

  my $indel_rate  = $self->indel_rate_param(0.1); # The *combined* total indel rate (half goes to insertions, half to deletions)
  my $indel_shape = $self->indel_shape_param('pow');

  my $reps   = $self->reps_param(20);
  my $length = $self->length_param(500);
  
  my $aln      = $self->aln_param('muscle');
  my $omegas   = $self->omega_param('lognormal_high');
  my $analysis = $self->analysis_param('slr');
  
  my $filter = $self->filter_param('none');
  my $filter_thresh = $self->filter_thresh_param(0);
  
  return $self->replace( $tree, $indel_rate, $indel_shape, $reps, $length, $aln, $omegas,
                         $analysis, $filter, $filter_thresh );
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

sub man_a {
  my $self = shift;

  my $tree_param = $self->tree_param('artificial');
  my $omega_param = $self->omega_param('anisimova_02_M3_hi');
  my $params = $self->replace($self->slrsim_defaults,$tree_param,$omega_param);
  return $params;
}

sub man_b {
  my $self = shift;

  my $tree_param = $self->tree_param('bglobin');
  my $omega_param = $self->omega_param('anisimova_02_M3');
  my $params = $self->replace($self->slrsim_defaults,$tree_param,$omega_param);
  return $params;
}

sub man_c {
  my $self = shift;

  my $tree_param = $self->tree_param('encode');
  my $omega_param = $self->omega_param('lognormal');
  my $params = $self->replace($self->slrsim_defaults,$tree_param,$omega_param);
  return $params;
}

sub fig_one {
  my $self = shift;

  my @array;
  foreach my $scheme ($self->man_a) {
#  foreach my $scheme ($self->man_a,$self->man_b,$self->man_c) {
#    foreach my $aln ('true','clustalw','muscle','mafft','prank','prank_codon','prank_codon_filter') {
#    foreach my $aln ('true', 'clustalw','mafft','prank','prank_filter') {
    foreach my $aln ('true','clustalw','mafft','prank') {
      my $p = {};
      if ($aln eq 'mafft_filter') {
        my $aln_p = $self->aln_param('mafft');
        my $filter_p = $self->filter_param('tcoffee');
        my $filter_thresh = $self->filter_thresh_param(9);
        $p = $self->replace($scheme,$aln_p,$filter_p,$filter_thresh);
        $p->{alignment_name} = 'mafft_filter';
      } elsif ($aln eq 'mafft_optimal') {
        my $aln_p = $self->aln_param('mafft');
        my $filter_p = $self->filter_param('oracle');
        my $filter_thresh = $self->filter_thresh_param(9);
        $p = $self->replace($scheme,$aln_p,$filter_p,$filter_thresh);
        $p->{alignment_name} = 'mafft_optimal';
      } else {
        my $aln_p = $self->aln_param($aln);
        $p = $self->replace($scheme,$aln_p);
      }

      push @array,@{$self->_fig_one_indel_sweep($p)};
    }
  }
  return \@array;
}

sub _fig_one_indel_sweep {
  my $self             = shift;
  my $scheme = shift;

  my $reps        = $self->reps_param(50);
  my $base_params = $self->replace( $scheme, $reps);

  my $tree = $base_params->{tree_name};
  my $aln = $base_params->{alignment_name};
  my $filt = $base_params->{filtering_name};
  if ($filt) {
    $base_params->{experiment_name} = "1_${tree}_${aln}_${filt}";
  } else {
    $base_params->{experiment_name} = "1_${tree}_${aln}";
  }

  my $tree_obj;
  if ( ref $tree && $tree->isa('Bio::EnsEMBL::Compara::NestedSet') ) {
    $tree_obj = $tree;
  } else {
    $tree_obj = $self->_get_tree_from_file( $base_params->{slrsim_tree_file} );
  }

  my @sets;
  my $params;
  my @lengths = map { $_ * 1 / 5 } 1 .. 10;
  my @indels  = map { $_ / 100 } 0 .. 10;
#  my @lengths = (1.5);
#  my @indels = (0.07);
  foreach my $path_length (@lengths) {
    $params = $self->replace( $base_params, $self->mean_path_param($path_length) );
    foreach my $indel (@indels) {
      $tree_obj = Bio::EnsEMBL::Compara::TreeUtils->scale_mean_to( $tree_obj, $path_length );
      my $total_length = Bio::EnsEMBL::Compara::TreeUtils->total_distance($tree_obj);

      $params = $self->replace(
        $params, {
          phylosim_insertrate => $indel,
          phylosim_deleterate => $indel,
          slrsim_label        => "${tree} ${path_length} ${indel}",
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

sub fig_two {
  my $self = shift;

  my @sets;
  foreach my $scheme ($self->man_b) {
#  foreach my $scheme ($self->man_a, $self->man_b, $self->man_c) {
    foreach my $aln ('mafft','prank') {
      foreach my $filter ('none','gblocks','prank','optimal','tcoffee') {
        foreach my $filter_threshold (9) {
          my $tree = $scheme->{tree_name};

          my $params;
          if ($filter eq 'true') {
            $params = $self->replace(
              $scheme,
              $self->aln_param('true'),
              $self->filter_param('none'),
              {slrsim_label => "${tree}_True Alignment_${filter}"}
              );
          } else {
            $params = $self->replace(
              $scheme,
              $self->aln_param($aln),
              $self->filter_param($filter),
              $self->filter_thresh_param($filter_threshold),
              {slrsim_label => "${tree}_${aln}_${filter}"}
              );
          }

          $params = $self->replace(
            $params,
            $self->reps_param(80),
            $self->indel_rate_param(0.14),
            $self->mean_path_param(1.5),
            {experiment_name => 'fig_two'}
            );
          push @sets, $params;

          $self->store_meta($params) if (scalar(@sets) == 1);
        }
      }
    }
  }
  return \@sets;
}

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

  return { slrsim_tree_file => $tree,
  tree_name => $tree_name};
}

sub omega_param {
  my $self       = shift;
  my $omega_name = shift;

  my $distr_params = $self->param('omega_distributions')->{$omega_name};
  $self->hash_print($distr_params);
  return $distr_params;
}

sub mean_path_param {
  my $self = shift;
  my $path = shift;

  return {slrsim_tree_mean_path => $path};
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

  if ($filter eq 'gblocks') {
    $f->{maximum_mask_fraction} = 1;
  } else {
    $f->{maximum_mask_fraction} = 0.6;
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
