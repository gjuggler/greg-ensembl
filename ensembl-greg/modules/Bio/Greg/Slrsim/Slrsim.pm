package Bio::Greg::Slrsim::Slrsim;

use strict;
use Bio::EnsEMBL::Hive::Process;
use Bio::Greg::Hive::Process;
use File::Path;
use FreezeThaw qw(freeze thaw cmpStr safeFreeze cmpStrHard);

use base (
  'Bio::Greg::Hive::Process',
  'Bio::Greg::Hive::PhyloSim',
  'Bio::Greg::Hive::Align',
  'Bio::Greg::Hive::AlignmentScores',
  'Bio::Greg::Hive::PhyloAnalysis',
  );

sub param_defaults {
  my $self = shift;

  my $phylosim = Bio::Greg::Hive::PhyloSim::param_defaults;
  my $align = Bio::Greg::Hive::Align::param_defaults;
  my $alignment_scores = Bio::Greg::Hive::AlignmentScores::param_defaults;
  my $phylo_analysis = Bio::Greg::Hive::PhyloAnalysis::param_defaults;


  return $self->replace($phylosim, $align, $alignment_scores, $phylo_analysis, {
    genes_table                        => 'genes',
    sites_table                        => 'sites',
                        });
}

sub fetch_input {
  my $self = shift;

  $self->load_all_params();

  $self->create_table_from_params( $self->compara_dba, $self->param('sites_table'),
                                   $self->_sites_table_structure );
  
  $self->create_table_from_params( $self->compara_dba, $self->param('genes_table'),
                                   $self->_genes_table_structure );
  
}

sub run {
  my $self = shift;

  my $tree = $self->get_tree;

  $self->param('phylosim_seq_length', '50');
  $self->param('analysis_action', 'paml');

  # First, simulate the true alignment.
  my $sim_obj = $self->_simulate_alignment($tree);
  my $true_aln = $sim_obj->{aln};
  my $true_sitewise_hash = $sim_obj->{omegas};
  my $true_pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($true_aln);
  
  # Now, align the sequences.
  my $inferred_aln = $self->_align($tree, $true_aln, $true_pep_aln);
  my $inferred_pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($inferred_aln);

  # Mask the alignment if needed.
  my $masked_aln = $self->_mask_alignment($tree, $inferred_aln);
  my $masked_pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($masked_aln);

  # Run the sitewise analysis.
  my $sitewise_hash = $self->_run_sitewise($tree, $masked_aln, $masked_pep_aln);

  # Collect and store results.
  $self->_collect_and_store_results($tree, 
                                    $masked_aln, $masked_pep_aln, $sitewise_hash, 
                                    $true_aln, $true_pep_aln, $true_sitewise_hash);

}

sub _simulate_alignment {
  my $self = shift;
  my $tree = shift;

  my $out_f = $self->_save_file('sim', 'perlobj');
  my $out_file = $out_f->{full_file};
  
  if (!-e $out_file) {
    print "  running simulation\n";
    my $obj = $self->simulate_alignment($tree);
    my $str = freeze($obj);
    open(OUT,">$out_file");
    print OUT freeze($obj);
    close(OUT);
  } else {
    print "  loading simulated aln from file\n";
  }

  open(IN, $out_file);
  my @lines = <IN>;
  close(IN);

  my ($obj) = thaw(join('',@lines));
  return $obj;
}

sub _align {
  my $self = shift;
  my $tree = shift;
  my $true_aln = shift;
  my $true_pep_aln = shift;

  my $out_f = $self->_save_file('inferred_aln', 'fasta');
  my $out_file = $out_f->{full_file};
  
  if (!-e $out_file) {
    print "  running alignment\n";
    # Calls Bio::Greg::Hive::Align->align
    my $inferred_aln = $self->align($tree, $true_aln, $true_pep_aln);
    Bio::EnsEMBL::Compara::AlignUtils->to_file($inferred_aln, $out_file);
  } else {
    print "  loading inferred aln from file\n";
  }

  my $inferred_aln = Bio::EnsEMBL::Compara::AlignUtils->from_file($out_file);
  return $inferred_aln;
}

sub _mask_alignment {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pep_aln = shift;

  my $out_f = $self->_save_file('masked_aln', 'fasta');
  my $out_file = $out_f->{full_file};
  
  if (!-e $out_file) {
    print "  masking alignment\n";
    # Calls Bio::Greg::Hive::AlignmentScores->mask_alignment
    my $masked_aln = $self->mask_alignment($tree, $aln, $pep_aln);
    Bio::EnsEMBL::Compara::AlignUtils->to_file($masked_aln, $out_file);
  } else {
    print "  loading masked aln from file\n";
  }

  my $masked_aln = Bio::EnsEMBL::Compara::AlignUtils->from_file($out_file);
  return $masked_aln;
}


sub _run_sitewise {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pep_aln = shift;

  my $out_f = $self->_save_file('sitewise', 'out');
  my $out_file = $out_f->{full_file};

  if (!-e $out_file) {
    print("  running sitewise analysis\n");
    my $output_lines = $self->run_sitewise_analysis($tree,$aln, $pep_aln);
    open(OUT, ">$out_file");
    print OUT join("", @{$output_lines});
    close(OUT);
  } else {
    print("  parsing sitewise results\n");
  }

  my $sitewise_results = $self->parse_sitewise_file($tree, $aln, $pep_aln, $out_file);
  return $sitewise_results;
}

sub _collect_and_store_results {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pep_aln = shift;
  my $sitewise_hash = shift;
  my $true_aln = shift;
  my $true_pep_aln = shift;
  my $true_sitewise_hash = shift;

  # Get the sequence to act as a reference in site-wise value comparisons.
  my $reference_id = $self->param('slrsim_ref') || '';
  my @seqs = $pep_aln->each_seq;
  my ($ref_seq) = grep { $_->id eq $reference_id } @seqs;
  $ref_seq = $seqs[0] unless ( defined $ref_seq );
  my $seq_str = $ref_seq->seq;
  $seq_str =~ s/-//g;
  my $id = $ref_seq->id;
  
  foreach my $seq_position ( 1 .. length($seq_str) ) {
    my $true_column = $true_pep_aln->column_from_residue_number( $id, $seq_position );
    my $aln_column = $pep_aln->column_from_residue_number( $id, $seq_position );
    my $true_dnds = $true_sitewise_hash->{$true_column}->{omega};
    my $aln_dnds = $sitewise_hash->{$aln_column}->{omega};
    
    if ( !( defined $aln_dnds && defined $true_dnds ) ) {
      if ( !defined $true_dnds ) {
        my $str = sprintf("No true dnds! aln:%s true:%s\n", $aln_column, $true_column);
        die($str);
      } elsif ( !defined $aln_dnds ) {
        sprintf("  no aln dnds! aln:%s true:%s\n", $aln_column, $true_column);
      }
    }

    my $obj;
    $obj->{true_type}      = $true_sitewise_hash->{$true_column}->{'type'}      || '';
    $obj->{true_ncod}      = $true_sitewise_hash->{$true_column}->{'ncod'}      || '';
    $obj->{true_dnds}      = $true_sitewise_hash->{$true_column}->{'omega'} || '';

    $obj->{aln_type}       = $sitewise_hash->{$aln_column}->{'type'}        || '';
    $obj->{aln_ncod}       = $sitewise_hash->{$aln_column}->{'ncod'}        || '';
    $obj->{aln_dnds}       = $sitewise_hash->{$aln_column}->{'omega'} || '';
    $obj->{aln_note}       = $sitewise_hash->{$aln_column}->{'note'}        || '';
    $obj->{aln_lrt}        = $sitewise_hash->{$aln_column}->{'lrt_stat'}    || '';
    $obj->{aln_dnds_lower} = $sitewise_hash->{$aln_column}->{'omega_lower'} || '';
    $obj->{aln_dnds_upper} = $sitewise_hash->{$aln_column}->{'omega_upper'} || '';

    $obj->{aln_position} = $aln_column;
    $obj->{seq_position} = $seq_position;
    
    my $params = $self->params;
    $params = $self->replace($params, $obj);

    $self->store_params_in_table( $self->dbc, $self->param('sites_table'), $params );
  }

  $self->store_params_in_table( $self->dbc, $self->param('genes_table'), $self->params);

  
}

sub _save_file {
  my $self = shift;
  my $filename_base = shift;
  my $ext = shift;

  my $rep = $self->param('slrsim_rep');
  my $id = $self->param('slrsim_label');

  my $filename = "${id}_${filename_base}_${rep}";

  my $file_params = {
    id => $id,
    filename => $filename,
    extension => $ext,
    subfolder => 'data',
  };

  my $file_obj = $self->save_file($file_params);
  return $file_obj;
}

sub _genes_table_structure {
  my $structure = {
    experiment_name => 'char64',

    slrsim_label                        => 'char64',
    alignment_score_threshold           => 'float',
    alignment_score_mask_character_cdna => 'char8',
    filtering_name                      => 'string',
    alignment_name                      => 'string',
    alignment_type                      => 'string',
    sitewise_action                     => 'string',
    sitewise_filter_order               => 'string',

    slrsim_rep         => 'int',
    slrsim_tree_file   => 'string',
    slrsim_tree_length => 'float',
    slrsim_ref         => 'string',

    phylosim_seq_length         => 'int',
    phylosim_omega_distribution => 'string',
    phylosim_meanlog            => 'float',
    phylosim_sdlog              => 'float',

    phylosim_insertrate  => 'float',
    phylosim_deleterate  => 'float',
    phylosim_insertmodel => 'string',
    phylosim_deletemodel => 'string',

    tree_length_slr  => 'float',
    tree_max_branch  => 'float',
    tree_mean_branch => 'float',

    sum_of_pairs_score => 'float',
    total_column_score => 'float',

    column_entropy_mean_true => 'float',
    column_entropy_mean_aln  => 'float',

    site_count               => 'float',
    unfiltered_site_count    => 'float',
    unfiltered_site_fraction => 'float',
    unique_keys              => 'data_id,experiment_name,slrsim_label'
  };
  return $structure;
}

sub _sites_table_structure {
  my $structure = {
    data_id          => 'int',
    parameter_set_id => 'int',
    node_id          => 'int',

    # Site-wise stuff.
    aln_position => 'int',
    seq_position => 'int',

    true_dnds    => 'float',
    true_type    => 'char16',

    aln_dnds       => 'float',
    aln_dnds_lower => 'float',
    aln_dnds_upper => 'float',
    aln_type       => 'char16',
    aln_note       => 'char16',
    aln_lrt => 'float',

    unique_keys => 'data_id,parameter_set_id,aln_position'
  };
  return $structure;
}


1;
