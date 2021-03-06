package Bio::Greg::Slrsim::CollectSlrsimStats;

use strict;
use Time::HiRes qw(sleep);

use Cwd;
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::Process;

use Bio::Greg::StatsCollectionUtils;

use base ('Bio::Greg::Hive::CollectSitewiseStats');

sub get_gene_table_structure {
  my $self = shift;

  my $added_structure = {
    experiment_name => 'char64',

    slrsim_label                        => 'string',
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

  # Add our structure on top of the base structure defined in CollectSitewiseStats.
  my $structure = $self->SUPER::get_gene_table_structure;
  $structure = $self->replace_params( $structure, $added_structure );

  delete $structure->{tree_newick};

  return $structure;
}

sub get_sites_table_structure {
  my $self = shift;

  my $added_structure = {
    data_id          => 'int',
    parameter_set_id => 'int',
    node_id          => 'int',

    # Site-wise stuff.
    aln_position => 'int',
    seq_position => 'int',

    true_dnds    => 'float',
    true_type    => 'char16',
    true_ncod    => 'int',

    aln_dnds       => 'float',
    aln_dnds_lower => 'float',
    aln_dnds_upper => 'float',
    aln_type       => 'char16',
    aln_note       => 'char16',
    aln_ncod       => 'int',
    aln_lrt => 'float',
    nongap_bl  => 'float',

    unique_keys => 'data_id,parameter_set_id,aln_position'
  };

# We're not using any of the base structure defined in CollectSitewiseStats, so just return this hash.
  return $added_structure;
}

sub data_for_gene {
  my $self = shift;

  # This is used by the StatsCollectionUtils methods.
  $self->param( 'get_all_sites', 1 );

  my $data = $self->SUPER::data_for_gene();
  return undef unless ( defined $data );
  print "Original data:\n";
  Bio::EnsEMBL::Compara::ComparaUtils->hash_print($data);

  my $sa_true;
  my $sa_aln;
  my $cdna_true;
  my $cdna_aln;
  my @true_entropies;
  my @aln_entropies;
  my $sum_of_pairs_score;
  my $total_column_score;
  my $tree;
  eval {
    my $cur_params = $self->params;
    my $true_aln_params =
      $self->replace_params( $cur_params,
      { alignment_table => 'protein_tree_member', alignment_score_filtering => 0 } );
    print "Getting TRUE alignment...\n";
    ( $tree, $sa_true, $cdna_true ) =
      Bio::EnsEMBL::Compara::ComparaUtils->tree_aln_cdna( $self->compara_dba, $true_aln_params );

    print "Getting INFERRED alignment...\n";
    ( $tree, $sa_aln, $cdna_aln ) =
      Bio::EnsEMBL::Compara::ComparaUtils->tree_aln_cdna( $self->compara_dba, $cur_params );

    Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $cdna_true, { length => 180 } );
    Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $cdna_aln,  { length => 180 } );
  };
  die( "Hold up: " . $@ ) if ($@);
  return if ( !$sa_true || !$sa_aln );

  print "Calculating stuff...\n";

#$data->{sum_of_pairs_score} = Bio::EnsEMBL::Compara::AlignUtils->sum_of_pairs_score( $sa_true, $sa_aln );
#$data->{total_column_score} = Bio::EnsEMBL::Compara::AlignUtils->total_column_score( $sa_true, $sa_aln );
#$data->{column_entropy_mean_true} =
#  Bio::EnsEMBL::Compara::AlignUtils->average_column_entropy($sa_true);
#$data->{column_entropy_mean_aln} =
#  Bio::EnsEMBL::Compara::AlignUtils->average_column_entropy($sa_aln);

  $data->{site_count}            = $self->site_count($sa_aln);
  $data->{unfiltered_site_count} = $self->unfiltered_site_count($sa_aln);
  $data->{unfiltered_site_fraction} =
    $self->unfiltered_site_count($sa_aln) / $self->site_count($sa_aln);

  # Get the SLR-inferred tree.
  my $newick   = $self->param('slr_tree');
  my $slr_tree = Bio::EnsEMBL::Compara::TreeUtils->from_newick($newick);
  $data->{'tree_newick_slr'} = $newick;
  $data->{'tree_length_slr'} = $self->tree_length($slr_tree);

  $data->{'parameter_set_id'} = 0;
  return $data;
}

sub data_for_site {
  my $self         = shift;
  my $ref_name     = shift;
  my $seq_position = shift;

  my $data = $self->params;

  my $cdna_true = $self->param('cdna_true');
  my $sa_true   = $self->param('sa_true');
  my $cdna_aln  = $self->param('cdna_aln');
  my $sa_aln    = $self->param('sa_aln');

  my $tree = $self->param('temp_tree');
  if ( !defined $tree ) {
    my $tree = $self->get_tree;
    $self->param( 'temp_tree', $tree );
  }

  if ( !defined $cdna_true ) {
    my $true_aln_params =
      $self->replace_params( $self->params,
      { alignment_table => 'protein_tree_member', alignment_score_filtering => 0 } );
    ( $tree, $sa_true, $cdna_true ) =
      Bio::EnsEMBL::Compara::ComparaUtils->tree_aln_cdna( $self->compara_dba, $true_aln_params );
    $self->param( 'cdna_true', $cdna_true );
    $self->param( 'sa_true',   $sa_true );

    ( $tree, $sa_aln, $cdna_aln ) =
      Bio::EnsEMBL::Compara::ComparaUtils->tree_aln_cdna( $self->compara_dba, $self->params );
    $self->param( 'cdna_aln', $cdna_aln );
    $self->param( 'sa_aln',   $sa_aln );
  }

  my $true_omegas = $self->param('true_omegas');
  my $aln_omegas  = $self->param('aln_omegas');
  if ( !defined $true_omegas ) {

    # Get all the site-wise data from the omega table.
    my $aln_table_name = $data->{'omega_table'};
    my $sth1           = $self->compara_dba->dbc->prepare(
      "SELECT aln_position,omega,type,note,ncod,lrt_stat FROM sitewise_omega WHERE node_id=?;");
    my $cmd =
      "SELECT aln_position,omega,omega_lower,omega_upper,type,note,ncod,lrt_stat FROM $aln_table_name WHERE data_id=?;";
    print "Getting true omegas...\n";
    my $sth2 = $self->compara_dba->dbc->prepare($cmd);
    $sth1->execute( $self->node_id );
    $sth2->execute( $self->data_id );
    $true_omegas = $sth1->fetchall_hashref('aln_position');
    $aln_omegas  = $sth2->fetchall_hashref('aln_position');

    $self->param( 'true_omegas', $true_omegas );
    $self->param( 'aln_omegas',  $aln_omegas );
  }
  print "node: " . $self->node_id . "\n";

  my @true_entropies;
  my @aln_entropies;

  #my @true_entropies = Bio::EnsEMBL::Compara::AlignUtils->column_entropies($cdna_true);
  #my @aln_entropies  = Bio::EnsEMBL::Compara::AlignUtils->column_entropies($cdna_aln);

  my $obj;
  my $true_col = $sa_true->column_from_residue_number( $ref_name, $seq_position );
  my $aln_col = $sa_aln->column_from_residue_number( $ref_name, $seq_position );

  $obj->{seq_position} = $seq_position;
  $obj->{aln_position} = $aln_col;
  $obj->{true_dnds}    = $true_omegas->{$true_col}->{'omega'};
  $obj->{aln_dnds}     = $aln_omegas->{$aln_col}->{'omega'};

  #printf("t:%.3f a:%.3f\n",$obj->{true_dnds},$obj->{aln_dnds});
  if ( !( defined $obj->{aln_dnds} && defined $obj->{true_dnds} ) ) {
    if ( $data->{'analysis_action'} eq '' || $data->{'analysis_action'} eq 'none' ) {
      # Do nothing.
      $obj->{aln_dnds}  = 0;
      $obj->{true_dnds} = 0;
    } else {
      if ( !defined $obj->{true_dnds} ) {
        my $str = sprintf("No true dnds! aln:%s  %s  true:%s\n", $aln_col, $obj->{aln}, $true_col);
        die($str);
        next;
      } elsif ( !defined $obj->{aln_dnds} ) {
        print "No aln dnds!\n";
        #next;
      }
    }
  }

  $obj->{aln_type}       = $aln_omegas->{$aln_col}->{'type'}        || '';
  $obj->{true_type}      = $true_omegas->{$true_col}->{'type'}      || '';
  $obj->{aln_note}       = $aln_omegas->{$aln_col}->{'note'}        || '';
  $obj->{true_ncod}      = $true_omegas->{$true_col}->{'ncod'}      || '';
  $obj->{aln_ncod}       = $aln_omegas->{$aln_col}->{'ncod'}        || '';
  $obj->{aln_lrt}        = $aln_omegas->{$aln_col}->{'lrt_stat'}    || '';
  $obj->{aln_dnds_lower} = $aln_omegas->{$aln_col}->{'omega_lower'} || '';
  $obj->{aln_dnds_upper} = $aln_omegas->{$aln_col}->{'omega_upper'} || '';

  #  $obj->{true_entropy} = sprintf( "%.3f", $true_entropies[$true_col] || 0 );
  #  $obj->{aln_entropy}  = sprintf( "%.3f", $aln_entropies[$aln_col]   || 0 );

  my $true_slice = $sa_true->slice( $true_col, $true_col );

  #Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $true_slice);
  printf "t -> %.3f %.3f %.3f\n", $obj->{true_dnds}, $obj->{true_ncod}, $obj->{true_entropy};
  my $aln_slice = $sa_aln->slice( $aln_col, $aln_col );

  #Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $aln_slice);
  printf "a -> %.3f %.3f %.3f\n", $obj->{aln_dnds}, $obj->{aln_ncod}, $obj->{aln_entropy};

  $obj->{parameter_set_id} = 0;

  # Store values in our output table.
  $data = $self->replace_params( $data, $obj );

  return $data;
}

# Need a custom site-wise loop to go through each residue in the reference sequence and get the corresponding site data.
sub get_sites_data {
  my $self = shift;

  my $tree = $self->get_tree;
  my $aln  = $tree->get_SimpleAlign;

  # Get the sequence to act as a reference in site-wise value comparisons.
  my $reference_id = $self->param('slrsim_ref') || '';
  my @seqs = $aln->each_seq;
  my ($ref_seq) = grep { $_->id eq $reference_id } @seqs;
  $ref_seq = $seqs[0] unless ( defined $ref_seq );

  my $seq_str = $ref_seq->seq;
  $seq_str =~ s/-//g;
  foreach my $sequence_position ( 1 .. length($seq_str) ) {
    my $ref_aa = substr( $seq_str, $sequence_position - 1, 1 );
    print "Position $sequence_position $ref_aa\n";
    my $site_data = $self->data_for_site( $ref_seq->id, $sequence_position );
    if ( defined $site_data ) {
      $self->store_params_in_table( $self->compara_dba, $self->param('sites_table'), $site_data );
    }
  }
}

1;
