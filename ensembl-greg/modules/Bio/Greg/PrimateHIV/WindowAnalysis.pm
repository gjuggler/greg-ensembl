package Bio::Greg::PrimateHIV::WindowAnalysis;

use strict;
use Bio::Greg::Codeml;
use Bio::Greg::Hive::PhyloAnalysis;
use File::Path;

use base ( 'Bio::Greg::Hive::Process', 'Bio::Greg::Hive::PhyloAnalysis', 'Bio::Greg::StatsCollectionUtils',
  'Bio::Greg::Hive::Align');

my $TREE = 'Bio::EnsEMBL::Compara::TreeUtils';

sub param_defaults {
  return {
    aln_type                    => 'genomic_mammals',
    output_table => 'stats_windows'
  };
}

sub fetch_input {
  my ($self) = @_;

  # Fetch parameters from all possible locations.
  $self->load_all_params();

  # Create tables if necessary.
  $self->create_table_from_params( $self->compara_dba, $self->param('output_table'),
                                   $self->get_table_structure );  
}

sub data_label {
  my $self = shift;
  
  return $self->param('parameter_set_shortname');
}

sub run {
  my $self = shift;

  my $params = $self->params;
  my $tree   = $self->get_tree;

  $params->{'fail_on_altered_tree'} = 0;

  $self->param( 'reference_species', 9606 );
  my $ref_species = $self->param('reference_species');
  my @members     = $tree->leaves;
  my ($ref_member) = grep { $_->taxon_id == $self->param('reference_species') } @members;
  if ($self->param('gene_id')) {
    my $gn = $self->param('gene_id');
    my ($gene_name_member) = grep { $_->get_Gene->external_name eq $gn || $_->stable_id eq $gn || $_->gene_member->stable_id eq $gn || $_->get_Transcript->stable_id eq $gn} @members;
    $ref_member = $gene_name_member if (defined $gene_name_member);
    print "Reference member found by gene ID [$gn]!\n" if ($self->debug);
  }
  $self->param('ref_member_id',$ref_member->stable_id);

  dir("No ref member!") unless (defined $ref_member);

  my $gene_name = $ref_member->get_Gene->external_name || $ref_member->gene_member->stable_id;
  $self->param('gene_name',$gene_name);

  my $c_dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
    -url => 'mysql://ensadmin:ensembl@ensdb-archive:5304/ensembl_compara_58' );
  my $params = $self->params;
  my $tree_aln_obj =
    Bio::EnsEMBL::Compara::ComparaUtils->get_compara_or_genomic_aln( $c_dba, $tree, $ref_member,
    $self->params );
  my $aln = $tree_aln_obj->{aln};
  $tree = $tree_aln_obj->{tree};
  my $extra  = $tree_aln_obj->{extra};
  $self->set_params($extra);

  Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $aln, { width => 150, full => 1 } ) if ($self->debug);

  if (!defined $tree) {
    $self->fail_and_die("Tree undefined: [$tree]\n");
  }
  
  if (scalar($tree->leaves) < 2) {
    $self->fail_and_die("Tree too small!".' '.$self->param('aln_type').' '.$tree->newick_format);
  }

  if ($aln->length < 50) {
    $self->fail_and_die("Alignment too short!".' '.$self->param('aln_type').' '.$tree->newick_format.' '.$aln->length);
  }

  $aln = Bio::EnsEMBL::Compara::AlignUtils->sort_by_tree($aln,$tree);

  Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $aln, { width => 150, full => 1 } ) if ($self->debug);
  my $pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($aln);
  Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $pep_aln, { width => 150, full => 1 } ) if ($self->debug);

  my $aln_type = $self->param('aln_type');
  my $p_set_id = $self->param('parameter_set_id');
  my $subfolder = "primate_hiv_alns/$p_set_id/";

  # Output the tree and aln.
  my $aln_filename = $gene_name . "_" . $self->data_label;
  my $aln_file_obj =
    $self->save_aln( $aln,
    { hash_subfolders => 1, subfolder => $subfolder, filename => $aln_filename } );
  $self->param( 'aln_file', $aln_file_obj->{rel_file} );

  my $tree_filename = $gene_name . "_" . $self->data_label;
  my $tree_file_obj =
    $self->save_file(
    { extension => 'nh', hash_subfolders => 1, subfolder => $subfolder, filename => $tree_filename } );
  Bio::EnsEMBL::Compara::TreeUtils->to_file($tree,$tree_file_obj->{full_file});
  $self->param( 'tree_file', $tree_file_obj->{rel_file} );

  # Get a file to save SLR results in.
  my $slr_filename = $gene_name . "_" . $self->data_label;
  my $slr_file_obj =
    $self->save_file({ extension => 'out', hash_subfolders => 1, subfolder => $subfolder, filename => $tree_filename } );
  $self->param( 'slr_file', $slr_file_obj->{rel_file} );

  # Run SLR.
  my $slr_params = Bio::Greg::Hive::PhyloAnalysis->param_defaults;
  my $all_params = $self->replace($slr_params,$self->params);
  $all_params->{output_to_file} = $slr_file_obj->{full_file};

  # If the SLR output already exists, don't run again.
  my $results;
  my $slr_full_file = $slr_file_obj->{full_file};
  if (-e $slr_full_file) {
    open(IN,"$slr_full_file");
    my @output = <IN>;
    close(IN);
    print join("",@output)."\n" if ($self->debug);
    $results = $self->parse_slr_output(\@output,$all_params);
  } else {
    $results = $self->run_sitewise_dNdS($tree,$aln,$all_params);
  }

  $self->param('slr_dnds',$results->{omega});
  $self->param('slr_kappa',$results->{kappa});

  my $hash = $self->results_to_psc_hash($results,$pep_aln);

  # Store windowed p-values.
  foreach my $size (10, 30, 50, 100, 9999) {
    $self->run_with_windows($size,$size/2,$aln,$tree,$ref_member,$hash);
  }

}

sub write_output {
  my $self = shift;
  
}

sub run_with_windows {
  my $self = shift;
  my $w_size = shift;
  my $w_step = shift;
  my $aln = shift;
  my $tree = shift;
  my $ref_member = shift;
  my $psc_hash = shift;

  my $params = $self->params;
  my $cur_params = $self->replace( $params, {} );

  print "REF: ".$ref_member->stable_id."\n" if ($self->debug);

  my @seqs = $aln->each_seq;
  my ($ref_seq) = grep { $_->id eq 'ens_'. $ref_member->taxon_id } @seqs;
  my $ref_tx = $ref_member->get_Transcript;

  my $len = $ref_member->seq_length;

  my $no_windows_yet = 1;
  for (my $i=1; $i < $len; $i += $w_step) {
    my $lo = $i;
    my $hi = $i + $w_size;
    $hi = $len if ($hi > $len);

    last if (!$no_windows_yet && $hi - $lo < $w_size - $w_step);
    $no_windows_yet = 0;
    printf ">>>> PEPTIDE WINDOW: %d %d\n",$lo,$hi if ($self->debug);

    my $lo_coords = $self->get_coords_from_pep_position($ref_member,$lo);
    my $hi_coords = $self->get_coords_from_pep_position($ref_member,$hi);
    my $aln_coord_lo = $aln->column_from_residue_number($ref_seq->id,$lo);
    my $aln_coord_hi = $aln->column_from_residue_number($ref_seq->id,$hi);

    my $cur_params = $self->replace($params,{
      stable_id_peptide => $ref_member->stable_id,
      stable_id_transcript => $ref_tx->stable_id,
      stable_id_gene => $ref_member->get_Gene->stable_id,
      
      peptide_window_start => $lo,
      peptide_window_end => $hi,
      peptide_window_width => $w_size,
      aln_window_start => $aln_coord_lo,
      aln_window_end => $aln_coord_hi,
      hg19_chr_name => $lo_coords->{hg19_name},
      hg18_chr_name => $lo_coords->{hg18_name},
      hg19_window_start => $lo_coords->{hg19_pos},
      hg18_window_start => $lo_coords->{hg18_pos},
      hg19_window_end => $hi_coords->{hg19_pos},
      hg18_window_end => $hi_coords->{hg18_pos},
      gene_name => $ref_member->get_Gene->external_name,
                                    });

    foreach my $sites ($psc_hash) {
      my @window_sites = map {$sites->{$_}} keys %$sites;
      @window_sites = grep {
        my $pos = $_->{aln_position};
        ($pos >= $aln_coord_lo && $pos < $aln_coord_hi);
      } @window_sites;
      my $pval = $self->combined_pval(\@window_sites,'fisher');
      print "window[$lo-$hi] pval[$pval]\n" if ($self->debug);
      my @pos_sites = grep {$_->{omega} > 1} @window_sites;
      my $added_params = {
        'n_leaves' => scalar($tree->leaves),
        'pval' => $pval,
        'n_sites' => scalar(@window_sites),
        'n_pos_sites' => scalar(@pos_sites),
      };
      $cur_params = $self->replace($cur_params,$added_params);
    }

    $cur_params->{data_id} = $cur_params->{node_id};

    $self->store_params_in_table($self->dbc,$self->param('output_table'),$cur_params);
  }

  $self->fail_and_die("No windows covered!") if ($no_windows_yet);

}

sub get_table_structure {
  my $self = shift;

  my $structure = {
    job_id => 'int',
    data_id   => 'int',
    parameter_set_id => 'int',

    gene_name => 'string',
    aln_type => 'char16',
    parameter_set_name => 'string',
    parameter_set_shortname => 'char16',

    stable_id_gene => 'string',
    stable_id_transcript => 'string',
    stable_id_peptide => 'string',

    peptide_window_start => 'int',
    peptide_window_end => 'int',
    peptide_window_width => 'int',

    aln_window_start => 'int',
    aln_window_end => 'int',

    'hg19_chr_name'     => 'string',
    'hg19_window_start' => 'int',
    'hg19_window_end'   => 'int',

    'hg18_chr_name'     => 'string',
    'hg18_window_start' => 'int',
    'hg18_window_end'   => 'int',

    filtered_Hsap => 'int',
    filtered_Ptro => 'int',
    filtered_Ggor => 'int',
    filtered_Ppyg => 'int',
    filtered_Mmul => 'int',

    pval => 'float',
    n_leaves => 'int',
    n_sites    => 'int',
    n_pos_sites    => 'int',

    slr_file => 'string',
    tree_file => 'string',
    aln_file => 'string',

    slr_dnds => 'float',
    slr_kappa => 'float',

    unique_keys => 'data_id,parameter_set_id,peptide_window_start,peptide_window_width',
  };

  return $structure;
}


1;
