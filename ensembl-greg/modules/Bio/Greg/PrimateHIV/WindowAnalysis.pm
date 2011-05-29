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
    output_table => 'stats_windows',
    window_sizes => '10,30,50,100,9999'
  };
}

sub fetch_input {
  my ($self) = @_;

  # Fetch parameters from all possible locations.
  $self->load_all_params();

  # Create tables if necessary.
  $self->create_table_from_params( $self->compara_dba, $self->param('output_table'),
                                   $self->get_table_structure );  

  Bio::EnsEMBL::Compara::ComparaUtils->load_registry();
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
    $ref_member = $gene_name_member if (defined $gene_name_member && $gene_name_member->taxon_id == $self->param('reference_species'));
    print "Reference member [".$ref_member->stable_id."] found by gene ID [$gn]!\n" if ($self->debug);
  }
  if (!defined $ref_member) {
    warn("Good reference member not found -- just using the first one");
    $ref_member = $members[0];
  }

  die("No ref member!") unless (defined $ref_member);
  $self->param('ref_member_id',$ref_member->stable_id);

  my $gene_name = $ref_member->get_Gene->external_name || $ref_member->gene_member->stable_id;
  $self->param('gene_name',$gene_name);
  $self->param('gene_id',$gene_name);

  my $c_dba = $self->compara_dba;
  my $params = $self->params;

  my ($tree, $aln) = $self->_get_aln($tree, $ref_member);

  my $out_f = $self->_save_file('tree', 'nh');
  my $out_full = $out_f->{full_file};
  Bio::EnsEMBL::Compara::TreeUtils->to_file($tree, $out_full);

  $self->param('tree_file', $out_f->{rel_file});

  my $p_set_id = $self->param('parameter_set_id');
  my $p_short = $self->data_label;

  my $sitewise_results = $self->_run_slr($tree, $aln);
  $self->param('slr_dnds',$sitewise_results->{omega});
  $self->param('slr_kappa',$sitewise_results->{kappa});

  my $pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($aln);


  # Extract only the site result objects from the sitewise hash.
  my $sites_hash;
  foreach my $i (1 .. $pep_aln->length) {
    $sites_hash->{$i} = $sitewise_results->{$i};
  }


  # Store windowed p-values.
  my $window_size_string = $self->param('window_sizes');
  my @window_sizes = split(/, ?/, $window_size_string);
  foreach my $size (@window_sizes) {
    $self->run_with_windows($size, $size/2, $aln, $pep_aln, $tree, $ref_member, $sites_hash);
  }
}

sub _run_slr {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  # Get a file to save SLR results in.
  my $out_f = $self->_save_file('slr', 'out');
  my $slr_full_file = $out_f->{full_file};

  $self->param('slr_file', $out_f->{rel_file});

  my $pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($aln);

  $self->param('analysis_action', 'slr');

  if (!-e $slr_full_file) {  
    print "  running SLR\n";
    my $output_lines = $self->run_sitewise_analysis($tree, $aln, $pep_aln);
    open(OUT, ">$slr_full_file");
    print OUT join("", @{$output_lines});
    close(OUT);
  } else {
    print "  parsing SLR\n";
  }

  my $results = $self->parse_sitewise_file($tree, $aln, $pep_aln, $slr_full_file);
  return $results;
}

sub _save_file {
  my $self = shift;
  my $filename_base = shift;
  my $ext = shift;

  my $id = $self->param('gene_id');
  my $pset = $self->param('parameter_set_shortname');

  my $filename = "${id}_${pset}_$filename_base";

  my $file_params = {
    id => "${id}_${pset}",
    filename => $filename,
    extension => $ext,
    subfolder => 'data',
  };

  my $file_obj = $self->save_file($file_params);
  return $file_obj;
}


sub _get_aln {
  my $self = shift;
  my $tree = shift;
  my $ref_member = shift;

  my $out_f = $self->_save_file('aln', 'fasta');
  my $out_file = $out_f->{full_file};
  $self->param('aln_file', $out_f->{rel_file});

  my $pep_f = $self->_save_file('pep_aln', 'fasta');
  my $pep_file = $pep_f->{full_file};
  $self->param('pep_aln_file', $pep_f->{rel_file});

  my $c_dba = $self->compara_dba;

  if (!-e $out_file) {
    print "  fetching aln...\n";
    my $tree_aln_obj =
      Bio::EnsEMBL::Compara::ComparaUtils->get_compara_or_genomic_aln( $c_dba, $tree, $ref_member,
                                                                       $self->params );
    if ($tree_aln_obj == -1) {
      $self->fail_and_die("Error getting tree or aln!");
    }

    my $aln = $tree_aln_obj->{aln};
    $tree = $tree_aln_obj->{tree};
    my $extra  = $tree_aln_obj->{extra};
    $self->set_params($extra);
    
    Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $aln, { width => 150, full => 1 } ) if ($self->debug);
    if (!defined $tree) {
      $self->fail_and_die("Tree undefined: [$tree]");
    }
    if (scalar($tree->leaves) < 2) {
      $self->fail_and_die("Tree too small!".' '.$self->param('aln_type').' '.$tree->newick_format);
    }
    if ($aln->length < 50) {
      $self->fail_and_die("Alignment too short!".' '.$self->param('aln_type').' '.$tree->newick_format.' '.$aln->length);
    }
    $aln = Bio::EnsEMBL::Compara::AlignUtils->sort_by_tree($aln,$tree);
    my $pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($aln);

    Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $aln, { width => 150, full => 1 } ) if ($self->debug);

    Bio::EnsEMBL::Compara::AlignUtils->to_file($aln, $out_file);
    Bio::EnsEMBL::Compara::AlignUtils->to_file($pep_aln, $pep_file);
  } else {
    print "  loading aln from file $out_file\n";
  }
  
  my $aln = Bio::EnsEMBL::Compara::AlignUtils->from_file($out_file);
  $tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_aln( $tree, $aln );
  return ($tree, $aln);
}

sub write_output {
  my $self = shift;  
}

sub run_with_windows {
  my $self = shift;
  my $w_size = shift;
  my $w_step = shift;
  my $aln = shift;
  my $pep_aln = shift;
  my $tree = shift;
  my $ref_member = shift;
  my $psc_hash = shift;

  my $params = $self->params;
  my $cur_params = $self->replace( $params, {} );

  my $ref_seq;
  my $ref_tx;
  my $len;

  if ($tree->isa('Bio::Tree::TreeI')) {
    $tree = Bio::EnsEMBL::Compara::TreeUtils->from_treeI($tree);
  }

  if ($ref_member->isa('Bio::EnsEMBL::Compara::Member')) {
    print "REF MEMBER: ".$ref_member->stable_id."\n" if ($self->debug);
    my @seqs = $pep_aln->each_seq;
    my $taxon_id = $ref_member->taxon_id;
    ($ref_seq) = grep { $_->id =~ m/$taxon_id/i } @seqs;

    if (!defined $ref_seq) {
      # Remember, we connect the genomic aln and member by genome_db name now.
      my $name = $ref_member->genome_db->name;
      ($ref_seq) = grep { $_->id =~ m/$name/i } @seqs;      
    }
    if (!defined $ref_seq) {
      # Remember, we connect the genomic aln and member by genome_db name now.
      my $name = $ref_member->name;
      ($ref_seq) = grep { $_->id =~ m/$name/i } @seqs;      
    }
    $ref_tx = $ref_member->get_Transcript;
  } else {
    $ref_seq = $ref_member;
    $ref_member = undef;
  }

  my $seq_str = $ref_seq->seq;
  my $seq_str_nogaps = $seq_str;
  $seq_str_nogaps =~ s/-//g;
  $len = length($seq_str_nogaps);

  print "$seq_str_nogaps\n";

  my @rows;

  my @all_sites;
  foreach my $i ( 1 .. $pep_aln->length) {
    push @all_sites, $psc_hash->{$i} if (defined $psc_hash->{$i});
  }
  print "Sites with values:". scalar(@all_sites)."\n";

  my $no_windows_yet = 1;
  for (my $i=1; $i < $len; $i += $w_step) {
    my $lo = $i;
    my $hi = $i + $w_size;
    $hi = $len if ($hi > $len);

    last if (!$no_windows_yet && $hi - $lo < $w_size - $w_step);
    $no_windows_yet = 0;
    printf ">>>> PEPTIDE WINDOW: %d %d\n",$lo,$hi if ($self->debug);

    my $cur_params = $params;

    if (defined $ref_member) {
      # Re-fetch a fresh reference member from the database.
      my $ref_member_id = $self->param('ref_member_id');
      print "ref_member_id: [$ref_member_id]\n" if ($self->debug);
      my $mba        = $self->compara_dba->get_MemberAdaptor;
      $ref_member = $mba->fetch_by_source_stable_id( undef, $ref_member_id );
      my $lo_coords = $self->get_coords_from_pep_position($ref_member,$lo);
      my $hi_coords = $self->get_coords_from_pep_position($ref_member,$hi);

      $cur_params = $self->replace($cur_params,{
        stable_id_peptide => $ref_member->stable_id,
        stable_id_transcript => $ref_member->get_Transcript->stable_id,
        stable_id_gene => $ref_member->get_Gene->stable_id,
        gene_name => $ref_member->get_Gene->external_name,
                                   }
        );
      if (defined $lo_coords && $lo_coords ne '') {
        $cur_params = $self->replace($cur_params, {
          hg19_chr_name => $lo_coords->{hg19_name},
          hg18_chr_name => $lo_coords->{hg18_name},
          hg19_window_start => $lo_coords->{hg19_pos},
          hg18_window_start => $lo_coords->{hg18_pos},
                                     }
          );
      }
      if (defined $hi_coords && $hi_coords ne '') {
        $cur_params = $self->replace($cur_params, {
          hg19_window_end => $hi_coords->{hg19_pos},
          hg18_window_end => $hi_coords->{hg18_pos},
                                     }
          );
      }
    }
    
    my $aln_coord_lo = $aln->column_from_residue_number($ref_seq->id,$lo);
    my $aln_coord_hi = $aln->column_from_residue_number($ref_seq->id,$hi);
    $cur_params = $self->replace($cur_params, {
      peptide_window_start => $lo,
      peptide_window_end => $hi,
      peptide_window_width => $w_size,
      aln_window_start => $aln_coord_lo,
      aln_window_end => $aln_coord_hi
                   });

    my $window_hash;
    map {
      my $pos = $_->{aln_position};
      if ($pos >= $aln_coord_lo && $pos < $aln_coord_hi) {
        $window_hash->{$pos} = $_;
      }
    } @all_sites;
    my @window_sites = values %$window_hash;

    my $pval = $self->combined_pval($window_hash,'fisher');
    
    # Get the mean dn/ds for the window
    my $omega_total = 0;
    map {$omega_total += $_->{omega} } @window_sites;
    my $mean_dnds = -1;
    if (scalar(@window_sites) > 0) {
      $mean_dnds = sprintf "%.3f", $omega_total / scalar(@window_sites);
    }
    
    print "window[$lo-$hi] pval[$pval]\n" if ($self->debug);
    my @pos_sites = grep {$_->{omega} > 1} @window_sites;
    my $added_params = {
      'n_leaves' => scalar($tree->leaves),
      'pval' => $pval,
      'n_sites' => scalar(@window_sites),
      'n_pos_sites' => scalar(@pos_sites),
      'mean_dnds' => $mean_dnds
    };
    $cur_params = $self->replace($cur_params,$added_params);
    
    $cur_params->{data_id} = $cur_params->{node_id};
    
    if ($self->within_hive) {
      $self->store_params_in_table($self->dbc,$self->param('output_table'),$cur_params);
    } else {
      $self->hash_print($cur_params);
      push @rows, $cur_params;
    }
  }
  
  if (scalar(@rows) > 0) {
    return @rows;
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
    pep_aln_file => 'string',

    slr_dnds => 'float',
    slr_kappa => 'float',

    unique_keys => 'data_id,parameter_set_id,peptide_window_start,peptide_window_width',
  };

  return $structure;
}


1;
