package Bio::Greg::Mammals::SitewiseMammals;

use strict;
use Bio::Greg::Codeml;
use File::Path;

use base (
  'Bio::Greg::Hive::Process', 'Bio::Greg::StatsCollectionUtils',
  'Bio::Greg::Hive::Align',
  'Bio::Greg::Hive::SitewiseMapper',
  'Bio::Greg::Hive::PhyloAnalysis'
);

sub param_defaults {
  return {
    aln_type                            => 'compara',
    genes_table                        => 'genes',
    sites_table                        => 'sites',
    quality_threshold                   => 30,
  };
}

sub fetch_input {
  my $self = shift;

  # Fetch parameters from all possible locations.
  $self->load_all_params();

  # Create tables if necessary.
  $self->create_table_from_params( $self->compara_dba, $self->param('sites_table'),
                                   $self->_sites_table_structure );

  $self->create_table_from_params( $self->compara_dba, $self->param('genes_table'),
                                   $self->_genes_table_structure );
}

sub run {
  my $self = shift;

  my $tree = $self->get_tree;
  my $ref_member = $self->_get_ref_member($tree);

  my $aln_f = $self->_save_file($tree,'_aln','fasta');
  $self->param('existing_alignment_file',$aln_f->{full_file});

  #print $tree->newick_format."\n";
  my $tree_aln_obj =
    Bio::EnsEMBL::Compara::ComparaUtils->get_compara_or_genomic_aln( $self->compara_dba, $tree, $ref_member,
                                                                     $self->params );
  my $aln = $tree_aln_obj->{aln};
  my $pep_aln = $tree_aln_obj->{pep_aln};
  $tree = $tree_aln_obj->{tree};
  my $extra  = $tree_aln_obj->{extra};
  $self->_output_files($tree,$aln,$pep_aln);
  $self->set_params($extra);  
  #$self->pretty_print($aln,{full=>1});

  $self->param('gene_name',$ref_member->display_label);
  $self->param('protein_id',$ref_member->stable_id);

  my $slr_hash;

  # Find reasons to skip sitewise analysis.
  my $skip_sitewise = 0;
  $skip_sitewise = 'Tree size' if (scalar($tree->leaves) <= 4);
  $skip_sitewise = 'Alignment length: '.$pep_aln->length if ($pep_aln->length > 10000);

  if (!$skip_sitewise) {
    #$self->_run_hyphy($tree,$aln);
    $slr_hash = $self->_run_and_store_slr($tree,$aln,$pep_aln);
  } else {
    $self->param('skipped_reason', $skip_sitewise);
    $slr_hash = {};
  }

  $self->_collect_gene_data($tree,$aln,$slr_hash);

}

sub _collect_gene_data {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $slr_hash = shift;

  my $ref_member = $self->_get_ref_member($tree);

  # Basic details.
  my $tx = $ref_member->get_Transcript;
  eval {
    $tx = $tx->transform('chromosome');
  };
  if (defined $tx) {
    $self->param('chr_name','chr'.$tx->slice->seq_region_name);
    $self->param('chr_start',$tx->coding_region_start);
    $self->param('chr_end',$tx->coding_region_end);
    $self->param('chr_strand',$tx->strand);
  }
  $self->param('gene_taxon_id',$ref_member->taxon_id);
  $self->param('gene_name',$ref_member->get_Gene->external_name);
  $self->param('gene_description',$ref_member->get_Gene->description);
  $self->param('gene_id',$ref_member->get_Gene->stable_id);
  $self->param('transcript_id',$ref_member->get_Transcript->stable_id);
  $self->param('protein_id',$ref_member->get_Transcript->translation->stable_id);
  $self->param('seq_length',$ref_member->seq_length);
  $self->param('aln_length',$aln->length / 3);
  
  # Calculate GC content.
  $self->param( 'gc_cds', sprintf( "%.3f", $self->gc_content($ref_member) ) );
  $self->param( 'gc_3', sprintf( "%.3f", $self->gc3_content($ref_member) ) );
  $self->param( 'gc_genomic', sprintf( "%.3f", $self->genomic_gc_content($ref_member) ) );

  # Tree properties.
  $self->param('tree_mean_path',$self->mean_path($tree));
  $self->param('leaf_count',scalar($tree->leaves));
  $self->param('duplication_count',$self->duplication_count($tree));

  $self->param('sitewise_value_count',$self->sitewise_count($slr_hash));
  $self->param('pval_fisher',$self->combined_pval($slr_hash,'fisher'));
  $self->calculate_fractions($slr_hash); # Calculates the f_pos, f_neg, and f_neut

  $self->store_params_in_table($self->dbc,$self->param('genes_table'),$self->params);
}

sub _output_files {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  my $aln_f = $self->_save_file($tree,'_aln','fasta');
  my $tree_f = $self->_save_file($tree,'_tree','nhx');
  Bio::EnsEMBL::Compara::TreeUtils->to_file($tree,$tree_f->{full_file});
  Bio::EnsEMBL::Compara::AlignUtils->to_file( $aln, $aln_f->{full_file} );

  $self->param('aln_file',$aln_f->{rel_file});
  $self->param('tree_file',$tree_f->{rel_file});
}

sub _run_hyphy {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  my $hyphy_params = Bio::Greg::Hive::PhyloAnalysis->param_defaults;
  $hyphy_params = $self->replace($hyphy_params,$self->params);

  $self->run_hyphy($tree,$aln,$hyphy_params);
}

sub _run_and_store_slr {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pep_aln = shift;

  my $slr_f = $self->_save_file($tree,'_slr','out');
  $self->param('slr_file',$slr_f->{rel_file});
  my $output_filename = $slr_f->{full_file};
  
  # Prepare parameters.
  my $slr_params = Bio::Greg::Hive::PhyloAnalysis->param_defaults;
  my $all_params = $self->replace($slr_params,$self->params);
  $all_params->{output_to_file} = $output_filename;

  # Run SLR on the tree & alignment, and returns the results as a hash 
  my $slr_results;
  if (-e $output_filename) {
    warn("Warning: SLR results file already exists -- using that instead of running again!");
    open(IN,"$output_filename");
    my @output = <IN>;
    close(IN);
    $slr_results = $self->parse_slr_output(\@output,$all_params);
  } else {
    $slr_results = $self->run_sitewise_dNdS($tree,$aln,$all_params);
  }

  # Turn results array into hashref data.
  my $slr_hash = $self->results_to_psc_hash($slr_results,$pep_aln);

  # Collect Pfam, exon, and filter data.
  my $pfam_sites = $self->collect_pfam($tree,$pep_aln);
  my $exon_sites = $self->collect_exons($tree,$pep_aln);
  my $filter_sites = $self->do_filter($tree,$pep_aln);

  my $ref_member = $self->_get_ref_member($tree);
  my $ref_seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id($pep_aln,$ref_member->stable_id);  

  # Add ungapped branch lengths to the slr hash.
  $self->add_ungapped_branch_lengths($tree, $pep_aln, $slr_hash);

  foreach my $key (keys %$slr_hash) {
    my $site_obj = $slr_hash->{$key};
    my $aln_position = $site_obj->{aln_position};
    my $cur_params = $self->params;

    $cur_params = $self->replace($cur_params,$site_obj);

    if ($ref_member->taxon_id == 9606) {
      my $ref_position = $ref_seq->location_from_column($aln_position);
      if (defined $ref_position && $ref_position->location_type() eq 'EXACT') {
        my $ref_coords = $self->get_coords_from_pep_position($ref_member,$ref_position->start);
        $cur_params = $self->replace($cur_params,$ref_coords);
      }
    }

    if (defined $pfam_sites->{$aln_position}) {
      $cur_params = $self->replace($cur_params, $pfam_sites->{$aln_position});
    }    
    if (defined $exon_sites->{$aln_position}) {
      $cur_params = $self->replace($cur_params, $exon_sites->{$aln_position});
    }
    if (defined $filter_sites->{$aln_position}) {
      $cur_params = $self->replace($cur_params, $filter_sites->{$aln_position});
    }
    
    $self->store_params_in_table($self->dbc,$self->param('sites_table'),$cur_params);
  }

  return $slr_hash;
}

sub _get_ref_member {
  my $self = shift;
  my $tree = shift;

  my $ref_member;
  my @members = $tree->leaves;

  # Try for human first.
  ($ref_member) = grep { $_->taxon_id == 9606 } @members;
  
  # Then mouse.
  ($ref_member) = grep { $_->taxon_id == 10090 } @members if (!defined $ref_member);

  # Then dog.
  ($ref_member) = grep { $_->taxon_id == 9615 } @members if (!defined $ref_member);

  # Then, take whatever we can get.
  ($ref_member) = $members[0] if (!defined $ref_member);
  
  $self->param('ref_member_id',$ref_member->stable_id);
  return $ref_member;
}

sub _save_file {
  my $self = shift;
  my $tree = shift;
  my $file_suffix = shift;
  my $file_ext = shift;

  my $ref_member = $self->_get_ref_member($tree);
  my $id = $ref_member->stable_id;

  my $shortname = $self->param('parameter_set_shortname');
  my $subfolder = 'data';
  
  my $params = {
    id => $id,
    filename => $id . '_' . $shortname . $file_suffix,
    extension => $file_ext,
    subfolder => $subfolder
  };
  return $self->save_file($params);
}


sub _sites_table_structure {
  my $self = shift;

  my $structure = {
    # IDs.
    node_id => 'int',
    parameter_set_id => 'tinyint',

    # Genomic locations on the reference sequence
    chr_name => 'char8',
    chr_start => 'int',
    chr_end => 'int',

    # Pfam annotations.
    pfam_domain => 'char16',
    pfam_position => 'int',
    pfam_score => 'int',
    
    # Exon annotation.
    exon_position => 'char8',
    splice_distance => 'int',

    # Basic SLR-derived attributes.
    aln_position => 'int',
    ncod => 'int',
    nongap_bl => 'float',           # This is calculated using StatsCollectionUtils' method add_ungapped_branch_lengths
    omega => 'float',
    omega_lower => 'float',
    omega_upper => 'float',
    lrt_stat => 'float',
    type => 'char16',
    note => 'char16',
    random => 'char16',

    unique_keys => 'node_id,parameter_set_id,aln_position'
  };

  return $structure;
}

sub _genes_table_structure {
  my $self = shift;

  my $structure = {
    # IDs.
    node_id => 'int',
    job_id => 'int',
    parameter_set_id => 'tinyint',
    parameter_set_shortname => 'char8',

    # Human gene attributes.
    chr_name => 'string',
    chr_start =>'int',
    chr_end => 'int',
    gene_name => 'char32',
    gene_description => 'string',
    gene_taxon_id => 'int',
    gene_id => 'char16',
    transcript_id => 'char16',
    protein_id => 'char16',

    # Alignment properties.
    aln_length => 'int',
    seq_length => 'int',
    filtered_site_total => 'int',

    # General gene / tree attributes.
    gc_cds => 'float',
    'gc_3' => 'float',
    gc_genomic => 'float',
    tree_mean_path => 'float',
    leaf_count => 'float',
    duplication_count => 'float',

    # Hyphy calculations.
    hyphy_dnds => 'float',
    hyphy_dnds_lo => 'float',
    hyphy_dnds_hi => 'float',

    # SLR calculations.
    slr_dnds => 'float',
    slr_kappa => 'float',

    # Sitewise stats.
    sitewise_value_count => 'float',
    pval_fisher => 'float',
    f_neg => 'float',
    f_pos => 'float',
    f_neut => 'float',

    n_neg => 'int',
    n_pos => 'int',
    n_neut => 'int',

    # Files output.
    slr_file => 'string',
    tree_file => 'string',
    aln_file => 'string',

    unique_keys => 'node_id,parameter_set_id'
  };

  return $structure;
}

1;
