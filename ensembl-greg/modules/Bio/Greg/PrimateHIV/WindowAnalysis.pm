package Bio::Greg::PrimateHIV::WindowAnalysis;

use strict;
use Bio::Greg::Codeml;
use Bio::Greg::Hive::PhyloAnalysis;
use File::Path;
use List::Util;

use Bio::Align::DNAStatistics;

use base (
  'Bio::Greg::Hive::Process', 
  'Bio::Greg::Hive::PhyloAnalysis', 
  'Bio::Greg::StatsCollectionUtils',
  'Bio::Greg::Hive::Align',
  'Bio::Greg::Hive::CountSubstitutions');

my $TREE = 'Bio::EnsEMBL::Compara::TreeUtils';

sub param_defaults {
  return {
    aln_type                    => 'compara',
    output_table => 'stats_windows',
    window_sizes => '9999'
  };
}

sub fetch_input {
  my ($self) = @_;

  # Fetch parameters from all possible locations.
  $self->load_all_params();

  $self->init_table($self->genes_table);

  # Create tables if necessary.
  $self->create_table_from_params( $self->compara_dba, 'windows',
                                   $self->windows_table );

  $self->create_table_from_params( $self->compara_dba, 'vars',
                                   $self->vars_table );

}

sub write_output {
  my $self = shift;

  $self->write_table('genes');
}

sub data_label {
  my $self = shift;
  
  return $self->param('parameter_set_shortname');
}

sub run {
  my $self = shift;

  #$self->param('force_recalc', 1);

  my $gene_id = $self->param('gene_id');
  my $mba = $self->compara_dba->get_MemberAdaptor;
  my $member = $mba->fetch_by_source_stable_id(undef, $gene_id);
  my $gene = $member->get_Gene;
  $self->param('gene_name', $gene->external_name);
  $self->param('data_id', $member->dbID);

  my $pta = $self->compara_dba->get_ProteinTreeAdaptor;
  my $protein_tree = $pta->fetch_by_Member_root_id($member);
  #print $protein_tree->ascii."\n";

  my $ortholog_tree = Bio::EnsEMBL::Compara::ComparaUtils->get_one_to_one_ortholog_tree($self->compara_dba, $member, 'ortholog.*');
  if (!defined $ortholog_tree) {
    $self->fail_and_die("no_orthologs", "Nothing left from one2one_orthologs: [$ortholog_tree]");
  }
  $ortholog_tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_clade($self->compara_dba, $ortholog_tree, 'Primates');
  if (!defined $ortholog_tree || scalar($ortholog_tree->leaves) < 2) {
    $self->fail_and_die("no_orthologs", "Nothing other than human peptide left after restricting to clade: [$ortholog_tree]");
  }

  my $tree = $ortholog_tree;
  my $params = $self->params;

  print $ortholog_tree->ascii."\n";

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

  $ref_member = $member;

  die("No ref member!") unless (defined $ref_member);
  my $ref_member = $ref_member->get_canonical_peptide_Member;
  $self->param('gene_id',$ref_member->gene_member->get_Gene->stable_id);
  $self->param('transcript_id',$ref_member->get_Transcript->stable_id);
  $self->param('protein_id',$ref_member->stable_id);
  $self->param('ref_member', $ref_member);
  my $gene_name = $ref_member->get_Gene->external_name || $ref_member->get_Gene->stable_id;
  $self->param('gene_name',$gene_name);

  $self->param('chr_name', ''.$ref_member->get_Transcript->slice->seq_region_name);
  $self->param('chr_start', $ref_member->get_Transcript->coding_region_start);
  $self->param('chr_end', $ref_member->get_Transcript->coding_region_end);
  $self->param('chr_strand', ''.$ref_member->get_Transcript->slice->strand);

  if ($self->param('chr_name') =~ m/MT/gi) {
    $self->param('bioperl_codontable_id', 2);
    $self->param('icode', 1);
  } else {
    $self->param('bioperl_codontable_id', 1);
    $self->param('icode', 0);
  }

  $self->_collect_vars($ref_member);

  my $c_dba = $self->compara_dba;
  my $params = $self->params;

  $ref_member->name($ref_member->stable_id);
  my ($tree, $aln) = $self->_get_aln($tree, $ref_member);
  $self->_out_aln($aln, 'orig');

  $aln = Bio::EnsEMBL::Compara::AlignUtils->remove_empty_seqs($aln);
  if (scalar($aln->each_seq) < 3) {
    $self->fail_and_die("seq_count", "Not enough sequences for PAML analysis!");
  }
  if ($aln->length > 15000) {
    $self->fail_and_die("aln_length", "Alignment too long - over 15k nucleotides!");
  }

  my ($tree, $aln) = $self->_remove_paralogs($tree, $aln, $ref_member);
  my $genome_tree = $self->_get_genome_tree($tree);
  $tree = $genome_tree;
#  $self->_out_aln($aln, 'paralogs_removed');

  if (scalar($aln->each_seq) < 3) {
    $self->fail_and_die("seq_count", "Not enough sequences for PAML analysis!");
  }

  $aln = $self->_realign($tree, $aln);
  $aln = Bio::EnsEMBL::Compara::AlignUtils->flatten_to_sequence( $aln, $ref_member->name);
#  $self->_out_aln($aln, 'realigned');

  $tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_aln($tree, $aln);

  $aln = $self->_mask_aln($tree, $aln);
#  $self->_out_aln($aln, 'masked');

  if (scalar($aln->each_seq) < 3) {
    $self->fail_and_die("seq_count", "Not enough sequences for PAML analysis!");
  }

  $self->_dnds($aln, 9606, 9598, 'chimp');
  $self->_dnds($aln, 9606, 9600, 'orang');
  $self->_dnds($aln, 9606, 9544, 'rhesus');

  $self->_run_paml($tree, $aln, $ref_member);

  my $out_f = $self->_save_file('tree', 'nh');
  my $out_full = $out_f->{full_file};
  Bio::EnsEMBL::Compara::TreeUtils->to_file($tree, $out_full);

  $self->param('tree_file', $out_f->{rel_file});

  my $p_set_id = $self->param('parameter_set_id');
  my $p_short = $self->data_label;

  # Run SLR
  my $sitewise_results = $self->_run_slr($tree, $aln);
  # Store SLR sitewise
  my $pep_aln = $self->_tx_aln($aln);
  $self->store_sitewise($tree, $pep_aln, $sitewise_results, {omega_table => 'slr_sites'});

  $self->param('slr_dnds',$sitewise_results->{omega});
  $self->param('slr_kappa',$sitewise_results->{kappa});

  $self->param('n_seqs', scalar($pep_aln->each_seq));
  $self->param('aln_length', $pep_aln->length);
  $self->param('seq_length', $ref_member->seq_length);

  # Do the windows stuff.
  my $pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($aln, $self->params);
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

sub _realign {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  
  my $aln_file = $self->_save_file('realigned_aln', 'fasta');
  my $aln_f = $aln_file->{full_file};
  if (!-e $aln_f || $self->param('force_recalc')) {
    print "  realigning with Prank...\n";
    $self->param('aligner', 'prank_codon');
    my $pep_aln = $self->_tx_aln($aln);
    $aln = $self->align(undef, $aln, $pep_aln);
    Bio::EnsEMBL::Compara::AlignUtils->to_file($aln, $aln_f);
  }
  print "  loading re-alignment from file\n";
  $aln = Bio::EnsEMBL::Compara::AlignUtils->from_file($aln_f);

  return $aln;
}

sub _mask_aln {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  my $m0 = $self->_save_file('m0_mask', 'perlobj');
  my $m0_f = $m0->{full_file};
  if ($self->param('force_recalc') == 1) {
    unlink($m0_f);
  }

  $aln = $self->mask_substitution_runs($tree, $aln, $m0_f);

  my $masked_file = $self->_save_file('masked_aln', 'fasta');
  my $masked_f = $masked_file->{full_file};
  if (!-e $masked_f || $self->param('force_recalc')) {
    Bio::EnsEMBL::Compara::AlignUtils->to_file($aln, $masked_f);
  }
  return $aln;
}

sub _collect_vars {
  my $self = shift;
  my $ref_member = shift;

  my $ref_tx = $ref_member->get_Transcript;
  my $ref_tx_id = $ref_tx->stable_id;

  my $registry = 'Bio::EnsEMBL::Registry';
  my $var_a = $registry->get_adaptor('Human','variation','variationfeature');

  my $gene = $ref_member->get_Gene;
  my $slice = $gene->feature_Slice;
  
  my $var_features = $var_a->fetch_all_by_Slice($slice);
  my $all_pn = 0;
  my $all_ps = 0;
  my $filt_pn = 0;
  my $filt_ps = 0;
  my $allkg_pn = 0;
  my $allkg_ps = 0;
  foreach my $var_f (@$var_features) {
    my $sets = $var_f->get_all_VariationSets;
    my $is_kg = 0;
    my $set_name = '';
    foreach my $set (@$sets) {
      if ($set->name =~ m/^1000 genomes/i) {
        $is_kg = 1;
      }
      $set_name = $set->name;
    }

    my $ts_vs = $var_f->get_all_TranscriptVariations;
    foreach my $ts_v (@$ts_vs) {
      next unless ($ts_v->transcript_stable_id eq $ref_tx_id);

      my $consequences = $ts_v->consequence_type;
      my $consequence = @{$consequences}[0];

      my $allele_freq = 0;
      # Keep only synonymous or non-synonymous consequences.
      if ($consequence =~ m/synonymous/gi) {
        $allele_freq = allele_freq($var_f);
        printf "  %s %s %s %.2f\n", $var_f->variation_name, $var_f->allele_string, $consequence, $allele_freq;

        if ($consequence =~ m/NON_SYNONYMOUS/i) {
          $all_pn++;
          $allkg_pn++ if ($is_kg);
          $filt_pn++ if ($is_kg && $allele_freq > 0.15);
        } else {
          $all_ps++;
          $allkg_ps++ if ($is_kg);
          $filt_ps++ if ($is_kg && $allele_freq > 0.15);
        }

        my $var_p = $self->replace($self->params, {
          consequence => $consequence,
          set_name => $set_name,
          allele_freq => $allele_freq,
          allele_string => $var_f->allele_string,
          seq_region_start => $var_f->seq_region_start,
          var_name => $var_f->variation_name
                                   });
        $self->store_params_in_table($self->dbc, 'vars', $var_p);
      }
    }
  }
  
  $self->param('all_pn', $all_pn);
  $self->param('all_ps', $all_ps);
  $self->param('filt_pn', $filt_pn);
  $self->param('filt_ps', $filt_ps);
  $self->param('allkg_pn', $allkg_pn);
  $self->param('allkg_ps', $allkg_ps);
}

sub _dnds {
  my $self = shift;
  my $aln = shift;
  my $tax_a = shift;
  my $tax_b = shift;
  my $lbl = shift;

  sub _get_seq {
    my $aln = shift;
    my $tax_id = shift;

    my $ptrn;
    $ptrn = "ENSP0" if ($tax_id == 9606);
    $ptrn = "ENSMMUP0" if ($tax_id == 9544);
    $ptrn = "ENSPPYP0" if ($tax_id == 9600);

    my @seqs = $aln->each_seq;
    my ($seq) = grep {$_->id =~ m/$ptrn/i} @seqs;
    return $seq;
  };

  my $seq_a = _get_seq($aln, $tax_a);
  my $seq_b = _get_seq($aln, $tax_b);

  return unless (defined $seq_a);
  return unless (defined $seq_b);

  my $dn = 0;
  my $ds = 0;

  for (my $i=1; $i < $seq_a->length / 3; $i++) {
    my $lo = ($i - 1) * 3 + 1;
    my $hi = $lo + 2;
    my $codon_a = $seq_a->subseq($lo, $hi);
    my $codon_b = $seq_b->subseq($lo, $hi);

    next if ($codon_a =~ m/-/ || $codon_b =~ m/-/);
    next if ($codon_a =~ m/n/i || $codon_b =~ m/n/i);

    my $diff = 0;
    $diff = 1 if ($codon_a ne $codon_b);

    my $codontable_id = $self->param('bioperl_codontable_id');
    $codontable_id = 1 unless (defined $codontable_id);
    my $aa_a = new Bio::PrimarySeq(-seq => $codon_a)->translate( -codontable_id => $codontable_id)->seq;
    my $aa_b = new Bio::PrimarySeq(-seq => $codon_b)->translate( -codontable_id => $codontable_id)->seq;

    my $nsyn = 0;
    $nsyn = 1 if ($aa_a ne $aa_b);

    #print "$aa_a $aa_b $codon_a $codon_b  non-synonymous\n" if ($nsyn);
    #print "$aa_a $aa_b $codon_a $codon_b  synon\n" if ($diff && !$nsyn);
    
    $dn++ if ($nsyn);
    $ds++ if ($diff && !$nsyn);
  }

  $self->param($lbl."_dn", $dn);
  $self->param($lbl."_ds", $ds);
}

sub _run_paml {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $ref_member = shift;

  my $m0 = $self->_save_file('m0', 'txt');
  my $m7 = $self->_save_file('m7', 'txt');
  my $m8 = $self->_save_file('m8', 'txt');

  $aln = Bio::EnsEMBL::Compara::AlignUtils->copy_aln($aln);

  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);

  my $res;
  my $lines;
  my $params;

  if (!-e $m0->{full_file} || $self->param('force_recalc')) {
    # M0
    my $params = {
      model => 0,
      fix_blength => 0,
      getSE => 1,
      fix_omega => 0,
      omega => 0.2,
      Small_Diff => 1e-7
    };
    print "  running m0...\n";
    $res = Bio::Greg::Codeml->branch_model_likelihood( $treeI, $aln, $self->worker_temp_directory, $self->replace($self->params, $params) );
    $lines = $res->{lines};
    _out($m0->{full_file}, $lines);
  }

  # Copy branch lengths from M0 results to our treeI.
  $lines = _in($m0->{full_file});
  my $m0_tree = Bio::Greg::Codeml->parse_codeml_results($lines);

  print $m0_tree->ascii."\n";
  print $treeI->ascii."\n";
  Bio::EnsEMBL::Compara::TreeUtils->transfer_branchlengths($m0_tree, $treeI);

  $self->store_subs($tree, $aln, $lines, $ref_member, 'subs');

  my ($dnds, $dnds_se) = Bio::Greg::Codeml->parse_m0_dnds($lines);
  $self->param('paml_dnds', $dnds);
  $self->param('paml_dnds_se', $dnds_se);

  my $m0_lnl = Bio::Greg::Codeml->extract_lnL($lines);
  $self->param('m0_lnl', $m0_lnl);

  if (!-e $m7->{full_file} || $self->param('force_recalc')) {
    # M7
    $params = {
      model => 0,
      NSsites => 7,
      fix_blength => 0,
      getSE => 1,
      Small_Diff => 1e-6
    };
    print "  running m7...\n";
    $res = Bio::Greg::Codeml->branch_model_likelihood( $treeI, $aln, $self->worker_temp_directory, $self->replace($self->params, $params) );
    $lines = $res->{lines};
    _out($m7->{full_file}, $lines);
  }

  $lines = _in($m7->{full_file});
  my $m7_lnl = Bio::Greg::Codeml->extract_lnL($lines);
  $self->param('m7_lnl', $m7_lnl);

  if (!-e $m8->{full_file} || $self->param('force_recalc')) {
    # M8
    $params = {
      model => 0,
      NSsites => 8,
      fix_blength => 0,
      omega => 0.3,
      fix_omega => 0,
      getSE => 0,
      Small_Diff => 1e-7
    };
    print "  running m8...\n";
    $res = Bio::Greg::Codeml->branch_model_likelihood( $treeI, $aln, $self->worker_temp_directory, $self->replace($self->params, $params) );
    $lines = $res->{lines};
    _out($m8->{full_file}, $lines);
  }

  $lines = _in($m8->{full_file});
  my $m8_lnl = Bio::Greg::Codeml->extract_lnL($lines);
  $self->param('m8_lnl', $m8_lnl);

  my $pep_aln = $self->_tx_aln($aln);
  my $sitewise_results = $self->parse_paml_output($tree, $aln, $pep_aln, $lines);
  $self->store_sitewise($tree, $pep_aln, $sitewise_results, {omega_table => 'paml_sites'});

}

sub _out {
  my $file = shift;
  my $lines = shift;

  open(OUT, ">$file");
  print OUT join("", @$lines) . "\n";
  close(OUT);
}

sub _in {
  my $file = shift;
  open(IN, $file);
  my @lines = <IN>;
  close(IN);
  return \@lines;
}

sub _get_genome_tree {
  my $self = shift;
  my $tree = shift;

  my $params = {
    keep_species => join(',', map {$_->taxon_id} $tree->leaves)
  };
  my $genome_tree = Bio::EnsEMBL::Compara::ComparaUtils->get_genome_tree_subset($self->compara_dba, $params);

  # Match up genome tree names to gene tree IDs.
  my $id_hash;
  map {$id_hash->{$_->taxon_id} = $_->name} $tree->leaves;
  foreach my $node ($genome_tree->nodes) {
    if ($node->is_leaf) {
      $node->name($id_hash->{$node->taxon_id});
    } else {
      $node->name('');
    }
  }

  my $taxid_to_name = {
    9606 => 'Human',
    9593 => 'Gorilla',
    9598 => 'Chimpanzee',
    9601 => 'Orangutan',
    9544 => 'Macaque',
    9483 => 'Marmoset',
    9478 => 'Tarsier',
    30608 => 'Mouse lemur',
    30611 => 'Bushbaby'
  };

  my @missing_names = ();
  foreach my $tx_id (keys %$taxid_to_name) {
    push @missing_names, $taxid_to_name->{$tx_id} if (! defined $id_hash->{$tx_id});
  }


  @missing_names = sort {$a cmp $b} @missing_names;
  my $missing_string = join(", ", @missing_names);
  $self->param('missing_species', $missing_string);

  return $genome_tree;
}

sub _get_paml_tree {
  my $self = shift;
  my $aln = shift;
  
  my $tree_str = qq^(((((((Human, Chimpanzee), Gorilla), Orangutan), Macaque), Marmoset), Tarsier), (MouseLemur, Bushbaby))^;
  my $tree = Bio::EnsEMBL::Compara::TreeUtils->from_newick($tree_str);
  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
  
  foreach my $leaf ($treeI->leaves) {
    foreach my $i ('', '_1', '_2', '_3') {
      my $new_child = new $leaf;
      $new_child->name($leaf->name.$i);
      
      $leaf->add_child($new_child);
    }
    $leaf->name('');
  }
  
  $tree = Bio::EnsEMBL::Compara::TreeUtils->from_treeI($treeI);
  $tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_aln($tree, $aln);
  $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
  return $treeI;
}


sub _run_slr {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  # Get a file to save SLR results in.
  my $out_f = $self->_save_file('slr', 'out');
  my $slr_full_file = $out_f->{full_file};

  $self->param('slr_file', $out_f->{rel_file});

  my $pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($aln, $self->params);

  $self->param('analysis_action', 'slr');

  if (!-e $slr_full_file || $self->param('force_recalc')) {  
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

  my $id = $self->param('gene_name');
  my $filename = "${id}_$filename_base";

  my $file_params = {
    id => "${id}",
    filename => $filename,
    extension => $ext,
    subfolder => 'data',
  };

  my $file_obj = $self->save_file($file_params);
  return $file_obj;
}

sub _out_aln {
  my $self = shift;
  my $aln = shift;
  my $name = shift;

  my $aln_f = $self->_save_file($name, 'fasta');
  my $pdf_f = $self->_save_file($name, 'pdf');

  if (!-e $aln_f->{full_file} || !-e $pdf_f->{full_file}) {
    my $tmp = $self->worker_temp_directory;
    my $aln_file = $aln_f->{full_file};
    my $pdf_file = $pdf_f->{full_file};

    Bio::EnsEMBL::Compara::AlignUtils->to_file($aln, $aln_file);  

    my $cmd = qq^
library(phylosim)
library(RColorBrewer)
library(plyr)
source("~/src/greg-ensembl/projects/phylosim/PhyloSimPlots.R")
source("~/src/greg-ensembl/projects/primate_hiv/hiv_manuscript_plots.R")

aln <- read.aln("${aln_file}")
rownames(aln) <- ensp.to.species(rownames(aln))
aln <- sort.aln.if(aln, sort.order())

pep.aln <- aln.tx(aln)
n.pages <- ceiling( ncol(pep.aln) / 125)
pdf(file="${pdf_file}", width=30, height=3*n.pages)
sim <- PhyloSim(); sim\$.alignment <- pep.aln
plotAlignment(sim, axis.text.size=12, aln.plot.chars=T, aln.char.text.size=3, num.pages=n.pages)
dev.off()
^;
  #Bio::Greg::EslrUtils->run_r($cmd);
  }
}

sub _get_aln {
  my $self = shift;
  my $tree = shift;
  my $ref_member = shift;

  my $out_f = $self->_save_file('aln', 'fasta');
  my $out_file = $out_f->{full_file};
  $self->param('aln_file', $out_f->{rel_file});

  my $c_dba = $self->compara_dba;

  if (!-e $out_file || $self->param('force_recalc')) {
    print "  fetching aln...\n";
    my $tree_aln_obj =
      Bio::EnsEMBL::Compara::ComparaUtils->get_compara_or_genomic_aln( $c_dba, $tree, $ref_member,
                                                                       $self->params );
    if ($tree_aln_obj == -1) {
      $self->fail_and_die("error_fetching", "Error getting tree or aln!");
    }

    my $aln = $tree_aln_obj->{aln};
    $tree = $tree_aln_obj->{tree};
    my $extra  = $tree_aln_obj->{extra};
    $self->set_params($extra);
    
    Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $aln, { width => 150, full => 1 } ) if ($self->debug);
    if (!defined $tree) {
      $self->fail_and_die("tree_undef", "Tree undefined: [$tree]");
    }
    if (scalar($tree->leaves) < 2) {
      $self->fail_and_die("small_tree", "Tree too small!".' '.$self->param('aln_type').' '.$tree->newick_format);
    }
    if ($aln->length < 50) {
      $self->fail_and_die("small_aln", "Alignment too short!".' '.$self->param('aln_type').' '.$tree->newick_format.' '.$aln->length);
    }
    $aln = Bio::EnsEMBL::Compara::AlignUtils->sort_by_tree($aln,$tree);
    my $pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($aln, $self->params);
    Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $aln, { width => 150, full => 1 } ) if ($self->debug);
    Bio::EnsEMBL::Compara::AlignUtils->to_file($aln, $out_file);
  } else {
    print "  loading aln from file $out_file\n";
  }
  
  my $aln = Bio::EnsEMBL::Compara::AlignUtils->from_file($out_file);
  $tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_aln( $tree, $aln );
  return ($tree, $aln);
}

sub _remove_paralogs {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $ref_member = shift;

  my $tax_id_hash;
  foreach my $leaf ($tree->leaves) {
    my $tx_id = $leaf->taxon_id;
    if (!defined $tax_id_hash->{$tx_id}) {
      $tax_id_hash->{$tx_id} = 0;
    } else {
      $tax_id_hash->{$tx_id} = 1;
    }
  }

  my $paralog_count = 0;
  map {$paralog_count += $tax_id_hash->{$_}} keys %$tax_id_hash;
  $self->store_param('species_with_paralogs', $paralog_count);

  # Calculate distances between each sequence. Turn Ns to gaps so they don't mess up
  # the calculations.
  my $map = {
    'N' => '-'
  };
  my $aln_copy = Bio::EnsEMBL::Compara::AlignUtils->translate_chars($aln, $map);
  my $stats = Bio::Align::DNAStatistics->new();
  my $jcmatrix = $stats->distance(-align => $aln_copy, 
                                  -method => 'D_JukesCantor');

  my $ref_id = $ref_member->stable_id;  
  my  ($ref_seq) = grep {$_->id eq $ref_id} $aln->each_seq;
  foreach my $taxon_id (keys %$tax_id_hash) {
    if ($tax_id_hash->{$taxon_id} > 0) {
      my @entries = grep {$_->taxon_id == $taxon_id} $tree->leaves;
      @entries = sort {$a->seq_length <=> $b->seq_length || 
                         $a->stable_id cmp $b->stable_id} @entries;

      my $min_dist = 9999;
      my $min_id = $entries[0];
      foreach my $entry (@entries) {
        my $other_id = $entry->stable_id;
        my $dist = $jcmatrix->get_entry($ref_id, $other_id);
        printf "%s %.3f\n", $other_id, $dist;
        my ($other_seq) = grep {$_->id eq $other_id} $aln->each_seq;

        $min_id = $other_id if ($dist <= $min_dist && $dist > 0);
        $min_dist = $dist if ($dist <= $min_dist && $dist > 0);
      }

      print "Min ID for $taxon_id: $min_id $min_dist\n";

      foreach my $entry (@entries) {
          # Remove all but the closest paralog from the alignment.
        if ($entry->stable_id ne $min_id) {
          my $cur_id = $entry->stable_id;
          print "Removing $cur_id\n";
          $aln = Bio::EnsEMBL::Compara::AlignUtils->remove_seq_from_aln($aln, $entry->stable_id);
        }
      }
    }
  }

  $tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_aln( $tree, $aln );
  return ($tree, $aln);
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
      my $ref_member_id = $ref_member->stable_id;
      print "ref_member_id: [$ref_member_id]\n" if ($self->debug);
      my $mba        = $self->compara_dba->get_MemberAdaptor;
      $ref_member = $mba->fetch_by_source_stable_id( undef, $ref_member_id );
      my $lo_coords = $self->get_coords_from_pep_position($ref_member, $lo);
      my $hi_coords = $self->get_coords_from_pep_position($ref_member, $hi);

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
    my @decent_omegas;
    my $omega_total = 0;
    foreach my $site (@window_sites) {
      my $cur_omega = $site->{omega};
      # We cap the max dN/dS at 3 for the purpose of calculating the mean.
      $cur_omega = 3 if ($cur_omega > 3);
      push @decent_omegas, $cur_omega;
      $omega_total += $cur_omega;
    }
    my $mean_dnds = -1;
    if (scalar(@window_sites) > 0) {
      $mean_dnds = sprintf "%.3f", $omega_total / scalar(@window_sites);
    }

    my $std_dev = $self->standard_deviation(\@decent_omegas);
    my $omega_lo = $mean_dnds - $std_dev;
    my $omega_hi = $mean_dnds + $std_dev;
    $omega_lo = 0 if ($omega_lo < 0);

    print "window[$lo-$hi] pval[$pval]\n" if ($self->debug);
    my @pos_sites = grep {$_->{omega} > 1} @window_sites;
    my $added_params = {
      'n_leaves' => scalar($tree->leaves),
      'pval' => $pval,
      'n_sites' => scalar(@window_sites),
      'n_pos_sites' => scalar(@pos_sites),
      'mean_omega' => $mean_dnds,
      'omega_lower' => $omega_lo,
      'omega_upper' => $omega_hi
    };
    $cur_params = $self->replace($cur_params,$added_params);
    
    if ($self->within_hive) {
      $self->store_params_in_table($self->dbc, 'windows', $cur_params);
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

sub genes_table {
  return {
    job_id => 'int',
    data_id => 'int',
    data_prefix => 'char4',
    
    gene_name => 'char32',
    gene_id => 'char32',
    transcript_id => 'char32',
    protein_id => 'char32',

    chr_name => 'char16',
    chr_start => 'int',
    chr_end => 'int',
    chr_strand => 'int',

    n_seqs => 'int',
    aln_length => 'int',
    seq_length => 'int',

    slr_dnds => 'float',
    slr_kappa => 'float',
    paml_dnds => 'float',
    paml_dnds_se => 'float',

    m0_lnl => 'float',
    m7_lnl => 'float',
    m8_lnl => 'float',

    all_ps => 'int',
    all_pn => 'int',
    filt_ps => 'int',
    filt_pn => 'int',
    allkg_pn => 'int',
    allkg_ps => 'int',
    
    chimp_ds => 'int',
    chimp_dn => 'int',
    orang_ds => 'int',
    orang_dn => 'int',
    rhesus_ds => 'int',
    rhesus_dn => 'int',

    masked_nucs => 'int',
    masked_ids => 'string',
    missing_species => 'string',
    species_with_paralogs => 'int',

    unique_keys => 'data_id',    
  };
}

sub vars_table {
  my $self = shift;

  my $structure = {
    job_id => 'int',
    data_id => 'int',
    gene_name => 'char32',
    seq_region_start => 'int',
    var_name => 'char32',

    set_name => 'char32',
    allele_string => 'char32',
    allele_freq => 'float',
    consequence => 'char32',

    unique_keys => 'data_id,var_name'
  };
  
  return $structure;
}

sub windows_table {
  my $self = shift;

  my $structure = {
    job_id => 'int',
    data_id   => 'int',
    data_prefix => 'char4',

    gene_name => 'char32',
    gene_id => 'char32',
    transcript_id => 'char32',
    protein_id => 'char32',

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

    pval => 'float',
    mean_omega => 'float',
    omega_lower => 'float',
    omega_upper => 'float',
    n_leaves => 'int',
    n_sites    => 'int',
    n_pos_sites    => 'int',

    unique_keys => 'data_id,peptide_window_start,peptide_window_width',
  };

  return $structure;
}

sub allele_freq {
	my $vf = shift;

	my $ref_allele = $vf->ref_allele_string;

	my %freqs;

	for my $allele (@{$vf->get_all_Alleles}) {
		
		if(!defined $allele->population){ 
			#print "fail: pop\n";
			next;
		}

		my $pop_name = $allele->population->name;

		if($pop_name !~ /^1000GENOMES/){ 
			#print "fail: not 1000G $pop_name\n";
			next 
		};

		if($allele->allele eq $ref_allele){
			#print "fail: ref_allele " . $allele->allele . "\n";
			next;
		}

		if(defined $allele->frequency){ 
			$freqs{$pop_name} += $allele->frequency;
		} else {
			#print "fail: not freq\n";
		}

	}

	my @tmp = values %freqs;

	my $nr_freqs = scalar @tmp;
	
	if($nr_freqs == 0){
		return 'NA';
	}

	my $average_freq = ( List::Util::sum(@tmp)/$nr_freqs);

	#print Dumper \%freqs;

	return $average_freq;

}


1;
