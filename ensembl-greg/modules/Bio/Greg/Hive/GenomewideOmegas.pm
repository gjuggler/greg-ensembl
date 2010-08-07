package Bio::Greg::Hive::GenomewideOmegas;

use strict;
use Bio::Greg::Codeml;
use File::Path;
use Bio::Greg::Gorilla::Utils;
use Bio::Align::Utilities qw(:all);

use Time::HiRes qw(sleep);
use DateTime;
use DateTime::Format::MySQL;

use base ('Bio::Greg::Hive::Process', 'Bio::Greg::Hive::Align');

my $TREE = 'Bio::EnsEMBL::Compara::TreeUtils';

sub table_def {
  my $self = shift;

  return {
    data_id => 'int',
 
    label => 'char32',

    method => 'char32',
    branch_model => 'int',
    Mgene => 'int',
    cleandata => 'int',
    
    work_started => 'timestamp',
    work_finished => 'timestamp',

    newick_t => 'string',
    newick_dnds => 'string',
    newick_ds => 'string',
    newick_species => 'string',
    foreground_species => 'string',

    dnds_Hsap => 'float',
    dnds_Ggor => 'float',
    dnds_Ppyg => 'float',
    dnds_Mmul => 'float',

    dnds_Ptro => 'float',
    dnds_Cfam => 'float',
    dnds_Mmus => 'float',
    dnds_Background => 'float',

    n_species => 'int',
    n_sites => 'int',
    n_genes => 'int',

    target_alignment_length => 'int',
    remove_funky_codons => 'int',
    unroot_species_tree => 'int',
    
    unique_keys => 'data_id,label'
  };
}

sub fetch_input {
  my ($self) = @_;

  my $default_params = 
    {
     job_role => '', # set to 'fan_jobs' to fan out alignment gathering jobs.
     
     branch_model => 2,
     Mgene => -1,
     method => 'paml',
     cleandata => 0,
     
     foreground_species => '9606,9593',
     keep_species => '9606,9593,9600,9544',
     
     table => 'stats_dnds',
     alignment_cache_dir => '/nfs/users/nfs_g/gj1/scratch/gorilla/dnds',
     always_fetch_new_alignments => 1,
     
     aln_export_job_size => 100,
     
     target_alignment_length => 6000000,
     remove_funky_codons => 1,
     remove_gappy_columns => 1,
     unroot_species_tree => 0,
     prank_realign => 0
    };
  
  # Fetch parameters from all possible locations.
  $self->load_all_params($default_params);
  $self->create_table_from_params($self->compara_dba, $self->param('table'), $self->table_def);
}


sub species_set_id {
  my $self = shift;

  my @keepers = $self->get_taxon_ids_from_keepers_list($self->param('keep_species'));
  
  my $i=0;
  map {$i += $_} @keepers;

  return $i;
}

sub run {
  my $self = shift;

  my $dt = DateTime->now(time_zone => 'GMT');
  my $f = DateTime::Format::MySQL->format_datetime($dt);
  $self->param('work_started',$f);
  $self->param('label',$self->data_id);

  my $cur_params = $self->params;
  my $codeml_params = $self->params;

  $self->store_params_in_table( $self->dbc, $self->param('table'), $cur_params);

  my @keepers = $self->get_taxon_ids_from_keepers_list($self->param('keep_species'));
  #my @keepers = split(',',$self->param('keep_species'));
#  @keepers = grep {$_ != 9598} @keepers;

  print "Keeping: @keepers\n";
  my $f_binom = sub {my $self = shift;$self->binomial};
  my $f_shortname = sub {
    my $self = shift;
    if ($self->is_leaf && $self->isa("Bio::EnsEMBL::Compara::NCBITaxon")) {
      return $self->short_name;
    } elsif ($self->is_leaf) {
      return $self->genome_db->short_name(@_);
    } else {
      return '';
    }
  };
  my $f_taxid = sub {my $self = shift;$self->is_leaf ? $self->taxon_id(@_) : ''};
  my $f_name = sub {my $self = shift;$self->is_leaf ? $self->name(@_) : ''};

  # Prepare the taxonomic tree.
  my $tax_tree = Bio::EnsEMBL::Compara::ComparaUtils->get_genome_taxonomy_below_level($self->compara_dba,'Chordata');
  $tax_tree = Bio::EnsEMBL::Compara::TreeUtils->keep_members_by_method_call($tax_tree,\@keepers,'taxon_id');

  Bio::EnsEMBL::Compara::ComparaUtils->fix_genome_polytomies($tax_tree);

  print "Taxonomic tree: ".$tax_tree->newick_format."\n";

  my $branch_model = $self->param('branch_model');

  # Always fall back to a one-rate model if we've just got two species.  
  $branch_model = 0 if (scalar $tax_tree->leaves == 2);

  my $labeled_tree = $tax_tree;
  my @foreground;
  print "BRANCH MODEl $branch_model\n";
  if ($branch_model == 0) {
    # One-ratio model.
    $codeml_params->{model} = 0;
  } elsif ($branch_model == 1) {
    # Free-ratios model. Nothing extra to do here.
    $codeml_params->{model} = 1;
  } elsif ($branch_model == 2) {
    # Foreground-background model. Create a categories map and label up our foreground branch(es).
    $codeml_params->{model} = 2;

    my $categories;
    my $i=1;
    @foreground = ('9999'); # Code for all the 'background' species.
    push @foreground, split(',',$self->param('foreground_species'));
    # Sort the taxon IDs before iterating so we have a stable ordering of omega categories.
    foreach my $fg_taxon_id (sort {$a <=> $b} @foreground) {
      $categories->{$fg_taxon_id} = $i++;
    }
    $labeled_tree = $self->categorize_nodes($tax_tree,$categories);
  }

  if ($self->param('unroot_species_tree') == 1) {
    $labeled_tree = Bio::EnsEMBL::Compara::TreeUtils->unroot($labeled_tree);
  }
  print "Taxonomic tree: ".$labeled_tree->newick_format."\n";

  # If we're given a node_id_list, just output aln files for those node_ids and exit.
  if ($self->param('job_role') eq 'fetch_alns') {
    $self->combine_all_alignments(\@keepers);
    return;
  } elsif ($self->param('job_role') eq 'fan_jobs') {
    $self->fan_jobs(\@keepers);
    return;
  }

  # Turn the labeled Ensembl ncbi taxon tree object into a BioPerl TreeI-compliant object.
  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($labeled_tree,$f_shortname);

  # Get an alignment file of all concatenated genes.
  my $base = $self->param('alignment_cache_dir');
  $self->get_output_folder($base);
  my $folder = $self->get_output_folder;
  my $file_aln = "$folder/aln_".$self->data_id.".fasta";
  my $file_genelist = "$folder/aln_".$self->data_id."_genes.txt";
  my $aln;
  print "ALIGNMENT FILE: $file_aln\n";
   # If the cache doesn't already exist.
  if (!-e $file_aln || $self->param('always_fetch_new_alignments') == 1) {
    $aln = $self->combine_all_alignments(\@keepers);
    my $out = new Bio::AlignIO(-format => 'fasta',
			       -file => '>'.$file_aln);
    $out->force_displayname_flat(1);
    $out->write_aln($aln);

    open(OUT,">$file_genelist");
    print OUT join("\n",@{$self->param('gene_id_list')});
    close(OUT);
  } else {
    print "LOADING BIG ALIGNMENT FROM FILE!\n";
    my $in = new Bio::AlignIO(-format => 'fasta',
			       -file => $file_aln);
    $aln = $in->next_aln;
  }
  my $tx;
  foreach my $taxon ($tax_tree->leaves) {
    $tx->{$taxon->taxon_id} = $taxon->short_name;
  }
  $aln = Bio::EnsEMBL::Compara::AlignUtils->translate_ids($aln, $tx); 
  $self->pretty_print($aln,{length => 200});

  # Run Codeml.
  my $tmpdir = $self->worker_temp_directory;
  $tmpdir = '/tmp/pamltest';
  mkpath([$tmpdir]);
  
  if ($self->param('Mgene') > -1) {
    # Provide per-gene codon counts to the Codeml program.
    $codeml_params->{gene_codon_counts} = $self->param('gene_codon_counts');
    $codeml_params->{Mgene} = $self->param('Mgene');
  } else {
    delete $codeml_params->{Mgene};
  }

  my $results = Bio::Greg::Codeml->branch_model_likelihood($treeI,$aln,$tmpdir,$codeml_params);

  if ($branch_model == 0) {
    my @omegas = @{$results->{omegas}};
    my $omega = $omegas[0];

    foreach my $short_name (values %$tx) {
      $cur_params->{'dnds_'.$short_name} = $omega;      
    }
  }

  if ($branch_model == 2) {
    my @omegas = @{$results->{omegas}};
    print "@omegas\n";

    for (my $i=0; $i < scalar @omegas; $i++) {
      my $omega = $omegas[$i];
      
      # match up this omega value with the foreground branch
      my $taxon_id = $foreground[$i];
      my $key = Bio::Greg::Gorilla::Utils->taxon_short($taxon_id);
      #$self->store_tag("dnds_${key}",$omega);
      $cur_params->{'dnds_'.$key} = $omega;
    }
  }

  if ($branch_model == 1) {
    my $dnds = $results->{dnds};

    foreach my $key (keys %$dnds) {
      print "$key:\n";
      my $obj = $dnds->{$key};
      
      $cur_params->{'dnds_'.$key} = $obj->{'dN/dS'};
      
      #$self->store($obj,$key);
      Bio::EnsEMBL::Compara::ComparaUtils->hash_print($dnds->{$key});
    }
  }

  $cur_params->{n_sites} = $aln->length;
  $cur_params->{branch_model} = $branch_model;

  my $T = 'Bio::EnsEMBL::Compara::TreeUtils';
  $cur_params->{newick_dnds} = $T->to_newick($results->{dnds_tree});
  $cur_params->{newick_t} = $T->to_newick($results->{t_tree});
  $cur_params->{newick_ds} = $T->to_newick($results->{ds_tree});
  $cur_params->{newick_species} = $T->to_newick($treeI);
  $cur_params->{n_species} = scalar($tax_tree->leaves);
  if (defined $self->param('gene_codon_counts')) {
    $cur_params->{n_genes} = scalar(@{$self->param('gene_codon_counts')});
  }

  my $dt_finish = DateTime->now(time_zone => 'GMT');
  my $finish = DateTime::Format::MySQL->format_datetime($dt_finish);
  $cur_params->{'work_finished'} = $finish;
  
  $self->store_params_in_table( $self->dbc, $self->param('table'), $cur_params);
}

sub get_prank_aln {
  my $self = shift;
  my $cdna_aln = shift;
  my $tree = shift;

  print " -> Aligning with Prank!\n";
  print $tree->newick_format."\n";
  my $n = $cdna_aln->length;
  print "Before: $n\n";

  # Align with Prank to try and de-align incorrectly called exons.
  my $pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($cdna_aln);
  my $prank_params = {
    alignment_prank_f => 0
  };
  my $new_pep = $self->align_with_prank($pep_aln,$tree,$prank_params);
  my $new_aln = Bio::EnsEMBL::Compara::AlignUtils->peptide_to_cdna_alignment($new_pep,$tree);

  $n = $new_aln->length;
  print "After $n:\n";
  $self->pretty_print($new_aln,{length => 200});

  return $new_aln;
}

sub store {
  my $self = shift;
  my $obj = shift;
  my $prefix = shift;

  return if ($prefix =~ m/^\d+$/);

  $self->param('node_id',1);
  foreach my $key (keys %$obj) {
    next unless ($key eq 'dN/dS');
    my $k = $prefix.'_'.$key;
    $self->store_tag($k,$obj->{$key});
    #$self->store_meta({$k => $obj->{$key}});
  }
}

sub get_node_file {
  my $self = shift;
  my $node_id = shift;

  my $species_set_id = $self->species_set_id;
  return $self->get_hashed_file("genomewide_alns","${species_set_id}_${node_id}.fasta");
}

sub fan_jobs {
  my $self = shift;
  my $keepers_arrayref = shift;

  my @node_ids;
    my $sth = $self->dbc->prepare("select node_id from protein_tree_tag where tag='cc_root_Primates';");
    $sth->execute;
    my $ref = $sth->fetchall_arrayref([0]);
    @node_ids = map {$_->[0]} @{$ref};
    $sth->finish;

    sub output_nodes {
      my $self = shift;
      my @nodes = @_;

      my $output_id = $self->string_to_hash($self->input_id);
      $output_id->{node_id_list} = join(',',@nodes);
      $output_id->{job_role} = 'fetch_alns';
      print "Fanning out nodes:\n";
      $self->hash_print($output_id);
      my ($output_job_id) = @{$self->dataflow_output_id($output_id,99)};
      print " --> Flowed out to job $output_job_id\n";
      sleep(0.1);
    }
   
    my $n = $self->param('aln_export_job_size');
    my @temp_array;
    for (my $i=0; $i < scalar(@node_ids)-1; $i++) {
      my $node_id = $node_ids[$i];
#      if (!-e $self->get_node_file($node_id)) {
	#print "Need to get file: $node_id\n";
	push @temp_array, $node_id;
#      }
      if (scalar(@temp_array) > $n) {
	$self->output_nodes(@temp_array);
	@temp_array = ();
      }
    }
    if (scalar(@temp_array) > 0) {
      $self->output_nodes(@temp_array);
    }
}

sub combine_all_alignments {
  my $self = shift;
  my $keepers_arrayref = shift;

  my $cur_params = $self->params;
  my $fetch_params = $self->params;

  my @taxon_ids = @{$keepers_arrayref};

  my @alns;

  my @node_ids;
  if (defined $self->param('node_id_list')) {
    @node_ids = split(',',$self->param('node_id_list'));
    #@node_ids = @node_ids[0..4];
  } else {
    my $sth = $self->dbc->prepare("select node_id from protein_tree_tag where tag='cc_root_Primates';");
    $sth->execute;
    my $ref = $sth->fetchall_arrayref([0]);
    @node_ids = map {$_->[0]} @{$ref};
    $sth->finish;
  }

  # Look at the stats_trees table for some idea if this tree's any good or not.
  my $sth = $self->dbc->prepare("select * from stats_trees where node_id=? limit 1");

  my $seq_hash = {};
  map {$seq_hash->{$_} = [];} @taxon_ids;

  my $tree;
  my $len;
  my $kept_gene_count = 0;
  my @gene_codon_counts = ();
  my @gene_id_list = ();
  my $stop=0;
  foreach my $node_id (@node_ids) {
    #last if ($kept_gene_count >= 1000);
    last if (defined $self->param('target_alignment_length') && 
	     $len > $self->param('target_alignment_length') && 
	     $self->param('target_alignment_length') > -1 && 
	     $self->param('job_role') ne 'fetch_alns');

    delete $fetch_params->{tree};
    $fetch_params->{node_id} = $node_id;
    $tree = $self->get_tree($fetch_params);

    if ($self->param('job_role') eq 'fetch_alns') {
      
      $sth->execute($node_id);
      my $stats = $sth->fetchrow_hashref;
      
      # RULES FOR INCLUDING A TREE IN THE ANALYSIS:
      my $why_not = undef;
      #  1) Has one-to-one orthology in our species of interest.
      $why_not = "Too few leaves" if (scalar($tree->leaves) < 2);
      $why_not = "Not one-to-one"  if (!Bio::EnsEMBL::Compara::TreeUtils->has_one_to_one_orthology($tree,\@taxon_ids));
      #  2) Has a reasonably-sized gene tree (when including all mammals in the nodeset)
      $why_not = "Tree too big or small!" if ($stats->{tree_length} < 0.5 || $stats->{tree_length} > 10);
      if ($why_not) {
	print "  -> Skipping $node_id [$why_not]\n";
	next;
      }
    }

    my $aln_file = $self->get_node_file($node_id);
    my $aln;
    if ($self->param('job_role') eq 'fetch_alns') {
      $aln = $self->load_and_save_aln($tree,$aln_file,$keepers_arrayref);
    } elsif (-e $aln_file) {
      print " --> File exists!\n";
      my $in = new Bio::AlignIO(-format => 'fasta',
			       -file => $aln_file);
      $aln = $in->next_aln;
      my $tx_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($aln);
      print " --> ".$aln->length()/3 . "\n";
      $in->close;
    } else {
	print "No alignment file found for $node_id! Skipping...\n";
	next;
    }

    next unless (defined $aln && ref $aln && $aln->length > 0);
    $self->pretty_print($aln);

    $len += $aln->length;
    push @gene_codon_counts, ($aln->length / 3);
    push @alns, $aln;

#    if ($self->param('job_role') eq 'fetch_alns') {
      my @leaves = $tree->leaves;
      my ($human_member) = grep {$_->taxon_id == 9606} @leaves;
      my $id = $node_id;
      if (defined $human_member) {
        push @gene_id_list, $human_member->stable_id;
        $id = $human_member->stable_id;
      }
      print "$id added! $kept_gene_count genes [$len nucs] so far\n";
#    }

    $kept_gene_count++;

    $cur_params->{n_sites} = $len;
    $cur_params->{n_genes} = $kept_gene_count;
    $self->store_params_in_table( $self->dbc, $self->param('table'), $cur_params);
  }

  $sth->finish;

  return undef if (scalar(@alns) == 0);

  # Make a final alignment.
  my $aln = cat(@alns);
#  my $aln = new Bio::SimpleAlign;
#  foreach my $id (keys %$seq_hash) {
#    my $seq = join('',@{$seq_hash->{$id}});
#    $aln->add_seq(Bio::LocatableSeq->new(-seq => $seq,
#					 -id => $id));
#  }

  print "Aln length: ".$aln->length."\n";
  $self->pretty_print($aln,{length => 200});

  $self->param('gene_codon_counts',\@gene_codon_counts);
  $self->param('gene_id_list',\@gene_id_list);
  return $aln;
}

sub load_and_save_aln {
  my $self = shift;
  my $tree = shift; 
  my $aln_file = shift;
  my $keepers_arrayref = shift;

  my @taxon_ids = @{$keepers_arrayref};

  my $fetch_params = $self->params;
  $fetch_params->{tree} = $tree;
  my @leaves = $fetch_params->{tree}->leaves;
  my $aln = $self->get_cdna_aln($fetch_params);
  
  my $id_to_taxon;
  map {$id_to_taxon->{$_->stable_id} = $_->taxon_id} @leaves;
  
  my ($orig_aln, $translated_aln, $filtered_aln);
  
  $orig_aln = $aln;

  if ($self->param('prank_realign') == 1) {
    $aln = $self->get_prank_aln($aln,$tree);
  }

  $translated_aln = Bio::EnsEMBL::Compara::AlignUtils->translate_ids($aln,$id_to_taxon);
  
  die unless ($aln->length % 3 == 0);
  
  #$self->pretty_print($aln,{length=>150});

  if ($self->param('remove_funky_codons') == 1) {
    my $first_keeper = '9606';
    my $second_keeper = '9598';
    my $third_keeper = '9593';
    $filtered_aln = Bio::EnsEMBL::Compara::AlignUtils->remove_funky_stretches($translated_aln,$first_keeper,$second_keeper,3);
    $filtered_aln = Bio::EnsEMBL::Compara::AlignUtils->remove_funky_stretches($translated_aln,$second_keeper,$third_keeper,3);
    #$filtered_aln = Bio::EnsEMBL::Compara::AlignUtils->remove_triple_mutated_codons($translated_aln,$first_keeper,$second_keeper);
  } else {
    $filtered_aln = $translated_aln;
  }

  if (my $obj = Bio::EnsEMBL::Compara::AlignUtils->has_stop_codon($filtered_aln)) {
    $self->hash_print($obj);
    die("Stop codon found BEFORE de-gapping! ".$tree->newick_format);
    return undef;
  }

  if ($self->param('remove_gappy_columns') == 1) {
    $self->pretty_print($filtered_aln,{length=>150});
    print " >> Removing gaps...\n";
    $filtered_aln = Bio::EnsEMBL::Compara::AlignUtils->remove_gappy_columns_in_threes($filtered_aln);
    #$self->pretty_print($filtered_aln,{length=>150});
  } else {
  }


  my $why_not;
#  if ($orig_aln->length - $filtered_aln->length > 30) {
#    $why_not = "Too many funky codons!";
#  }
  if ($filtered_aln->length < 30) {
    $why_not = "Too short filtered aln!";
  }

  $aln = $filtered_aln;

  if ($why_not) {
    print "  -> Skipping! [$why_not]\n";
    next;
  }

  my $tx_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($filtered_aln);
  print "Aln length: "+$tx_aln->length."\n";
  $self->pretty_print($tx_aln,{full=>1});


  # Fill in each sequence in the hash with the current alignment, or gaps if that species is currently missing.
  my $seq_hash = {};
  map {$seq_hash->{$_} = [];} @taxon_ids;
  foreach my $id (keys %$seq_hash) {
    my @chars = @{$seq_hash->{$id}};
    my $seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id($aln,$id);
    if (defined $seq) {
      push @chars, split(//,$seq->seq);
    } else {
      print "Missing seq for $id from alignment!\n";
      push @chars, split(//,'-' x $aln->length);
    }
    $seq_hash->{$id} = \@chars;
  }

  # Combine the seq_hash sequences into a new final simplealign.
  my $aln = new Bio::SimpleAlign;
  foreach my $id (keys %$seq_hash) {
    my $seq = join('',@{$seq_hash->{$id}});
    $aln->add_seq(Bio::LocatableSeq->new(-seq => $seq,
					 -id => $id));
  }

  if ($aln->length %3 != 0) {
    die("Alignment length not divisible by 3! ".$tree->newick_format);
    return undef;
  }

  if (Bio::EnsEMBL::Compara::AlignUtils->has_stop_codon($aln)) {
    die("Stop codon found after de-gapping! ".$tree->newick_format);
    return undef;
  }

  
  ## Write aln to file!
  my $out = new Bio::AlignIO(-format => 'fasta',
			     -file => '>'.$aln_file);
  $out->force_displayname_flat(1);
  $out->write_aln($aln);

  return $aln;
}

# Uses a mapping from ncbi_taxon_id => rate_category to create a PAML-formatted
# branch model tree for terminal branches of various species.
sub categorize_nodes {
  my $self = shift;
  my $tree = shift;
  my $taxon_categories = shift;
  my $name_f = shift;

  # Create a copy of the tree.
  my $copy = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree);

  foreach my $leaf ($copy->nodes) {
    # If this leaf's taxon_id exists in the mapping, append that value (along with '#') to the label.
    # Otherwise, leave it alone.
    my $category = $taxon_categories->{$leaf->taxon_id} || '';
    $leaf->name($leaf->short_name);
    $leaf->name($leaf->short_name."#$category") if (defined $category && $category ne '');
  }

  return $copy;
}


1;
