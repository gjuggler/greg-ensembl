package Bio::Greg::Hive::GenomewideOmegas;

use strict;
use Bio::Greg::Codeml;
use File::Path;
use Bio::Greg::Gorilla::Utils;

use base ('Bio::Greg::Hive::Process');

my $TREE = 'Bio::EnsEMBL::Compara::TreeUtils';

sub table_def {
  my $self = shift;

  return {
    data_id => 'int',
 
    method => 'char32',
    branch_model => 'int',
    Mgene => 'int',
    
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
    dnds_Mmul => 'float',

    n_species => 'int',
    n_sites => 'int',
    n_genes => 'int',
    n_gappy_sites => 'int',
    n_nongappy_sites => 'int',
    
    unique_keys => 'data_id,method'
  };
}

sub fetch_input {
  my ($self) = @_;

  my $default_params = {
    branch_model => 2,
    Mgene => -1,
    foreground_species => '9606,9593',
    keep_species => '9606,9593,9600,9544',
    method => 'paml',
    table => 'stats_dnds',
    alignment_cache_dir => '/nfs/users/nfs_g/gj1/scratch/gorilla/dnds',
    always_fetch_new_alignments => 1,
    target_alignment_length => '100000'
  };

  # Fetch parameters from all possible locations.
  $self->load_all_params($default_params);
  $self->create_table_from_params($self->compara_dba, $self->param('table'), $self->table_def);
}


sub run {
  my $self = shift;

  my $cur_params = $self->params;
  my $codeml_params = $self->params;

  my @keepers = split(',',$self->param('keep_species'));
#  @keepers = grep {$_ != 9598} @keepers;

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
  my $tax_tree = Bio::EnsEMBL::Compara::ComparaUtils->get_genome_taxonomy_below_level($self->compara_dba,'Mammalia');
  $tax_tree = Bio::EnsEMBL::Compara::TreeUtils->keep_members_by_method_call($tax_tree,\@keepers,'taxon_id');

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
    @foreground = split(',',$self->param('foreground_species'));
    # Sort the taxon IDs before iterating so we have a stable ordering of omega categories.
    foreach my $fg_taxon_id (sort {$a <=> $b} @foreground) {
      $categories->{$fg_taxon_id} = $i++;
    }
    $labeled_tree = $self->categorize_nodes($tax_tree,$categories);
  }

  print "Labeled tree: \n";
  print "  ". $labeled_tree->newick_format($f_shortname)."\n";

  # Turn the labeled Ensembl ncbi taxon tree object into a BioPerl TreeI-compliant object.
  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($labeled_tree,$f_shortname);

  # Get an alignment file of all concatenated genes.
  my $base = $self->param('alignment_cache_dir');
  $self->get_output_folder($base);
  my $folder = $self->get_output_folder;
  my $file = "$folder/aln_".$self->data_id.".fasta";
  my $aln;
  # If the cache doesn't already exist.
  if (!-e $file || $self->param('always_fetch_new_alignments')) {
    $aln = $self->combine_all_alignments(\@keepers);
    my $out = new Bio::AlignIO(-format => 'fasta',
			       -file => '>'.$file);
    $out->force_displayname_flat(1);
    $out->write_aln($aln);
  } else {
    my $in = new Bio::AlignIO(-format => 'fasta',
			       -file => $file);
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
  }

  $codeml_params->{cleandata} = 1; # Don't use sites with any gaps.
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
  $cur_params->{n_genes} = scalar(@{$self->param('gene_codon_counts')});
  
  $self->store_params_in_table( $self->db_handle, $self->param('table'), $cur_params);
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

sub combine_all_alignments {
  my $self = shift;
  my $keepers_arrayref = shift;

  my @taxon_ids = @{$keepers_arrayref};

  my $sth = $self->db_handle->prepare("select node_id from protein_tree_tag where tag='cc_root_Primates';");
  $sth->execute;
  my $ref = $sth->fetchall_arrayref([0]);
  my @node_ids = map {$_->[0]} @{$ref};
  $sth->finish;

  # Look at the stats_trees table for some idea if this tree's any good or not.
  $sth = $self->db_handle->prepare("select * from stats_trees where node_id=? limit 1");

  my $seq_hash = {};
  my $fetch_params = $self->params;

  my $len;
  my $kept_gene_count = 0;
  my @gene_codon_counts = ();
  foreach my $node_id (@node_ids) {
    #last if ($kept_gene_count >= 1000);
    last if (defined $self->param('target_alignment_length') && $len > $self->param('target_alignment_length'));

    print "node: $node_id\n";
    $fetch_params->{node_id} = $node_id;
    my $tree = $self->get_tree($fetch_params);

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

    $fetch_params->{tree} = $tree;
    my @leaves = $fetch_params->{tree}->leaves;
    my $aln = $self->get_cdna_aln($fetch_params);
    my $id_to_taxon;
    map {$id_to_taxon->{$_->stable_id} = $_->taxon_id} @leaves;
    $aln = Bio::EnsEMBL::Compara::AlignUtils->translate_ids($aln,$id_to_taxon);

    die unless ($aln->length % 3 == 0);
    my $orig_aln = $aln;
    my $first_keeper = $taxon_ids[0];
    my $second_keeper = $taxon_ids[1];
    $aln = Bio::EnsEMBL::Compara::AlignUtils->remove_triple_mutated_codons($aln,$first_keeper,$second_keeper);

    #  3) Has more than 5 'funky' codons (codons with all 3 nucleotides mutated b/t human and gorilla)
    $why_not = "Too many funky codons!" if ($orig_aln->length - $aln->length > 5);

    if ($why_not) {
      print "  -> Skipping $node_id [$why_not]\n";
      next;
    }

    print "$node_id ok! $kept_gene_count $len\n";
    $kept_gene_count++;

    # Add any new species from the alignment to the hash.
    map {$seq_hash->{$_->id} = [] if (!defined $seq_hash->{$_->id});} $aln->each_seq;

    # Fill in each sequence in the hash with the current alignment, or gaps if that species is currently missing.
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
    print $aln->length."\n";

    $len += $aln->length;
    push @gene_codon_counts, ($aln->length / 3);
  }

  $sth->finish;

  # Make a final alignment.
  my $aln = new Bio::SimpleAlign;
  foreach my $id (keys %$seq_hash) {
    my $seq = join('',@{$seq_hash->{$id}});
    $aln->add_seq(Bio::LocatableSeq->new(-seq => $seq,
					 -id => $id));
  }

  print "Aln length: ".$aln->length."\n";
  $self->pretty_print($aln,{length => 200});
  #$aln = Bio::EnsEMBL::Compara::AlignUtils->remove_gappy_columns_in_threes($aln);
  #print "Aln length: ".$aln->length."\n";
  #$self->pretty_print($aln,{length => 200});

  $self->param('gene_codon_counts',\@gene_codon_counts);
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
  my $copy = $tree->copy;

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
