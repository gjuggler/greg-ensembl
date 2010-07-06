package Bio::Greg::Hive::GenomewideOmegas;

use strict;
use Bio::Greg::Codeml;
use File::Path;

use base ('Bio::Greg::Hive::Process');

my $TREE = 'Bio::EnsEMBL::Compara::TreeUtils';

sub table_def {
  my $self = shift;

  return {
    data_id => 'int',
 
    method => 'char32',
    
    newick_t => 'string',
    newick_dnds => 'string',
    newick_ds => 'string',
    newick_species => 'string',

    dnds_Hsap => 'float',
    dnds_Ggor => 'float',
    dnds_Ppyg => 'float',
    dnds_Mmul => 'float',

    n_species => 'int',
    n_sites => 'int',
    n_gappy_sites => 'int',
    n_nongappy_sites => 'int',
    
    unique_keys => 'data_id,method'
  };
}

sub fetch_input {
  my ($self) = @_;

  my $default_params = {
    foreground_species => '9606,9593',
    keep_species => '9606,9593,9600,9544',
    method => 'paml',
    table => 'stats_dnds'
  };

  # Fetch parameters from all possible locations.
  $self->load_all_params($default_params);
  $self->create_table_from_params($self->compara_dba, $self->param('table'), $self->table_def);
}


sub run {
  my $self = shift;

  my $cur_params = $self->params;

  my @keepers = split(',',$self->param('keep_species'));
  my @foreground = split(',',$self->param('foreground_species'));

  @keepers = grep {$_ != 9598} @keepers;

  my $f_binom = sub {my $self = shift;$self->binomial};
  my $f_taxid = sub {my $self = shift;$self->is_leaf ? $self->short_name : ''};

  # Prepare the tree:
  my $tax_tree = Bio::EnsEMBL::Compara::ComparaUtils->get_genome_taxonomy_below_level($self->compara_dba,'Mammalia');
  $tax_tree = Bio::EnsEMBL::Compara::TreeUtils->keep_members_by_method_call($tax_tree,\@keepers,'taxon_id');
  my $categories;
  my $i=1;
  foreach my $fg_taxon_id (@foreground) {
    $categories->{$fg_taxon_id} = $i++;
  }
  my $labeled_tree = $self->categorize_nodes($tax_tree,$categories);
  print $labeled_tree->newick_format($f_taxid)."\n";
  # Turn the labeled Ensembl tree object into a BioPerl TreeI-compliant object.
  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($labeled_tree,$f_taxid);

  # Get an alignment file of all concatenated genes.
  my $base = '/nfs/users/nfs_g/gj1/scratch/gorilla/dnds';
  $self->get_output_folder($base);
  my $folder = $self->get_output_folder;
  my $file = "$folder/aln_".$self->data_id.".fasta";
  my $aln;
  if (!-e $file) {
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
 print "Final length: ".$aln->length."\n";
  $self->pretty_print($aln,{length => 200});

  # Use the Codeml.pm helper function to run Codeml.
  my $tmpdir = $self->worker_temp_directory;
  $tmpdir = '/tmp/pamltest';
  mkpath([$tmpdir]);
  my $params = $self->params;
  $params->{model} = 2; # Use pre-specified foreground / background branches.
  $params->{cleandata} = 1; # Don't use sites with any gaps.
  my $results = Bio::Greg::Codeml->branch_model_likelihood($treeI,$aln,$tmpdir,$params);

  my @omegas = @{$model->{omegas}};
  for (my $i=0; $i < scalar @omegas; $i++) {
    my $omega = $omegas[$i];
    
    $self->store_tag("${key}omega_${i}",$omega);
  }

#  my $dnds = $results->{dnds};
#
#  foreach my $key (keys %$dnds) {
#    print "$key:\n";
#    my $obj = $dnds->{$key};
#
#    $cur_params->{'dnds_'.$key} = $obj->{'dN/dS'};
#    
#    $self->store($obj,$key);
#    Bio::EnsEMBL::Compara::ComparaUtils->hash_print($dnds->{$key});
#  }

  $cur_params->{n_sites} = $aln->length;

  my $T = 'Bio::EnsEMBL::Compara::TreeUtils';
  $cur_params->{newick_dnds} = $T->to_newick($results->{dnds_tree});
  $cur_params->{newick_t} = $T->to_newick($results->{t_tree});
  $cur_params->{newick_ds} = $T->to_newick($results->{ds_tree});
  $cur_params->{newick_species} = $T->to_newick($treeI);
  $cur_params->{n_species} = scalar($tax_tree->leaves);
  
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

  my $seq_hash = {};
  my $fetch_params = $self->params;

  my $len;
  foreach my $node_id (@node_ids[0..50]) {
    print "node: $node_id\n";
    $fetch_params->{node_id} = $node_id;
    my $tree = $self->get_tree($fetch_params);
    next if (scalar($tree->leaves) < 2);
    next if (!Bio::EnsEMBL::Compara::TreeUtils->has_one_to_one_orthology($tree,\@taxon_ids));
    print " ok!\n";

    $fetch_params->{tree} = $tree;
    my @leaves = $fetch_params->{tree}->leaves;
    
    my $aln = $self->get_cdna_aln($fetch_params);
    my $id_to_taxon;
    map {$id_to_taxon->{$_->stable_id} = $_->taxon_id} @leaves;

    $aln = Bio::EnsEMBL::Compara::AlignUtils->translate_ids($aln,$id_to_taxon);

    #$self->pretty_print($aln,{length => 200});
    die unless ($aln->length % 3 == 0);

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
  }

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

  return $aln;
}


# Uses a mapping from ncbi_taxon_id => rate_category to create a PAML-formatted
# branch model tree for terminal branches of various species.
sub categorize_nodes {
  my $self = shift;
  my $tree = shift;
  my $taxon_categories = shift;

  # Create a copy of the tree.
  my $copy = $tree->copy;

  foreach my $leaf ($copy->nodes) {
    # If this leaf's taxon_id exists in the mapping, append that value (along with '#') to the label.
    # Otherwise, leave it alone.
    my $category = $taxon_categories->{$leaf->taxon_id} || '';
    $leaf->taxon_id($leaf->taxon_id."#$category") if (defined $category && $category ne '');
  }

  return $copy;
}


1;
