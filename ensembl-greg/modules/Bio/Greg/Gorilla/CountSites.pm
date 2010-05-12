package Bio::Greg::Gorilla::CountSites;

use strict;
use Time::HiRes qw(time gettimeofday tv_interval);
use Cwd;
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::NestedSet;
use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Hive::Process;
use Time::HiRes qw(sleep);
use Bio::EnsEMBL::Registry;
use Bio::Greg::StatsCollectionUtils;

use base ('Bio::Greg::Hive::Process','Bio::Greg::StatsCollectionUtils');

#
# Some global-ish variables.
#

# OUTPUT FILES / OBJECTS / STATES.
my %node_set_hash;

my $TREE = "Bio::EnsEMBL::Compara::TreeUtils";
my $ALN = "Bio::EnsEMBL::Compara::AlignUtils";
my $COMPARA = "Bio::EnsEMBL::Compara::ComparaUtils";

my $counts_genes_def = {
  'node_id' => 'int',

  human_gene => 'string',
  human_chr => 'string',
  human_start => 'int',
  human_end => 'int',
  human_strand => 'int',

  'CGH' => 'int',
  'C.G.H' => 'int',
  'H.CG' => 'int',
  'C.GH' => 'int',
  'G.CH' => 'int',
  'syn_CGH' => 'int',
  'syn_C.G.H' => 'int',
  'syn_H.CG' => 'int',
  'syn_C.GH' => 'int',
  'syn_G.CH' => 'int',

  most_frequent_pattern => 'string',
  most_frequent_pattern_excess => 'float',

  synon => 'int',
  nonsynon => 'int',
  constant => 'int',
  gap => 'int',
  
  tree_newick => 'string',
  tree_pattern => 'string',
  tree_length => 'float',
  tree_max_path => 'float',
  seq_length_mean => 'float',
  orig_leaf_count => 'int',
  };

my $counts_sites_def = {
  node_id => 'int',
  aln_position => 'int',
  chr_name => 'string',
  chr_start => 'int',
  chr_end => 'int',
  chr_strand => 'string',

  type => 'string',
  gap_type => 'string',

  pattern    => 'string',
  codon_a    => 'string',
  codon_b    => 'string',
  has_cpg    => 'string',

  codon_genome => 'string',
  };


sub fetch_input {
  my ($self) = @_;

  ### DEFAULT PARAMETERS ###
  my $params = {
    gorilla_count_species => '9593,9598,9606',
    gorilla_map_taxon => 9593,
    counts_sites_table => 'counts_sites',
    counts_genes_table => 'counts_genes'
  };
  ##########################

  # Fetch parameters from all possible locations.
  $self->load_all_params($params);

  # Create table if necessary.
  $self->create_table_from_params($self->compara_dba,$self->params->{counts_sites_table},$counts_sites_def);
  $self->create_table_from_params($self->compara_dba,$self->params->{counts_genes_table},$counts_genes_def);

}

sub run {
  my $self = shift;

  my $tree = $self->get_tree;

  print $tree->newick_format."\n";

  my $species_str = $self->param('gorilla_count_species');
  my @species_list = split(',',$species_str);

  my $taxon_to_letter = {
    9606 => 'H',
    9598 => 'C',
    9593 => 'G',
    9600 => 'O'
  };

  my @keeper_leaves = $TREE->get_leaves_for_species($tree,\@species_list);

  print map {$_->stable_id." "} @keeper_leaves;
  print "\n";

  # Test whether it's a good 1-1-1 orthology.
  my $is_good_tree = 1;
  my %keeper_hash;
  map {$keeper_hash{$_->taxon_id}=1} @keeper_leaves;
  map {$is_good_tree = 0 if (!defined $keeper_hash{$_})} @species_list;
  $is_good_tree = 0 if ($#keeper_leaves != $#species_list);
  if (!$is_good_tree) {
    print "Not one-to-one orthology: Doing nothing!\n";
    $self->autoflow_inputjob(0);
    return;
  }

  my @keeper_ids = map {$_->node_id} @keeper_leaves;
  #print "@keeper_ids\n";

  $tree = $TREE->extract_subtree_from_leaves($tree,\@keeper_ids);
  print $tree->newick_format()."\n";

  my $aln = $tree->get_SimpleAlign(-cdna => 1);
  $aln = $ALN->sort_by_tree($aln,$tree);

  $ALN->pretty_print($aln,{length=>200});
  
  my $gene_data = {
    'CGH' => 0,
    'C.G.H' => 0,
    'H.CG' => 0,
    'C.GH' => 0,
    'G.CH' => 0,
    'syn_CGH' => 0,
    'syn_C.G.H' => 0,
    'syn_H.CG' => 0,
    'syn_C.GH' => 0,
    'syn_G.CH' => 0,
    synon => 0,
    nonsynon => 0,
    constant => 0,
    gap => 0,
    most_frequent_pattern => '',
    most_frequent_pattern_excess => 0,
    tree_newick => $tree->newick_format(),
    tree_length => $self->tree_length($tree),
    tree_max_path => $self->max_path($tree),
    tree_pattern => $self->get_tree_pattern($tree,$taxon_to_letter),
    seq_length_mean => $self->seq_length_mean($tree)
  };

  $gene_data->{orig_leaf_count} = $self->root_node_gene_count($tree);
  # Collect human protein.
  my @human_proteins = grep { $_->taxon_id == 9606 } $tree->leaves;
  my @human_genes    = map  { $_->get_Gene } @human_proteins;
  if ( scalar @human_proteins > 0 ) {
    my $member = $human_proteins[0];
    $gene_data->{'human_gene'}    = $member->get_Gene->stable_id;
  }
  # Collect protein coords.
  if ( scalar @human_proteins > 0) {
    my $member = $human_proteins[0];
    my $tscr_orig = $member->get_Transcript;
    my $tscr      = $tscr_orig->transform("chromosome");
    if ( defined $tscr ) {
      my $chr    = "chr" . $tscr->slice->seq_region_name;
      my $strand = $tscr->strand;
      my $start  = $tscr->coding_region_start;
      my $end    = $tscr->coding_region_end;
      $gene_data->{human_chr}    = $chr;
      $gene_data->{human_start}  = $start;
      $gene_data->{human_end}    = $end;
      $gene_data->{human_strand} = $strand;
    }
  }


  for (my $i=1; $i < $aln->length; $i+= 3) {
    my $slice = $aln->slice($i,$i+2);

    # Skip all-gap alignment slices.
    next if (scalar $slice->each_seq == 0);
    
    my $codon_hashref;
    my $aa_hashref;
    my $taxon_id_hash;
    foreach my $member ($tree->leaves) {
      #foreach my $seq ($slice->each_seq) {
      #my $member = $tree->find_leaf_by_name($seq->id);
      $taxon_id_hash->{$member->taxon_id} = 1;
      my ($seq) = $slice->each_seq_with_id($member->stable_id);
      my $seq_aa;
      my $seq_codon;
      if (!defined $seq) {
	#print "gap!\n";
	$seq_codon = "---";
	$seq_aa = "-";
      } else {
	$seq_codon = $seq->seq;
	$seq_aa = $seq->translate->seq;
      }
      $codon_hashref->{$seq_codon} = [] if (!defined $codon_hashref->{$seq_codon});
      $aa_hashref->{$seq_aa}++;
      push @{$codon_hashref->{$seq_codon}},$member->taxon_id;
    }
    
    my $type = 'NULL';

    $type = 'synonymous' if (scalar keys %$aa_hashref == 1);
    $type = 'nonsynonymous' if (scalar keys %$aa_hashref > 1);
    $type = 'constant' if (scalar keys %$codon_hashref == 1);
    $type = 'gap' if (defined $aa_hashref->{'-'});

    my @codon_string_arr;
    my $codon_value_hash;
    my $j=0;
    my $has_cpg = 0;
    foreach my $key (keys %$codon_hashref) {
      $j++;

      $has_cpg = 1 if ($key =~ m/cg/i);

      my @taxon_ids = @{$codon_hashref->{$key}};
      my @chars = map {$taxon_to_letter->{$_}} @taxon_ids;
      my $taxa_string = join("",sort @chars);
      $codon_value_hash->{$taxa_string} = $key;
      push @codon_string_arr,$taxa_string;
    }
    my @codon_arr = sort {length $a <=> length $b || $a cmp $b} @codon_string_arr;
    my $codon_string = join(".",@codon_arr);

    my $codon_a = $codon_value_hash->{$codon_arr[0]};
    my $codon_b = '';
    my $codon_b = $codon_value_hash->{$codon_arr[1]} if (scalar @codon_arr > 1);

    my @species_letters = map {$taxon_to_letter->{$_}} keys %$taxon_id_hash;
    my $gap_type = join("", map {$_ if ($codon_value_hash->{$_} eq '---');} @species_letters);

    my $coord_data = $self->get_genomic_coord($tree,$aln,$i,$self->param('gorilla_map_taxon'));
    #die("No coords!") unless ($coord_data);

    $gene_data->{$codon_string}++ if ($type ne 'gap');
    $gene_data->{'syn_'.$codon_string}++ if ($type ne 'gap' && $type ne 'nonsynonymous');
    $gene_data->{nonsynon}++ if ($type eq 'nonsynonymous');
    $gene_data->{synon}++ if ($type eq 'synonymous');
    $gene_data->{constant}++ if ($type eq 'constant');
    $gene_data->{gap}++ if ($type eq 'gap');

    my $codon_genome = $coord_data->{char};
    if ($codon_genome ne '') {
      if (! grep {$codon_genome eq $_} values %$codon_value_hash) {
	print ("Genomic codon $codon_genome not found in codon hash!\n");
      }
    }

    my $site_data = {
      aln_position => $i, # Give alignment coordinates with the gaps NOT removed.
      parameter_set_id => $self->parameter_set_id,

      codon_a => $codon_a,
      codon_b => $codon_b,
      has_cpg => $has_cpg,

      codon_genome => $codon_genome,
      
      type => $type,
      gap_type => $gap_type,
      pattern => $codon_string
    };
    
    $site_data = $self->replace_params($self->params,$site_data,$coord_data);

    #$COMPARA->hash_print($site_data);    
    $self->store_params_in_table($self->compara_dba,$self->param('counts_sites_table'),$site_data);
    #$ALN->pretty_print($slice,{length=>200});
  }

  # Get the most frequent pattern and excess.
  my $most_frequent_pattern = '';
  my $most_frequent_count = 0;
  foreach my $key ('H.CG','C.GH','G.CH') {
    if ($gene_data->{$key} >= $most_frequent_count) {
      $most_frequent_count = $gene_data->{$key};
      $most_frequent_pattern = $key;
    }
  }
  $gene_data->{most_frequent_pattern} = $most_frequent_pattern;
  my $other_count = 0;
  foreach my $key ('H.CG','C.GH','G.CH') {
    $other_count = $gene_data->{$key} if ($key ne $most_frequent_pattern && $gene_data->{$key} > $other_count);
  }
  $gene_data->{most_frequent_pattern_excess} = ($most_frequent_count+1) / ($other_count+1);

  # Apply the gene-wide data and store a row in the gene-count table.
  $gene_data = $self->replace_params($self->params,$gene_data);
  $self->store_params_in_table($self->compara_dba,$self->param('counts_genes_table'),$gene_data);
}

sub get_tree_pattern {
  my $self = shift;
  my $tree = shift;
  my $taxon_to_letter = shift;

  $tree = $tree->copy;
  foreach my $leaf ($tree->leaves) {
    $leaf->stable_id($leaf->taxon_id);
  }
  my $taxon_tree = $tree->newick_format('no_bl');

  foreach my $leaf ($tree->leaves) {
    my $taxon_id = $leaf->taxon_id;
    my $char = $taxon_to_letter->{$leaf->taxon_id};
    $taxon_tree =~ s/$taxon_id/$char/ig;
  }
  $taxon_tree =~ s/[\(\)]/\./ig;
  $taxon_tree =~ s/[,;]//g;
  $taxon_tree =~ s/\.+/\./g;
  $taxon_tree =~ s/^\.+//g;
  $taxon_tree =~ s/\.+$//g;

  # Split, re-sort and re-join the string.
  my @toks = split("\\.",$taxon_tree);
  @toks = map {join "", sort(split("",$_))} @toks;
  my $final_pattern = join(".",sort {length $a <=> length $b || $a cmp $b} @toks);
  print "$final_pattern\n";
  return $final_pattern;
}

sub get_genomic_coord {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $aln_position = shift;
  my $taxon_id = shift;

  foreach my $leaf ($tree->leaves) {
    next unless ($leaf->taxon_id == $taxon_id);
    
    my ($seq) = $aln->each_seq_with_id( $leaf->stable_id );
    my $seq_str = $seq->seq;
    
    my $tscr = $leaf->get_Transcript;
    #print STDERR " ->" . $tscr->stable_id . "\n";
    $tscr = $tscr->transform("chromosome");
    next unless ( defined $tscr );
    my $chr = "chr" . $tscr->slice->seq_region_name;
    
    my $char = substr($seq_str,$aln_position-1,3);
    my $loc = $seq->location_from_column($aln_position);
    next if ( !defined $loc || $loc->location_type() eq 'IN-BETWEEN' );
    #print "HEY!" unless ($loc->start == $aln_position);
    my $cdna_start = $loc->start - 1;
    my $cdna_end = $loc->end - 1;
    my @gc_arr = $tscr->cdna2genomic( $cdna_start + $tscr->cdna_coding_start, $cdna_end + $tscr->cdna_coding_start+2 );

    #print "Start: ".$tscr->start."\n";

    my $gc = $gc_arr[0];
    next unless ( $gc && $gc->isa("Bio::EnsEMBL::Mapper::Coordinate") );

    my $strand = "+";
    $strand = "-" if ( $gc->strand == -1 );

    my $start = $gc->start;
    my $end   = $gc->end;

    #print " --> $start\n";
    my $obj = {
      chr_name          => $chr,
      chr_start        => $start,
      chr_end          => $end,
      aln_position => $aln_position,
      char         => $char,
      stable_id    => $leaf->stable_id,
      node_id      => $self->param('node_id'),
      member_id    => $leaf->dbID,
      chr_strand       => $strand
      };
    return $obj;    
  }
  return undef;
}

sub write_output {

}


1
