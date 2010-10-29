package Bio::Greg::Hive::CountSubstitutions;

use strict;
use Bio::Greg::StatsCollectionUtils;

use base ( 'Bio::Greg::Hive::Process', 'Bio::Greg::StatsCollectionUtils' );

my $TREE    = "Bio::EnsEMBL::Compara::TreeUtils";
my $ALN     = "Bio::EnsEMBL::Compara::AlignUtils";
my $COMPARA = "Bio::EnsEMBL::Compara::ComparaUtils";

sub codon_subs_table_def {
  my $self = shift;
  return {
    parameter_set_id => 'int',
    node_id          => 'int',    # Node ID of the Compara tree.

    codeml_node_number => 'int',
    codeml_node_id => 'char32',   # Node ID retrieved from Codeml.
    leaves_beneath => 'string',   # Sorted list of all leaves beneath this node.
    member_id      => 'int',      # Member ID of the sequence (if it's a terminal branch).
    transcript_id  => 'char32',
    taxon_id       => 'int',
    chr            => 'char8',    # Genomic coordinates.
    chr_strand     => 'int',
    chr_start      => 'int',
    seq_pos        => 'int',      # Position on the sequence in AA coordinates.
    aln_pos        => 'int',
    cdna_pos       => 'int',

    codon_from_context => 'char8',
    codon_from         => 'char4',
    codon_to           => 'char4',
    codon_genomic      => 'char4',
    aa_from            => 'char1',
    aa_to              => 'char1',
    confidence         => 'float',    # Codeml's confidence score.

    codon_cpg => 'tinyint',           # 1 if 'codon_from' overlapped *any* CpG sites.

    # These fields describe the single nucleotide mutation.
    nuc_from  => 'char1',
    nuc_to    => 'char1',
    mut_count => 'tinyint',           # Number of mutations within the codon.
    mut_pos   => 'tinyint',           # position in the codon: 1, 2, or 3
    mut_cpg   => 'tinyint',           # 0 = no cpg overlap; 1 = C-position; 2 = G-position.
    mut_rev_cpg =>
      'tinyint',    # Same as above, but gives the mutation's cpg position on the opposite strand.
    mut_syn  => 'tinyint',    # synonymous
    mut_nsyn => 'tinyint',    # nonsynonymous
    mut_sw   => 'tinyint',    # strong-to-weak
    mut_ws   => 'tinyint',    # weak-to-strong
    mut_ts   => 'tinyint',    # transition
    mut_tv   => 'tinyint',    # transversion

    unique_keys => 'data_id,codeml_node_number,aln_pos',
    extra_keys => 'codeml_node_id,taxon_id'
  };
}

sub store_substitution_in_table {
  my $self    = shift;
  my $table   = shift;
  my $sub_obj = shift;

}

sub extract_substitution_info {
  my $self     = shift;
  my $node = shift;
  my $sub_obj  = shift;
  my $tree     = shift;
  my $aln      = shift;
  my $cdna_aln = shift;

  # sub_obj has the following: pos, codon_a, codon_b, aa_a, aa_b, confidence
  my $id = $node->id;
  my $sub_obj_id           = $sub_obj->{id};
  my $aln_pos      = $sub_obj->{pos};
  my $cdna_aln_pos = ( $aln_pos - 1 ) * 3 + 1;
  my $codon_a      = $sub_obj->{codon_a};
  my $codon_b      = $sub_obj->{codon_b};
  my $aa_a         = $sub_obj->{aa_a};
  my $aa_b         = $sub_obj->{aa_b};
  my $confidence   = $sub_obj->{confidence};

  next if ( $aa_a eq '*' || $aa_b eq '*' );
  next if ( $codon_a =~ m/[nx]/i || $codon_b =~ m/[nx]/i );

  my $nuc_left        = '';
  my $nuc_right       = '';
  my $codon_context_a = '';
  my $codon_context_b = '';

  my $final_params = {};

  # Get the list of leaves beneath this node.
  my $leaves_beneath = '';
  if ($node->is_Leaf) {
    $leaves_beneath = $node->id;
  } else {
    my @leaves;
    foreach my $sub_node ($node->get_all_Descendents) {
      push @leaves,$sub_node if ($sub_node->is_Leaf);
    }
    @leaves = sort {$a->id cmp $b->id} @leaves;
    $leaves_beneath = join(',',map {$_->id} @leaves);
  }

  # Only extract sequence-based info if we're a leaf node.
  if ($node->is_Leaf) {    
    
    # Use location_from_column to get the sequence position.
    my $seq = $aln->get_seq_by_id($id);
    $self->throw("No seq for id[$id]!") unless ($seq);
    my $seq_str        = $seq->seq;
    my $seq_str_nogaps = $seq_str;
    $seq_str_nogaps =~ s/-//g;    # Remove gaps
    my $location = $seq->location_from_column($aln_pos);
    $self->throw("Undefined seq position!") unless ($location);
    $self->throw("In between!") if ( $location->location_type eq 'IN-BETWEEN' );
    my $seq_pos = $location->start;
    
    my $aa_b_from_aln = substr( $seq_str_nogaps, $seq_pos - 1, 1 );
    $self->throw( sprintf( "AA doesn't match: expected[%s] got[%s]", $aa_b, $aa_b_from_aln ) )
      unless ( $aa_b eq $aa_b_from_aln );
    
    # Dummy-check that the Codeml codon matches what we get from the CDNA align.
    my $cdna_seq = $cdna_aln->get_seq_by_id($id);
    my $cdna_str = $cdna_seq->seq;
    $cdna_str =~ s/-//g;    # Remove gaps.
    $self->throw("No CDNA seq!") unless ($cdna_seq);
    
    # 1-based amino acid to 1-based dna: 1 -> 1, 2 -> 4, etc.
    my $cdna_pos = ( $seq_pos - 1 ) * 3 + 1;
    my $codon_b_from_aln = substr( $cdna_str, $cdna_pos - 1, 3 );
    if ( $codon_b ne $codon_b_from_aln ) {
      print
        "codon_b[$codon_b] from_aln[$codon_b_from_aln] aln_pos[$aln_pos] cdna_aln_pos[$cdna_aln_pos] seq_pos[$seq_pos] cdna_seq_pos[$cdna_pos] \n";
      warn("Codon doesn't match!");
    }
    
    if ( $cdna_pos > 2 && $cdna_pos < length($cdna_str) - 5 ) {
      $nuc_left  = substr( $cdna_str, $cdna_pos - 1 - 1, 1 );
      $nuc_right = substr( $cdna_str, $cdna_pos - 1 + 3, 1 );
      $codon_context_a = $nuc_left . $codon_a . $nuc_right;
      $codon_context_b = $nuc_left . $codon_b . $nuc_right;
    }
    
    my ($member) = grep { $_->name eq $id } $tree->nodes;
    
    my $taxon_id = $member->taxon_id;
    my $member_id;
    my $transcript_id;
    my $name = $member->taxon->short_name;
    my $cdna_start;
    my $chr;
    my $chr_codon_start;
    my $chr_codon_strand;
    my $genomic_codon = '';
    my $genomic_context = '';
    
    my $alias = $member->taxon->ensembl_alias;
    my $slice_a = Bio::EnsEMBL::Registry->get_adaptor( $alias, 'core', 'slice' );
    
    if ( $self->param('aln_type') =~ m/compara/i ) {
      
      # Get genomic coordinates from the Compara member.
      
      my ($member) = grep { $_->name eq $id } $tree->nodes;
      $member_id = $member->member_id;
      my $tx = $member->get_Transcript;
      $transcript_id = $tx->stable_id;
      my $slice = $tx->slice;
      my $db    = $slice->adaptor->db;
      my @css   = @{ $db->get_CoordSystemAdaptor->fetch_all };
            
      my $name = $member->taxon->short_name;
      $chr = $tx->slice->seq_region_name;
      
      foreach my $i ( -2, -1, 0, 1, 2, 3, 4) {
        my $cdna_start        = $cdna_pos + $i - 1;
        my @genomic_coord_arr = $tx->cdna2genomic( $cdna_start + $tx->cdna_coding_start,
                                                   $cdna_start + $tx->cdna_coding_start );
        @genomic_coord_arr = grep { $_->isa("Bio::EnsEMBL::Mapper::Coordinate") } @genomic_coord_arr;
        if ( scalar(@genomic_coord_arr) == 0 ) {
          warn("No coords! $name $chr $cdna_start\n");
        }
        if ( scalar(@genomic_coord_arr) > 1 ) {
          warn("Multiple coords! $name $chr $cdna_start\n");
        }
        foreach my $gc (@genomic_coord_arr) {
          next unless ( $gc->isa("Bio::EnsEMBL::Mapper::Coordinate") );
          my $chr_strand = $gc->strand;
          my $chr_start  = $gc->start;
          my $slice = $slice_a->fetch_by_region( undef, $chr, $chr_start, $chr_start, $chr_strand );

          $genomic_codon .= $slice->seq if ($i >= 0 && $i <= 2);

          $genomic_context .= $slice->seq;
          $chr_codon_start  = $chr_start  if ( $i == 0 );
          $chr_codon_strand = $chr_strand if ( $i == 0 );          
        }
      }
      
      return undef unless ( defined $chr_codon_start );
    } else {
      
      # Get the genomic coordinates that were stored while we were collecting the genomic alignment.
      my $gc = $cdna_aln->annotation->{_genomic_coords};
      
      foreach my $i ( 0, 1, 2) {
        my $index      = $cdna_aln_pos + $i;
        my $key        = $member->taxon->short_name . '_' . $index;
        my $coord      = $gc->{$key};
        my $chr_start  = $coord->start;
        my $chr_strand = $coord->strand;
        $chr_codon_start  = $chr_start  if ( $i == 0 );
        $chr_codon_strand = $chr_strand if ( $i == 0 );
        
        $chr = $coord->id;
        my $chr_cs = $coord->coord_system;
        my $slice = $slice_a->fetch_by_region( undef, $chr, $chr_start, $chr_start, $chr_strand );
        $genomic_codon .= $slice->seq;
      }
    }
    
    if ( $genomic_codon ne $codon_b ) {
      #Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $aln,      { full => 1,width=>150 } );
      #Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $cdna_aln, { full => 1,width=>150 } );
      $self->hash_print($sub_obj);
      
      #print "aln :". $cdna_seq->seq."\n";
      #print "str :". $cdna_str."\n";
      #print "cds :". $member->sequence_cds."\n";
      print "gen_context $genomic_context\n";
      print
        "$name aln_pos[$cdna_aln_pos] chr[$chr $chr_codon_start $chr_codon_strand] aln[$codon_b] genomic[$genomic_codon]\n";

      warn("Genomic codon doesn't match alignment!");

      if ($genomic_context =~ m/($codon_b)/ig) {
        # Continue on, the nearby context has the codon so we're not too far off.
        my $index = pos($genomic_context) - length($codon_b);
        print "Index: $index\n";
        my $shift = ($index-2);
        $chr_codon_start += $shift;
        $genomic_codon = substr($genomic_context,$index,3);
        print "Fixed it: $codon_b $genomic_codon\n";
        if ($codon_b ne $genomic_codon) {
          die("Even after fixing, genomic codon doesn't match alignment for $name cdna_aln_pos[$cdna_aln_pos] aln_codon[$codon_b] genomic_codon[$genomic_context]!");
        }
      } else {
        die("Genomic codon doesn't match alignment for $name cdna_aln_pos[$cdna_aln_pos] aln_codon[$codon_b] genomic_codon[$genomic_codon]!");
      }
    }
    

    my $leaf_params = {
      member_id     => $member_id,
      transcript_id => $transcript_id,
      taxon_id      => $taxon_id,
      seq_pos    => $seq_pos,
      cdna_pos   => $cdna_pos,
      chr        => $chr,
      chr_strand => $chr_codon_strand,
      chr_start  => $chr_codon_start,
      codon_from_context => $codon_context_a,
      codon_genomic      => $genomic_codon,
    };
    
    $final_params = $self->replace($final_params,$leaf_params);
  }

  # Synon/nonsynon.
  my $mut_syn  = 0;
  my $mut_nsyn = 0;
  $mut_syn = 1 if ( $aa_a eq $aa_b );
  $mut_nsyn = 1 if ( $aa_a ne $aa_b );
  
  # Isolate the nucleotide difference.
  my $mut_count = 0;
  my $mut_pos   = -1;
  my $nuc_a;
  my $nuc_b;
  foreach my $i ( 1, 2, 3 ) {
    my $tmp_nuc_a = substr( $codon_a, $i - 1, 1 );
    my $tmp_nuc_b = substr( $codon_b, $i - 1, 1 );
    if ( $tmp_nuc_a ne $tmp_nuc_b ) {
      $mut_pos = $i;
      $mut_count++;
      $nuc_a = $tmp_nuc_a;
      $nuc_b = $tmp_nuc_b;
    }
  }
  my $both_nucs = $nuc_a . $nuc_b;
  if ( $mut_count > 1 ) {

    # If we have multiple substitutions, we'll only use the latter substitution
    # to collect the rest of the stats.
  }

  # Weak-strong, transition-transversion.
  my $mut_ws = 0;
  my $mut_sw = 0;
  my $mut_ts = 0;
  my $mut_tv = 0;
  $mut_ws = 1 if ( $nuc_a =~ m/[at]/i && $nuc_b =~ m/[gc]/i );
  $mut_sw = 1 if ( $nuc_a =~ m/[gc]/i && $nuc_b =~ m/[at]/i );
  $mut_ts = 1 if ( $both_nucs =~ m/(ag|ga|ct|tc)/i );
  $mut_tv = 1 if ( $both_nucs =~ m/(ac|ca|at|ta|gc|cg|gt|tg)/i );

  # CpG overlaps.
  my $mut_cpg = 0;
  $mut_cpg = 1
    if ( ( $mut_pos == 1 && $codon_a =~ m/cg./i )
    || ( $mut_pos == 2 && $codon_a =~ m/.cg/i ) );
  $mut_cpg = 2
    if ( ( $mut_pos == 2 && $codon_a =~ m/cg./i )
    || ( $mut_pos == 3 && $codon_a =~ m/.cg/i ) );
  if ( $nuc_left ne '' ) {

    # Look for CpG overlaps on the boundary, too.
    $mut_cpg = 2 if ( ( $mut_pos == 1 && $codon_context_a =~ m/cg.../i ) );
  }
  if ( $nuc_right ne '' ) {
    $mut_cpg = 1 if ( ( $mut_pos == 3 && $codon_context_a =~ m/...cg/i ) );
  }
  my $codon_cpg = 0;
  if ( $codon_context_a ne '' ) {
    $codon_cpg = 1 if ( $codon_context_a =~ m/cg/i );
  } else {
    $codon_cpg = 1 if ( $codon_a =~ m/cg/i );
  }

  # Look on opposite strand to assign mut_rev_cpg and codon_rev_cpg
  my $mut_rev_cpg = 0;

  $mut_rev_cpg = 2
    if ( ( $mut_pos == 1 && $codon_a =~ m/cg./i )
    || ( $mut_pos == 2 && $codon_a =~ m/.cg/i ) );
  $mut_rev_cpg = 1
    if ( ( $mut_pos == 2 && $codon_a =~ m/cg./i )
    || ( $mut_pos == 3 && $codon_a =~ m/.cg/i ) );

  if ( $nuc_right ne '' ) {
    $mut_rev_cpg = 2 if ( $mut_pos == 3 && $codon_context_a =~ m/...cg/i );
  }

  if ( $nuc_left ne '' ) {
    $mut_rev_cpg = 1 if ( $mut_pos == 1 && $codon_context_a =~ m/cg.../i );
  }

  my $base_params = {
    codeml_node_id => $id,
    codeml_node_number => $sub_obj_id,
    leaves_beneath => $leaves_beneath,
    aln_pos    => $aln_pos,
    
    codon_from         => $sub_obj->{codon_a},
    codon_to           => $sub_obj->{codon_b},
    aa_from            => $sub_obj->{aa_a},
    aa_to              => $sub_obj->{aa_b},
    confidence         => $sub_obj->{confidence},

    codon_cpg => $codon_cpg,

    nuc_from    => $nuc_a,
    nuc_to      => $nuc_b,
    mut_count   => $mut_count,
    mut_pos     => $mut_pos,
    mut_ws      => $mut_ws,
    mut_sw      => $mut_sw,
    mut_ts      => $mut_ts,
    mut_tv      => $mut_tv,
    mut_cpg     => $mut_cpg,
    mut_rev_cpg => $mut_rev_cpg,
    mut_syn     => $mut_syn,
    mut_nsyn    => $mut_nsyn,
  };

  $final_params = $self->replace($final_params,$base_params);
  return $final_params;
}

1;
