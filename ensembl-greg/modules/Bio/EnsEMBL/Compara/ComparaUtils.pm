package Bio::EnsEMBL::Compara::ComparaUtils;

use strict;

use Bio::TreeIO;
use Bio::EnsEMBL::Compara::LocalMember;
use Bio::EnsEMBL::Compara::ProteinTree;
use Bio::EnsEMBL::Compara::NestedSet;
use Bio::EnsEMBL::Utils::Exception;

use Bio::EnsEMBL::Compara::AlignUtils;
use Bio::EnsEMBL::Compara::TreeUtils;

use Bio::EnsEMBL::Registry;

use Bio::Greg::IndexedFasta;
use Bio::DB::Fasta;
use Bio::DB::Qual;

use File::Path;

#
# A grab bag of useful methods for Compara modules.
#

my $TREE    = "Bio::EnsEMBL::Compara::TreeUtils";
my $ALN     = "Bio::EnsEMBL::Compara::AlignUtils";
my $COMPARA = "Bio::EnsEMBL::Compara::ComparaUtils";

*Bio::EnsEMBL::Compara::NestedSet::leaves = sub {
  my $self = shift;
  return @{$self->get_all_leaves};
};
*Bio::EnsEMBL::Compara::NestedSet::nodes = sub {
  my $self = shift;
  return @{$self->get_all_nodes};
};
*Bio::EnsEMBL::Compara::NestedSet::sorted_children = sub {
  my $self = shift;

  my @sortedkids = 
     sort { scalar($b->leaves) <=> scalar($a->leaves)
              ||
            $a->is_leaf <=> $b->is_leaf
                     ||
            $a->get_child_count <=> $b->get_child_count
                     ||
            $a->distance_to_parent <=> $b->distance_to_parent
          }  @{$self->children;};
  return \@sortedkids;
};
*Bio::EnsEMBL::Compara::NestedSet::get_all_leaves = sub {
  my $self = shift;
  
  my $leaves = {};
  $self->_recursive_get_all_leaves($leaves);
  my @leaf_list = values %{$leaves};
#  my @leaf_list = sort {$a->node_id <=> $b->node_id} values(%{$leaves});
  return \@leaf_list;
};
*Bio::EnsEMBL::Compara::NestedSet::ascii = sub {
  my $self = shift;
  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($self);
  return $treeI->ascii(@_);
};
*Bio::EnsEMBL::Compara::NestedSet::enclosed_leaves_string = sub {
  my $self = shift;

  if ($self->is_leaf) {
    return $self->id;
  }

  my @leaves_beneath = map {$_->name} $self->leaves;
  my $leaf_string = join("|", sort {$a cmp $b} @leaves_beneath);
  return $leaf_string;
};

*Bio::EnsEMBL::Compara::NestedSet::taxon_id = sub {
  my $self = shift;
  my $taxon_id = shift;

  if (defined $taxon_id) {
    $self->set_tagvalue('taxon_id', $taxon_id);
  }
  return $self->get_tagvalue('taxon_id') || undef;
};

*Bio::EnsEMBL::Compara::NestedSet::branch_length = sub {
  my $self = shift;
  my $length = shift;

  if (defined $length) {
    $self->distance_to_parent($length);
  }
  return $self->distance_to_parent;
};
*Bio::EnsEMBL::Compara::NestedSet::id = sub {
  my $self = shift;
  return $self->name;
};

*Bio::EnsEMBL::Compara::Member::name = sub {
  my $self = shift;
  my $value = shift;
  if(defined($value)) { $self->add_tag('name', $value); }
  else { $value = $self->get_tagvalue('name'); }
  if (!defined $value) {
    $self->add_tag('name', $self->stable_id);
    $value = $self->stable_id;
  }
  return $value;
};

*Bio::EnsEMBL::Compara::NestedSet::_internal_nhx_format = sub {
  my $self = shift;
  my $format_mode = shift;
  my $nhx = "";

  if($self->get_child_count() > 0) {
    $nhx .= "(";
    my $first_child=1;
    foreach my $child (@{$self->sorted_children}) {  
      $nhx .= "," unless($first_child);
      $nhx .= $child->_internal_nhx_format($format_mode);
      $first_child = 0;
    }
    $nhx .= ")";
  }
  
  $nhx .= $self->name;
  $nhx .= sprintf(":%1.4f", $self->distance_to_parent);

  $nhx .= "[&&NHX";
  my $tags = $self->get_tagvalue_hash;
  foreach my $key (keys %$tags) {
    next if ($key eq 'name');
    $nhx .= ":$key=".$tags->{$key};
  }
  $nhx .= "]";
  return $nhx;
};

*Bio::SimpleAlign::tx = sub {
  my $self = shift;
  return Bio::EnsEMBL::Compara::AlignUtils->translate($self, {});
};


sub load_registry {
  my $class = shift;
  if ( $ENV{'USER'} =~ /gj1/ ) {
    Bio::EnsEMBL::Registry->load_registry_from_multiple_dbs(
      {
        -host => 'ens-livemirror',
        -user => 'ensadmin',
        -pass => 'ensembl',
        #-verbose => 1,
      },
#      {
#        -host => 'ensdb-archive',
#        -port => 5304,
#        -user => 'ensro',
#        -user => 'ensadmin',
#        -pass => 'ensembl',
#        #-verbose => 1
#      }
      );
    Bio::EnsEMBL::Registry->set_disconnect_when_inactive(1);
    my $compara_dba = Bio::EnsEMBL::Registry->get_DBAdaptor( 'multi', 'compara' );
    if ($compara_dba) {
      printf " >> Using Compara DB at [%s/%s]\n",$compara_dba->dbc->host,$compara_dba->dbc->dbname;
    }
  } else {
    #Bio::EnsEMBL::Registry->no_version_check(1);
  }
}

sub get_one_to_one_ortholog {
  my $class          = shift;
  my $compara_dba    = shift;
  my $member         = shift;
  my $other_taxon_id = shift;

  my $ma   = $compara_dba->get_MemberAdaptor;
  my $gdba = $compara_dba->get_GenomeDBAdaptor;
  my $ha   = $compara_dba->get_HomologyAdaptor;

  my $other_gdb   = $gdba->fetch_by_taxon_id($other_taxon_id);
  my $gene_member = $member->gene_member;

  # Need to search with a gene member.
  my @homologies = @{ $ha->fetch_all_by_Member_paired_species( $gene_member, $other_gdb->name ) };

  foreach my $homology (@homologies) {
    my $desc = $homology->description;
    if ( $desc eq 'ortholog_one2one' ) {
      my @genes = @{ $homology->gene_list };
      my ($other_member) = grep { $_->taxon_id == $other_taxon_id } @genes;

      # Want to return the canonical peptide.
      $other_member = $other_member->get_canonical_peptide_Member;
      return ( $homology, $other_member ) if ( defined $other_member );
    }
  }
  return undef;
}

sub protein_tree_taxonomy_distance {
  my $class     = shift;
  my $orig_tree = shift;

  my $tree = $orig_tree->copy;

  # First, we need to create a newick string with the protein tree IDs normalized into taxon_id_[n]
  my %taxon_hash;
  foreach my $member ( @{ $tree->get_all_leaves } ) {
    my $taxon_id = $member->taxon_id;
    $taxon_hash{$taxon_id} = 0 if ( !defined $taxon_hash{$taxon_id} );
    $taxon_hash{$taxon_id} = $taxon_hash{$taxon_id} + 1;
    $member->stable_id( "'" . $taxon_id . "X" . $taxon_hash{$taxon_id} . "'" );
  }

# Now, take the NCBI taxonomy tree and normalize in a similar way. Insert sister nodes where we have more than one
# sequence per species in the above tree.
  my $ncbi_tree = $class->get_genome_taxonomy_below_level( $orig_tree->adaptor->db );

  # Remove taxa where we have no genes.
  my @remove_me = ();
  foreach my $member ( @{ $ncbi_tree->get_all_leaves } ) {
    my $taxon_id = $member->taxon_id;
    push @remove_me, $member if ( !defined $taxon_hash{$taxon_id} );
  }
  foreach my $node (@remove_me) {
    $TREE->delete_lineage( $tree, $node );
    $tree->minimize_tree;
  }
  $tree->minimize_tree;

  # Duplicate nodes where we have more than one gene per taxon.
  foreach my $member ( @{ $ncbi_tree->get_all_leaves } ) {
    my $taxon_id             = $member->taxon_id;
    my $gene_count_for_taxon = $taxon_hash{$taxon_id};
    $member->name( "'" . $taxon_id . "X1'" );
    my $dist = $member->distance_to_parent;
    for ( my $i = 2 ; $i <= $gene_count_for_taxon ; $i++ ) {
      my $new_leaf = new Bio::EnsEMBL::Compara::NCBITaxon;
      $member->parent->add_child($new_leaf);
      $new_leaf->name( "'" . $taxon_id . "X" . $i . "'" );

      #$new_leaf->distance_to_parent($dist);
      $new_leaf->distance_to_parent(0);
    }
  }

  # Format names a bit.
  foreach my $member ( @{ $ncbi_tree->get_all_nodes } ) {
    if ( !$member->is_leaf ) {
      $member->name('');
    }
  }
  $ncbi_tree = $TREE->remove_elbows($ncbi_tree);

  #print $tree->newick_format . "\n";
  #print $ncbi_tree->newick_format . "\n";
  $TREE->robinson_foulds_dist( $tree, $ncbi_tree );
  $TREE->k_tree_dist( $tree, $ncbi_tree );
}

sub cigar_line {
  my ( $class, $str ) = @_;
  my $cigar_line = "";

  $_ = $str;
  my @groups = /(-+|[^-]+)/ig;
  foreach my $token (@groups) {
    my $len  = length($token);
    my $char = 'M';
    if ( index( $token, '-' ) > -1 ) {
      $char = 'D';
    }
    if ( $len > 1 ) {
      $cigar_line .= $len;
    }
    $cigar_line .= $char;
  }
  return $cigar_line;
}

# Stores a given SimpleAlign (and optional cdna align) into $output_table.
sub store_SimpleAlign_into_table {
  my $class        = shift;
  my $output_table = shift;
  my $tree         = shift;
  my $aa_aln       = shift;
  my $cdna_aln     = shift;

  my $pta = $tree->adaptor;
  my $mba = $pta->db->get_MemberAdaptor;

  #$pta->protein_tree_member($output_table);

  # Make sure everything's the same length.
  my $length;
  foreach my $seq ($aa_aln->each_seq) {
    $length = $seq->length unless (defined $length);
    die("Aa_aln not all the same length!") unless ($length == $seq->length);
  }

  #
  # Align cigar_lines to members and store
  #
  foreach my $node ( @{ $tree->get_all_leaves } ) {

    #print $node." ". $node->node_id." ".$node->stable_id."\n";
    my $name = $node->stable_id;

    if ( defined $aa_aln ) {

      print "membr: ".$node->member_id."\n";
      print "seq id: ".$node->sequence_id."\n";
      $pta->dbc->do(
        "UPDATE member set sequence_id=NULL where member_id=" . $node->member_id . ";" );

      if ($node->sequence_id) {
        $pta->dbc->do("DELETE FROM sequence WHERE sequence_id=".$node->sequence_id.";");
      }
      
      my $seq_obj = $ALN->get_seq_with_id( $aa_aln, $name );
      if ( !defined $seq_obj ) {
        throw("No sequence with ID $name found in the alignment!");
      }
      my $aa_seq     = $seq_obj->seq;
      my $cigar_line = $class->cigar_line($aa_seq);
      $node->cigar_line($cigar_line);
      $aa_seq =~ s/-//g;
      $node->sequence($aa_seq);
    }

    if ( defined $cdna_aln ) {
      my $cdna_seq = $ALN->get_seq_with_id( $cdna_aln, $name )->seq;
      $cdna_seq =~ s/-//g;
      $node->sequence_cds($cdna_seq);
    }
    
    $node->sequence_id(0);
    
    # Call the ProteinTreeAdaptor method to do the actual database dirty work.
    $mba->store($node);
    $pta->store($node);
  }
}

# Reads the MCoffee scores from $mcoffee_scores file and stores them into $output_table.
sub store_MCoffee_scores_into_table {
  my $class          = shift;
  my $tree           = shift;
  my $mcoffee_scores = shift;    # File where the mcoffee scores are held.
  my $output_table   = shift;

  my $pta = $tree->adaptor;

  #
  # Read in the scores file manually.
  #
  my %score_hash;
  my %overall_score_hash;
  my $FH = IO::File->new();
  $FH->open($mcoffee_scores) || throw("Could not open alignment scores file [$mcoffee_scores]");
  <$FH>;                         #skip header
  my $i = 0;
  while (<$FH>) {
    $i++;
    next if ( $i < 7 );          # skip first 7 lines.
    next if ( $_ =~ /^\s+/ );    #skip lines that start with space
    if ( $_ =~ /:/ ) {
      my ( $id, $overall_score ) = split( /:/, $_ );
      $id            =~ s/^\s+|\s+$//g;
      $overall_score =~ s/^\s+|\s+$//g;
      next;
    }
    chomp;
    my ( $id, $align ) = split;
    $score_hash{$id} ||= '';
    $score_hash{$id} .= $align;
  }
  $FH->close;

  foreach my $member ( @{ $tree->get_all_leaves } ) {

    #
    # Do a manual insert into the given table.
    #

    my $table_name = $output_table;
    my $sth        = $pta->prepare(
      "INSERT ignore INTO $table_name 
                               (node_id,member_id,method_link_species_set_id,cigar_line)  VALUES (?,?,?,?)"
    );
    my $score_string = $score_hash{ $member->sequence_id };
    $score_string =~ s/[^\d-]/9/g
      ; # Convert non-digits and non-dashes into 9s. This is necessary because t_coffee leaves some leftover letters.

    $sth->execute( $member->node_id, $member->member_id, $member->method_link_species_set_id,
      $score_string );
    $sth->finish;
  }
}

# Finds empty sequences in $aln (SimpleAlign) and removes the according sequence from $tree (NestedSet).
sub remove_empty_seq_leaves {
  my $class = shift;
  my $aln   = shift;
  my $tree  = shift;

  foreach my $seq ( $aln->each_seq ) {
    my $seq_str       = $seq->seq;
    my $collapsed_str = $seq_str;
    $collapsed_str =~ s/[-NX]//gi;

    # REMOVE THE SEQUENCE IF IT ONLY HAS GAPS OR N's
    if ( length($collapsed_str) == 0 ) {
      print "Removing empty seq from tree/aln: " . $seq->id . "\n";
      $aln->remove_seq($seq);
      my $node = $tree->find_leaf_by_name( $seq->id );
      if ($node) {
        $tree->remove_nodes( [$node] );
      }
    }
  }
}

sub aligned_members_equal_length {
  my $class = shift;
  my $tree = shift;

  my $length;
  foreach my $leaf ($tree->leaves) {
    my $aln_str = $leaf->alignment_string;
    $length = length($aln_str) if (!defined $length);
    if ($length != length($aln_str)) {
      return 0;
    }
  }  
  return 1;
}

sub fetch_score_aln {
  my $class = shift;
  my $aa_aln = shift;
  my $tree = shift;
  my $params = shift;

  my $table = $params->{'alignment_score_table'};
  $table = 'protein_tree_member_score' unless (defined $table);
  my $pta = $tree->adaptor;

  my $new_aln = new $aa_aln;
  foreach my $leaf ($tree->leaves) {

    # Grab the score line for each leaf node.
    my $id = $leaf->stable_id;
    my $member_id = $leaf->member_id;
    my $cmd = "SELECT cigar_line FROM $table where member_id=$member_id;";
    my $sth = $pta->prepare($cmd);
    $sth->execute();
    my $data = $sth->fetchrow_hashref();
    $sth->finish();
    my $cigar = $data->{'cigar_line'};
    
    $cigar = $leaf->alignment_string unless ($cigar);
    my $new_seq = new Bio::LocatableSeq(-seq => $cigar, -id => $id);
    $new_aln->add_seq(-seq => $new_seq, -id => $id);
  }

  return $new_aln;
}

# Returns a masked alignment, either amino acid or DNA-level.
sub fetch_masked_alignment {
  my $class       = shift;
  my $aa_aln      = shift;
  my $cdna_aln    = shift;
  my $tree        = shift;
  my $params      = shift;

  my $aln = $cdna_aln;

  my $default_params = {
    sequence_quality_filtering           => 0,
    sequence_quality_threshold           => 3,
    sequence_quality_mask_character_aa   => 'X',
    sequence_quality_mask_character_cdna => 'N',
   
    cdna => 1
  };
  $params = $class->replace_params( $default_params, $params );

  #
  # SEQUENCE QUALITY MASKING
  #
  my $ann = $aln->annotation;
  if ( $params->{quality_threshold} && $params->{quality_threshold} > 0 ) {
    print "Filtering by sequence quality...\n";
    my @members = $tree->leaves;
    @members = sort {$a->stable_id cmp $b->stable_id} @members;
    foreach my $member ( @members ) {
      #print $member->stable_id . "\n";
      my @quals = $class->get_quality_array_for_member($member);

      print $member->stable_id."\n";
      my $stable_id = $member->stable_id;
      my $taxon_id = $member->taxon_id;
      my $length_nucs = $member->seq_length * 3;
      print "  $length_nucs\n";

      my $filtered_nucs = 0;

      if ( scalar(@quals) > 0 ) {
        # Get aln seq.
        my $seq = $aln->get_seq_by_id( $member->name );

        my $seq_nogaps = $seq->seq;
        $seq_nogaps =~ s/-//g;

        #printf "%s %.100s\n", $stable_id, $seq_nogaps;

        # Thread the quality scores into the aligned sequence.
        my $aligned_qualref = $class->thread_quality_scores_into_aligned_seq($seq,\@quals);

        # Filter this compara alignment by seq quality.
        my $quality_threshold = $params->{quality_threshold};

        my $filtered = $class->filter_alignment_seq_in_threes($seq, $aligned_qualref, $quality_threshold);
        my $n_filtered = $class->count_filtered_sites( $seq->seq, $filtered );
        $filtered_nucs = $n_filtered;

        $seq->seq($filtered);
      }

      if ($params->{store_filtered_codons}) {
        my $process = $params->{process};
        my $p = $process->params;
        $p = $process->replace($p, {
          taxon_id => $taxon_id,
          stable_id => $stable_id,
          filtered_nucs => $filtered_nucs,
          length_nucs => $length_nucs
                               });
        $process->create_table_from_params( $process->dbc, 'qual_filt',
                                         $class->filtered_codon_table );
        $process->store_params_in_table($process->dbc, 'qual_filt', $p);
        
      }
    }
  }

  $aln = $ALN->sort_by_tree( $aln, $tree );
  $aln->annotation($ann);
  return $aln;
}

sub filtered_codon_table {
  my $class = shift;

  return {
    data_id => 'int',
    taxon_id => 'int',
    stable_id => 'char32',
    filtered_nucs => 'int',
    length_nucs => 'int',
    unique_keys => 'data_id,taxon_id,stable_id'
  };
}

sub thread_quality_scores_into_aligned_seq {
  my $self = shift;
  my $aligned_seq = shift;
  my $qual_arrayref = shift;

  my $aln_str = $aligned_seq->seq;
  my @quals = @{$qual_arrayref};

  my @aligned_quals;
  for (my $i=0; $i < length($aln_str); $i++) {
    my $char = substr($aln_str,$i,1);
    if ($char eq '-') {
      push @aligned_quals, 0;
    } else {
      my $next_qual = shift @quals;
      #print "$i $next_qual\n";
      push @aligned_quals, $next_qual;
    }
  }
  #print "$aln_str\n";
  #print join(' ',@aligned_quals)."\n";
  return \@aligned_quals;
}

sub filter_alignment_seq_by_qual_array {
  my $class             = shift;
  my $seq               = shift;
  my $qual_arrayref     = shift;
  my $quality_threshold = shift;

  # Sequence and quality arrayref should both be in alignment coordinates, with gaps included.

  my $seq_str = $seq->seq;
  my @quals   = @$qual_arrayref;

  if ( length($seq_str) != scalar(@quals) ) {
    my $sl = length($seq_str);
    my $ql = scalar(@quals);

    # We often see the cds alignment missing the last stop codon, so a diff of 3 is OK.
    if ( $ql - $sl != 3 ) {
      print "seq : $seq_str\n";
      warn( "Warning: aln seq and qual array not same length (seq=$sl qual=$ql) for id "
          . $seq->id
          . "!" );
    }
  }

  #print join('-',@quals)."\n";

  my $filter_count = 0;
  for ( my $i = 0 ; $i < scalar(@quals) ; $i++ ) {
    next if ( substr( $seq_str, $i, 1 ) =~ m/[N\-\.]/i );

    my $q = $quals[$i];
    if ( $q < $quality_threshold ) {
      #print $q."\n";
      # Replace that nucleotide with an N.
      substr $seq_str, ($i), 1, 'N';
      $filter_count++;
    }
  }
  print STDERR "Filtered [$filter_count] from " . $seq->id . "\n";
  #print STDERR $seq_str."\n";
  return $seq_str;
}

# IMPORTANT: Always give the amino acid alignment
sub mask_aln_by_trimal {
  my $class  = shift;
  my $tree   = shift;
  my $aln    = shift;
  my $aa_aln = shift;
  my $params = shift;

  # Write temporary alignment.
  my $dir = '/tmp/eslr_temp';
  mkpath( [$dir] );
  my $node_id = $tree->node_id;
  my $aln_f   = $dir . "/aln{$node_id}.fa";
  Bio::EnsEMBL::Compara::AlignUtils->to_file( $aa_aln, $aln_f );

  # Build a command for TrimAl.
  my $cmd    = "trimal -in $aln_f -out $aln_f -colnumbering -cons 30 -gt 0.33 -gw 3";
  my $output = `$cmd`;
  print "OUTPUT:" . $output . "\n";

  chomp $output;
  my @cons_cols = split( /[\s,]+/, $output );
  print "@cons_cols\n";

  if ( $params->{'cdna'} ) {

    # use the $aa_aln for getting the columns, and multiply by 3.
    my @cdna_cons_cols;
    foreach my $i (@cons_cols) {
      my $cdna_i = int($i) * 3;
      push @cdna_cons_cols, ($cdna_i);
      push @cdna_cons_cols, ( $cdna_i + 1 );
      push @cdna_cons_cols, ( $cdna_i + 2 );
    }

    #print "@cons_cols\n";
    #print "@cdna_cons_cols\n";
    $aln =
      Bio::EnsEMBL::Compara::AlignUtils->mask_columns( $aln, \@cdna_cons_cols,
      $params->{'trimal_mask_character'} );
  } else {
    $aln =
      Bio::EnsEMBL::Compara::AlignUtils->mask_columns( $aln, \@cons_cols,
      $params->{'trimal_mask_character'} );
  }

  rmtree($dir);
  return $aln;
}

sub mask_aln_by_sequence_quality {
  my $class  = shift;
  my $tree   = shift;
  my $aln    = shift;
  my $params = shift;

  my $threshold   = $params->{'sequence_quality_threshold'};
  my $cdna_option = $params->{'cdna'};
  my $mask_char   = $params->{'sequence_quality_mask_character_aa'};
  $mask_char = $params->{'sequence_quality_mask_character_cdna'} if ($cdna_option);

  my $pta      = $tree->adaptor;
  my @twox_ids = (
    9978, 9371,   9739,  9478,  42254, 30538, 10020, 9365, 59463, 9358, 9813, 37347,
    9593, 132908, 30611, 30608, 9785,  43179, 9986,  9685, 9361
  );

  my $qual_hash_ref;
  my $sequence_quality;
  foreach my $leaf ( @{ $tree->get_all_leaves } ) {
    my $id          = $leaf->stable_id;
    my $member_id   = $leaf->member_id;
    my $sequence_id = $leaf->sequence_id;

    my $cmd = "SELECT sequence FROM sequence_quality where sequence_id = $sequence_id;";
    my $sth = $pta->prepare($cmd);
    $sth->execute();
    my $data = $sth->fetchrow_hashref();
    $sth->finish();
    $sequence_quality = $data->{'sequence'};

    if ( !defined $sequence_quality || $sequence_quality eq '' ) {

      #print "NO seq quality: $id\n";
      $sequence_quality = '9' x length( $leaf->sequence );
    }

    #    }

    my $cigar_line = $leaf->cigar_line;

    #print $leaf->stable_id."\n";
    my $qual_cigar_line = $ALN->combine_seq_cigar( $sequence_quality, $cigar_line );

    # Convert the protein mask into a CDNA mask by repeating each char 3 times.
    my @arr = split( //, $qual_cigar_line );
    if ($cdna_option) {
      @arr = map { ( $_ . '' ) x 3 } @arr;
    }
    $qual_cigar_line = join( "", @arr );
    $qual_hash_ref->{$id} = $qual_cigar_line;
  }

  $aln = $ALN->mask_below_score( $aln, $threshold, $qual_hash_ref, $mask_char );
  return $aln;
}

# Fetches an alternate set of ProteinTree alignments from $alt_member_table.
sub fetch_alternate_tree_aln {
  my $class            = shift;
  my $tree             = shift;
  my $alt_member_table = shift;

  foreach my $leaf ( @{ $tree->get_all_leaves } ) {

    # "Release" the stored / cached values for the alignment strings.
    undef $leaf->{'cdna_alignment_string'};
    undef $leaf->{'alignment_string'};

    # Grab the correct cigar line for each leaf node.
    my $id  = $leaf->member_id;
    my $cmd = "SELECT cigar_line FROM $alt_member_table where member_id=$id;";
    my $sth = $tree->adaptor->prepare($cmd);
    $sth->execute();
    my $data = $sth->fetchrow_hashref();
    $sth->finish();
    my $cigar = $data->{'cigar_line'};

    die "No cigar line for member $id!\n" unless ($cigar);
    $leaf->cigar_line($cigar);
  }
  return $tree;
}

######
###### DUMPING / LOADING TREES & ALIGNMENTS
######

# Writes the given ProteinTree to the given directory. Returns the filename.
sub dump_ProteinTree_tree {
  my $class     = shift;
  my $tree      = shift;
  my $fastafile = shift;

  my $treeI   = $class->to_TreeI($tree);
  my $treeout = Bio::TreeIO->new(
    '-format' => 'newick',
    '-file'   => ">$fastafile"
  );
  $treeout->write_tree($tree);
  $treeout->close();

  return $fastafile;
}

# Dumps a ProteinTree's sequences to a Fasta file.
sub dump_ProteinTree_seqs {
  my $class                   = shift;
  my $tree                    = shift;
  my $fastafile               = shift;
  my $include_exon_boundaries = shift if (@_);

  $fastafile =~ s/\/\//\//g;    # converts any // in path to /
  return $fastafile if ( -e $fastafile );    # Just return filename if it exists already.

  my $sa = $class->get_ProteinTree_seqs( $tree, $include_exon_boundaries );

  open( OUTSEQ, ">$fastafile" )
    or throw("Error opening $fastafile for write!");
  foreach my $seq ( $sa->each_seq ) {
    print OUTSEQ ">" . $seq->id . "\n" . $seq->seq . "\n";
  }
  close OUTSEQ;
  return $fastafile;
}

# Returns a ProteinTree's sequences as a SimpleAlign object (with NO gaps).
sub get_ProteinTree_seqs {
  my $class                   = shift;
  my $tree                    = shift;
  my $include_exon_boundaries = shift if (@_);

  my $sa = Bio::SimpleAlign->new();

  my $member_list = $tree->get_all_leaves;
  foreach my $member ( @{$member_list} ) {
    my $seq = '';
    if ($include_exon_boundaries) {
      $seq = $member->get_exon_bounded_sequence;
    } else {
      $seq = $member->sequence;
    }

    #$seq =~ s/(.{72})/$1/g;
    #chomp $seq;
    my $lseq = new Bio::LocatableSeq( -seq => $seq, -id => $member->stable_id );
    $sa->add_seq($lseq);
  }

  return $sa;
}

sub tree_aln_cdna {
  my $class  = shift;
  my $dba    = shift;
  my $params = shift;

  my $tree = $class->get_tree_for_comparative_analysis( $dba, $params );
  $tree->minimize_tree;

  my $aa = $tree->get_SimpleAlign;
  my $cdna = $tree->get_SimpleAlign( -cdna => 1 );

  my $filtered_cdna =
    Bio::EnsEMBL::Compara::ComparaUtils->fetch_masked_alignment( $aa, $cdna, $tree, $params, 1 );

  my $filtered_aa = Bio::EnsEMBL::Compara::AlignUtils->translate($filtered_cdna, $params);

  my $score_aln = Bio::EnsEMBL::Compara::ComparaUtils->fetch_score_aln($aa,$tree,$params);

  #Bio::EnsEMBL::Compara::AlignUtils->pretty_print($score_aln,{full=>1,length=>150});

  return ($tree,$filtered_aa,$filtered_cdna,$score_aln);
}

sub get_taxon_ids_from_keepers_list {
  my $self        = shift;
  my $dba         = shift;
  my $list_string = shift;

  # Get all genome_dbs, and map common names to taxon_ids.
  my $name_to_taxon_id;
  my @gdbs = Bio::EnsEMBL::Compara::ComparaUtils->get_all_genomes($dba);
  foreach my $gdb (@gdbs) {
    my $taxon_id = $gdb->taxon_id;
    my $taxon    = $gdb->taxon;

    $name_to_taxon_id->{ lc $gdb->short_name }      = $taxon_id;
    $name_to_taxon_id->{ lc $gdb->name }            = $taxon_id;
    $name_to_taxon_id->{ lc $taxon->binomial }      = $taxon_id;
    $name_to_taxon_id->{ lc $taxon->common_name }   = $taxon_id;
    $name_to_taxon_id->{ lc $taxon->ensembl_alias } = $taxon_id;

    #$name_to_taxon_id->{lc $taxon->ensembl_alias_name} = $taxon_id;
    $name_to_taxon_id->{ lc $taxon->species } = $taxon_id;

  }

  my @tokens = split( /, ?/, $list_string );
  my @results;
  foreach my $token (@tokens) {
    if ( defined $name_to_taxon_id->{ lc $token } ) {
      push @results, $name_to_taxon_id->{ lc $token };
    } else {
      push @results, $token;
    }
  }
  return @results;
}

sub restrict_aln_to_tree {
  my $class = shift;
  my $aln   = shift;
  my $tree  = shift;

  my $ok_id_hash;

  foreach my $leaf ( $tree->leaves ) {
    $ok_id_hash->{ $leaf->name } = 1;
    if ($leaf->isa("Bio::EnsEMBL::Compara::Member")) {
      #print "MEMBER\n";
      $ok_id_hash->{ $leaf->stable_id } = 1;
      my $taxon = $leaf->taxon;
      $ok_id_hash->{ $leaf->genome_db->name }    = 1;
      $ok_id_hash->{ $taxon->taxon_id }    = 1;
      $ok_id_hash->{ $taxon->common_name } = 1;
      $ok_id_hash->{ $taxon->short_name }  = 1;
      $ok_id_hash->{ $taxon->binomial }    = 1;
    }
  }

  my $new_aln = $aln->new;
  foreach my $seq ( $aln->each_seq ) {
    my $new_seq = Bio::LocatableSeq->new( -seq => $seq->seq, -id => $seq->id );
    if ( $ok_id_hash->{ $seq->id } ) {
      $new_aln->add_seq($new_seq);
    }
  }
  return $new_aln;
}

sub get_species_subtree {
  my $class        = shift;
  my $dba          = shift;
  my $species_tree = shift;
  my $params       = shift;

  $species_tree = $species_tree->copy();

  my $taxon_a = $dba->get_NCBITaxonAdaptor;

  my %keep_hash;
  my %remove_hash;
  my @keep_species = $class->get_taxon_ids_from_keepers_list( $dba, $params->{keep_species} );
  my @keep_names = map {
    my $taxon    = $taxon_a->fetch_node_by_taxon_id($_);
    my $binomial = $taxon->binomial;
    $binomial =~ s/\s/_/g;
    $binomial;
  } @keep_species;
  map { $keep_hash{$_} = 1 } @keep_names;

  foreach my $leaf ( @{ $species_tree->get_all_leaves } ) {

    #print "label: ".$leaf->name."\n";
    my $name = $leaf->name;
    if ( !exists $keep_hash{$name} ) {
      $remove_hash{$name} = 1;
    }
  }
  my @remove_names = keys %remove_hash;

  # Add fake node_ids for use by the TreeUtils methods.
  my $count = 1;
  map { $_->node_id( $count++ ) } @{ $species_tree->get_all_nodes };
  if ( scalar @remove_names > 0 ) {
    $species_tree = $TREE->remove_members_by_method_call( $species_tree, \@remove_names, 'name' );
  }

  return $species_tree;
}

sub get_one_to_one_ortholog_tree {
  my $class = shift;
  my $dba = shift;
  my $member = shift;
  my $keep_types = shift;

  $keep_types = 'ortholog_one2one' unless (defined $keep_types);

  my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor('human', 'core', 'gene');
  my $pta = $dba->get_ProteinTreeAdaptor;
  my $member_adaptor = $dba->get_MemberAdaptor;
  my $homology_adaptor = $dba->get_HomologyAdaptor;

  if (defined $member->gene_member) {
    $member = $member->gene_member;
  }

  my $tree = $pta->fetch_by_gene_Member_root_id($member);

  return undef if (scalar($tree->leaves) < 2);

  my $orth_types; 
 my $homs = $homology_adaptor->fetch_all_by_Member($member);
  #print " Fetching homs\n";
  foreach my $hom (@$homs) {
    my $desc = $hom->description;
    next unless ($desc =~ m/$keep_types/i );
    my @members = @{$hom->gene_list};
    foreach my $m (@members) {
      #print $m->stable_id." ".$m->member_id."\n";
      my $pep_m = $m->get_canonical_peptide_Member;
      $orth_types->{$pep_m->member_id} = $desc;
    }
  }

  if (defined $tree) {
    #print "  Collecting keepers\n";
    my $keep_members;
    foreach my $leaf ($tree->leaves) {
      #print $leaf->stable_id." ". $leaf->gene_member_id."\n";
      if ($orth_types->{$leaf->gene_member_id}) {
        $keep_members->{$leaf} = $leaf;
      }
      if ($orth_types->{$leaf->member_id}) {
        $keep_members->{$leaf} = $leaf;
      }
    }
    
    #print "  Extracting subtree\n";
    my @keeps = values %$keep_members;
    $tree = $TREE->extract_subtree_from_leaf_objects($tree, \@keeps);
    
    if (defined $tree) {
      foreach my $member ($tree->leaves) {
        my $type = $orth_types->{$member->member_id};
        $member->add_tag('orth_type',$type);
      }
    }
  }

  warn("Nothing left after extracting one2one!") unless ( defined $tree );
  return $tree;
}

sub restrict_tree_to_clade {
  my $class = shift;
  my $dba = shift;
  my $tree = shift;
  my $clade_name = shift;
  my $genome_tree_cache = shift;

  my $genome_tree = $genome_tree_cache;
  if (!defined $genome_tree_cache) {
    $genome_tree = $class->get_genome_tree($dba);
  }
  
  my $clade_node;
  foreach my $node ($genome_tree->nodes) {
    if ($node->name eq $clade_name) {
      $clade_node = $node;
      last;
    }
  }
  
  die("No clade found to restrict tree!") unless (defined $clade_node);

  my $tax_id_hash = {};
  map {$tax_id_hash->{$_->taxon_id} = 1} $clade_node->leaves;
  
  my @keep_names = ();
  foreach my $leaf ($tree->leaves) {
    if ($tax_id_hash->{$leaf->taxon_id}) {
      push @keep_names, $leaf->name;
    }
  }

  $tree = Bio::EnsEMBL::Compara::TreeUtils->extract_subtree_from_names($tree, \@keep_names, 0);
  return $tree;
}

sub get_tree_for_comparative_analysis {
  my $class  = shift;
  my $dba    = shift;
  my $params = shift;

  my $default_params = {
    keep_species   => '',
    remove_species => '',
    flow_node_set => ''
  };
  $params = $class->replace_params( $default_params, $params );

  #my $ps_params = $class->load_params_from_param_set($dba->dbc,$params->{'parameter_set_id'});
  #$params = $class->replace_params($params,$ps_params);

  # Fetch the tree.
  my $tree;
  if ( !defined $params->{tree} ) {
    my $node_id = $params->{'node_id'};
    $node_id = $params->{'protein_tree_id'} unless ( defined $node_id );
    die("Need a node ID for fetching a tree!") unless ($node_id);
    my $pta = $dba->get_ProteinTreeAdaptor;

    # GJ 2009-09-18
    if ( $params->{'alignment_table'} ) {
      print "ALIGNMENT TABLE:" . $params->{alignment_table} . "\n";
      #$pta->protein_tree_member( $params->{'alignment_table'} );
    }
    if ( $params->{'alignment_score_table'} ) {
      #$pta->protein_tree_score( $params->{'alignment_score_table'} );
    }

    print "Fetching tree from node $node_id...\n" if ($params->{debug});
    $tree = $pta->fetch_node_by_node_id($node_id);
  } else {
    $tree = $params->{tree};
  }

  
  # Remove members not identified as one2one orthologs of the gene_id member.
  if ($params->{flow_node_set} eq 'one2one') {
    die ("Gene ID not set w/ method one2one!") unless (defined $params->{gene_id});

    print "  keeping only one2one orthologs...\n";
    #print "  original tree:" . $tree->newick_format."\n";

    my $id = $params->{gene_id};
    
    my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor('human', 'core', 'gene');
    my $member_adaptor = $dba->get_MemberAdaptor;
    my $homology_adaptor = $dba->get_HomologyAdaptor;

    my @genes;
    push @genes, @{ $gene_adaptor->fetch_all_by_external_name($id)};
    if (scalar(@genes) == 0) {
      my $stable_id_gene = $gene_adaptor->fetch_by_stable_id($id);
      push @genes, $stable_id_gene if (defined $stable_id_gene);
    }
    my $member;
    foreach my $gene (@genes) {
      my $s_id = $gene->stable_id;
      $member = $member_adaptor->fetch_by_source_stable_id(undef, $s_id);
      last if (defined $member);
    }

    $member = $member_adaptor->fetch_by_source_stable_id(undef, $id) if (!defined $member);
    die("Cannot find gene for id $id!") unless (defined $member);

    my $gene_member = $member->gene_member;
    #print "MEMBER: ".$member->stable_id."\n";

    my $one2one_stable_ids;
    my $homs = $homology_adaptor->fetch_all_by_Member($gene_member);
    foreach my $hom (@$homs) {
      #print "$hom\n";
      next unless ($hom->description eq 'ortholog_one2one');
      my @members = @{$hom->gene_list};
      foreach my $m (@members) {
        #print $m->stable_id."\n";
        $one2one_stable_ids->{$m->stable_id} = 1;
      }
    }

    my $keep_members;
    foreach my $leaf ($tree->leaves) {
      if ($one2one_stable_ids->{$leaf->gene_member->stable_id}) {
        $keep_members->{$leaf} = $leaf;
      }
    }
    #print join(",",values %$keep_members)."\n";
    my @keeps = values %$keep_members;
    $tree = $TREE->extract_subtree_from_leaf_objects($tree, \@keeps);

    warn("Nothing left after extracting one2one!") unless ( defined $tree );

    return -1 if (!defined $tree);
  }


  my $keep_species   = $params->{'keep_species'};
  my $remove_species = $params->{'remove_species'};
#  print "KEEP{$keep_species} REMOVE{$remove_species}\n";
  my @all_ids;
  my %keep_hash;
  my %remove_hash;

  # Remove nodes not within our desired taxonomic subtree.
  if ( $keep_species ne '' ) {
    
    # Collect all taxon IDs within the tree.
    @all_ids = $TREE->get_species_in_tree($tree);
    
    
    # Find a list of all species in the tree, and add any species NOT in the keepers list to the remove list.
    my @ks = $class->get_taxon_ids_from_keepers_list( $dba, $keep_species );
    map { $keep_hash{$_} = 1 } @ks;
    foreach my $tax_id (@all_ids) {
      if ( !exists $keep_hash{$tax_id} ) {
        $remove_hash{$tax_id} = 1;
      }
    }
  }

  # Explicitly add anything that shows up in the 'remove_species' parameter to the remove list.
  if ( $remove_species ne '' ) {
    my @rs = split( ",", $remove_species );
    map { $remove_hash{$_} = 1 } @rs;
  }

  # Do the actual removal of everything in the remove_hash.
  my @taxon_ids = keys %remove_hash;
  if ( scalar @taxon_ids > 0 ) {
    $tree = $TREE->remove_members_by_taxon_id( $tree, \@taxon_ids );
    throw("Nothing left after removing species!") unless ( defined $tree );
  }

  # IMPORTANT: Re-root the tree so we get rid of parents above this one.
  $tree->re_root;
  return $tree;
}

my $qual_base = "/nfs/users/nfs_g/gj1/scratch/2x_quality/assemblies/";

my $ens_to_index = {};

sub get_quality_array_for_member {
  my $class  = shift;
  my $member = shift;

  # GJ 2009-07-01
  my $gdb  = $member->genome_db;
  my $name = $gdb->short_name;
  $name =~ s/ /_/g;

  my $tx = $member->get_Transcript;

  print $tx->stable_id."\n";

  my $cache_object = {};

  my @qual_array      = ();
  my $fake_char_count = 0;
  foreach my $exon ( @{ $tx->get_all_translateable_Exons } ) {
    my $slice = $exon->slice;
    print $slice->name."\n"; 
   my $exon_slice = $slice->sub_Slice( $exon->start, $exon->end, $exon->strand );
    if (!defined $exon_slice) {
      my $ens_seq = $exon->seq->seq;
      my $len = length($ens_seq);
      my @extra_qs = (0) x $len;
      push @qual_array, @extra_qs;
      next;
    }
    #print $exon_slice->name."\n";

    my ($dna,$qual_str) = $class->get_bases_for_slice( $exon_slice, $gdb, $cache_object );

    my $ens_seq = $exon->seq->seq;
    if ( $dna ne '' ) {

      #print "$dna\n";
      #print "$qual_str\n";
      #print "$ens_seq\n";
      my @q_array = split( ' ', $qual_str );
      if ( $dna ne $ens_seq ) {
        warn("Ensembl and local dna not the same!");
        print "$dna\n";
        print "$qual_str\n";
        print "$ens_seq\n";
        my $diff = length($ens_seq) - scalar(@q_array);
        if ( $diff > 0 ) {
          my @extra_qs = (0) x $diff;
          push @q_array, @extra_qs;
          $dna .= 'N' x $diff;
          warn("Padding dna sequence with Ns!");
        } else {
          die("Qual array is *longer* than ens seq!");
        }
      }
      push @qual_array, @q_array;
      die("Seq and qual length not the same!") unless ( length($dna) == scalar(@q_array) );
    } else {

      #print "No scores: ens_seq[$ens_seq]\n";
      # Test for an ensembl sequence of all N's -- if this is the case, we should pad the qual
      # array with zeros.
      my $tmp = $ens_seq;
      $tmp =~ s/n//gi;
      if ( length($tmp) == 0 ) {
        #warn("Ensembl sequence for exon is all Ns!") if (debug);
        my @qs = (0) x length($ens_seq);
        push @qual_array, @qs;
        $fake_char_count += length($ens_seq);
      }
    }
  }
  if ( scalar(@qual_array) == $fake_char_count ) {

    # Return nothing if all we have is fake chars.
    return ();
  }
  return @qual_array;
}

# Returns the NCBI taxnomy of Ensembl genomes below a given taxonomic clade.
sub get_genome_taxonomy_below_level {
  my $class         = shift;
  my $dba           = shift;
  my $root_taxon_id = shift || 33154;
  my $params = shift || {};

  my @gdbs = $class->get_all_genomes($dba);

  my @ncbi_ids = map { $_->taxon_id } @gdbs;

  my $taxon_a = $dba->get_NCBITaxonAdaptor;
  $taxon_a->clear_cache;

  # Try first with a taxon label.
  my $root = $taxon_a->fetch_node_by_name($root_taxon_id);
  if ( !defined $root ) {
    my $quoted_id = $dba->dbc->db_handle->quote($root_taxon_id);
    $root = $taxon_a->fetch_node_by_taxon_id($quoted_id);
  }

  die("Cannot find genome taxonomy for root: [$root_taxon_id]!\n") unless (defined $root);

  # Collect all genome_db leaves, plus their internal lineages, into an array.
  my %keepers;
  $keepers{ $root->node_id } = $root;
  foreach my $gdb (@gdbs) {

    # THIS IS EXTREMELY IMPORTANT!!!
    # If you run things on the hive and don't remove this cached variable, you'll
    # get invisible and impossible to debug errors where the gdb's taxon inexplicably
    # has no parental node. This is a hack, but does the trick.
    delete $gdb->{'_taxon'};
    my $tx    = $gdb->taxon;
    my $tx_id = $tx->taxon_id;

    if ( !$TREE->has_ancestor_node_id( $tx, $root ) ) {
      #warn( "taxon not below tax level: [" . $tx->node_id . "]" );
      next;
    } else {
      print "Okay: " . $tx->name . "\n" if ($params->{debug});
    }

    my $node = $tx;
    while ( defined $node ) {
      $keepers{ $node->node_id } = $node;
      last if ( !$TREE->has_ancestor_node_id( $node, $root ) );
      $node = $node->parent;
    }
    print "\n" if ($params->{debug});
  }

  my @nodes = values %keepers;
  print "GenomeDB tree size: " . scalar(@nodes) . "\n" if ($params->{debug});

  my $new_tree = $taxon_a->_build_tree_from_nodes( \@nodes );
  $new_tree = $new_tree->minimize_tree;

  return $new_tree;
}

sub get_genomes_within_clade {
  my $class = shift;
  my $dba   = shift;
  my $clade = shift || 1;

  my @gdbs = $class->get_all_genomes($dba);
  my $species_tree = $class->get_genome_taxonomy_below_level( $dba, $clade );

  my @genomes;
  foreach my $gdb (@gdbs) {
    my $leaf = $species_tree->find_node_by_node_id( $gdb->taxon->taxon_id );
    push @genomes, $gdb if ($leaf);
  }

  return @genomes;
}

sub fix_genome_polytomies {
  my $class = shift;
  my $tree  = shift;

  # Fix up any multifurcations in the tree. For now we'll follow:
  # http://mbe.oxfordjournals.org/cgi/content/full/26/6/1259/FIG6

  my $node;
  my $new_node;

  my $taxid_hash = {
    taxon_id => 1000000
  };

  # Fix homo/pan/gor.
  $class->_fix_multifurcation($tree, 'Homininae', [ 'Homo sapiens', 'Pan troglodytes' ], $taxid_hash );

  # Fix squirrel and rodents.
  $class->_fix_multifurcation( $tree, 'Sciurognathi', [ 'Murinae', 'Dipodomys ordii' ], $taxid_hash );  

  ### Fix the laurasiatheria mess.
  # Fix cow & dolphin.
  $new_node = $class->_fix_multifurcation( $tree, 'Cetartiodactyla', [ 'Bos taurus', 'Tursiops truncatus' ], $taxid_hash );  
  $new_node->name('cow-dolphin');
  # Bring horse out next to the cetartiodactyla.
  $new_node = $class->_fix_multifurcation( $tree, 'Laurasiatheria', [ 'Cetartiodactyla', 'Equus caballus' ], $taxid_hash );  
  $new_node->name('horse-et-al');
  $new_node = $class->_fix_multifurcation( $tree, 'Laurasiatheria', [ 'Chiroptera', 'Carnivora' ], $taxid_hash );  
  $new_node->name('bats-carnivores');
  $new_node = $class->_fix_multifurcation( $tree, 'Laurasiatheria', [ 'horse-et-al', 'bats-carnivores' ], $taxid_hash );    
  # Fix cow & dolphin.
  $class->_fix_multifurcation( $tree, 'Cetartiodactyla', [ 'cow-dolphin', 'Vicugna pacos' ], $taxid_hash );  

  # Fix Afrotheria.
  $class->_fix_multifurcation( $tree, 'Afrotheria', [ 'Procavia capensis', 'Loxodonta africana' ], $taxid_hash );

  # Atlantogenata following http://www.pnas.org/content/106/13/5235.full
  $class->_fix_multifurcation( $tree, 'Eutheria', [ 'Afrotheria', 'Xenarthra' ], $taxid_hash );
  $tree->find_node_by_name('Afrotheria')->parent->name('Atlantogenata');

  # Fix Euarchontoglires, Laurasiatheria, Afrotheria.
  $class->_fix_multifurcation( $tree, 'Eutheria', [ 'Laurasiatheria', 'Euarchontoglires' ], $taxid_hash );
  $class->_fix_multifurcation( $tree, 'Euarchontoglires', [ 'Primates', 'Tupaia belangeri'], $taxid_hash );

  my $n = $tree->find_node_by_name('Fungi/Metazoa group');
  $n->name('Eukaryota') if ($n);
  $n = $tree->find_node_by_name('Opisthokonta');
  $n->name('Eukaryota') if ($n);

  #print $tree->ascii;

  return $tree;
}

sub _fix_multifurcation {
  my $class                   = shift;
  my $tree = shift;
  my $polytomy_name                = shift;
  my $new_group_name_arrayref = shift;
  my $taxid_hash = shift;

  my $polytomy = $tree->find_node_by_name($polytomy_name);
  
  if (!defined $polytomy || !$TREE->is_polytomy($polytomy)) {
    warn("Tree node not found for polytomy [$polytomy_name]\n");
    return undef;
  }

  #print "Polytomy: ". $polytomy->newick_format . "\n";
  #print "Fixing multifurcation: " . $polytomy->newick_format . "\n";

  my @new_group_names = @$new_group_name_arrayref;

  my $new_group = new $polytomy;
  $new_group->taxon_id($taxid_hash->{taxon_id});
  $taxid_hash->{taxon_id} = $taxid_hash->{taxon_id} + 1;

  $new_group->distance_to_parent(1);
  $polytomy->add_child($new_group);

  my @children = @{ $polytomy->sorted_children };

  foreach my $name (@new_group_names) {
    my ($node) = grep { $_->find_node_by_name($name) } @children;
    die("No node found for $name!") unless ($node);
    $new_group->add_child($node);
  }

  #print $new_group->name."  ". $new_group->taxon_id."\n";
  #print "New group: ". $new_group->newick_format . "\n";
  return $new_group;
}

# Gets the genome tree representing the current species being analyzed.
# Uses the keep_species parameter to restrict the tree.
sub get_genome_tree_subset {
  my $class = shift;
  my $dba = shift;
  my $params = shift;
  my $gdb_tree_hash = shift;

  my $full_tree = $gdb_tree_hash;
  if (!defined $full_tree) {
    $full_tree = $class->get_genome_tree($dba);
  }
  #print $full_tree->newick_format."\n";
  my @keepers = $class->get_taxon_ids_from_keepers_list( $dba, $params->{keep_species} );

  my $tree = $full_tree;
  if (scalar(@keepers) > 0) {
    $tree = Bio::EnsEMBL::Compara::TreeUtils->keep_members_by_method_call($tree,\@keepers,'taxon_id');
  }
  return $tree;
}

sub get_genome_tree {
  my $class = shift;
  my $dba   = shift;

  my @gdbs = $class->get_all_genomes($dba);
  my $species_tree = $class->get_genome_taxonomy_below_level( $dba );
  $species_tree = $class->fix_genome_polytomies($species_tree);
  foreach my $gdb (@gdbs) {
    my $tx    = $gdb->taxon;
    my $tx_id = $tx->taxon_id;
    my $leaf  = $species_tree->find_node_by_node_id($tx_id);
    $leaf->name( $tx->binomial ) if ( defined $leaf );
  }

  return $species_tree;
}

# Returns the NCBI taxonomy tree of Ensembl genomes in NHX format.
sub get_genome_tree_with_extras {
  my $class  = shift;
  my $dba    = shift;
  my $params = shift;

  my @gdbs = $class->get_all_genomes($dba);
#  my $species_tree = $class->get_genome_taxonomy_below_level( $dba, 1 );
  my $species_tree = $class->get_genome_taxonomy_below_level( $dba, 33154 );
  $species_tree = $class->fix_genome_polytomies($species_tree);

  my $labels_option = $params->{'labels'};
  my $include_imgs  = $params->{'images'};

  @gdbs = sort { $a->taxon->binomial cmp $b->taxon->binomial } @gdbs;
  foreach my $gdb (@gdbs) {
    my $tx    = $gdb->taxon;
    my $tx_id = $tx->taxon_id;

    #print $gdb->name . "\n";
    my $gdb_dba  = $gdb->db_adaptor;
    my $coverage = 'high';
    if ($gdb_dba) {
      my $meta = $gdb_dba->get_MetaContainer;
      $coverage = @{ $meta->list_value_by_key('assembly.coverage_depth') }[0];
    }

    my $leaf = $species_tree->find_node_by_node_id($tx_id);
    next unless $leaf;

    $leaf->add_tag('coverage', $coverage);

    if ( $coverage eq 'low' ) {
      #print "  -> 2x!\n";
      $leaf->add_tag("low_coverage", 1);
    } else {
      #$leaf->add_tag("coverage", 1);
    }

#    $leaf->add_tag( "ncol", "red" ) if ( $coverage eq 'low' );
#    $leaf->add_tag( "bcol", "gray" );

    if ($include_imgs) {
      my $underbar_species = $tx->binomial;
      $underbar_species =~ s/ /_/g;
      $leaf->add_tag( "img", "http://www.ensembl.org/img/species/pic_${underbar_species}.png" );
    }
    if ( $labels_option eq 'mnemonics' ) {
      my $sql =
        "select stable_id from member where taxon_id=$tx_id and source_name='ENSEMBLPEP' limit 1;";
      my $sh = $dba->dbc->prepare($sql);
      $sh->execute();
      my @vals = $sh->fetchrow_array();
      my $ensp = $vals[0];
      @vals = ( $ensp, $tx_id, $tx->ensembl_alias, $tx->binomial, $tx->short_name );
      my $pretty_str = sprintf( "%-20s  %6s  %-22s  %-30s  %-8s", @vals );
      $leaf->name($pretty_str);
    } elsif ( $labels_option eq 'binomial' ) {
      $leaf->name( $tx->binomial );
    } else {
      $leaf->name( $tx->ensembl_alias );
    }
  }
  return $species_tree;
}

sub get_genome_tree_nhx {
  my $class  = shift;
  my $dba    = shift;
  my $params = shift;

  return $class->get_genome_tree_with_extras( $dba, $params )->nhx_format;
}

sub get_all_genomes {
  my $class = shift;
  my $dba   = shift;

  my $gda     = $dba->get_GenomeDBAdaptor();
  my $all_dbs = $gda->fetch_all();
  my @all_genomes;
  foreach my $db (@$all_dbs) {
    push @all_genomes, $db if ( $db->taxon );
  }
  return @all_genomes;
}


# Return a hashref object from a parameter string.
sub load_params_from_string {
  my $class        = shift;
  my $param_string = shift;
  my $debug        = shift;

  my $param_object;
  my $tmp_params = eval($param_string);
  foreach my $key ( keys %$tmp_params ) {
    print( "  $key\t=>\t", $tmp_params->{$key}, "\n" ) if ($debug);
    $param_object->{$key} = $tmp_params->{$key};
  }

  return $param_object;
}

# Create a string from a hashref.
sub hash_to_string {
  my $class   = shift;
  my $hashref = shift;

  my %hash = %{$hashref};

  my $str     = "{";
  my @keyvals = ();
  foreach my $key ( sort keys %hash ) {
    push( @keyvals, "'$key'" . "=>" . "'" . $hash{$key} . "'" );
  }

  $str .= join( ",", @keyvals );
  $str .= "}";
  return $str;
}

sub hash_print {
  my $class   = shift;
  my $hashref = shift;

  print "{\n";
  foreach my $key ( sort keys %{$hashref} ) {
    printf( "    %-40.40s => %-40s\n", $key, $hashref->{$key} );
  }
  print "}\n";
}

# Eval a hashref from a string.
sub string_to_hash {
  my $class  = shift;
  my $string = shift;

  my $hash = eval $string;
  return $hash if ( defined $hash );

  return {};
}

sub load_params_from_tree_tags {
  my $class   = shift;
  my $dba     = shift;
  my $node_id = shift;

  my $pta  = $dba->get_ProteinTreeAdaptor;
  my $tree = $pta->fetch_node_by_node_id($node_id);

  my $tags = $tree->get_tagvalue_hash;

  my @param_objs;
  my $simple_tags;
  foreach my $tag ( keys %$tags ) {
    if ( $tag =~ /params/i ) {
      push @param_objs, $class->load_params_from_tag( $tree, $tag );
    } else {
      $simple_tags->{$tag} = $tags->{$tag};
    }
  }
  push @param_objs, $simple_tags;

  return $class->replace_params(@param_objs);
}

sub load_params_from_tag {
  my $class = shift;
  my $tree  = shift;
  my $tag   = shift;

  my $tag_string = $tree->get_tagvalue($tag);
  return $class->string_to_hash($tag_string);
}

sub load_params_from_param_set {
  my $class        = shift;
  my $dbc          = shift;
  my $param_set_id = shift;

  return {} unless ( defined $param_set_id );

  my $params;
  my $cmd =
    qq^SELECT parameter_value FROM parameter_set WHERE parameter_set_id=$param_set_id AND parameter_name="params";  ^;
  my $sth = $dbc->prepare($cmd);
  eval { $sth->execute(); };
  if ($@) {
    return {};
  }
  my @row;
  while ( @row = $sth->fetchrow_array ) {
    $params = eval( $row[0] );
  }
  $sth->finish;

  my $cmd =
    qq^SELECT parameter_value FROM parameter_set WHERE parameter_set_id=$param_set_id AND parameter_name="name";  ^;
  $sth = $dbc->prepare($cmd);
  $sth->execute();
  my @row;
  while ( @row = $sth->fetchrow_array ) {
    $params->{parameter_set_name} = $row[0];
  }
  $sth->finish;

  return $params;
}

sub replace_params {
  my $class  = shift;
  my @params = @_;

  my $final_params = shift @params;
  while (@params) {
    my $p = shift @params;
    if ( !ref $p ) {
      $p = eval($p);
    }

    my $new_params = {};
    foreach my $key ( keys %{$final_params} ) {
      $new_params->{$key} = $final_params->{$key};
    }
    foreach my $key ( keys %{$p} ) {
      $new_params->{$key} = $p->{$key};
    }
    $final_params = $new_params;
  }

  return $final_params;
}

# Loads a given tree and alignment (with optional CDNA alignment too) into a ProteinTree object, ready to be stored in the Compara database.
# @created GJ 2009-02-17
sub load_ProteinTree_from_files {
  my $class    = shift;
  my $treeF    = shift;    # Newick formatted string or file.
  my $alnF     = shift;    # Fasta formatted string or file, main or protein alignment
  my $cdnaAlnF = shift;    # (optional) Fasta formatted string or file, secondary or CDNA alignment.

  my $treeI;
  my $sa;
  my $cdna_sa;

  $treeI   = $TREE->to_treeI($treeF);
  $sa      = $ALN->to_aln($alnF);
  $cdna_sa = $ALN->to_aln($cdnaAlnF) if ( defined $cdnaAlnF );

  throw("Error loading tree $treeF")     unless ( defined $treeI );
  throw("Error loading alignment $alnF") unless ( defined $sa );

  my $tree = $TREE->from_treeI($treeI);

  # Re-dress the NestedSet as a ProteinTree.
  bless( $tree, 'Bio::EnsEMBL::Compara::ProteinTree' );

  foreach my $leaf ( @{ $tree->get_all_leaves } ) {
    my $id       = $leaf->name;
    my $seq      = $ALN->get_seq_with_id( $sa, $id );
    my $cdna_seq = $ALN->get_seq_with_id( $cdna_sa, $id ) if ( defined $cdna_sa );

    if ( defined $seq ) {

      # Calculate and set the cigar line.
      my $cigar = $class->cigar_line( $seq->seq );
      $leaf->cigar_line($cigar);

      # Remove gaps from sequence and set it.
      my $seqstr = $seq->seq;
      $seqstr =~ s/-//g;
      $leaf->sequence($seqstr);
    }

    if ( defined $cdna_seq ) {
      bless( $leaf, "Bio::EnsEMBL::Compara::LocalMember" );

      # Set the CDNA sequence if necessary.
      my $cdna_seqstr = $cdna_seq->seq;
      $cdna_seqstr =~ s/-//g;
      $leaf->cdna_sequence($cdna_seqstr);
      $leaf->cdna_sequence_id(0);
    }
  }
  return $tree;
}

sub restrict_tree_to_aln {
  my $class = shift;
  my $tree  = shift;
  my $aln   = shift;

  my $used_seqs;

  #Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $aln, { width => 150,full => 1});
  
  my @keepers;
  foreach my $leaf ( $tree->leaves ) {
    # Try binomial, taxon_id, and alias.
    my $seq;

    $seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id( $aln, $leaf->name );

    if ($leaf->can('stable_id') && $leaf->stable_id) {
      $seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id( $aln, $leaf->stable_id );
    }
    if ($leaf->can('genome_db') && $leaf->genome_db) {
      $seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id( $aln, $leaf->genome_db->name )
        if ( !defined $seq );
    }
    if ($leaf->can('taxon') && $leaf->taxon) {
      $seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id( $aln, $leaf->taxon->binomial )
        if ( !defined $seq );
      $seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id( $aln, $leaf->taxon->ensembl_alias )
        if ( !defined $seq );
      $seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id( $aln, $leaf->taxon->taxon_id )
        if ( !defined $seq );
    }
    next if ($used_seqs->{$seq});
    if (defined $seq) {
      $used_seqs->{$seq} = $seq;
      push @keepers, $leaf;
    }
  }
  my $pruned_tree =
    Bio::EnsEMBL::Compara::TreeUtils->extract_subtree_from_leaf_objects( $tree, \@keepers );
  return $pruned_tree;
}

sub get_species_tree_for_aln {
  my $class       = shift;
  my $compara_dba = shift;
  my $aln         = shift;

  my $species_tree = $class->get_genome_tree($compara_dba);

  die("Species tree is too small!") unless ( scalar( $species_tree->leaves ) > 20 );

  my @keepers;
  foreach my $leaf ( $species_tree->leaves ) {

    # Try binomial, taxon_id, and alias.
    my $seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id( $aln, $leaf->binomial );
    $seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id( $aln, $leaf->ensembl_alias )
      if ( !defined $seq );
    $seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id( $aln, $leaf->taxon_id )
      if ( !defined $seq );
    push @keepers, $leaf if ( defined $seq );
  }
  my $pruned_tree =
    Bio::EnsEMBL::Compara::TreeUtils->extract_subtree_from_leaf_objects( $species_tree, \@keepers );
  return $pruned_tree;
}

sub get_compara_or_genomic_aln {
  my $class      = shift;
  my $c_dba      = shift;
  my $tree       = shift;
  my $ref_member = shift;
  my $params     = shift;

  my $num_species = Bio::EnsEMBL::Compara::TreeUtils->species_count($tree);
  my $taxon_id_hashref;
  map {
    $taxon_id_hashref->{$_->taxon_id} = 1;
  } $tree->leaves;

  my $quality_threshold = $params->{quality_threshold};
  $quality_threshold = $class->default_sequence_quality_threshold
    unless ( defined $quality_threshold );

  my $aln;
  my $extra_info;

  # Use an existing alignment file if one is defined.
  my $aln_f = $params->{existing_alignment_file};
  if ($aln_f && -e $aln_f) {
    warn("Trying to load existing alignment file at [$aln_f]...");
    my $aln = Bio::EnsEMBL::Compara::AlignUtils->from_file($aln_f);
    $aln = Bio::EnsEMBL::Compara::ComparaUtils->restrict_aln_to_tree( $aln, $tree );
    $tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_aln( $tree, $aln );
    $aln = Bio::EnsEMBL::Compara::AlignUtils->sort_by_tree( $aln, $tree );
    my $pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($aln);
    return {
      tree    => $tree,
      aln     => $aln,
      pep_aln => $pep_aln,
      extra   => {}
    };
  }

  my $aln_type = $params->{aln_type};

  if ( $aln_type =~ m/genomic/i ) {
    my ( $cdna, $aa );
    if ( $aln_type eq 'genomic_primates' ) {
      my $params = {
        mlss_type => 'epo',
        species_set => 'primates',
        quality_threshold => $quality_threshold,
        restrict_to_species => $taxon_id_hashref
      };
      ( $cdna, $aa, $extra_info ) = Bio::EnsEMBL::Compara::ComparaUtils->genomic_aln_for_member(
        $c_dba,
        $ref_member,
        $params
      );
    } elsif ( $aln_type eq 'genomic_mammals' ) {
      my $params = {
        mlss_type => 'epo',
        species_set => 'mammals',
        quality_threshold => $quality_threshold,
        restrict_to_species => $taxon_id_hashref
      };
      ( $cdna, $aa, $extra_info ) =
        Bio::EnsEMBL::Compara::ComparaUtils->genomic_aln_for_member( $c_dba, $ref_member, $params);
    } elsif ( $aln_type eq 'genomic_all' ) {
      my $params = {
        mlss_type => 'epo_low_coverage',
        quality_threshold => $quality_threshold,
        restrict_to_species => $taxon_id_hashref
      };
      ( $cdna, $aa, $extra_info ) =
        Bio::EnsEMBL::Compara::ComparaUtils->genomic_aln_for_member( $c_dba, $ref_member, $params);
    }
    $aln = $cdna;
    $aln = Bio::EnsEMBL::Compara::ComparaUtils->restrict_aln_to_tree( $aln, $tree );

    $tree = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree);
    $tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_aln( $tree, $aln );

    #map {print "   [".$_."]\n";} $tree->leaves;
  } else {
    $params->{alignment_score_filtering}  = 0;
    $params->{sequence_quality_filtering} = 0;
    $params->{quality_threshold}          = $quality_threshold;

    my $cdna = $tree->get_SimpleAlign( -cdna => 1 );
    my $aa = $tree->get_SimpleAlign();

    $aln = $cdna;
    $aln = $class->fetch_masked_alignment( $aa, $cdna, $tree, $params);

    # Collect alignment annotations right after getting the masked aln.
    my $annotation = $aln->annotation;
    foreach my $key ( $annotation->get_all_annotation_keys ) {
      my @values = $annotation->get_Annotations($key);
      my $val    = $values[0];
      $extra_info->{$key} = $val->value;
    }
    
    unless ($params->{keep_blank_columns}) {
      $aln = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns_in_threes($aln);
    }
    #Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $aln, { width => 150, full => 1 } );

    $tree = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree);
    $tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_aln( $tree, $aln );
  }

  if (!defined $tree || scalar($tree->leaves) < 2) {
    return -1;
  }

  if ($params->{rename_sequences_by_taxon}) {
    my $seen_nodes;
    foreach my $node ($tree->nodes) {
      my $seen_count = 0;
      if ($node->is_leaf) {
        if (defined $seen_nodes->{$node->taxon_id}) {
          $seen_count = $seen_nodes->{$node->taxon_id};
          $seen_nodes->{$node->taxon_id}++;
        } else {
          $seen_count = 0;
          $seen_nodes->{$node->taxon_id} = 1;
        }
      }
      
      my $ens_string;
      if ($node->is_leaf) {
        $ens_string = 'ens_' . $node->taxon_id . '_' . $seen_count;
      } else {
        $ens_string = 'ens_' . $node->node_id . '_' . $seen_count;
      }
#      print "$ens_string\n";
      
      if ($node->is_leaf) {
        my $aln_map;
        $aln_map->{ $node->genome_db->name } = $ens_string;
        $aln_map->{ $node->taxon->binomial } = $ens_string;
        $aln_map->{ $node->stable_id } = $ens_string;
        $aln = Bio::EnsEMBL::Compara::AlignUtils->translate_ids( $aln, $aln_map );
        $node->name($ens_string);
      }
      
      my $tree_map;
      $tree_map->{$node->node_id} = $ens_string;
      $tree = Bio::EnsEMBL::Compara::TreeUtils->translate_ids( $tree, $tree_map );
    }
  } elsif ($aln_type =~ m/genomic/i) {
    # If we collected a genomic alignment and don't want to rename the IDs, then currently
    # we have an alignment with IDs like 'homo_sapiens' and a tree with IDs like 'ENSP000001'.
    # Let's try to match alignment sequences with tree sequences, and rename them accordingly.

    my $taxon_seen_count;
    foreach my $leaf ($tree->leaves) {
      # AlignSlice.pm (which does the multi-align output) adds numbers to track the 2nd,
      # 3rd, etc. sequence per species. So we need to try to match up the 2nd, 3rd etc. sequence
      # in the tree with the alignment.
      my $gdb_key = $leaf->genome_db->name;
      my $name = $leaf->name;
      my $aln_map;
      $aln_map->{ $leaf->taxon->binomial } = $name;
      $aln_map->{ $leaf->stable_id } = $name;

      if (defined $taxon_seen_count->{$gdb_key}) {
        # Tell the alignment to rename homo_sapiens2 to 'ENSP00000xyz'
        $aln_map->{ $gdb_key.$taxon_seen_count->{$gdb_key} } = $name;
        $taxon_seen_count->{$gdb_key}++;
      } else {
        # Tell the alignment to rename homo_sapiens to 'ENSP0000xyz'
        $aln_map->{ $gdb_key } = $name;
        $taxon_seen_count->{$gdb_key} = 2;
      }
      $aln = Bio::EnsEMBL::Compara::AlignUtils->translate_ids( $aln, $aln_map );
    }
  }

  if ( scalar( $tree->leaves ) < 2 ) {
    return;
  }

  my $trimmed_num_species = Bio::EnsEMBL::Compara::TreeUtils->species_count($tree);
  if ( $trimmed_num_species != $num_species ) {
    warn("Tree didn't stay the same!! before: $num_species  after: $trimmed_num_species");
    
    if ($params->{'fail_on_altered_tree'}) {
      print "FAILING!!!\n";
      return undef;
    }
  }

  #Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $aln, { width => 150, full => 1 } );

  # Compare the alignment pep sequence to the original transcript (sanity check).
  if ($aln_type =~ m/(genomic)/i && 
        !Bio::EnsEMBL::Compara::AlignUtils->contains_sequence( $aln, $ref_member->sequence_cds, $params )) {

    my $pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($aln, $params);
    
    Bio::EnsEMBL::Compara::AlignUtils->pretty_print($aln,{full => 1});
    Bio::EnsEMBL::Compara::AlignUtils->pretty_print($pep_aln,{full => 1});
    print $ref_member->sequence_cds." ".$ref_member->stable_id."\n";
    print $ref_member->sequence." ".$ref_member->stable_id."\n";
    #warn("Alignment doesn't contain the exact ref member CDS... checking for amino acid identity.");
    
    die("Alignment doesn't contain ref member cds!") unless ($extra_info->{off_phase_start});
  }

  # Flatten and filter genomic aligns.
  if ( $aln_type =~ m/(genomic)/i || $aln_type =~ m/(compara)/i && $params->{flatten_aln}) {
    print "Before flattening: " . $aln->length . "\n";
    print $ref_member->stable_id."\n";
    print $ref_member->name."\n";
    $aln = Bio::EnsEMBL::Compara::AlignUtils->flatten_to_sequence( $aln, $ref_member->name);
    $aln = Bio::EnsEMBL::Compara::AlignUtils->filter_stop_codons($aln, $params);
    $aln = Bio::EnsEMBL::Compara::AlignUtils->ensure_multiple_of_three($aln);
    if ( Bio::EnsEMBL::Compara::AlignUtils->has_stop_codon($aln, $params) ) {
      die("STOP CODON!!!");
    }
    print "After flattening: " . $aln->length . "\n";
  }

  $aln = Bio::EnsEMBL::Compara::AlignUtils->sort_by_tree( $aln, $tree );
  my $pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($aln, $params);

  return {
    tree    => $tree,
    aln     => $aln,
    pep_aln => $pep_aln,
    extra   => $extra_info
  };
}

sub genomic_aln_for_tree {
  my $class        = shift;
  my $tree         = shift;
  my $ref_taxon_id = shift;

  my (@ref_members) = grep { $_->taxon_id == $ref_taxon_id } @{ $tree->get_all_leaves };
  die "No member for taxon_id $ref_taxon_id found!" unless ( scalar @ref_members > 0 );
  my $ref_member = $ref_members[0];

  #print $ref_member->stable_id . "\n";
  #return $class->genomic_aln_for_member($ref_member);
}

sub get_bases_for_slice {
  my $class        = shift;
  my $slice        = shift;
  my $gdb          = shift;
  my $cache_object = shift;

  my $taxon_id = $gdb->taxon_id;
  my $name     = $gdb->taxon->ensembl_alias;

  my $qual_base = "/nfs/users/nfs_g/gj1/scratch/2x_quality/assemblies";

  my $seq_name   = $slice->seq_region_name;
  my $start      = $slice->start;
  my $end        = $slice->end;
  my $strand     = $slice->strand;
  my $short_name = $gdb->taxon->short_name;

  my $bases_file = "${qual_base}/${short_name}.bases";
  my $quals_file = "${qual_base}/${short_name}.quals";
  my $seq        = '';
  my $quals      = '';
  my $qseq       = ();
  my $cur_slice;
  if ( -e $bases_file ) {
    my $slice_a = Bio::EnsEMBL::Registry->get_adaptor( $name, 'core', 'slice' );
    $slice->adaptor($slice_a);
    my $csa = Bio::EnsEMBL::Registry->get_adaptor( $name, 'core', 'coordsystem' );
    my @coord_systems;

    # Cache the indexedfastas within the gdb object.
    if (!defined $cache_object->{_qdb}) {
      $cache_object->{_qdb} = Bio::Greg::IndexedFasta->new($quals_file);
      $cache_object->{_db} = Bio::Greg::IndexedFasta->new($bases_file);
    }
    my $qdb = $cache_object->{_qdb};
    my $db = $cache_object->{_db};

    my $qobj;
    my $id = '';

    # First, we create a 'template' sequence string.
    my $length        = $end - $start + 1;
    my $seq_template  = 'N' x $length;
    my @qual_template = (-1) x $length;

    # Find the right coordinate system.
    my $coordsystem;
    my $tried_ids = '';
    foreach my $cs_name ( 'chromosome', 'contig', 'scaffold', 'supercontig' ) {
      $coordsystem = $cs_name;
      my $cs = $csa->fetch_by_name($cs_name);
      next unless ( defined $cs );
      my $proj = $slice->project($cs_name);
      foreach my $seg (@$proj) {
        my $projected_slice = $seg->to_Slice();
        $id = $projected_slice->seq_region_name;
        $tried_ids .= ' ' . $id;
        if ( $coordsystem eq 'chromosome' ) {
          $id = 'chr' . $id;
        }
        last if ( $db->has_key($id) );
      }
      last if ( $db->has_key($id) );
    }

    if (!$id || !$db->has_key($id)) {
      print "Quals file exists, but no db key found for $tried_ids\n";
    }

    # Return nothing if none of the coordinate systems worked.
    return '' unless ( $id && $db->has_key($id) );

    # Now, go through the projections and layer whatever seq / qual we find onto the template.
    my $projections = $slice->project($coordsystem);
    foreach my $seg (@$projections) {
      my $from_start  = $seg->from_start;
      my $from_end    = $seg->from_end;
      my $from_length = $from_end - $from_start + 1;
      my $to_slice    = $seg->to_Slice();

      my $id = $to_slice->seq_region_name;
      $id = 'chr'.$id if ($coordsystem eq 'chromosome');

      # Get the projected seq.
      $seq = $db->get_sequence_region( $id, $to_slice->start, $to_slice->end );
      $seq = uc($seq);
      #print ("Contig[".$to_slice->strand."] slice[".$slice->strand."]\n");
      if ($to_slice->strand == -1) {
        my $bioseq = Bio::PrimarySeq->new( -seq => $seq );
        $seq = $bioseq->revcom()->seq;
      }

      if ( length($seq) != $from_length ) {
        print "$seq\n";
        print "$from_length ".length($seq)."\n";
#        die("Length(seq) != from_length!");
      }

      substr( $seq_template, $from_start-1, $from_length, $seq );

      # Get the projected quals.
      $quals = $qdb->get_qual_region( $id, $to_slice->start, $to_slice->end );
      my @quals = split( ' ', $quals );

      if ($to_slice->strand == -1) {
        @quals = reverse @quals;
      }
      splice( @qual_template, $from_start-1, $from_length, @quals );
    }

    $quals = join(' ',@qual_template);
    $seq = $seq_template;
  } else {
    #print "No bases or quals file: " . $quals_file."\n";
  }

  if ($taxon_id == 59463) {
    print "$seq\n";
    print "$quals\n";
  }

  return ($seq,$quals);
}

sub default_sequence_quality_threshold {
  my $class = shift;
  return 30;
}

sub combine_slice_positions_into_ranges {
  my $class              = shift;
  my $slice_pos_arrayref = shift;

  my @slice_range_output = ();    # array of hashrefs.

  # Create an arrayref of positions for each slice.
  my $cur_slice;
  my $cur_pos;
  my $cur_start;
  my $cur_end;
  foreach my $slice_pos (@$slice_pos_arrayref) {
    my ( $slice, $pos ) = @$slice_pos;
    if (
         defined $cur_slice
      && $slice->strand == $cur_slice->strand
      && $slice->seq_region_name eq $cur_slice->seq_region_name

      # Never combine gap positions.
      && $slice->seq_region_name !~ m/gap/i
      ) {
      if ( $pos == $cur_pos + 1 || $pos == $cur_pos - 1 ) {
        $cur_end = $pos;
        $cur_pos = $pos;
        next;
      }
    }
    if ( defined $cur_slice ) {
      push @slice_range_output, {
        slice => $cur_slice,
        start => $cur_start,
        end   => $cur_end
        };
    }
    $cur_slice = $slice;
    $cur_start = $pos;
    $cur_end   = $pos;
    $cur_pos   = $pos;
  }

  # Push out the last slice range.
  if ( defined $cur_slice ) {
    push @slice_range_output, {
      slice => $cur_slice,
      start => $cur_start,
      end   => $cur_end
      };
  }
  return \@slice_range_output;
}

sub genomic_aln_for_member {
  my $class       = shift;
  my $compara_dba = shift;
  my $ref_member  = shift;
  my $params      = shift;

  my $mirror_dba = $compara_dba;
  my $mba        = $mirror_dba->get_MemberAdaptor;

  my $tx;
  $ref_member = $mba->fetch_by_source_stable_id( undef, $ref_member->stable_id );
  $tx = $ref_member->get_Transcript;

  return $class->genomic_aln_for_transcript($compara_dba,$tx,$params);
}

sub genomic_aln_for_transcript {
  my $class       = shift;
  my $compara_dba = shift;
  my $tx  = shift;
  my $params      = shift;

  my $debug = $params->{debug};
  $debug = 1 unless ( defined $debug );

  my $quality_threshold = $params->{quality_threshold};
  $quality_threshold = $class->default_sequence_quality_threshold
    unless ( defined $quality_threshold );

  my $taxon_id_hashref = $params->{restrict_to_species};

  if (defined $taxon_id_hashref) {
    print "  restricting to species: " . join(',', keys %$taxon_id_hashref) . "\n";
  }

  my $mirror_dba = $compara_dba;
  my $as_a       = $mirror_dba->get_AlignSliceAdaptor;
  my $mlss_a     = $mirror_dba->get_MethodLinkSpeciesSetAdaptor;

  # Get the EPO MLSS object.
  my $mlss_type = $params->{mlss_type} || 'epo';
  my $mlss;
  if ( $mlss_type eq 'epo_low_coverage' ) {
    my $mlss_list = $mlss_a->fetch_all_by_method_link_type('EPO_LOW_COVERAGE');
    $mlss = @{$mlss_list}[0];
  } else {
    my $species_set = 'mammals';
    $species_set = $params->{species_set} if ( defined $params->{species_set} );
    if ($mlss_a->can('fetch_by_method_link_type_species_set_name')) {
      $mlss = $mlss_a->fetch_by_method_link_type_species_set_name( 'EPO', $species_set );
    } else {      
      my $mlss_list = $mlss_a->fetch_all_by_method_link_type('EPO');
      foreach my $cur_mlss (@$mlss_list) {
        my $name = $cur_mlss->name;
        $mlss = $cur_mlss if ($name =~ m/$species_set/gi);
      }
    }
  }
  my $type      = $mlss->method_link_type;
  my $mlss_name = $mlss->name;
  $params->{debug} = 1;
  print "Fetching genomic alignments [$mlss_name]...\n" if ( $params->{debug} );

  my $slice_a = $tx->adaptor->db->get_SliceAdaptor;

  my @exons = @{ $tx->get_all_translateable_Exons };
  my $off_phase_start = 0;
  my @alns;
  my @aa_alns;

  my $seq_hash;
  my $qual_hash;
  my $gdb_hash;

  my $cur_aln_col = 0;

  # Store all the species showing up in this alignment.
  my $all_species;
  foreach my $exon (@exons) {
    my $sub_slice = $exon->slice->sub_Slice( $exon->start, $exon->end, $exon->strand );
    my $align_slice = $as_a->fetch_by_Slice_MethodLinkSpeciesSet( $sub_slice, $mlss, 0 );
    foreach my $a_s_slice ( @{ $align_slice->get_all_Slices } ) {
      my $gdb  = $a_s_slice->genome_db;
      my $name = $gdb->name;

      # Skip if this genome is outside our taxon_id list.
      next if (defined $taxon_id_hashref && !defined $taxon_id_hashref->{$gdb->taxon_id});

      if ( $name !~ m/ancestral/gi ) {
        #print "  $name\n" if ($debug);
        $all_species->{$name} = 1;
        $seq_hash->{$name}    = '' if ( !defined $seq_hash->{$name} );
        $qual_hash->{$name}   = [] if ( !defined $qual_hash->{$name} );
      }
    }
  }

  # Now go through each exon and get the alignment (& quality scores) for each slice.
  my $exon_i = 0;
  foreach my $exon (@exons) {
    $exon_i++;
    my $slice = $exon->slice;
    print $exon_i."  ". $slice->seq_region_name . " " . $exon->stable_id." ".$exon->start . " " . $exon->end . "\n" if ($debug);

    my $start     = $exon->coding_region_start($tx);
    my $end       = $exon->coding_region_end($tx);
    my $strand    = $exon->strand;
    my $frame     = $exon->frame;
    my $phase     = $exon->phase;
    my $end_phase = $exon->end_phase;

    if ($exon_i == 1 && $exon->phase > 0) {
      print "OFF_PHASE_START!!!\n";
      $off_phase_start = 1;
      print $exon->stable_id." ".$exon->phase."\n";
      if ($exon->strand == -1) {
        $end = $end + $exon->phase;
      } else {
        $start = $start - $exon->phase;
      }
    }

    my $sub_slice = $slice->sub_Slice( $start, $end, $strand );

    # We fetch the reference-flattened sequence (expand=0) so we can
    # accurately reconstruct the gnomic coords later on.
    my $align_slice = $as_a->fetch_by_Slice_MethodLinkSpeciesSet( $sub_slice, $mlss, 0 );
    my @a_s_slices = @{ $align_slice->get_all_Slices };

    foreach my $gdb_name ( keys %$all_species ) {
      my ($a_s_slice) = grep { $_->genome_db->name eq $gdb_name } @a_s_slices;
      if ( !defined $a_s_slice ) {
        print "$gdb_name ALLGAPS\n";

        # No slices for this species found; let's fill it up with gaps and -1s.
        $seq_hash->{$gdb_name} .= '-' x $sub_slice->length;
        push @{ $qual_hash->{$gdb_name} }, ( (-1) x $sub_slice->length );
      } else {

        #    foreach my $a_s_slice ( @{ $align_slice->get_all_Slices } ) {
        my $composed_a_s_seq = '';

        my $aln_start = $a_s_slice->start;
        my $aln_end   = $a_s_slice->end;
        my $gdb       = $a_s_slice->genome_db;
        my $name      = $gdb->name;
        $gdb_hash->{$name} = $gdb;
        #printf "%-20s %-6s %-6s\n",$name,$aln_start,$aln_end if ($debug);

        # Skip if it's outside our taxon ID list.
        next if (defined $taxon_id_hashref && !defined $taxon_id_hashref->{$gdb->taxon_id});

        my $cache_object = {};

        if ( $a_s_slice->isa("Bio::EnsEMBL::Compara::AlignSlice::Slice")
          && $name !~ m/ancestral/gi ) {
          my $slice_a =
            Bio::EnsEMBL::Registry->get_adaptor( $gdb->taxon->ensembl_alias, 'core', 'slice' );
          #print "  orig position\n" if ($debug);

          my @orig_positions;
          for ( my $i = 1 ; $i <= $a_s_slice->length ; $i++ ) {
            my ( $slice, $orig_position ) = $a_s_slice->get_original_seq_region_position($i);

            #print " $i " if ($debug);
            #print "$name $i $orig_position\n" if ($name eq 'Macaca mulatta');
            push @orig_positions, [ $slice, $orig_position ];
          }
          #print "\n" if ($debug);
          my $slice_ranges = $class->combine_slice_positions_into_ranges( \@orig_positions );


          #print "  slice ranges\n" if ($debug);
          my @slices = ();
          foreach my $slice_range (@$slice_ranges) {
            my $start = $slice_range->{start};
            my $end   = $slice_range->{end};
            my $slice = $slice_range->{slice};

            if ( $end < $start ) {
              my $tmp = $end;
              $end   = $start;
              $start = $tmp;
            }
            if ( $slice->seq_region_name =~ m/gap/i ) {
              push @slices, $slice;
            } else {
              my $sub_slice =
                $slice_a->fetch_by_region( undef, $slice->seq_region_name, $start, $end,
                $slice->strand );
              push @slices, $sub_slice;
            }
          }

          my $i = 0;
          foreach my $p_slice (@slices) {
            $i++;

            my $slice_seq = $p_slice->seq;

            #$slice_seq =~ s/\./-/g; # dots to dashes.
            my $dna;
            my $qual;
            my @q_arr;
            if ( $p_slice->seq_region_name =~ m/gap/i ) {
              $dna   = $p_slice->seq;
              $qual  = join( ' ', '99' x $p_slice->length );
              @q_arr = (99) x $p_slice->length;

              #print "$name $dna $qual\n";
            } else {
              #print "  getting quality...\n" if ($debug);
              ($dna,$qual) = $class->get_bases_for_slice( $p_slice, $gdb, $cache_object );
              @q_arr = split( ' ', $qual );
            }

            $seq_hash->{$name}  = '' if ( !defined $seq_hash->{$name} );
            $qual_hash->{$name} = [] if ( !defined $qual_hash->{$name} );

            if ( $dna ne '' ) {
              die("Not the same lengths!") unless ( length($dna) == scalar(@q_arr) );

              if ( $slice_seq ne $dna ) {
                print "Seqs not the same:\n";
                print "Slice: $slice_seq " . $p_slice->name . "\n";
                print "DNA  : $dna\n";

                # Case 1: padded in beginning with Ns.
                if ( $slice_seq =~ m/^n/i && length($slice_seq) > length($dna) ) {
                  my $no_ns = $slice_seq;
                  $no_ns =~ s/^n+//i;
                  my $diff = length($slice_seq) - length($dna);
                  $dna = ( 'N' x $diff ) . $dna;
                  push @q_arr, (0) x $diff;
                  warn("Padded BEG of seq for $name with Ns!");
                }

                # Case 2: padded at end with Ns.
                if ( $slice_seq =~ m/n$/i && length($slice_seq) > length($dna) ) {
                  my $no_ns = $slice_seq;
                  $no_ns =~ s/n+$//i;
                  my $diff = length($slice_seq) - length($dna);
                  $dna = $dna . ( 'N' x $diff );
                  unshift @q_arr, (0) x $diff;
                  warn("Padded END of seq for $name with Ns!");
                }
              }
              #die("Still not the same seqs for $name (even after padding Ns)!\n")
              #  unless ( $slice_seq eq $dna );

              # We're all OK now. Add the sequence to the hash.
              $seq_hash->{$name} .= $dna;
              push @{ $qual_hash->{$name} }, @q_arr;
              $composed_a_s_seq .= $dna;
            } else {

              # If we didn't get DNA from the sequence&quality files, use the slice's dna.
              $seq_hash->{$name} .= $p_slice->seq;
              my @q = (99) x $p_slice->length;
              push @{ $qual_hash->{$name} }, @q;
              $composed_a_s_seq .= $p_slice->seq;
              
              #warn( "No DNA for $name " . $p_slice->name );
            }
          }

          # Translate dots to dashes.
          my $a_s_seq = $a_s_slice->seq;
          $composed_a_s_seq =~ s/\./-/g;
          $a_s_seq          =~ s/\./-/g;

          #        print "$composed_a_s_seq $name\n";

          # If the composed sequence is all gaps, remove it from the seq and qual hash.
          #        if ( $composed_a_s_seq ne '' && $composed_a_s_seq !~ m/[^-]/ ) {
          #          warn("Only gaps in seq for $name!");
          #          delete $seq_hash->{$name};
          #          delete $qual_hash->{$name};
          #          $composed_a_s_seq = '';
          #        }

          if ( $composed_a_s_seq ne '' && $a_s_seq ne $composed_a_s_seq ) {
            print "a_s : " . $a_s_slice->seq . "\n";
            print "cmp : " . $composed_a_s_seq . "\n";
            #die("Not the same for $name!");
          }
        }

        if ( $name =~ m/ancestral/gi ) {

          # Try to get the tree which this ancestral sequence represents...
          # See ensembl-webcode/modules/EnsEMBL/Web/Component/Compara_Alignments.pm
          # and documentation for AlignSlice.pm
          my @slice_mapper_objs = @{ $a_s_slice->get_all_Slice_Mapper_pairs };
          foreach my $obj (@slice_mapper_objs) {
            my $tree_newick = $obj->{slice}->{_tree};

            #print "Newick: $tree_newick\n";
            my $ancestral_tree = Bio::EnsEMBL::Compara::TreeUtils->from_newick($tree_newick);
            my $num_leaves     = scalar( $ancestral_tree->leaves );

            # Generate a unique name string from the subtree covered by this ancestral seq.
            my @leaf_initials = sort map { substr( $_->name, 0, 1 ) } $ancestral_tree->leaves;

            # Gotta create a dummy copy of the genome db object too.
            my $gdb = $a_s_slice->genome_db;
#            $gdb = new $gdb;
#            $gdb->name( 'Anc.' . join( '', @leaf_initials ) );
#            $a_s_slice->genome_db($gdb);
          }
        }
      }
    }

    my $sa = $align_slice->get_SimpleAlign();
    my @remove_seqs = grep {$_->id =~ m/ancestral/i} $sa->each_seq;
    foreach my $seq (@remove_seqs) {
      $sa = Bio::EnsEMBL::Compara::AlignUtils->remove_seq_from_aln($sa,$seq->id);
    }
    #Bio::EnsEMBL::Compara::AlignUtils->pretty_print($sa,{full => 1}) if ($debug);

    my $ann = $sa->annotation;
    $ann->{_genomic_coords} = {};

    # Store the genomic coordinates for each species at each position (ya this is wasteful,
    # but I think it's the only way...)
    #print "  storing coordinates...\n" if ($debug);

    my $should_store_coordinates = 0;
    if ($should_store_coordinates) {

      my $genome_db_name_counter;
      foreach my $a_s_slice ( @{ $align_slice->get_all_Slices } ) {
        next unless ( defined $a_s_slice->genome_db->taxon_id );
        my $gdb = $a_s_slice->genome_db;
        next if ( $gdb->name =~ m/ancestral/gi );
        
        my $name_in_aln = $a_s_slice->genome_db->name.($genome_db_name_counter->{$a_s_slice->genome_db->name} or "");
        if (!defined($genome_db_name_counter->{$a_s_slice->genome_db->name})) {
          $genome_db_name_counter->{$a_s_slice->genome_db->name} = 2;
        } else {
          $genome_db_name_counter->{$a_s_slice->genome_db->name}++;
        }
  
        foreach my $i ( 1 .. $sa->length ) {
          my ( $slice, $pos ) = $a_s_slice->get_original_seq_region_position($i);
  
          my $global_aln_pos = $cur_aln_col + $i;
          my $coordinate     = Bio::EnsEMBL::Mapper::Coordinate->new();
          my $char = $a_s_slice->sub_Slice($i,$i)->seq;
          #print "$global_aln_pos ".$slice->name." char[$char]\n" if ($gdb->taxon_id==9606);
          $coordinate->start($pos);
          $coordinate->end($pos);
          $coordinate->coord_system( $slice->coord_system );
          $coordinate->id( $slice->seq_region_name );
          $coordinate->strand( $slice->strand );
  
          #print $global_aln_pos." ". $coordinate->start."\n";
  
          if (defined $ann->{_genomic_coords}->{ $name_in_aln . '_' . $global_aln_pos }) {
            my $orig = $ann->{_genomic_coords}->{ $name_in_aln . '_' . $global_aln_pos };
            my $new = $coordinate;
  
            if ($new->id =~ m/gap/i) {
              # Do nothing.
            } elsif ($new->id !~ m/gap/i && $orig->id =~ m/gap/i) {
              # Replace the original with the current one.
              $ann->{_genomic_coords}->{ $name_in_aln . '_' . $global_aln_pos } = $new;
            } elsif ($new->id !~ m/gap/i && $orig->id !~ m/gap/i) {
              warn("Overlapping duplicate genomic coordinates for [$name_in_aln $global_aln_pos]!");
              my $cur = $orig;
              printf "ORIG: %s %s %d-%d\n",$cur->id,$cur->coord_system->name,$cur->start,$cur->end;
              my $cur = $new;
              printf "CURR: %s %s %d-%d\n",$cur->id,$cur->coord_system->name,$cur->start,$cur->end;
            }
          }  else {
            $ann->{_genomic_coords}->{ $name_in_aln . '_' . $global_aln_pos } = $coordinate;
          }
        }
      }
    } else {
      #warn("Skipping storing coordinates -- you won't get nuc-level coordinates now!");
    }

    #Bio::EnsEMBL::Compara::AlignUtils->pretty_print($sa,{full=>1,length=>150});
    push @alns, $sa;
    $cur_aln_col += $sa->length;
  }

  # Combine all the exonic alignments together.
  my $cdna_aln = Bio::EnsEMBL::Compara::AlignUtils->combine_alns(@alns);

  print "Storing " . scalar( keys %{ $cdna_aln->annotation->{_genomic_coords} } ) . " coords\n";

  my $extra_info;
  $extra_info->{off_phase_start} = $off_phase_start;

  my $filtered_site_total = 0;
  # Filter available quality-scored sequences at the given threshold.
  if ( $quality_threshold > 0 ) {
    print ("Filtering genomic alignment by sequence quality...\n");
    foreach my $name ( keys %$seq_hash ) {
      
      my $seq = $cdna_aln->get_seq_by_id($name);

      my $seq_str = $seq->seq;
      my $dna_seq = $seq_hash->{$name};
      my @quals   = @{ $qual_hash->{$name} };

      printf "%s %.20s\n",$name,$dna_seq;
      print join(" ", @quals)."\n";

      if ( length($dna_seq) != scalar(@quals) ) {
        print "dna[" . length($dna_seq) . " qual[" . scalar(@quals) . "]\n";
        print "dna : $dna_seq\n";
        print join( "-", @quals ) . "\n";
        die("DNA seq and qual length not the same!");
      }
      if ( length($dna_seq) != length($seq_str) ) {
        print "$seq_str $name\n";
        print "$dna_seq $name\n";
        die("DNA seq and aln seq length not the same!");
      }
      if ( length($seq_str) != scalar(@quals) ) {
        print "$seq_str $name\n";
        die("Aln seq and qual length not the same!");
      }

      # Filter this genomic alignment by seq quality.
      my $filtered =
        $class->filter_alignment_seq_by_qual_array( $seq, \@quals, $quality_threshold );
      my $n_filtered = $class->count_filtered_sites( $seq_str, $filtered );
      $filtered_site_total += $n_filtered;

      my $gdb        = $gdb_hash->{$name};
      my $short_name = $gdb->taxon->short_name;
      $extra_info->{ 'filtered_' . $short_name } = $n_filtered;
      $seq->seq($filtered);
    }
  }

  $cdna_aln = $class->dots_to_gaps($cdna_aln);

  $extra_info->{filtered_site_total} = $filtered_site_total;
  my $aa_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($cdna_aln);

  Bio::EnsEMBL::Compara::AlignUtils->pretty_print($aa_aln, {full => 1});
  Bio::EnsEMBL::Compara::AlignUtils->pretty_print($cdna_aln, {full => 1});

  return ( $cdna_aln, $aa_aln, $extra_info );
}

sub filter_alignment_seq_in_threes {
  my $class = shift;
  my $seq = shift;
  my $qual_arrayref = shift;
  my $threshold = shift;

  my $seq_str = $seq->seq;
  my @quals   = @$qual_arrayref;

  if ( length($seq_str) != scalar(@quals) ) {
    my $sl = length($seq_str);
    my $ql = scalar(@quals);

    # We often see the cds alignment missing the last stop codon, so a diff of 3 is OK.
    if ( $ql - $sl != 3 ) {
      print "seq : $seq_str\n";
      warn( "Warning: aln seq and qual array not same length (seq=$sl qual=$ql) for id "
          . $seq->id
          . "!" );
    }
  }

  my $filter_count = 0;
  for (my $i=0; $i < length($seq_str)-2; $i+= 3) {
    my $qual = @quals[$i];
    my $qual2 = @quals[$i+1];
    my $qual3 = @quals[$i+2];
    if ($qual < $threshold || $qual2 < $threshold || $qual3 < $threshold) {
      my $existing_codon = substr $seq_str, ($i), 3;
      next if ( $existing_codon =~ m/[N\-\.]/i );
      substr $seq_str, ($i), 3, 'NNN';
      $filter_count += 3;
    }
  }

  print STDERR "Filtered [$filter_count] from " . $seq->id . "\n";
  #print STDERR $seq_str."\n";
  return $seq_str;
}

sub dots_to_gaps {
  my $class = shift;
  my $aln = shift;

  foreach my $seq ($aln->each_seq) {
    my $seq_str = $seq->seq;
    $seq_str =~ s/\./-/g;
    $seq->seq($seq_str);
  }
  return $aln;
}

sub count_filtered_sites {
  my $class        = shift;
  my $orig_str     = shift;
  my $filtered_str = shift;

  my $tmp = $orig_str;
  $tmp =~ s/[^n]//gi;
  my $n_orig = length($tmp);
  $tmp = $filtered_str;
  $tmp =~ s/[^n]//gi;
  my $n_filtered = length($tmp);
  return ( $n_filtered - $n_orig );
}

sub remove_from_tree_and_aln {
  my $class = shift;
  my $tree = shift;
  my $aln = shift;
  my $id_to_remove = shift;

  $tree = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree);

  my ($leaf) = grep {$_->name eq $id_to_remove} $tree->leaves;
  die("Leaf not found") unless (defined $leaf);
  $tree = $tree->remove_nodes([$leaf]);
  $tree = $tree->minimize_tree;

  my ($seq) = grep {$_->id eq $id_to_remove} $aln->each_seq;
  $aln = Bio::EnsEMBL::Compara::AlignUtils->remove_seq_from_aln($aln, $id_to_remove);
  
  return ($tree, $aln);
}

1;
