package Bio::EnsEMBL::Compara::ComparaUtils;

use Bio::TreeIO;
use Bio::EnsEMBL::Compara::LocalMember;
use Bio::EnsEMBL::Compara::ProteinTree;
use Bio::EnsEMBL::Utils::Exception;

use Bio::EnsEMBL::Compara::AlignUtils;
use Bio::EnsEMBL::Compara::TreeUtils;

use File::Path;

#
# A grab bag of useful methods for Compara modules.
#

my $TREE = "Bio::EnsEMBL::Compara::TreeUtils";
my $ALN = "Bio::EnsEMBL::Compara::AlignUtils";
my $COMPARA = "Bio::EnsEMBL::Compara::ComparaUtils";

sub cigar_line {
  my ($class,$str) = @_;
  my $cigar_line="";
  
  $_ = $str;
  my @groups = /(-+|[^-]+)/ig;
  foreach my $token (@groups) {
    my $len = length($token);
    my $char = 'M';
    if (index($token,'-') > -1) {
      $char = 'D';
    }
    if ($len > 1) {
      $cigar_line .= $len;
    }
    $cigar_line .= $char;
  }
  return $cigar_line;
}


# Stores a given SimpleAlign (and optional cdna align) into $output_table.
sub store_SimpleAlign_into_table {
  my $class = shift;
  my $output_table = shift;
  my $tree = shift;
  my $aa_aln = shift;
  my $cdna_aln = shift;
  
  my $pta = $tree->adaptor;
  my $mba = $pta->db->get_MemberAdaptor;

  $pta->protein_tree_member($output_table);
  
  #
  # Align cigar_lines to members and store
  #
  foreach my $node (@{$tree->get_all_leaves}) {
    #print $node." ". $node->node_id." ".$node->stable_id."\n";    
    my $name =  $node->stable_id;

    if (defined $aa_aln) {
      #print "membr: ".$node->member_id."\n";
      #print "seq id: ".$node->sequence_id."\n";
      $pta->dbc->do("UPDATE member set sequence_id=NULL where member_id=".$node->member_id.";");
      #if ($node->sequence_id) {
      #  $pta->dbc->do("DELETE FROM sequence WHERE sequence_id=".$node->sequence_id.";");
      #}
      my $seq_obj = $ALN->get_seq_with_id($aa_aln,$name);
      if (!defined $seq_obj) {
	warn("No sequence with ID $name found in the alignment!") ;
	next;
      }
      my $aa_seq = $seq_obj->seq;      
      my $cigar_line = $class->cigar_line($aa_seq);
      $node->cigar_line($cigar_line);
      $aa_seq =~ s/-//g;
      $node->sequence($aa_seq);
    }
    if (defined $cdna_aln) {
      #print "CDNA ID:".$node->cdna_sequence_id."\n";
        $pta->dbc->do("UPDATE member set cdna_sequence_id=NULL where member_id=".$node->member_id.";");
      #if ($node->cdna_sequence_id) {
      #  my $cmd = "DELETE FROM sequence WHERE sequence_id=".$node->cdna_sequence_id.";";
      #  $pta->dbc->do($cmd);
      #}
      my $cdna_seq = $ALN->get_seq_with_id($cdna_aln,$name)->seq;
      bless($node,"Bio::EnsEMBL::Compara::LocalMember");
      $cdna_seq =~ s/-//g;
      $node->cdna_sequence($cdna_seq);
    }

    $node->sequence_id(0);
    $node->cdna_sequence_id(0);

    # Call the ProteinTreeAdaptor method to do the actual database dirty work.
    $pta->store($node);
    $mba->store($node);
  }

  # Clean up 'hanging' sequences not referenced by any member:
#  $pta->dbc->do(qq^delete s.* FROM sequence s LEFT JOIN member m
#    ON (m.sequence_id=s.sequence_id OR m.cdna_sequence_id=s.sequence_id) 
#    WHERE (m.sequence_id IS NULL OR m.cdna_sequence_id IS NULL);
#    ^);
}


# Reads the MCoffee scores from $mcoffee_scores file and stores them into $output_table.
sub store_MCoffee_scores_into_table {
  my $class = shift;
  my $tree = shift;
  my $mcoffee_scores = shift;  # File where the mcoffee scores are held.
  my $output_table = shift;

  my $pta = $tree->adaptor;

  #
  # Read in the scores file manually.
  #
  my %score_hash;
  my %overall_score_hash;
  my $FH = IO::File->new();
  $FH->open($mcoffee_scores) || throw("Could not open alignment scores file [$mcoffee_scores]");
  <$FH>; #skip header
  my $i=0;
  while(<$FH>) {
    $i++;
    next if ($i < 7); # skip first 7 lines.
    next if($_ =~ /^\s+/);  #skip lines that start with space
    if ($_ =~ /:/) {
      my ($id,$overall_score) = split(/:/,$_);
      $id =~ s/^\s+|\s+$//g;
      $overall_score =~ s/^\s+|\s+$//g;
      next;
    }
    chomp;
    my ($id, $align) = split;
    $score_hash{$id} ||= '';
    $score_hash{$id} .= $align;
  }
  $FH->close;
  
  foreach my $member (@{$tree->get_all_leaves}) {
    #
    # Do a manual insert into the given table.
    #
    
    my $table_name = $output_table;
    my $sth = $pta->prepare("INSERT ignore INTO $table_name 
                               (node_id,member_id,method_link_species_set_id,cigar_line)  VALUES (?,?,?,?)");
    my $score_string = $score_hash{$member->sequence_id};
    $score_string =~ s/[^\d-]/9/g;   # Convert non-digits and non-dashes into 9s. This is necessary because t_coffee leaves some leftover letters.

    $sth->execute($member->node_id,$member->member_id,$member->method_link_species_set_id,$score_string);
    $sth->finish;
  }
}


# Finds empty sequences in $aln (SimpleAlign) and removes the according sequence from $tree (NestedSet).
sub remove_empty_seq_leaves {
  my $class = shift;
  my $aln = shift;
  my $tree = shift;

  foreach my $seq ($aln->each_seq) {
    my $seq_str = $seq->seq;
    my $collapsed_str = $seq_str;
    $collapsed_str =~ s/[-NX]//gi;

    # REMOVE THE SEQUENCE IF IT ONLY HAS GAPS OR N's
    if (length($collapsed_str) == 0) {
      print "Removing empty seq from tree/aln: ".$seq->id."\n";
      $aln->remove_seq($seq);
      my $node = $tree->find_leaf_by_name($seq->id);
      if ($node) {
	$tree->remove_nodes([$node]);
      }
    }
  }
}


sub fetch_masked_aa {
  my $class = shift;
  my $aa_aln = shift;
  my $tree = shift;
  my $params = shift;
  
  return $class->fetch_masked_alignment($aa_aln,undef,$tree,$params,0);
}

# Returns a masked alignment, either amino acid or DNA-level.
sub fetch_masked_alignment
{
  my $class = shift;
  my $aa_aln = shift;
  my $cdna_aln = shift;
  my $tree = shift;
  my $params = shift;
  my $cdna_option = shift;

#  print "$tree\n";
  my $dbc = $tree->adaptor->dbc;
  my $ps_params = $class->load_params_from_param_set($dbc,$params->{'parameter_set_id'});
  $params = $class->replace_params($params,$ps_params);

  #$class->hash_print($params);

  my $aln;
  $aln = $cdna_aln if ($cdna_option);
  $aln = $aa_aln if (!$cdna_option);

  my $default_params = {
    sequence_quality_filtering => 0,
    sequence_quality_threshold => 2,
    sequence_quality_mask_character => 'X',

    alignment_score_filtering => 0,
    alignment_score_threshold => 'auto',
    alignment_score_table => 'protein_tree_member_score',
    alignment_score_mask_character => 'X',
    
    trimal_filtering => 0,
    trimal_filtering_params => '',
    trimal_mask_character => 'X',

    cdna => $cdna_option
  };
  $params = $class->replace_params($default_params,$params);


  # Mask out bits of alignment.
  if ($params->{'alignment_score_filtering'}) {
    # Load up all the site-wise alignment quality scores.
    my $pta = $tree->adaptor;
    my $use_alignment_scores = 1;
    my $hash_ref;
    foreach my $leaf (@{$tree->get_all_leaves}) {
      # Grab the score line for each leaf node.
      my $id = $leaf->stable_id; # Must be stable_id to match the aln object.
      my $member_id = $leaf->member_id;
      my $table = $params->{'alignment_score_table'};
      my $cmd = "SELECT cigar_line FROM $table where member_id=$member_id;";
      my $sth = $pta->prepare($cmd);
      $sth->execute();
      my $data = $sth->fetchrow_hashref();
      $sth->finish();
      my $cigar = $data->{'cigar_line'};
      
      my $score_sequence;
      my $aln_score_string;
      
      # Bailout early if we don't find a score.
      if (!defined $cigar || $cigar eq "") {
	if ($use_alignment_scores) {
	  print "No alignment scores for member ".$leaf->stable_id." ($member_id) found! Not using alignment scores...\n";
	  $use_alignment_scores = 0;
	}
	my $cigar_line = $leaf->cigar_line;
	my $len = length($leaf->sequence);
	my $score_string = '9' x $len;
	$aln_score_string = $ALN->combine_seq_cigar($score_string,$cigar_line);
      } else {
	$aln_score_string = $cigar;
      }
      
      my @arr = split(//,$aln_score_string);
      if ($cdna_option) {
	@arr = map { ($_ . '') x 3 } @arr;
      }
      $aln_score_string = join("",@arr);
      $hash_ref->{$id} = $aln_score_string;
    }

    # Find a reasonable threshold based on the average scores in the alignment.
    my $threshold_param = $params->{'alignment_score_threshold'};
    
    my $threshold;
    if ($threshold_param eq 'auto') {
      my %id_to_avg;
      my $score_sum = 0;
      my $score_residues = 0;
      foreach my $key (keys %{$hash_ref}) {
        my $cigar = $hash_ref->{$key};
        $cigar =~ s/\D//g;
        my @chars = split(//,$cigar);
        my $nchars = scalar @chars ;
        last if ($nchars == 0); 	# Early exit if we don't have any chars
        my $sum = 0;
        map {$score_sum += $_; $sum += $_} @chars;
        
        $score_residues += $nchars;
        my $avg = $sum / $nchars;
        $id_to_avg{$key} = $avg;
      }
      $threshold = 6;
      my $total_avg = 0;
      if (!$use_alignment_scores) {
        $threshold = 6;
      } elsif ($score_residues != 0) {
        $total_avg = $score_sum / $score_residues;
        $threshold = int($total_avg - 1.5);   # takes the integer value of the average, and subtracts 2.
        $threshold = 3 if ($threshold < 3);   # Minimum threshold of 3 seems OK.
      }
    } else {
      $threshold = $params->{'alignment_score_threshold'};
    }
    printf " -> Filtering table: %s  threshold: %d  avg: %.3f)\n",$table,$threshold,$total_avg;
    $aln = $ALN->mask_below_score($aln,$threshold,$hash_ref,$params->{'alignment_score_mask_character'});
  }

  #
  # SEQUENCE QUALITY MASKING
  #
  if ($params->{'trimal_filtering'}) {
    printf " -> Masking sequences with trimal!\n";
    $aln = $class->mask_aln_by_trimal($tree,$aln,$aa_aln,$params);
  }

  #
  # SEQUENCE QUALITY MASKING
  #
  if ($params->{'sequence_quality_filtering'}) {
    printf " -> Masking sequences at quality threshold: %d\n",$params->{'sequence_quality_threshold'};
    $aln = $class->mask_aln_by_sequence_quality($tree,$aln,$params);
  }

  $aln = $ALN->sort_by_tree($aln,$tree);
  return $aln;
}

# IMPORTANT: Always give the amino acid alignment
sub mask_aln_by_trimal {
  my $class = shift;
  my $tree = shift;
  my $aln = shift;
  my $aa_aln = shift;
  my $params = shift;

  # Write temporary alignment.
  my $dir = '/tmp/eslr_temp';
  mkpath([$dir]);
  my $node_id = $tree->node_id;
  my $aln_f = $dir."/aln{$node_id}.fa";
  Bio::EnsEMBL::Compara::AlignUtils->to_file($aa_aln,$aln_f);

  # Build a command for TrimAl.
  my $cmd = "trimal -in $aln_f -out $aln_f -colnumbering -cons 30 -gt 0.33 -gw 3";
  my $output = `$cmd`;
  print "OUTPUT:". $output."\n";

  chomp $output;
  my @cons_cols = split(/[\s,]+/,$output);
  print "@cons_cols\n";

  if ($params->{'cdna'}) {
    # use the $aa_aln for getting the columns, and multiply by 3.
    my @cdna_cons_cols;
    foreach my $i (@cons_cols) {
      my $cdna_i=int($i)*3;
      push @cdna_cons_cols,($cdna_i);
      push @cdna_cons_cols,($cdna_i+1);
      push @cdna_cons_cols,($cdna_i+2);
    }
    #print "@cons_cols\n";
    #print "@cdna_cons_cols\n";
    $aln = Bio::EnsEMBL::Compara::AlignUtils->mask_columns($aln,\@cdna_cons_cols,$params->{'trimal_mask_character'});
  } else {
    $aln = Bio::EnsEMBL::Compara::AlignUtils->mask_columns($aln,\@cons_cols,$params->{'trimal_mask_character'});
  }

  #Bio::EnsEMBL::Compara::AlignUtils->pretty_print($aln,{length=>200});
  rmtree($dir);
  return $aln;
}


sub mask_aln_by_sequence_quality {
  my $class = shift;
  my $tree = shift;
  my $aln = shift;
  my $params = shift;

  my $threshold = $params->{'sequence_quality_threshold'};
  my $cdna_option = $params->{'cdna'};
  my $mask_char = $params->{'sequence_quality_mask_character'};

  my $pta = $tree->adaptor;
  my @twox_ids = (9978,9371,9739,9478,42254,30538,
		  10020,9365,59463,9358,9813,37347,
		  9593,132908,30611,30608,9785,43179,9986,9685,9361);

  my $qual_hash_ref;
  my $sequence_quality;
  foreach my $leaf ($tree->leaves) {
    my $id = $leaf->stable_id;
    my $member_id = $leaf->member_id;
    my $sequence_id = $leaf->sequence_id;

#    if (!grep {$_ eq $leaf->taxon_id} @twox_ids) {
#      # Only mask with leaves that are 2x genomes.
#      #print $leaf->stable_id." (NOT 2x)\n";
#      
#      # Put a dummy score string into the hash.
#      my $len = length($leaf->sequence);
#      $sequence_quality = '9' x $len;
#    } else {
      #print $leaf->stable_id." (IS 2x!)\n";
      
      my $cmd = "SELECT sequence FROM sequence_quality where sequence_id = $sequence_id;";
      my $sth = $pta->prepare($cmd);
      $sth->execute();
      my $data = $sth->fetchrow_hashref();
      $sth->finish();
      $sequence_quality = $data->{'sequence'};

      if (!defined $sequence_quality || $sequence_quality eq '') {
	#print "NO seq quality: $id\n";
	$sequence_quality = '9' x length($leaf->sequence);
      }
#    }
    
    my $cigar_line = $leaf->cigar_line;
    #print $leaf->stable_id."\n";
    my $qual_cigar_line = $ALN->combine_seq_cigar($sequence_quality,$cigar_line);
    
    # Convert the protein mask into a CDNA mask by repeating each char 3 times.
    my @arr = split(//,$qual_cigar_line);
    if ($cdna_option) {
	@arr = map { ($_ . '') x 3 } @arr;
      }
    $qual_cigar_line = join("",@arr);
    $qual_hash_ref->{$id} = $qual_cigar_line;
  }

  $aln = $ALN->mask_below_score($aln,$threshold,$qual_hash_ref,$params->{'sequence_quality_mask_character'});
  return $aln;
}


# Fetches an alternate set of ProteinTree alignments from $alt_member_table.
sub fetch_alternate_tree_aln {
  my $class = shift;
  my $tree = shift;
  my $alt_member_table = shift;

  foreach my $leaf (@{$tree->get_all_leaves}) {
    # "Release" the stored / cached values for the alignment strings.
    undef $leaf->{'cdna_alignment_string'};
    undef $leaf->{'alignment_string'};

    # Grab the correct cigar line for each leaf node.
    my $id = $leaf->member_id;
    my $cmd = "SELECT cigar_line FROM $alt_member_table where member_id=$id;";
    my $sth = $pta->prepare($cmd);
    $sth->execute();
    my $data = $sth->fetchrow_hashref();
    $sth->finish();
    my $cigar = $data->{'cigar_line'};

    die "No cigar line for member $id!\n" unless ($cigar);
    $leaf->cigar_line($cigar);
  }
  return $tree;
}

# Extracts "Human Gene Subtrees" from the given ProteinTree.
sub get_human_gene_subtrees {
  my $class = shift;
  my $tree = shift;

# How to find the human gene subtrees? Go through a tree and:
# 1. Always split a tree on a duplication node IF the resulting two subtrees both contain:
#   a) Greater than 4 leaves, and
#   b) at least one human gene on each side.
# 2. After that splitting process, only include the subtrees that contain a human gene.

  sub is_duplication {
    my $node = shift;
    my $val = $node->get_tagvalue("Duplication");
    return defined($val) && $val ne '' && $val > 0;
  };

  sub does_tree_have_human_genes {
    my $node = shift;
    my @leaves = $node->leaves;
    my @hum_genes = grep {$_->taxon_id == 9606} @leaves;
    return 1 if (scalar(@hum_genes) > 0);
    return 0;
  }
  
  sub is_tree_big_enough {
    my $node = shift;
    my $node_id = $node->node_id;
    
    return (scalar($node->leaves) > 4);
    
    #my $cmd = "SELECT leaf_count($node_id);";
    #my @arr = $dbh->selectrow_array($cmd);
    #return ($arr[0] > 4);
  }

  my @root_nodes = (); 		# Root nodes which we'll add to the analysis.
  
  my @leaves = @{$tree->get_all_leaves};
  EACH_LEAF: foreach my $leaf (@leaves) {
    #print $leaf->node_id."\n";
    next unless ($leaf->taxon_id==9606);
    # Go up the tree,
    my $node = $leaf;
    my $root = $leaf->subroot;
    my $parent = undef;
    CLIMB_UP: while ($parent = $node->parent) {
      if (is_duplication($parent)) {
	# Identify the "other child."
	my $other_child;
	foreach my $child (@{$parent->children}) {
	  $other_child = $child if ($child != $node);
	}
	# If (1) other child has human genes, (2) other child is big enough, and (3) cur node is big enough, then we're done!
	#if (is_tree_big_enough($node) && is_tree_big_enough($other_child)) {
	if (does_tree_have_human_genes($other_child) && is_tree_big_enough($other_child) && is_tree_big_enough($node)) {
	  my @tmp = grep {$_->taxon_id == 9606} $node->leaves;
	  my $human_genes_subroot = scalar(@tmp);
	  printf("Root for human gene %s is %s [%s leaves out of %s for the full tree, contains %s human gene(s)]\n",
		 $leaf->stable_id,
		 $node->node_id,
		 scalar($node->leaves),
		 scalar($root->leaves),
		 $human_genes_subroot
		 );
	  push @root_nodes,$node;
	  next EACH_LEAF;
	}
      }
      
      # Climb up the tree, stop at the top.
      last CLIMB_UP if ($parent == $root);
      $node = $parent;
    }
    # If we reach here, we should be at the root node.
    printf("Using root node for human gene %s\n",
	   $leaf->stable_id
	   );
    push @root_nodes,$root;
  }

  $tree->release_tree;
  return \@root_nodes;
}


# Calculates and returns the alignment tags from a given ProteinTree.
sub create_aln_tags_from_ProteinTree {
  my $class = shift;
  my $tree = shift;

  my $tags;
  my $pta = $tree->adaptor;
  $aln_method = "[generic alignment]" unless ($aln_method);

  my $sa = $tree->get_SimpleAlign;

  # Alignment percent identity.
  my $aln_pi = $sa->average_percentage_identity;
  $tags->{'aln_percent_identity'} = $aln_pi;

  # Alignment length.
  my $aln_length = $sa->length;
  $tags->{'aln_length'} = $aln_length;

  # Alignment residue count.
  my $aln_num_residues = $sa->no_residues;
  $tags->{'aln_num_residues'} = $aln_num_residues;

  # Alignment "gappiness".
  my ($gaps,$non_gaps);
  foreach my $seq ($sa->each_seq) {
    my $str = $seq->seq;
    my $str2 = $str;

    $str =~ s/[-]//gi; # Hold all non-gaps.
    $str2 =~ s/[^-]//gi; # Only hold gaps.
    $gaps += length($str2);
    $non_gaps += length($str);
  }
  my $gappiness = sprintf("%.3f",($gaps / ($non_gaps+1)));
  $tags->{'aln_gappiness'} = $gappiness;

  return $tags;
}

# Loads up and stores the alignment tags for the ProteinTree $tree.
sub store_aln_tags {
  my $class = shift;
  my $tree = shift;

  my $tags = $class->create_aln_tags_from_ProteinTree($tree);
  $class->store_tags($tree,$tags);
}

# Stores the $tags key/val pairs into the $tree NestedSet.
sub store_tags {
  my $class = shift;
  my $tree = shift;
  my $tags = shift;

  foreach my $tag (keys %{$tags}) {
    $tree->store_tag($tag,$tags->{$tag});
  }
}


######
###### DUMPING / LOADING TREES & ALIGNMENTS
######

# Writes the given ProteinTree to the given directory. Returns the filename.
sub dump_ProteinTree_tree {
  my $class = shift;
  my $tree = shift;
  my $fastafile = shift;

  my $treeI = $class->to_TreeI($tree);
  my $treeout = Bio::TreeIO->new('-format' => 'newick',
                                 '-file'     => ">$file");
  $treeout->write_tree($tree);
  $treeout->close();

  return $file;
}

# Dumps a ProteinTree's sequences to a Fasta file.
sub dump_ProteinTree_seqs {
  my $class = shift;
  my $tree = shift;
  my $fastafile = shift;
  my $include_exon_boundaries = shift if (@_);

  $fastafile =~ s/\/\//\//g;  # converts any // in path to /
  return $fastafile if(-e $fastafile);  # Just return filename if it exists already.
  
  my $sa = $class->get_ProteinTree_seqs($tree,$include_exon_boundaries);

  open(OUTSEQ, ">$fastafile")
    or throw("Error opening $fastafile for write!");  
  foreach my $seq ($sa->each_seq) {
    print OUTSEQ ">". $seq->id. "\n" . $seq->seq . "\n";
  }
  close OUTSEQ;
  return $fastafile;
}

# Returns a ProteinTree's sequences as a SimpleAlign object (with NO gaps).
sub get_ProteinTree_seqs {
  my $class = shift;
  my $tree = shift;
  my $include_exon_boundaries = shift if (@_);

  my $sa = Bio::SimpleAlign->new();

  my $member_list = $tree->get_all_leaves;
  foreach my $member (@{$member_list}) {
    my $seq = '';
    if ($include_exon_boundaries) {
      $seq = $member->get_exon_bounded_sequence;
    } else {
      $seq = $member->sequence;
    }
    #$seq =~ s/(.{72})/$1/g;
    #chomp $seq;
    my $lseq = new Bio::LocatableSeq(-seq => $seq, -id => $member->stable_id);
    $sa->add_seq($lseq);
  }
  
  return $sa;
}


sub get_tree_and_alignment {
  my $class = shift;
  my $dba = shift;
  my $params = shift;

  my $p2 = $class->replace_params($params,{remove_subtree => 1});
  my $subtree_only;
  $subtree_only = $class->get_tree_for_comparative_analysis($dba,$p2);

  my $tree = $class->get_tree_for_comparative_analysis($dba,$params);

  print "Getting simple align...\n";
  my $sa = $tree->get_SimpleAlign(-exon_cased => $params->{'exon_cased'});
  print "  -> Done!\n";
  $sa = $class->fetch_masked_aa($sa,$tree,$params);
  
  # Now, mask out non-subtree residues if appropriate.
  if ($params->{'mask_outside_subtree'}) {
    my %leaf_ids;
    foreach my $leaf ($subtree_only->leaves) {
      next unless ($leaf->isa("Bio::EnsEMBL::Compara::Member"));
      $leaf_ids{$leaf->stable_id} = 1;
    }
    my @not_in_tree = grep {!defined $leaf_ids{$_->id}} $sa->each_seq;
    foreach my $seq (@not_in_tree) {
      #print $seq->seq."\n";
      $ALN->mask_seq_from_aln($sa,$seq->id,$params,$params->{subtree_mask_character});
    }
  }

  $sa = $ALN->sort_by_tree($sa,$tree);

  return ($tree,$sa);
}

# GJ 2009-01-15
sub get_tree_for_comparative_analysis {
  my $class = shift;
  my $dba = shift;
  my $params = shift;

  my $default_params = {
    keep_species => '',
    remove_species => ''
  };
  $params = $class->replace_params($default_params,$params);

  my $ps_params = $class->load_params_from_param_set($dba->dbc,$params->{'parameter_set_id'});
  $params = $class->replace_params($params,$ps_params);

  # Fetch the tree.
  my $node_id = $params->{'node_id'};
  $node_id = $params->{'protein_tree_id'} unless (defined $node_id);
  die ("Need a node ID for fetching a tree!") unless ($node_id);
  my $pta = $dba->get_ProteinTreeAdaptor;

  # GJ 2009-09-18
  if ($params->{'input_table'}) {
    $pta->protein_tree_member($params->{'input_table'});
    $pta->protein_tree_score($params->{'input_table'}.'_score');
  }

  my $tree = $pta->fetch_node_by_node_id($node_id);

  my $keep_species = $params->{'keep_species'};
  my $remove_species = $params->{'remove_species'};
  my @all_ids;
  my %keep_hash;
  my %remove_hash;

  # Collect all taxon IDs within the tree.
  @all_ids = $TREE->get_species_in_tree($tree);

  # Remove nodes not within our desired taxonomic subtree.
  if ($keep_species ne '') {
    # Find a list of all species in the tree, and add any species NOT in the keepers list to the remove list.
    my @ks = split(",",$keep_species);
    map {$keep_hash{$_}=1} @ks;
    foreach my $tax_id (@all_ids) {
      if (!exists $keep_hash{$tax_id}) {
	$remove_hash{$tax_id} = 1;
      }
    }
  }

  if ($remove_species ne '') {
    my @rs = split(",",$remove_species);
    map {$remove_hash{$_}=1} @rs;
  }

  # Do the actual removal of everything in the remove_hash.
  my @taxon_ids = keys %remove_hash;
  $tree = $TREE->remove_members_by_taxon_id($tree,\@taxon_ids);

  # IMPORTANT: Re-root the tree so we get rid of parents above this one.
  $tree->re_root;
  return $tree;
}


my $qual_base = "/nfs/users/nfs_g/gj1/scratch/2x_quality/assemblies/";

my $ens_to_index = {};

sub get_quality_string_for_member {
  my $class = shift;
  my $member = shift;

  # GJ 2009-07-01
  my $gdb = $member->genome_db;
  my $name = $gdb->name;
  $name =~ s/ /_/g;
  my $index_file = $qual_base . $name . ".quals.fa";

  print "Quals index file: ".$index_file."\n";
  if (!-e $index_file) {
    return undef;
  }

  my $tx = $member->transcript;
  print $tx->stable_id."\n";
  print "STRAND: ".$member->gene->strand."\n";
  return if (!defined $tx || !defined $tx->adaptor);
  my $adaptor = $tx->adaptor;
  my $db = $adaptor->db;

  my $if_qual = $ens_to_index->{$index_file};
  if (!defined $if_qual) {
    use Bio::Greg::IndexedFasta;
    $if_qual = Bio::Greg::IndexedFasta->new();
    $if_qual->load_indexed_fasta($index_file);
    $ens_to_index->{$index_file} = $if_qual;
  }

  my $base_file = $qual_base.$name.".bases.fa";
  my $if_base = undef;
  if (-e $base_file) {
    print "Bases index file: ".$base_file."\n";
    $if_base = Bio::Greg::IndexedFasta->new();
    $if_base->load_indexed_fasta($base_file);
  }
    
  my $ens_dna = "";
  my $genome_dna = "";
  my @whole_qual_array = ();

  my $obj_seq = $member->sequence_exon_bounded;
  my @pep_exons = split(/[obj]/,$obj_seq);
  
  foreach my $exon (@{$tx->get_all_translateable_Exons}) {
    my $pep_exon = shift @pep_exons;

    my $length = $exon->length;
    my $slice = $exon->slice;
    $exon = $exon->transfer($tx->slice);
    if (!defined $exon) {
      print "Error transferring exon!\n";
      $ens_dna .= "N" x $length;
      next;
    }
    my $dna_seq = $exon->seq->seq;
    $ens_dna .= $dna_seq;
    
    $exon = $exon->transform("contig");
    if (defined $exon) {
      my $contig_name = $exon->slice->seq_region_name;
      my $start = $exon->start-1;
      my $end = $exon->end;
      my $len = $end - $start;

      print "e: ".$exon->seq->seq."\n";

      #print "$contig_name $start-$end\n";
      #my $qual_str = $if_qual->get_sequence($contig_name);
      #my @quals = split(" ",$qual_str);
      #my @qual_slice = @quals[$start .. $end-1];
      #push @whole_qual_array,@qual_slice;

      my $qual_str = $if_qual->get_sequence_region($contig_name,$start,$end);
      if (!defined $qual_str) {
	return undef;
      }
      my @quals = split(" ",$qual_str);
      @quals = reverse(@quals) if ($tx->strand == -1);
      push @whole_qual_array,@quals;

      if (defined $if_base) {
	my $base_str = $if_base->get_sequence($contig_name);
	my @bases = split("",$base_str);
	my @bases_slice = @bases[$start .. $end-1];
	my $bases_seq = join("",@bases_slice);
	if ($tx->strand == -1) {
	  $bases_seq = new Bio::PrimarySeq(-seq=>$bases_seq)->revcom()->seq;
	}
	$genome_dna .= $bases_seq;
	print "f: ".$bases_seq."\n";
      }
    } else {
      $genome_dna .= "N" x $length;
      push @whole_qual_array,(0) x $length;
    }
  }

  if (defined $if_base) {
    my $ens_seq = new Bio::PrimarySeq(-seq => $ens_dna);
    $genome_seq = new Bio::PrimarySeq(-seq => $genome_dna);
    printf "%s \n",$ens_seq->seq;
    printf "%s \n",$genome_seq->seq;
    printf "%s \n",join(" ",@whole_qual_array);
    printf "%s \n",$ens_seq->translate->seq;
    printf "%s \n",$genome_seq->translate->seq;
  }
  
  my @avg_array = ();
  for (my $i=0; $i < scalar(@whole_qual_array)-2; $i+= 3) {
    my $sum = $whole_qual_array[$i] + $whole_qual_array[$i+1] + $whole_qual_array[$i+2];
    
    my $min = 100;
    map {$min = $_ if ($_ < $min);} @whole_qual_array[$i .. $i+2];
    #push @avg_array,$sum/3;
    push @avg_array,$min;
  }
  my @mapped_array = map {int($_/10)} @avg_array;
  my $qual_str = join("",@mapped_array);
  print $qual_str."\n";
  return $qual_str;
}



# Returns the NCBI taxnomy of Ensembl genomes below a given taxonomic clade.
sub get_genome_taxonomy_below_level {
  my $class = shift;
  my $dba = shift;
  my $root_taxon_id = shift;

  my @gdbs = $class->get_all_genomes($dba);
  my @ncbi_ids = map {$_->taxon->taxon_id} @gdbs;

  my $taxon_a = $dba->get_NCBITaxonAdaptor;

  # Try first with a taxon label.
  my $root = $taxon_a->fetch_node_by_name($root_taxon_id);
  if (!defined $root) {
    $root = $taxon_a->fetch_node_by_taxon_id($root_taxon_id);
  }
  
  # Collect all genome_db leaves, plus their internal lineages, into an array.
  my %keepers;
  $keepers{$root->node_id} = $root;
  foreach my $gdb (@gdbs) {
    my $tx = $gdb->taxon;
    my $tx_id = $tx->taxon_id;
    
    if (!$TREE->has_ancestor_node_id($tx,$root)) {
      #print "not below tax level!!!\n";
      next;
    }
    
    my $node = $tx;
    while (defined $node) {
      $keepers{$node->node_id} = $node;
      $node = $node->parent;
    }
  }

  my @nodes = values %keepers;
  #print "Size: ".scalar(@nodes)."\n";

  my $new_tree = $taxon_a->_build_tree_from_nodes(\@nodes);
  $new_tree = $new_tree->copy()->minimize_tree;
  return $new_tree;
}


# Returns the NCBI taxonomy tree of Ensembl genomes in NHX format.
sub get_genome_tree_nhx {
  my $class = shift;
  my $dba = shift;
  my $params = shift;

  my @gdbs = Bio::EnsEMBL::Compara::ComparaUtils->get_all_genomes($dba);
  my $species_tree = $class->get_genome_taxonomy_below_level($dba,1);

  my $labels_option = $params->{'labels'};
  my $include_imgs = $params->{'images'};

  @gdbs = sort {$a->taxon->binomial cmp $b->taxon->binomial} @gdbs;
  foreach my $gdb (@gdbs) {
    my $tx = $gdb->taxon;
    my $tx_id = $tx->taxon_id;

    my $sql = "select stable_id from member where taxon_id=$tx_id and source_name='ENSEMBLPEP' limit 1;";
    my $sh = $dba->dbc->prepare($sql);
    $sh->execute();
    my @vals = $sh->fetchrow_array();
    my $ensp = $vals[0];

    my $leaf = $species_tree->find_node_by_node_id($tx_id);
    next unless $leaf;

    $leaf->add_tag("bcol","gray");
    if ($include_imgs) {
      my $underbar_species = $tx->binomial;
      $underbar_species =~ s/ /_/g;
      $leaf->add_tag("img","http://www.ensembl.org/img/species/pic_${underbar_species}.png");
    }
    if ($labels_option eq 'mnemonics') {
      @vals = ($ensp,$tx_id,$tx->ensembl_alias,$tx->binomial,$tx->short_name);
      my $pretty_str = sprintf("%-20s  %6s  %-22s  %-30s  %-8s",@vals);
      $leaf->name($pretty_str);
    } else {
      $leaf->name($tx->ensembl_alias);
    }
  }
  return $species_tree->nhx_format;
}

sub get_all_genomes {
  my $class = shift;
  my $dba = shift;

  my $gda = $dba->get_GenomeDBAdaptor();
  my $all_dbs = $gda->fetch_all();
  return @{$all_dbs};
}

sub get_2x_genomes {
  my $class = shift;
  my $dba = shift;

  # From http://listserver.ebi.ac.uk/mailing-lists-archives/ensembl-dev/msg04319.html:
#    The 2X genomes to be included are:
#- elephant (Loxondonta africana)
#- armadillo (Dasypus novemcinctus
#- tenrec (Echinops telfairi)
#- rabbit (Oryctolagus cuniculus)
#- guinea pig (Cavia porcelus)
#- hedgehog (Erinaceus europaeus)
#- shrew (Sorex araneus)
#- microbat (Myotis lucifugus)
#- tree shrew (Tupaia belangeri)
#- squirrel (Spermophilus tridecemlineatus)
#- bushbaby (Otolemur garnetii)
#- pika (Ochotona princeps)
#- mouse lemur (Microcebus murinus)
#- cat (felis catus)
#- megabat (Pteropus vampyrus)
#- dolphin (Tursiops truncatus)
#- alpaca (Vicugna pacos)
#- kangaroo rat (Dipodomys ordii)
#- rock hyrax (Procavia capensis)
#- tarsier (Tarsius syrichta)
#- gorilla (Gorilla gorilla)
#- sloth (Choloepus hoffmanni)

  my @taxon_names = ("loxodonta africana", 
		     "dasypus novemcinctus", 
		     "echinops telfairi",
		     "oryctolagus cuniculus",
		     "erinaceus europaeus",
		     "sorex araneus",
		     "myotis lucifugus",
		     "tupaia belangeri",
		     "spermophilus tridecemlineatus",
		     "otolemur garnettii",
		     "ochotona princeps",
		     "microcebus murinus",
		     "felis catus",
		     "pteropus vampyrus",
		     "tursiops truncatus",
		     "vicugna pacos",
		     "dipodomys ordii",
		     "procavia capensis",
		     "tarsius syrichta",
		     "gorilla gorilla",
		     "choloepus hoffmanni"
		     );
  
  my $gda = $dba->get_GenomeDBAdaptor();

  my $all_dbs = $gda->fetch_all();
  my @twox_dbs;
  print "Fetching 2x genomes. Total dbs: " . scalar @{$all_dbs} . "\n";
  foreach my $gdb (@{$all_dbs}) {
    my $name = $gdb->name;

    #print "DB: $name\n";
    if (grep(/$name/i,@taxon_names)) {
      #print "2x found!! $name\n";
      @taxon_names = grep(! /$name/i,@taxon_names);
      push @twox_dbs, $gdb;
    }
  }
  
  print "2x DBs missing from Compara:\n  ";
  print join("\n  ",@taxon_names)."\n";
  return @twox_dbs;
}

sub remove_2x_offlimits_from_tree {
  my $class = shift;
  my $tree = shift;
  my $dba = shift;

  print "Removing 2x off-limits genomes from tree.\n";
  print "  Before:" . scalar(@{$tree->get_all_leaves}) . " leaves\n";
  
  my @off_limits_tax_ids = (9593, # Gorilla
			    9600  # Orangutan
			    );
  # Early exit if it's a single-leaf tree:
  return $tree if (scalar $tree->leaves == 1);

  my @remove_nodes = ();
  foreach my $leaf ($tree->leaves) {
    if (grep {$leaf->taxon_id == $_} @off_limits_tax_ids) {
      printf("  -> Removing off-limits leaf %s, taxon %s\n",
	     $leaf->stable_id,
	     $leaf->taxon_id
	     );
      push @remove_nodes,$leaf;
    }
  }

  foreach my $node (@remove_nodes) {
    $TREE->delete_lineage($tree,$node);
    $tree->minimize_tree();
  }
  $tree->minimize_tree();

  print "  After:" . scalar(@{$tree->get_all_leaves}) . " leaves\n";
  return $tree;
}


# Removes all 2x genomes from the given ProteinTree.
sub remove_2x_genomes_from_tree {
  my $class = shift;
  my $tree = shift;
  my $dba = shift;

  print "Removing 2x genomes from tree.\n";
  print "Before: " . scalar(@{$tree->get_all_leaves}) . " leaves\n";

  my @twox_gdbs = $class->get_2x_genomes($dba);
  my @twox_names = map {$_->name} @twox_gdbs;
  my @remove_nodes = ();
  foreach my $leaf (@{$tree->get_all_leaves}) {
    my $binomial = $leaf->taxon->binomial;
    if (grep {/$binomial/i} @twox_names) {
      print "  -> Deleting leaf: ".$leaf->stable_id."\n";
      push @remove_nodes, $leaf;
    }
  }
  foreach my $node (@remove_nodes) {
    $TREE->delete_lineage($tree,$node);
    $tree->minimize_tree();
  }
  $tree->minimize_tree();

  print "After: " . scalar(@{$tree->get_all_leaves}) . " leaves\n";
  return $tree;
}

# Return a hashref object from a parameter string.
sub load_params_from_string {
  my $class        = shift;
  my $param_object = shift;
  my $param_string = shift;
  my $debug = shift;
  
  return $param_object unless(defined $param_string);
  print("parsing parameter string : $param_string\n") if($debug);
  
  my $tmp_params = eval($param_string);
  return $param_object unless($tmp_params);
  
  foreach my $key (keys %$tmp_params) {
    print("  $key\t=>\t", $tmp_params->{$key}, "\n") if ($debug);
    $param_object->{$key} = $tmp_params->{$key};
  }

  return $param_object;
}

# Create a string from a hashref.
sub hash_to_string {
  my $class = shift;
  my $hashref = shift;

  my %hash = %{$hashref};

  my $str = "{";
  my @keyvals = ();
  foreach my $key (sort keys %hash) {
    push(@keyvals, "'$key'" . "=>" . "'" . $hash{$key} . "'");
  }

  $str .= join(",",@keyvals);
  $str .= "}";
  return $str;
}

sub hash_print {
  my $class = shift;
  my $hashref = shift;

  print "{\n";
  foreach my $key (sort keys %{$hashref}) {
    printf("    %-40.40s => %-40.40s\n",$key,$hashref->{$key});
  }
  print "}\n";
}

# Eval a hashref from a string.
sub string_to_hash {
  my $class = shift;
  my $string = shift;

  my $hash = eval $string;
  return $hash if (defined $hash);

  return {};
}

sub load_params_from_tag {
  my $class = shift;
  my $tree = shift;
  my $tag = shift;

  my $tag_string = $tree->get_tagvalue($tag);
  return $class->string_to_hash($tag_string);
}

sub load_params_from_param_set {
  my $class = shift;
  my $dbc = shift;
  my $param_set_id = shift;

  return {} unless (defined $param_set_id);

  my $params;

  my $cmd = qq^SELECT parameter_value FROM parameter_set WHERE parameter_set_id=$param_set_id AND parameter_name="params";  ^;
  my $sth = $dbc->prepare($cmd);
  $sth->execute();
  my $row = "";
  while (@row = $sth->fetchrow_array) {
    $params = eval($row[0]);
  }
  $sth->finish;

  my $cmd = qq^SELECT parameter_value FROM parameter_set WHERE parameter_set_id=$param_set_id AND parameter_name="name";  ^;
  $sth = $dbc->prepare($cmd);
  $sth->execute();
  my $row = "";
  while (@row = $sth->fetchrow_array) {
    $params->{parameter_set_name} = $row[0];
  }
  $sth->finish;

  return $params;
}

sub replace_params {
  my $class = shift;
  my $paramsA = shift;
  my $paramsB = shift;

  my $new_params = {};

  foreach my $key (keys %{$paramsA}) {
    $new_params->{$key} = $paramsA->{$key};
  }
  foreach my $key (keys %{$paramsB}) {
    $new_params->{$key} = $paramsB->{$key};
  }

  return $new_params;
}

# Loads a given tree and alignment (with optional CDNA alignment too) into a ProteinTree object, ready to be stored in the Compara database.
# @created GJ 2009-02-17
sub load_ProteinTree_from_files {
  my $class = shift;
  my $treeF = shift;		# Newick formatted string or file.
  my $alnF = shift;		# Fasta formatted string or file, main or protein alignment
  my $cdnaAlnF = shift;		# (optional) Fasta formatted string or file, secondary or CDNA alignment.

  my $treeI;
  my $sa;
  my $cdna_sa;

  $treeI = $TREE->to_treeI($treeF);
  $sa = $ALN->to_aln($alnF);
  $cdna_sa = $ALN->to_aln($cdnaAlnF) if (defined $cdnaAlnF);


  throw("Error loading tree $treeF") unless (defined $treeI);
  throw("Error loading alignment $alnF") unless (defined $sa);
  
  my $tree = $TREE->from_treeI($treeI);
  # Re-dress the NestedSet as a ProteinTree.
  bless($tree,'Bio::EnsEMBL::Compara::ProteinTree');

  foreach my $leaf ($tree->leaves) {
    my $id = $leaf->name;
    my $seq = $ALN->get_seq_with_id($sa,$id);
    my $cdna_seq = $ALN->get_seq_with_id($cdna_sa,$id) if (defined $cdna_sa);

    if (defined $seq) {
      # Calculate and set the cigar line.
      my $cigar = $class->cigar_line($seq->seq);
      $leaf->cigar_line($cigar);
      
      # Remove gaps from sequence and set it.
      my $seqstr = $seq->seq;
      $seqstr =~ s/-//g;
      $leaf->sequence($seqstr);
    }
    
    if (defined $cdna_seq) {
      bless($leaf,"Bio::EnsEMBL::Compara::LocalMember");
      # Set the CDNA sequence if necessary.
      my $cdna_seqstr = $cdna_seq->seq;
      $cdna_seqstr =~ s/-//g;
      $leaf->cdna_sequence($cdna_seqstr);
      $leaf->cdna_sequence_id(0);
    }
  }
  return $tree;
}

1;
