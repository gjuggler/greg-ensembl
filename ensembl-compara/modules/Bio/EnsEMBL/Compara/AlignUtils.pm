package Bio::EnsEMBL::Compara::AlignUtils;

use Bio::AlignIO;
use Bio::EnsEMBL::Compara::LocalMember;
use File::Path;

#
# A grab bag of useful methods for tree manipulations.
#

my $TREE = "Bio::EnsEMBL::Compara::TreeUtils";
my $ALN = "Bio::EnsEMBL::Compara::AlignUtils";
my $COMPARA = "Bio::EnsEMBL::Compara::ComparaUtils";


sub to_aln {
  my $class = shift;
  my $alnF = shift;

  if (ref $alnF) {
    return $alnF;
  } elsif (-e $alnF) {
    return $class->from_file($alnF);
  } else {
    return $class->from_string($alnF);
  }
}

sub from_file {
  my $class = shift;
  my $alnF = shift;
  my $in = Bio::AlignIO->new(-file => $alnF,
			     -format => "fasta");
  my $sa = $in->next_aln();
  $in->close();
  return $sa;
}

# Returns a SimpleAlign from a String.
sub from_string {
  my $class = shift;
  my $alnF = shift;

  open(my $fake_fh, "+<", \$alnF);
  my $in = Bio::AlignIO->new(-fh => $fake_fh,
			  -format => "fasta");
  my $sa = $in->next_aln();
  $in->close();
  return $sa;
}

sub to_file {
  my $class = shift;
  my $aln = shift;
  my $alnF = shift;

  open(OUT,">$alnF");
  foreach my $seq ($aln->each_seq) {
    print OUT ">".$seq->id."\n";
    print OUT $seq->seq."\n";
  }
  close(OUT);
  return $alnF;
}

# Counts the number of sequences that are neither gapped nor masked at the given column.
# @created GJ 2009-01-09
# @bugfix GJ 2009-03-20
sub get_nongaps_at_column {
  my $class = shift;
  my $aln = shift;
  my $pos = shift;
  
  my $nongap_count = 0;
  foreach my $seq ($aln->each_seq) {
    my $residue = $seq->subseq($pos,$pos);
    #print $residue;
    $nongap_count++ if ($residue !~ /[-x]/i);
  }
  return $nongap_count;
}

sub get_column_string {
  my $class = shift;
  my $aln = shift;
  my $pos = shift;

  my $str = "";
  foreach my $seq ($aln->each_seq) {
    my $residue = $seq->subseq($pos,$pos);
    $str .= $residue;
  }
  return $str;
}


# Retrieves a sequence from a Bio::AlignI object based on the given ID.
sub get_seq_with_id {
  my $class = shift;
  my $aln = shift;
  my $id = shift;
  
  foreach my $seq ($aln->each_seq) {
    return $seq if ($seq->id eq $id);
  }
  
  print("aln_get_seq_with_id: No sequence with id $id found!\n");
  return undef;
}

sub translate {
  my $class = shift;
  my $aln = shift;

  my $sa = new Bio::SimpleAlign;
  foreach my $seq ($aln->each_seq) {
    $sa->add_seq($seq->translate());
  }
  return $sa;
}


# Generates an Ensembl cigar line (MD-line) for each sequence in the alignment.
# Returns a hashref of sequences indexed by IDs.
sub cigar_lines {
  my $class = shift;
  my $aln = shift;
  
  my $lines;
  foreach my $pos (1..$aln->no_sequences) {   # $aln->no_sequences... what an awful method name!
    my $seq = $aln->get_seq_by_pos($pos);
    $lines->{$seq->display_id} = $COMPARA->cigar_line($seq->seq);
  }
  return $lines;
}


=head2 flatten_to_sequence
 Title     : flatten_to_sequence
 Usage     : $aln->flatten_to_sequence($my_favorite_sequence_index)
 Function  : Flattens the alignment (removes all insertions) relative to the sequence at the given index.
 Returns   : a new Bio::SimpleAlign object, flattened appropriately.
 Argument  : The index of the sequence against which to flatten the alignment. (Use get_pos_by_id to find a sequence index by name)
=cut
  
sub flatten_to_sequence {
  my $class = shift;
  my $aln = shift;
  my $pos = shift;
  
  my $ref_seq = $aln->get_seq_by_pos($pos);
  my $display_id = $ref_seq->display_id;
  my $seq_str = $ref_seq->seq;
  
  # Create a duplicate of this alignment.
  #my $new_aln = new $aln;
  #foreach my $seq ($aln->each_seq) {
  #	my $new_seq = new Bio::LocatableSeq(-seq => $seq->seq, -id =>$seq->display_id);
  #	$new_aln->add_seq($new_seq);
  #}
  
  # Go through the reference sequence and remove gapped columns in the new alignment.
  my $gap_char = $aln->gap_char();
  
  #FIXME: This is less efficient than it would be if we grouped the deletions into gaps. Sue me.
  my @remove_cols;
  while ($seq_str =~ m/[$gap_char]/g) {
    print pos($seq_str)."\n";
    my @start_end = (pos($seq_str)-1,pos($seq_str)-1);
    push @remove_cols, \@start_end;
  }
  my $new_aln = $aln->_remove_columns_by_num(\@remove_cols);
  return $new_aln;
}


# GJ 2009-03-17
# Returns a map with key = (new column #) and val = (old column #)
sub remove_blank_columns {
  my $class = shift;
  my $aln = shift;

  my $new_to_old = {};
  my $old_to_new = {};
  my @cols_to_remove = ();
  my $j = 1;
  for (my $i=1; $i <= $aln->length; $i++) {
    my $col = $class->get_column_string($aln,$i);
    if ($col =~ m/^[-]+$/) { # Regex for a blank line.
      # I guess the remove_columns_by_num method is zero-based? GJ 2009-03-19
      push @cols_to_remove, [$i-1,$i-1];
      $old_to_new->{$i} = -1;
    } else {
      $new_to_old->{$j} = $i;
      $old_to_new->{$i} = $j;
      $j++;
    }
  }
  $aln = $aln->_remove_columns_by_num(\@cols_to_remove);
  return [$aln,$new_to_old,$old_to_new];
}

# Removes blank columns in threes, and returns a codon-wise column mapping.
sub remove_blank_columns_in_threes {
  my $class = shift;
  my $aln = shift;

  my @cols_to_remove = ();
  for (my $i=1; $i <= $aln->length-2; $i+= 3) {
    my $col = $class->get_column_string($aln,$i);
    my $col2 = $class->get_column_string($aln,$i+1);
    my $col3 = $class->get_column_string($aln,$i+2);
    if ($col =~ m/^[-]+$/ && $col2 =~ m/^[-]+$/ && $col3 =~ m/^[-]+$/ ) { # Regex for a blank line.
      push @cols_to_remove, [$i-1,$i-1+2];
    }
  }
  $aln = $aln->_remove_columns_by_num(\@cols_to_remove);
  return $aln;
}

=head2 sort_by_tree
  
  Title     : sort_by_tree
  Usage     : Bio::EnsEMBL::Compara::AlnUtils->sort_by_tree($aln,$tree);
  Function  : Given a Bio::AlignI object and a Bio::Tree::TreeI object, creates a new AlignI sorted by the tree.
  Returns   : a new Bio::SimpleAlign object, sorted appropriately.
  Argument  : ($aln,$tree) - a Bio::AlignI object, and a Bio::TreeI or Bio::EnsEMBL::Compara::ProteinTree object

=cut
sub sort_by_tree {
  my $class = shift;
  my $aln = shift;
  my $tree = shift;
  
  if (! ref $tree) {
    # Check that the file exists, and try to load it.
    my $infile = $tree;
    $aln->throw("Tree file not found: $infile") unless (-e $infile);
    my $temp = Bio::TreeIO->new('-file' => $infile,
				'-format' => 'newick');
    $tree = $temp->next_tree;
    $temp->close();
  } else {
    # If we're given an object hashref, we assume it's a TreeI object.
    if ($tree->isa('Bio::EnsEMBL::Compara::ProteinTree')) {
      #$tree = $TREE->to_treeI($tree);
      my $new_aln = $aln->new;
      foreach my $leaf ($tree->leaves) {
        my $name = $leaf->stable_id;
        my $seq = $class->get_seq_with_id($aln,$name);
        $new_aln->add_seq($seq);
      }
      return $new_aln;
    } elsif (! $tree->isa('Bio::Tree::TreeI')) {
      # Try converting ProteinTree to TreeI.
      warn "sort_by_tree: given a hashref that is not a TreeI object!\n";
    }
  }
  
  my $new_aln = $aln->new;
  foreach my $leaf ($tree->get_leaf_nodes) {
    my $name = $leaf->id;
    my $seq = $class->get_seq_with_id($aln,$name);
    $new_aln->add_seq($seq);
  }
  return $new_aln;
}

# Takes zero-based column coordinates to keep.
sub mask_columns {
  my $class = shift;
  my $aln = shift;
  my $columns_ref = shift;
  my $mask_char = shift;

  my %keeper_sites;
  map {$keeper_sites{int($_)+1}=1} @{$columns_ref};

  foreach my $i (1 .. $aln->length) {
    if (!exists $keeper_sites{$i}) {
      foreach my $seq ($aln->each_seq) {
	my $seq_str = $seq->seq;
	substr($seq_str,$i-1,1) = $mask_char unless (substr($seq_str,$i-1,1) eq '-');
        $seq->seq($seq_str);
      }
    }
  }
  return $aln;
}

sub prank_filter {
  my $class = shift;
  my $aln = shift;
  my $tree = shift;
  my $params = shift;

  my $threshold = $params->{'prank_filtering_threshold'} || 5;
  my $char = $params->{'prank_mask_character'} || 'X';

  my ($scores,$blocks) = $class->get_prank_filter_matrices($aln,$tree,$params);
  my $filtered_aln = $class->mask_below_score($aln,$threshold,$scores,$char);

  return $filtered_aln;
}

sub get_prank_filter_matrices {
  my $class = shift;
  my $aln = shift;
  my $tree = shift;
  my $params = shift;

  my $dna_aln = $tree->get_SimpleAlign(-cdna => 1);

  my $node_id = $tree->node_id;
  my $mask_char = $params->{'prank_mask_character'};

  my $dir = "/tmp/prank_temp";
  mkpath([$dir]);
  my $aln_f = $dir."/aln_${node_id}.fasta";
  my $tree_f = $dir."/tree_${node_id}.nh";
  my $out_f = $dir."/aln_filtered_${node_id}";
  my $xml_f = $out_f.".0.xml";

  # Output tree and alignment.
  $class->to_file($aln,$aln_f);
  Bio::EnsEMBL::Compara::TreeUtils->to_file($tree,$tree_f);

  my $cmd = qq^prank_fix -d=$aln_f -t=$tree_f -e -o=$out_f^;
#  if (!-e $xml_f) {
    system($cmd);
#  }

  use XML::LibXML;
  use Bio::Greg::Node;

  # Grab information from Prank's XML output.
  my $parser = XML::LibXML->new();
  $xml_tree = $parser->parse_file($xml_f);
  my $root = $xml_tree->getDocumentElement;

  my $newick = ${$root->getElementsByTagName('newick')}[0]->getFirstChild->getData;
  print $newick."\n";
  my $rootNode = Bio::Greg::Node->new();
  $rootNode = $rootNode->parseTree($newick);

  my %nameToId;
  my %idToName;
  my %seqsByName;
  foreach my $lid (@{$root->getElementsByTagName('leaf')}) {
    my $tree_name = $lid->getAttribute('id');
    $nameToId{$lid->getAttribute('name')} = 
    $idToName{$lid->getAttribute('id')} = $lid->getAttribute('name');
    my $seq = $lid->findvalue('sequence');
    $seq =~ s/\s//g;
    $seqsByName{$lid->getAttribute('name')} = $seq;
  }
  my %postprob;
  foreach my $nid (@{$root->getElementsByTagName('node')}) {
    my $node = $nid->getAttribute('id');
    foreach my $pid (@{$nid->getElementsByTagName('probability')}) {
      my $prob = $pid->getAttribute('id');
      my $data = $pid->getFirstChild->getData;
      $data =~ s/\s//g;
      $postprob{$node}{$prob} = $data;
    }
  }
  my %nameToState;
  foreach my $sid (${$root->getElementsByTagName('model')}[0]->getElementsByTagName('probability')) {
    $nameToState{$sid->getAttribute('name')} = $sid->getAttribute('id');
  }

  # Create a stored 'leaf name string' for each XML node.
  my $leaf_names_to_xml;
  my @nodes = $rootNode->nodes();
  foreach my $node (@nodes) {
    next if ($node->isLeaf);
    my @nms = $node->leafNames;
    @nms = map {$idToName{$_}} @nms;
    @nms = sort @nms;
    my $leaf_names = join(" ",@nms);
    $leaf_names_to_xml->{$leaf_names} = $node->name;
  }

  my $leaf_names_to_node;
  my $node_to_xml;
  foreach my $node ($tree->nodes) {
    next if ($node->is_leaf);

    my @nms = map {$_->name} $node->leaves;
    @nms = sort @nms;
    my $leaf_names = join(" ",@nms);
    $leaf_names_to_node->{$leaf_names} = $node->name;
    if ($leaf_names_to_xml->{$leaf_names}) {
      $node_to_xml->{$node} = $leaf_names_to_xml->{$leaf_names};
    }
  }
  
  my %pp_hash;
  my @nodes = $rootNode->nodes();
  foreach my $node (@nodes) {
    next if ($node->isLeaf);
    
    my @score = split(/,/,$postprob{$node->name}{$nameToState{'postprob'}});
    $pp_hash->{$node->name} = ();
    for (my $i=0; $i < scalar(@score); $i++) {
      $pp_hash->{$node->name}[$i] = $score[$i];
    }
  }
  
  my $leaf_scores;
  foreach my $leaf ($tree->leaves) {
    my $total_dist = $leaf->distance_to_root;
    my $aln_len = $aln->length;
    
    my $aln_string = $leaf->alignment_string;
    my @aln_arr = split("",$aln_string);

    sub get_other_child {
      my $parent = shift;
      my $node = shift;

      foreach my $child (@{$parent->children}) {
        return $child if ($node != $child);
      }
      return undef;
    }

    my @scores;
    my $score_string = "";
    for (my $i=0; $i < $aln_len; $i++) {
      if ($aln_arr[$i] eq '-') {
        $score_string .= '-';
        next;
      }

      my $total_prob = 0;
      my $node = $leaf;
      while (my $parent = $node->parent) {
        my $bl = $node->distance_to_parent;
        my $xml_node = $node_to_xml->{$parent};
        my $post_prob = $pp_hash->{$xml_node}[$i];
        if ($post_prob != -1) {
          # We've got nucleotides aligned here. Sum up according to our rules.
          
          # Find the fraction of branch length encompassed by our node and the
          # other node being aligned. If we have more branch length, adjust the weights.
          my $total_bl = Bio::EnsEMBL::Compara::TreeUtils->total_distance($node);
          my $other_node = get_other_child($parent,$node);
          my $other_total_bl = Bio::EnsEMBL::Compara::TreeUtils->total_distance($other_node);
          
          my $bl_fraction = $other_total_bl / $total_bl;
          $bl_fraction = 1 if ($bl_fraction > 1);

          # Add a value proportional to post_prob, bl, and 'balance' fraction.
          $total_prob += $post_prob/100 * ($bl / $total_dist) * $bl_fraction;
          # Add a value to make up for a 'balance' fraction less than 1.
          $total_prob += 1 * ($bl / $total_dist) * (1 - $bl_fraction);

          #$total_prob += $post_prob/100 * ($bl / $total_dist);
        } else {
          # If the alignment has a gap in the other node, don't penalize it.
          $total_prob += 100/100 * ($bl / $total_dist);
        }
        $node = $parent;
      }
      my $score = $total_prob*10 - 1;
      $score = 0 if ($score < 0);
      $score_string .= sprintf("%1d",$score);
    }
    print $score_string."\n";
    print $aln_string."\n";
    $leaf_scores->{$leaf->name} = $score_string;
  }
  return ($leaf_scores,[]);
}

sub filter_sites {
  my $class = shift;
  my $aln = shift;
  my $seq_name = shift;
  my $sites_ref = shift;
  my $mask_char = shift;

  my @sites = keys %{$sites_ref};

  foreach my $seq ($aln->each_seq) {
    my $label = $seq->display_id;
    if ($label eq $seq_name) {
      # We're at the right sequences now. Filter out the sites!
      my $seq_str = $seq->seq;
      foreach my $site (@sites) {
        substr($seq_str,$site-1,1) = $mask_char unless (substr($seq_str,$site-1,1) eq '-');
        $seq->seq($seq_str);
      }
      last;
    }
  }
  return $aln;
}

sub gblocks_filter {
  my $class = shift;
  my $aln = shift;
  my $tree = shift;
  my $params = shift;

  my $gblocks_char = $params->{'gblocks_mask_character'} || "U";

  my $node_id = $tree->node_id;
  my $dbc = $tree->adaptor->db->dbc;
  my $sth = $dbc->prepare("select aln_start,aln_end FROM gblocks_aln WHERE node_id=$node_id;");
  $sth->execute();

  my %keeper_sites;
  while (my @row = $sth->fetchrow_array) {
    my $start = $row[0];
    my $end = $row[1];

    foreach my $i ($start .. $end) {
      $keeper_sites{$i} = 1;
    }
  }

  foreach my $i (1 .. $aln->length) {
    if (!defined $keeper_sites{$i}) {
      # Mask out this column with "U"s
      foreach my $seq ($aln->each_seq) {
        my $seq_str = $seq->seq;
        substr($seq_str,$i-1,1) = $gblocks_char unless (substr($seq_str,$i-1,1) eq '-');
        $seq->seq($seq_str);
      }
    }
  }
  return $aln;
}

=head2 mask_below_score
  Title     : mask_below_score
  Usage     : $aln->mask_below_score($threshold, $score_hashref)
  Function  : Masks out nucleotides (with N) or amino acids (with X) in the alignment that are LESS THAN
             the given threshold score. So, with a threshold of 3, residues with score 3 are kept, but score 2 is turned to X or N.
             This was mainly created to allow thresholding based on T_Coffee's score_ascii output.
 Returns   : a new Bio::SimpleAlign object with the masking applied. Xs are applied to protein alphabets, and Ns to nucleotides.
 Argument  : $score_hashref should be a hash_ref with keys as the sequence labels, and values
             as alignment-shaped strings with quality scores instead of residues.
=cut

sub mask_below_score {
  my $class = shift;
  my $aln = shift;
  my $threshold = shift;
  my $score_hashref = shift;
  my $mask_char = shift;

  if ($threshold > 9) {
    warn "WARNING: Bio::SimpleAlign::mask_below_score: threshold $threshold capped to the maximum of 9.\n";
    $threshold = 9;
  }
  
  if (!ref $score_hashref && -e $score_hashref) {
    # If we're just given a string, interpret is as a filename pointing to a T_Coffee score_aln output.
    my $filename = $score_hashref;
    $score_hashref = {};
    
    my $FH = IO::File->new();
    $FH->open($filename) || throw("Could not open alignment scores file [$filename]!");
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
      $score_hashref->{$id} ||= '';
      $score_hashref->{$id} .= $align;
    }
    $FH->close;
  } else {
    # Nothing to do here... the user should have set this up as a hash_ref with keys as the display_id and values as score_strings.
    # For example:
    #  Alignment string
    #    ABC---DEF
    #  Score string
    #    995---339
  }

  my $alphabet = 'protein'; # Default to a protein alphabet.

  # Create a new shell SimpleAlign object.
  my $new_aln = new $aln;

  # Go through each sequence in the alignment and mask accordingly.
  foreach my $seq ( $aln->each_seq ) {
    $alphabet = $seq->alphabet;

    my $label = $seq->display_id;
#    print $label."\n";
    my $score_string = $score_hashref->{$label};
#    print $score_string."\n";
    $score_string =~ s/[^\d-]/9/g;   # Convert non-digits and non-dashes into 9s.
    # (The above is necessary because t_coffee leaves leftover letters in the sequence score strings)
    $score_string =~ s/[-]/9/g;   # Convert dashes to 9s, because we will never mask those out.
    #print $score_string."\n";
    my @score_array = split(//,$score_string);
    
    my @mask_array = map {$_ >= $threshold ? 1 : 0 } @score_array;  # Map the score array to an array of 0s or 1s
    #print join "", @mask_array,"\n";
    my @seq_array = split(//,$seq->seq);
    foreach my $i (0..scalar(@seq_array)-1) {
      if ($mask_array[$i] == 0 && $alphabet eq 'protein') {
        my $char = $seq_array[$i];
        next if ($char eq 'X' || $char eq 'x');
        if ($char =~ /[a-z]/) {
          $seq_array[$i] = lc($mask_char);
        } else {
          $seq_array[$i] = uc($mask_char);
        }
      } elsif ($mask_array[$i] == 0) {
        if ($char =~ /[a-z]/) {
          $seq_array[$i] = lc($mask_char);
        } else {
          $seq_array[$i] = uc($mask_char);
        }
      }
    }
    my $new_str = join "", @seq_array;

    # We've got the new alignment sequence, just put it into the new SimpleAlign.
    my $new_seq = new Bio::LocatableSeq(-seq => $new_str, -id => $label);
    $new_aln->add_seq($new_seq);
  }

  return $new_aln;
}


# GJ 2009-03-01 : Combines a sequence string and cigar line to return an alignment string.
sub combine_seq_cigar {
  my $class = shift;
  my $sequence = shift;		# Sequence string. (not bioseq!)
  my $cigar_line = shift;	# CIGAR line, in Ensembl format.

#  print $sequence."\n";
#  print $cigar_line . "\n";

  $cigar_line =~ s/([MD])/$1 /g;
  
  my @cigar_segments = split " ",$cigar_line;
  my $alignment_string = "";
  my $seq_start = 0;
  foreach my $segment (@cigar_segments) {
    if ($segment =~ /^(\d*)D$/) {
      my $length = $1;
      $length = 1 if ($length eq "");
      $alignment_string .= "-" x $length;
    } elsif ($segment =~ /^(\d*)M$/) {
      my $length = $1;
      $length = 1 if ($length eq "");
      $alignment_string .= substr($sequence,$seq_start,$length);
      $seq_start += $length;
    }
  }
#  print $alignment_string."\n";
  return $alignment_string;
}


# Masks out the entire sequence with ID $seq_id from the AlignI object $aln.
sub mask_seq_from_aln {
  my $class = shift;
  my $aln = shift;
  my $seq_id = shift;
  my $params = shift; # Optional params hashref.
  my $mask_char = shift;

  my $mask_char_from_params = $params->{alignment_quality_mask_character} || "X";

  $mask_char = $mask_char_from_params unless ($mask_char);

  my $seq = $class->get_seq_with_id($aln,$seq_id);
  my $seq_str = $seq->seq;
  
  printf "MASKING ENTIRE SEQUENCE %s: %.120s\n",$seq_id,$seq_str;
  #print "MASKING ENTIRE SEQUENCE $seq_id: ".$seq_str."\n";

  if ($seq->alphabet eq 'protein') {
    $seq_str =~ s/[^-]/$mask_char/g;
  } else {
    $seq_str =~ s/[^-]/N/g;
  }
  #print "new str: $seq_str\n";
  $seq->seq($seq_str);
  return $aln;
}

sub pretty_print {
  my $class = shift;
  my $aln = shift;
  my $params = shift;

  my $full = $params->{'full'} || 0;
  my $length = $params->{'length'} || 50;

  $full = 1 if (!$full && $length > $aln->length + 10);

  if ($full) {
    my $num_slices = int($aln->length / $length);
    
    for (my $i=0; $i <= $num_slices; $i++) {
      my $start = $i*$length+1;
      my $end = ($i+1)*$length+1;

      my $top_string = "";
      for (my $j=0; $j < $length && ($start+$j) < $aln->length; $j++) {
	if ($j % 10 == 0) {
	  my $pos_string = "".($start+$j);
	  $top_string .= $pos_string;
	  $j += length($pos_string)-1;
	} else {
	  $top_string .= " ";
	}
      }
      printf(" %-20.20s  %s\n","",$top_string);

      foreach my $seq ($aln->each_seq) {
        $end = $seq->length if ($end > $seq->length);
	my $seq_str = $seq->subseq($start,$end);
	printf(">%-20.20s  %s\n",$seq->id,$seq_str);
      }
    }
    return;
  }
  
  $length = ($length - 20)/2;

  my $aln_length = ($aln->length - 2*$length)." cols";
  my $empty_str = " "x$length;
  printf(" %-20.20s  %s ... %10s ... %s\n",
	 "",
	 $empty_str,
	 $aln_length,
	 $empty_str
	 );
	 
  foreach my $seq ($aln->each_seq) {
    my $seq_str = $seq->seq;
    my $end = substr($seq_str,-$length);
    my $start = substr($seq_str,0,$length);
    my $mid_str = "";
    if (length($seq_str) > 2*$length) {
      my $mid = substr($seq_str,$length,-$length);
      my $bases = $mid;
      $bases =~ s/-//g;
      my $gaps = $mid;
      $gaps =~ s/[^-]//g;
      my $percent = length($gaps)/length($mid) * 100;
      $mid_str = sprintf("(%1.1d%% gaps)",
			 $percent
			 );
    }

    printf(">%-20.20s  %s ... %10s ... %s\n",$seq->id,$start,$mid_str,$end);
  }
}


sub output {
  my $class = shift;
  my $aln = shift;
  my $file = shift;

  open(OUT,">$file");
  foreach my $seq ($aln->each_seq) {
    my $seq_str = $seq->seq;
    printf OUT ">%s\n%s\n",$seq->id,$seq_str;
  }
  close(OUT);
}


1;# Keep perl happy
