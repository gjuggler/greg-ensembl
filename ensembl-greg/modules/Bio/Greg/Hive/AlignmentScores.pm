package Bio::Greg::Hive::AlignmentScores;

use strict;
use Cwd;
use POSIX qw(ceil floor);
use Time::HiRes qw(sleep);
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::NestedSet;
use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Hive::Process;

use File::Path qw(mkpath);

use base ('Bio::Greg::Hive::Process');

my $ALN = 'Bio::EnsEMBL::Compara::AlignUtils';

sub param_defaults {
  my $params = {
    alignment_table       => 'protein_tree_member',
    alignment_score_table => 'protein_tree_member_score',
    filter => 'prank'    # Options: 'gblocks', 'prank', 'trimal', 'indelign'
  };
  return $params;
}

sub fetch_input {
  my ($self) = @_;

  $self->load_all_params();

  my $no_filter_param = $self->replace_params( $self->params, { alignment_score_filtering => 0 } );

  my ( $tree, $aln ) =
    Bio::EnsEMBL::Compara::ComparaUtils->get_tree_and_alignment( $self->compara_dba,
    $no_filter_param );

  $self->param( 'tree', $tree );
  $self->param( 'aln',  $aln );

  Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $aln, { length => 200 } );
}

sub run {
  my $self = shift;

  my $tree = $self->param('tree');
  my $aln  = $self->param('aln');
  my $pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($aln);
  
  my $params = $self->params;

  my $input_table = $params->{'alignment_table'};
  my $action      = $params->{'filter'};
  my $table       = $params->{'alignment_score_table'};
  
  # Score hash is a hash of strings, keyed by sequence ID.
  my $score_hash = $self->_get_alignment_scores($tree, $aln, $pep_aln);

  $self->store_scores( $tree, $score_hash, $table );
}

# Masks the CDNA alignment.
sub mask_alignment {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pep_aln = shift;
  my $score_hash = shift;

  #my $score_hash = $self->_get_alignment_scores($tree, $aln, $pep_aln);

  my $threshold = 9;
  my $mask_character = 'N';

  # Triplicate the alignment scores to match the CDNA alignment.
  foreach my $seq ($aln->each_seq) {
    my $str = $score_hash->{$seq->id};
    
    my @arr = split( //, $str );
    @arr = map { ( $_ . '' ) x 3 } @arr;
    $str = join("",@arr);      
    $score_hash->{$seq->id} = $str;
  }
  
  printf " -> Masking sequences at alignment score threshold: >= %d\n",$threshold;
  $self->param('aln_mask_threshold', $threshold);
  $aln = $ALN->mask_below_score($aln, $threshold, $score_hash, $mask_character);
  return $aln;
}

sub _get_alignment_scores {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pep_aln = shift;

  my $action = $self->param('filter');
  my $params = $self->params;

  my $score_hash;
  if ( $action =~ m/gblocks/i ) {
    print " -> RUN GBLOCKS [$action]\n";
    $self->run_gblocks( $tree, $aln, $pep_aln, $params );
  } elsif ( $action =~ 'prank' ) {
    print " -> RUN PRANK [$action]\n";
    $params->{prank_filtering_scheme} = $action;
    $self->run_prank( $tree, $aln, $pep_aln, $params );
  } elsif ( $action =~ m/trimal/i ) {
    print " -> RUN TRIMAL [$action]\n";
    $self->run_trimal( $tree, $aln, $pep_aln, $params );
  } elsif ( $action =~ m/indelign/i ) {
    print " -> RUN INDELIGN [$action]\n";
    eval { $self->run_indelign( $tree, $aln, $pep_aln, $params ); };
    if ($@) {
      print "Indelign error for action [$action]: $@\n";
    }
  } elsif ( $action =~ m/(coffee|score)/i ) {
    print " -> RUN TCOFFEE [$action]\n";
    $self->run_tcoffee( $tree, $aln, $pep_aln, $params );
  } elsif ( $action =~ m/(oracle|optimal)/i) {
    $self->run_oracle($tree,$aln,$pep_aln, $params);
  } elsif ($action =~ m/columns/i) {
    $score_hash = $self->run_columns($tree,$aln,$pep_aln, $params);
  } elsif ($action =~ m/none/i) {
    $self->run_none($tree, $aln, $pep_aln, $params);
  } elsif ($action =~ m/true/i) {
    $self->run_none($tree, $aln, $pep_aln, $params);
  } elsif ($action =~ m/branchlength/i) {    
    $self->run_branchlength($tree, $aln, $pep_aln, $params);
  } elsif ($action =~ m/guidance/i) {
    $self->run_guidance($tree, $aln, $pep_aln, $params);
  } elsif ($action =~ m/psar/i) {
    $self->run_psar($tree, $aln, $pep_aln, $params);
  } else {
    $self->throw("Alignment score action not recognized: [$action]!");
  }

  # TODO: Get score deciles, and convert to 0-9 score strings.
  my $score_hash = $self->get_score_hash($pep_aln);

  return $score_hash;
}

sub get_score_hash {
  my $self = shift;
  my $pep_aln = shift;

  my $filter = $self->param('filter');
  my $max_mask_f = $self->param('maximum_mask_fraction');

  print "  MAX MASK F: $max_mask_f\n";

  my @scores = $self->sorted_scores($pep_aln);
  print "@scores\n";
  my $max_mask_index = int($max_mask_f * (scalar(@scores)-1) );
  my $max_mask_t = $scores[$max_mask_index];
  print "  MAX MASK Threshold: $max_mask_t\n";

  # Do a hard threshold based on the different filters.
  my $t = 0.9;
  $t = 5 if ($filter eq 'tcoffee');
  $t = 5 if ($filter eq 'gblocks');
  $t = 0.5 if ($filter eq 'guidance');
  $t = 0.5 if ($filter eq 'optimal');
  $t = 0.5 if ($filter eq 'branchlength');
  $t = 50 if ($filter eq 'prank');
  $t = 50 if ($filter eq 'prank_mean');
  $t = 0.9 if ($filter eq 'psar');
  
  # Reduce the threshold to the max_mask_t if the default threshold will
  # mask out *more* than [mask_mask_fraction] sites
  if ($t > $max_mask_t) {
    print "  USING MAX MASK T!\n";
    $t = $max_mask_t 
  }

  my $score_hash;
  foreach my $seq ($pep_aln->each_seq) {
    my $score_string = '';
    foreach my $i (1 .. $pep_aln->length) {
      my $score = $self->get_score($pep_aln, $seq->id, $i);

      if ($score eq '-') {
        $score_string .= '-';
      } else {        

        if ($score > $t) {
          $score = 9;
        } else {
          $score = 0;
        }
        $score_string .= $score;
      }
    }
#    print "$score_string\n";
    $score_hash->{$seq->id} = $score_string;
  }
  return $score_hash;  
}

sub sorted_scores {
  my $self = shift;
  my $pep_aln = shift;

  my @score_bin = ();
  foreach my $seq ($pep_aln->each_seq) {
#    print $seq->id."\n";
    foreach my $i (1 .. $pep_aln->length) {
      my $n_nongap = Bio::EnsEMBL::Compara::AlignUtils->get_nongaps_at_column($pep_aln, $i);
      next if ($n_nongap <= 1);
      my $score = $self->get_score($pep_aln, $seq->id, $i);
#      print "  $i $score\n";
      push @score_bin, $score if ($score ne '-');
    }
  }

  @score_bin = sort {$a <=> $b} @score_bin;
  return @score_bin;
}


sub get_score {
  my $self = shift;
  my $aln = shift;
  my $id = shift;
  my $pos = shift;

  my $scores = $self->{_scores};
  my $score = $scores->{$id}->{$pos};
  if (!defined $score) {
    die("No score defined for $id $pos");
  } else {
    return $score;
  }
}

sub set_score {
  my $self = shift;
  my $aln = shift;
  my $id = shift;
  my $pos = shift;
  my $score = shift;

  if (!defined $self->{_scores}) {
    $self->{_scores} = {};
  }
  my $scores = $self->{_scores};
  if (!defined $scores->{$id}) {
    $scores->{$id} = {};
  }

  my $residue = Bio::EnsEMBL::Compara::AlignUtils->get_residue($aln, $id, $pos);
  if ($residue eq '-') {
    $scores->{$id}->{$pos} = '-';
  } else {
    $scores->{$id}->{$pos} = $score;
    #print " $id $pos $score\n";
  }
}

sub store_scores {
  my $self         = shift;
  my $tree         = shift;
  my $score_hash   = shift;
  my $output_table = shift;

#$self->compara_dba->dbc->do("CREATE TABLE IF NOT EXISTS $output_table LIKE protein_tree_member_score");
  my $sth = $tree->adaptor->prepare(
    "REPLACE INTO $output_table (node_id,member_id,cigar_line) VALUES (?,?,?)");
  foreach my $leaf ( $tree->leaves ) {
    my $score_string = $score_hash->{ $leaf->stable_id };

    die("No score string found when saving!")
      unless ( defined $score_string && $score_string ne '' );
    $sth->execute( $leaf->node_id, $leaf->member_id, $score_string );
    printf "%20s %10s %s\n", $leaf->stable_id, $leaf->member_id, $score_string;
    sleep(0.1);
  }
  $sth->finish;
}

sub run_none {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pep_aln = shift;

  foreach my $seq ($pep_aln->each_seq) {
    foreach my $i (1 .. $pep_aln->length) {
      $self->set_score($pep_aln, $seq->id, $i, 1);
    }
  }
}

sub run_branchlength {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pep_aln = shift;
  my $params = shift;

  my $action = $self->param('filter');
  my $window_flank = 0;

  if ($action eq 'window_branchlength') {
    $window_flank = 2;
  }

  my $total_bl = Bio::EnsEMBL::Compara::TreeUtils->total_distance($tree);
  my @id_list = map { $_->id } $pep_aln->each_seq;
  my $seq_arrayref = Bio::EnsEMBL::Compara::AlignUtils->to_arrayrefs($pep_aln);

  my $tree_bl_hash = {};
  my $scores;
  foreach my $i ( 1 .. $pep_aln->length ) {
    my @nongap_ids_at_pos = grep { $seq_arrayref->{$_}->[$i] ne '-' } @id_list;
    my $nongap_bl = $self->get_subtree_bl( $tree, \@nongap_ids_at_pos, $tree_bl_hash );

    foreach my $seq ($aln->each_seq) {
      my $score = $nongap_bl / $total_bl;
      die unless (defined $score);
      $self->set_score($pep_aln, $seq->id, $i, $score);
    }
  }
}

sub run_columns {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pep_aln = shift;
  my $params = shift;

  my $num_sequences = scalar($tree->leaves);

  my @column_scores;
  foreach my $i (1 .. $pep_aln->length) {
    my $col_string = Bio::EnsEMBL::Compara::AlignUtils->get_column_string($pep_aln,$i);
    $col_string =~ s/[^-]//g;
    my $num_nongaps = length($col_string);
    my $gap_fraction = $num_nongaps / $num_sequences;
    my $score = $gap_fraction * 9;
    $score = 9 if ($score > 9);
    $score = 0 if ($score < 0);
    $score = sprintf("%1d",$score);

    push @column_scores,$score;
  }

  my %scores_hash;
  
  foreach my $leaf ( $tree->leaves ) {
    my $seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id($pep_aln, $leaf->name);
    my $aln_string = _apply_columns_to_leaf( \@column_scores, $seq->seq );
    $scores_hash{ $leaf->stable_id } = $aln_string;
  }

  return \%scores_hash;
}

sub run_oracle {
  my $self   = shift;
  my $tree   = shift;
  my $aln    = shift;
  my $pep_aln    = shift;
  my $params = shift;

  my $true_aln = $self->param('true_aln');
  my $true_pep_aln = $self->param('true_pep_aln');

  print $tree->newick_format."\n";

  # Goal: score each site by the fraction of the tree against which
  #   it is correctly aligned.

  # Plan:
  # 1) For each site, calculate the size of the subtree comprised of all
  #    sequences against which it's correctly aligned. Score as the fraction
  #    compared to the tree's total branch length.

  # Set-up:
  # Index the alignments by sequence residue #.
  my $ALNU     = 'Bio::EnsEMBL::Compara::AlignUtils';
  my $true_obj = $ALNU->to_arrayrefs($true_pep_aln);
  my $test_obj = $ALNU->to_arrayrefs($pep_aln);

  # Add all aligned pairs to the index hashtable.
  my $true_pairs = $ALNU->store_pairs( $true_pep_aln, $true_obj );
  my $test_pairs = $ALNU->store_pairs( $pep_aln,  $test_obj );

  my @id_list = map { $_->id } $pep_aln->each_seq;

  my $tree_bl_hash = {};
  my $scores_hash;
  map { $scores_hash->{$_} = '' } @id_list;
  foreach my $i ( 1 .. $pep_aln->length ) {
    my @nongap_ids_at_pos = grep { $test_obj->{$_}->[$i] ne '-' } @id_list;
    my $total_nongap_bl = $self->get_subtree_bl( $tree, \@nongap_ids_at_pos, $tree_bl_hash );

    print $ALNU->get_column_string( $pep_aln, $i ) . "\n";
    my $column_score_string = '';
    foreach my $this_seq_id (@id_list) {

      # Get a list of all other IDs to which this seq-residue is correctly aligned.
      my $this_res_num = $test_obj->{$this_seq_id}->[$i];

      if ( $this_res_num eq '-' ) {
        $self->set_score($pep_aln, $this_seq_id, $i, '-');
        $scores_hash->{$this_seq_id} .= '-';
        $column_score_string .= '-';
        next;
      }

      my @correctly_aligned_ids = ($this_seq_id);
      foreach my $other_seq_id (@id_list) {
        next if ( $this_seq_id eq $other_seq_id );
        my $other_res_num = $test_obj->{$other_seq_id}->[$i];

        #        next if ($other_res_num eq '-');

        my $pair_string = join( '_', $this_seq_id, $this_res_num, $other_seq_id, $other_res_num );
        if ( $true_pairs->{$pair_string} == 1 ) {

          #print "$pair_string yes!!\n";
          push @correctly_aligned_ids, $other_seq_id;
        }
      }
      my $correct_bl = $self->get_subtree_bl( $tree, \@correctly_aligned_ids, $tree_bl_hash );

      #printf "%s - %s  %.3f %.3f\n", $this_seq_id, join(',',@correctly_aligned_ids), $total_nongap_bl, $correct_bl;

      # Spread out a little bit towards lower scores...
      my $score = 0;
      if ($total_nongap_bl > 0) {
        $score = 0 + ($correct_bl/$total_nongap_bl);
      }
      $self->set_score($pep_aln, $this_seq_id, $i, $score);

      if ($total_nongap_bl > 0) {
        $score = 0 + ($correct_bl/$total_nongap_bl * 9.0);
      }
      $score = 9 if ($score > 9);
      $score = 0 if ($score < 0);
      $score = sprintf("%1d",$score);
      $scores_hash->{$this_seq_id} .= $score;
      $column_score_string .= $score;

      #printf "%d %-20s:%.3f\n",$i,$this_seq_id,$score;
    }
    print $column_score_string. "\n";
  }
  
#  return $scores_hash;
}

sub get_subtree_bl {
  my $self    = shift;
  my $tree    = shift;
  my $seq_ids = shift;
  my $bl_hash = shift;

  my @id_array = @$seq_ids;
  @id_array = sort {$a cmp $b} @id_array;
  my $key = join('_',@id_array);
  my $existing_value = $bl_hash->{$key};
  return $existing_value if (defined $existing_value);
    
  #print "No existing value! $key\n";
  my $subtree = Bio::EnsEMBL::Compara::TreeUtils->extract_subtree_from_leaves($tree,\@id_array, 1);

  my $total = Bio::EnsEMBL::Compara::TreeUtils->total_distance($subtree);
  if (scalar($subtree->leaves) == 1) {
    $total = 0;
  }
  $bl_hash->{$key} = $total;
  return $total;
}

sub run_indelign {
  my $self   = shift;
  my $tree   = shift;
  my $aln    = shift;
  my $pep_aln    = shift;
  my $params = shift;

  my ( $i_obj, $d_obj, $ins_rate, $del_rate ) =
    Bio::EnsEMBL::Compara::AlignUtils->indelign( $pep_aln, $tree, $params,
    $self->worker_temp_directory );
  my @ins = @$i_obj;
  my @del = @$d_obj;

  my $num_leaves         = scalar( $tree->leaves );
  my $tree_length        = Bio::EnsEMBL::Compara::TreeUtils->total_distance($tree);
  my $rate_sum           = $ins_rate + $del_rate;
  my $max_allowed_indels = $rate_sum * $tree_length * 1;
  $max_allowed_indels = floor($max_allowed_indels);
  $max_allowed_indels = 1 if ( $max_allowed_indels < 1 );

  my $blocks_string = "1" x $pep_aln->length;
  for ( my $i = 0 ; $i < $pep_aln->length ; $i++ ) {
    my $indel_sum = $ins[$i] + $del[$i];

    if ( $indel_sum > $max_allowed_indels ) {
      my $indel_excess = ( $indel_sum / $max_allowed_indels ) * 2;
      my $indel_score  = 9 - $indel_excess;
      $indel_score = 0 if ( $indel_score < 0 );
      substr $blocks_string, $i, 1, floor($indel_score);
    } else {
      substr $blocks_string, $i, 1, '9';
    }
  }

  my @column_scores = split( "", $blocks_string );

  #  print "@column_scores\n";
  my %scores_hash;
  foreach my $leaf ( $tree->leaves ) {
    my $seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id($pep_aln, $leaf->name);
    my $aln_string = _apply_columns_to_leaf( \@column_scores, $seq->seq );
    print "$aln_string\n";
    $scores_hash{ $leaf->stable_id } = $aln_string;
  }
  return \%scores_hash;
}

sub run_tcoffee {
  my $self   = shift;
  my $tree   = shift;
  my $aln    = shift;
  my $pep_aln    = shift;
  my $params = shift;

  Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $pep_aln, { length => 200 } );

  my $tmpdir = $self->worker_temp_directory;
  print "$tmpdir\n";
  system("rm -rf ${tmpdir}*.*");

  my $filename = "$tmpdir" . "tcoffee_aln.fasta";
  my $tmpfile  = Bio::AlignIO->new(
    -file   => ">$filename",
    -format => 'fasta'
  );
  $tmpfile->write_aln($pep_aln);
  $tmpfile->close;

  my $tmp = $tmpdir;
  $tmp = substr($tmp,0,-1);
  my $prefix = "export HOME_4_TCOFFEE=\"${tmp}\";";
  $prefix = "export DIR_4_TCOFFEE=\"${tmp}\";";
  $prefix .= "export METHODS_4_TCOFFEE=\"${tmp}\";";
  $prefix .= "export MCOFFEE_4_TCOFFEE=\"${tmp}\";";
  $prefix .= "export TMP_4_TCOFFEE=\"${tmp}\";";
  $prefix .= "export CACHE_4_TCOFFEE=\"${tmp}\";";
  $prefix .= "export NO_ERROR_REPORT_4_TCOFFEE=1;";
  $prefix .= "export NUMBER_OF_PROCESSORS_4_TCOFFEE=1;";
    ;    # GJ 2008-11-04. What a hack!

  my $outfile = $filename . ".score_ascii";

  my $bin = 't_coffee';
  my $cmd =
    qq^$bin -mode=evaluate -evaluate_mode t_coffee_slow -infile=$filename -outfile=$outfile -output=score_ascii -n_core=1 -multi_core=no^;
  print $cmd. "\n";
  my $rc = system( $prefix. $cmd );

  $rc == 0 or die("TCoffee failed: $?");

  my $scores_file = $outfile;
  my %score_hash;
  if ( -e $scores_file ) {
    my $FH = IO::File->new();
    $FH->open($scores_file) || die("Could not open tcoffee scores file!");
    <$FH>;    #skip header
    my $i = 0;
    while (<$FH>) {
      $i++;

      #      next if ($i < 7); # skip first 7 lines.
      next if ( $_ =~ /^\s+/ );    #skip lines that start with space
      if ( $_ =~ /:/ ) {

        #my ($id,$overall_score) = split(/:/,$_);
        #$id =~ s/^\s+|\s+$//g;
        #$overall_score =~ s/^\s+|\s+$//g;
        #print "___".$id."___".$overall_score."___\n";
        next;
      }
      chomp;
      my ( $id, $align ) = split;
      $_ = $id;
      $id =~ s^/.*^^;
      $score_hash{$id} = '' if ( !exists $score_hash{$id} );
      $score_hash{$id} .= $align;

      #print $id." ". $align."\n";
    }
    $FH->close;
  } else {
    die("No tcoffee scores file at $scores_file!!\n");
  }

  foreach my $seq ( $pep_aln->each_seq ) {
    my $id     = $seq->id;
    my $string = $score_hash{$id};
    die("No score string found for seq $id!") unless ( defined $string );

    $string =~ s/[^\d-]/9/g;
    # Convert non-digits and non-dashes into 9s. This is necessary because t_coffee leaves some leftover letters.
    my @chars = split('', $string);

    foreach my $i (1 .. $pep_aln->length) {
      $self->set_score($pep_aln, $seq->id, $i, $chars[$i-1]);
    }
  }
}

sub run_gblocks {
  my $self   = shift;
  my $tree   = shift;
  my $aln    = shift;
  my $pep_aln    = shift;
  my $params = shift;

  my $defaults = {
    t  => 'p',    # Type of sequence (p=protein,c=codon,d=dna)
    b3 => '8',    # Max # of contiguous nonconserved positions
    b4 => '3',    # Minimum length of a block
    b5 => 'a',    # Allow gap positions (n=none, h=with half,a=all)
    b6 => 'y',
  };
  my $use_params = $self->replace_params( $defaults, $params );

  printf("Sitewise_dNdS::run_gblocks\n") if ( $self->debug );

  my $aln_length = $pep_aln->length;
  my $tmpdir     = $self->worker_temp_directory;
  my $filename   = "$tmpdir" . "gblocks_aln.fasta";
  my $tmpfile    = Bio::AlignIO->new(
    -file   => ">$filename",
    -format => 'fasta'
  );
  $tmpfile->write_aln($pep_aln);
  $tmpfile->close;

  my @leaves             = $tree->leaves;
  my $num_leaves         = scalar(@leaves);
  my $b1 = int( ( $num_leaves + 1) / 2 + 0.5 );
  my $b2 = int( ($num_leaves * 0.86) + 0.5);

  #Example command: Gblocks 2138.fasta -t=p -b2=20 -b3=50 -b4=3 -b5=a -p=s
  my $cmd = sprintf(
    "Gblocks %s -t=%s -b1=%s -b2=%s -b3=%s -b4=%s -b5=%s -b6=%s -p=s\n",
    $filename,           $use_params->{'t'},  $b1, $b2,
    $use_params->{'b3'}, $use_params->{'b4'}, $use_params->{'b5'},
    $use_params->{'b6'}
  );
  print "GBLOCKS: $cmd\n";
  my $ret = system("$cmd");
  open FLANKS, "$filename-gb.txts" or die "$!\n";
  my $segments_string;
  while (<FLANKS>) {
    chomp $_;
    next unless ( $_ =~ /Flanks:/ );
    $segments_string = $_;
    last;
  }
  close FLANKS;
  $segments_string =~ s/Flanks\: //g;
  $segments_string =~ s/\s+$//g;

  print "SEGS: " . $segments_string . "\n";
  $_ = $segments_string;
  my @bs = /\[(.+?)\]/g;

  my $blocks_string = "0" x $pep_aln->length;
  foreach my $b (@bs) {
    my ( $start, $end ) = split( " ", $b );
    for ( my $i = $start ; $i <= $end ; $i++ ) {
      substr $blocks_string, $i - 1, 1, '9';
    }
  }

  my @column_scores = split( "", $blocks_string );

  foreach my $seq ($pep_aln->each_seq) {
    foreach my $i (1 .. $pep_aln->length) {
      my $col_score = $column_scores[$i-1];
      $self->set_score($pep_aln, $seq->id, $i, $column_scores[$i-1]);
    }
  }
}

sub _apply_columns_to_leaf {
  my $columns_ref = shift;
  my $aln_str = shift;

  my @column_scores = @{$columns_ref};

  my @aln_chars = split( "", $aln_str );

  for ( my $i = 0 ; $i < scalar(@aln_chars) ; $i++ ) {
    if ( $aln_chars[$i] ne '-' ) {
      $aln_chars[$i] = $column_scores[$i] || '0';
    }
  }

  return join( "", @aln_chars );
}

sub run_trimal {
  my $self   = shift;
  my $tree   = shift;
  my $aln    = shift;
  my $pep_aln    = shift;
  my $params = shift;

  # Write temporary alignment.
  my $dir   = $self->worker_temp_directory;
  my $aln_f = $dir . "/aln_" . $self->data_id . ".fa";
  Bio::EnsEMBL::Compara::AlignUtils->to_file( $pep_aln, $aln_f );

  # Build a command for TrimAl.
  my $trim_params = '';
  $trim_params = $params->{'trimal_filtering_params'}
    if ( $params->{'trimal_filtering_params'} );    # Custom parameters if desired.
  my $cmd = "trimal -gt 0.3 -cons 20 -in $aln_f -out $aln_f -colnumbering $trim_params";
  print "CMD: $cmd\n";
  my $output = `$cmd`;

  #print "OUTPUT:". $output."\n";

  chomp $output;
  my @cons_cols = split( /[,]/, $output );
  print join( ",", @cons_cols ) . "\n";

  my @column_scores = 0 x $pep_aln->length;
  foreach my $column (@cons_cols) {
    $column_scores[$column] = 9;
  }

  my %scores_hash;
  foreach my $leaf ( $tree->leaves ) {
    my $aln_string = _apply_columns_to_leaf( \@column_scores, $leaf );
    $scores_hash{ $leaf->stable_id } = $aln_string;
  }

  return \%scores_hash;
}

sub run_guidance {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pep_aln = shift;
  my $params = shift;
  
  $params->{temp_dir} = $self->worker_temp_directory;

  my $dir = $self->worker_temp_directory;

  system("rm -rf ${dir}*");
  mkpath($dir);

  my $node_id = $tree->node_id;
  my $seqs_f = $dir."seqs_${node_id}.fasta";
  my $aln_f = $dir."aln_${node_id}.fasta";
  my $tree_f = $dir."tree_${node_id}.nh";
  Bio::EnsEMBL::Compara::AlignUtils->dump_ungapped_seqs($pep_aln,$seqs_f);
  Bio::EnsEMBL::Compara::AlignUtils->to_file($pep_aln,$aln_f);
  Bio::EnsEMBL::Compara::TreeUtils->to_file($tree,$tree_f);  

  my $guidance_dir = $self->base . "/projects/slrsim/guidance.v1.01/";
  my $cwd = cwd();
  chdir $guidance_dir;

  my $reps = 100;
  my $msa_program = "MAFFT";
  if ($self->param('aligner') =~ m/prank/i) {
    $msa_program = "PRANK";
    $reps = 30;
  }

  my $cmd = qq^perl www/Guidance/run_calc.pl --seqFile ${seqs_f} --msaFile ${aln_f} --msaProgram ${msa_program} --seqType aa --bootstraps ${reps} --seqCutoff 0 --colCutoff 0 --outDir ${dir}^;
  my $rc = system($cmd);
  $rc == 0 or die("GUIDANCE failed: $? ($!)");

  # Parse the GUIDANCE output into residue scores.
  my $seq_codes_f = $dir."Seqs.Codes";
  my $residue_scores_f = $dir."MSA.${msa_program}.Guidance_res_pair_res.scr";

  my $num_to_seq_id;
  open(IN, "$seq_codes_f");
  my @lines = <IN>;
  close(IN);
  foreach my $line (@lines) {
    chomp $line;
    my @toks = split("\\s+", $line);
    print "[".$toks[1]."]\n";
    $num_to_seq_id->{$toks[1]} = $toks[0];
  }

  my $row_col_scores;
  open(IN, "$residue_scores_f");
  my @lines = <IN>;
  close(IN);
  print shift @lines; # Remove header line.
  foreach my $line (@lines) {
    chomp $line;
    $line =~ s/^\s+//;
    my @toks = split("\\s+", $line);
    my $key = $num_to_seq_id->{$toks[1]} . $toks[0];
    $row_col_scores->{$key} = $toks[2];
  }

  my @id_list = map { $_->id } $pep_aln->each_seq;
  my $scores_hash;
  foreach my $seq_id (@id_list) {
    my $score_string = '';
    foreach my $i ( 1 .. $pep_aln->length) {
      my $score = $row_col_scores->{$seq_id.$i};

      if ($score) {
        $self->set_score($pep_aln, $seq_id, $i, $score);
        $score = 0 + ($score * 9.0);
        $score = 9 if ($score > 9);
        $score = 0 if ($score < 0);
        $score = sprintf("%1d",$score);
        $score_string .= $score;
      } else {
        $self->set_score($pep_aln, $seq_id, $i, '-');
        $score_string .= '-';
      }

    }
    $scores_hash->{$seq_id} = $score_string;
  }

  chdir $cwd;

#  return $scores_hash;
}

sub run_psar {
  my $self   = shift;
  my $tree   = shift;
  my $aln    = shift;
  my $pep_aln    = shift;
  my $params = shift;

  $params->{temp_dir} = $self->worker_temp_directory;

  my $dir = $self->worker_temp_directory;

  print join(" ", map {$_->id} $aln->each_seq)."\n";

  my $aln_f = $dir."aln.fasta";
  Bio::EnsEMBL::Compara::AlignUtils->to_file($aln,$aln_f);

  my $params_f = $dir."params.txt";
  my $samples_dir = $dir."samples";
  my $params_s = qq^# [REQUIRED] Parameter for the nucleotide substitutions
SUBSTPAR=0.6
# [REQUIRED] Nucleotide background probabilities
PA=0.3
PC=0.2
PG=0.2
PT=0.3
# [OPTIONAL]State transition probabilities for states, M, IS, and IA
#MtoIS=0.02
#MtoIA=0.02
#IStoIS=0.8
#IAtoIA=0.8
#IAtoIS=0.00001
#IStoIA=0.00001^;
  open(OUT, ">$params_f");
  print OUT $params_s."\n";
  close(OUT);

  my $cwd = cwd();
  chdir $dir;
  
  my $cmd = qq^psar ${aln_f} ${params_f} ${samples_dir}^;
  my $rc = system($cmd);
  ($rc == 0) or die("PSAR failed: $? ($!)");

  chdir $cwd;

  # Parse the PSAR output into residue scores -- we need to take an average pair-score
  # for each residue.
  my $pair_scores_f = $dir."PSAR_pair.txt";
  die unless (-e $pair_scores_f);
  open(IN, $pair_scores_f);
  my @lines = <IN>;
  close(IN);

  my $pairs_hash;
  shift @lines; # Remove the header line.
  foreach my $line (@lines) {
    chomp $line;
    my @toks = split("\\s+", $line);
    my ($seq1, $seq2, $pos1, $pos2, $score) = @toks;
    my $pair_key = join('_',$seq1,$seq2,$pos1,$pos2);
    $pairs_hash->{$pair_key} = $score;
    $pair_key = join('_',$seq2,$seq1,$pos2,$pos1);
    $pairs_hash->{$pair_key} = $score;
    #print $pair_key . " $score\n";
  }

  my @seqs = $aln->each_seq;
  my @pep_seqs = $pep_aln->each_seq;
  my $seq_x;
  my $seq_y;
  foreach my $i (1 .. $pep_aln->length) {
    foreach my $x (0 .. scalar(@seqs)-1) {
      $seq_x = $seqs[$x];
      my $pep_seq = $pep_seqs[$x];
      my $pep_location = $pep_seq->location_from_column($i);
      my $score_sum = 0;
      my $score_denominator = 0;
      if (defined $pep_location && $pep_location->location_type() eq 'EXACT') {
        my $start_j = ($i - 1) * 3 + 1;
        foreach my $j ( $start_j,  $start_j + 1, $start_j + 2 ) {
          my $location = $seq_x->location_from_column($j);
          foreach my $y (0 .. scalar(@seqs)-1) {
            next if ($y == $x);
            $seq_y = $seqs[$y];
            my $location_b = $seq_y->location_from_column($j);
            if (defined $location && $location->location_type() eq 'EXACT' &&
                defined $location_b && $location_b->location_type() eq 'EXACT') {
              my $loc_x = $location->start;
              my $loc_y = $location_b->start;

              my $pair_id = join('_', $x+1, $y+1, $loc_x, $loc_y);
              my $pair_score = $pairs_hash->{$pair_id};
              #print $pair_id."\n";
              #print "$i $j $x $y $pair_score\n";
              die unless (defined $pair_score);

              $score_sum += $pair_score;
              $score_denominator += 1;
            } else {
              # Nothing - gap in other sequence!
            }
          }
        }
      } else {
        $self->set_score($pep_aln, $seq_x->id, $i, '-');
      }
      
      # Take the mean pairwise score across all codon positions.
      my $mean_score = 1;
      if ($score_denominator > 0) {
        $mean_score = $score_sum / $score_denominator;
      }
      #print $seq_x->id." ".$i." ".$mean_score."\n";
      $self->set_score($pep_aln, $seq_x->id, $i, $mean_score);
    }
  }
}

sub run_prank {
  my $self   = shift;
  my $tree   = shift;
  my $aln    = shift;
  my $pep_aln    = shift;
  my $params = shift;

  $params->{temp_dir} = $self->worker_temp_directory;

  $self->_get_prank_filter_matrices( $tree, $aln, $pep_aln, $params );
}

sub _get_prank_filter_matrices {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pep_aln = shift;
  my $params = shift;

  my $node_id = $tree->node_id;

  my $dir = $self->worker_temp_directory . 'prank_tmp';
  mkpath([$dir]);
  my $aln_f = $dir."/aln_${node_id}.fasta";
  my $tree_f = $dir."/tree_${node_id}.nh";
  my $out_f = $dir."/aln_filtered_${node_id}";
  my $xml_f = $out_f.".0.xml";

  # Output tree and alignment.
  $ALN->to_file($pep_aln,$aln_f);
  Bio::EnsEMBL::Compara::TreeUtils->to_file($tree,$tree_f);

  my $prank_bin = "prank";

  my $cmd = qq^${prank_bin} -d=$aln_f -t=$tree_f -e -o=$out_f^;
  system($cmd);

  my $module = 'XML::LibXML';
  eval "use $module";
  use Bio::Greg::Node;

  # Grab information from Prank's XML output.
  my $parser;
  eval {
    $parser = XML::LibXML->new();
  };
  if ($@) {
    print "$@\n";
    return;
  }
  my $xml_tree = $parser->parse_file($xml_f);
  my $root = $xml_tree->getDocumentElement;

  my $newick = ${$root->getElementsByTagName('newick')}[0]->getFirstChild->getData;
  #print $newick."\n";
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
  
  my $pp_hash;
  my $leaf_scores;
  my @nodes = $rootNode->nodes();
  foreach my $node (@nodes) {
    next if ($node->isLeaf);
    my @score = split(/,/,$postprob{$node->name}{$nameToState{'postprob'}});
    $pp_hash->{$node->name} = ();
    for (my $i=0; $i < scalar(@score); $i++) {
      $pp_hash->{$node->name}[$i] = $score[$i];
    }
  }

  my $aln_len = $pep_aln->length;

  if ($params->{'prank_filtering_scheme'} =~ m/prank_column/i) {
    my @score_array;

    for (my $i=0; $i < $aln_len; $i++) {
      my $score_sum = 0;
      my $bl_sum = 0;
      foreach my $node ($tree->nodes) {
        my $bl = $node->distance_to_parent;

        my $xml_node = $node_to_xml->{$node};
        my $post_prob = $pp_hash->{$xml_node}[$i];
        if (defined $xml_node && defined $post_prob && $post_prob != -1) {
          $score_sum += $post_prob * $bl;
          $bl_sum += 1 * $bl;
        } else {

        }
      }
      $score_array[$i] = 0;
      if ($bl_sum > 0) {
        my $weighted_score = $score_sum / $bl_sum;
        $score_array[$i] = $weighted_score;
      }
    }

    foreach my $leaf ($tree->leaves) {
      my $seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id($pep_aln, $leaf->name);
      my $score_string = "";
      my $aln_string = $seq->seq;
      my @aln_arr = split("",$aln_string);
      for (my $i=0; $i < $aln_len; $i++) {
        if ($aln_arr[$i] eq '-') {
          $score_string .= '-';
        } else {
          my $sitewise_score = $score_array[$i] / 10;
          $sitewise_score = 0 if ($sitewise_score < 0);
          $sitewise_score = 9 if ($sitewise_score > 9);
          $score_string .= sprintf("%1d",$sitewise_score);
        }
      }
      $leaf_scores->{$leaf->name} = $score_string;
    }
  } elsif ($params->{'prank_filtering_scheme'} eq 'prank' || $params->{'prank_filtering_scheme'} =~ m/prank_mean/i) {
    foreach my $leaf ($tree->leaves) {
      my $seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id($pep_aln, $leaf->name);
      my $score_string = "";
      my $aln_string = $seq->seq;
      my @aln_arr = split("",$aln_string);

      for (my $i=0; $i < $aln_len; $i++) {
        if ($aln_arr[$i] eq '-') {
          $score_string .= '-';

          # Set the score.
          $self->set_score($pep_aln, $seq->id, $i+1, '-');
        } else {
          my $sitewise_score = 9;

          my $bl_sum = 0;
          my $pp_sum = 0;
          my $node = $leaf;
          while (my $parent = $node->parent) {
            my $bl = $node->distance_to_parent;

            my $xml_node = $node_to_xml->{$parent};
            my $post_prob = $pp_hash->{$xml_node}[$i];
            if (defined $xml_node && defined $post_prob && $post_prob != -1) {
              $bl_sum += 1 * $bl;
              $pp_sum += $post_prob * $bl;
            } else {
            }
            $node = $parent;
          }
          $bl_sum = 0.01 if ($bl_sum == 0);
          my $sitewise_score = $pp_sum / $bl_sum;

          # Set the score.
          $self->set_score($pep_aln, $seq->id, $i+1, $sitewise_score);

          $sitewise_score = 0 if ($sitewise_score < 0);
          $sitewise_score = 9 if ($sitewise_score > 9);
          $score_string .= sprintf("%1d",$sitewise_score);
        }
      }
      $leaf_scores->{$leaf->name} = $score_string;
    }

  } elsif ($params->{'prank_filtering_scheme'} =~ m/prank_min/i) {
    foreach my $leaf ($tree->leaves) {
      my $score_string = "";
      my $seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id($pep_aln, $leaf->name);
      my $aln_string = $seq->seq;
      my @aln_arr = split("",$aln_string);

      for (my $i=0; $i < $aln_len; $i++) {
        if ($aln_arr[$i] eq '-') {
          $score_string .= '-';

          # Set the score.
          $self->set_score($pep_aln, $seq->id, $i+1, '-');
        } else {
          
          my $min_pp = 100;
          my $node = $leaf;
          while (my $parent = $node->parent) {

            my $xml_node = $node_to_xml->{$parent};
            my $post_prob = $pp_hash->{$xml_node}[$i];
            if (defined $xml_node && defined $post_prob && $post_prob != -1) {
              
              $min_pp = $post_prob if ($post_prob < $min_pp);
            }
            $node = $parent;
          }

          # Set the score.
          $self->set_score($pep_aln, $seq->id, $i+1, $min_pp);

          my $sitewise_score = $min_pp / 10;
          $sitewise_score = 0 if ($sitewise_score < 0);
          $sitewise_score = 9 if ($sitewise_score > 9);
          $score_string .= sprintf("%1d",$sitewise_score);
        }
      }
      $leaf_scores->{$leaf->name} = $score_string;
    }

  } elsif ($params->{'prank_filtering_scheme'} =~ m/prank_treewise/i) {
    
    foreach my $leaf ($tree->leaves) {
      my $total_dist = $leaf->distance_to_root;
      my $aln_len = $pep_aln->length;
      my $seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id($pep_aln, $leaf->name);
      my $aln_string = $seq->seq;
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
      $leaf_scores->{$leaf->name} = $score_string;
    }
    
  }
  return $leaf_scores;
}


1;
