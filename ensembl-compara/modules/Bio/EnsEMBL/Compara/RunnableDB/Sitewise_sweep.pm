#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Compara::RunnableDB::Sitewise_sweep

=cut

=head1 DESCRIPTION

A quick hacky module to run SLR to grab the estimated dS per branch,
filtering over a range of values from the t_coffee alignment scores.

input_id/parameters format eg: "{'protein_tree_id'=>1234}"
    protein_tree_id : use 'id' to fetch a cluster from the ProteinTree

=cut

package Bio::EnsEMBL::Compara::RunnableDB::Sitewise_sweep;

use strict;
use Getopt::Long;
use IO::File;
use File::Basename;
use Time::HiRes qw(time gettimeofday tv_interval);

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

use Bio::AlignIO;
use Bio::TreeIO;
use Bio::SimpleAlign;

use Cwd;

use Bio::EnsEMBL::Hive;
our @ISA = qw(Bio::EnsEMBL::Hive::Process);

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data from the database.
    Returns :   none
    Args    :   none

=cut


sub fetch_input {
  my( $self) = @_;

  $self->check_job_fail_options;
  $self->throw("No input_id") unless defined($self->input_id);

  #create a Compara::DBAdaptor which shares the same DBI handle
  #with the Pipeline::DBAdaptor that is based into this runnable
  $self->{'comparaDBA'} = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new
    (
     -DBCONN=>$self->db->dbc
    );

  ### DEFAULT PARAMETERS ###
  my $p = '';

  $p = 'input_table';
  $self->{$p} = 'protein_tree_member';

  $p = 'output_table';
  $self->{$p} = 'sitewise_sweep';

  $p = 'mask_steps';
  my @mask_arr = (1,2,3,4,5,6,7,8,9);
  $self->{$p} = \@mask_arr;

  $p = 'saturated';
  $self->{$p} = 0.01;
  #########################

  # Load up the global analysis parameters
  $self->get_params($self->parameters);
  
  # Load up the per-job input parameters. Note that any duplicate keys will override the analysis-wide parameter from the above call.
  $self->get_params($self->input_id);
  
  $self->check_if_exit_cleanly;
  
  ### LOAD UP THE TREE AND EXIT IF IT AINT RIGHT ###
  if (defined($self->{'protein_tree_id'})) {
      my $pta = $self->{'comparaDBA'}->get_ProteinTreeAdaptor();
      $self->{'protein_tree'} = $pta->fetch_node_by_node_id($self->{'protein_tree_id'});
      
      # Allow for a different input table.
      if (! ($self->{'input_table'} eq 'protein_tree_member')) {
	  # We've got an alternate member table, so let's grab the cigar line from there.
	  my $tree = $self->{'protein_tree'};
	  my $table = $self->{'input_table'};
	  
	  foreach my $leaf (@{$tree->get_all_leaves}) {
	      # "Release" the stored / cached values for the alignment strings.
	      undef $leaf->{'cdna_alignment_string'};
	      undef $leaf->{'alignment_string'};
	      
	      # Grab the correct cigar line for each leaf node.
	      my $id = $leaf->member_id;
	      my $cmd = "SELECT cigar_line FROM $table where member_id=$id;";
	      my $sth = $pta->prepare($cmd);
	      $sth->execute();
	      my $data = $sth->fetchrow_hashref();
	      $sth->finish();
	      my $cigar = $data->{'cigar_line'};
	      
	      throw "No cigar line for member $id!\n" unless ($cigar);
	      $leaf->cigar_line($cigar);
	  }
      }
  }

  throw "undefined ProteinTree as input!\n" unless (defined $self->{'protein_tree'});
  my $num_leaves = $self->{'protein_tree'}->num_leaves;
  #if ($num_leaves < 4) {
  #    $self->{'protein_tree'}->release_tree;
  #    $self->{'protein_tree'} = undef;
  #    throw "Protein Tree leaf count below 4, FAIL THE JOB!!!";
  #}
  #########################

  ### LOAD UP THE CODON ALIGNMENT AND NEWICK FORMAT TREE ###
  $p = 'cdna';
  $self->{$p} = $self->{'protein_tree'}->get_SimpleAlign(-cdna => 1);
  $p = 'newick';
  $self->{$p} = $self->{'protein_tree'}->newick_format("int_node_id");
  #########################

  $self->{'protein_tree'}->print_tree(10) if ($self->debug);
  return 1;
}


=head2 run

    Title   :   run
    Usage   :   $self->run
    Function:   runs the analysis.
    Returns :   none
    Args    :   none

=cut


sub run {
  my $self = shift;
  $self->check_if_exit_cleanly;
  $self->run_sitewise_dNdS;
}


=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   stores SLR annotations
    Returns :   none
    Args    :   none

=cut


sub write_output {
  my $self = shift;

  $self->store_sitewise_dNdS;

  $self->check_if_exit_cleanly;

  $self->{'protein_tree'}->release_tree;
}


sub DESTROY {
  my $self = shift;

  if ($self->{'protein_tree'}) {
    printf("Sitewise_dNdS::DESTROY  releasing tree\n") if($self->debug);
    $self->{'protein_tree'}->release_tree;
    $self->{'protein_tree'} = undef;
  }

  $self->SUPER::DESTROY if $self->can("SUPER::DESTROY");
}


##########################################
#
# internal methods
#
##########################################

sub get_params {
  my $self         = shift;
  my $param_string = shift;

  return unless($param_string);
  print("parsing parameter string : ",$param_string,"\n") if($self->debug);

  # Turn the given param string into a hashref.
  my $params = eval($param_string);
  return unless($params);

  my $p = '';

  # PROTEIN_TREE_ID
  $p = 'protein_tree_id';
  $self->{$p} = $params->{$p} if (defined($params->{$p}));

  # INPUT_TABLE: the table from which to retrieve alignments.
  $p = 'input_table';
  $self->{$p} = $params->{$p} if (defined($params->{$p}));

  # OUTPUT_TABLE: the table to send SLR output to.
  $p = 'output_table';
  $self->{$p} = $params->{$p} if (defined($params->{$p}));

  # GENCODE: the genetic code to use for SLR.
  # Options: 'universal' or 'mammalian' mitochondrial.
  $p = 'gencode';
  $self->{$p} = $params->{$p} if (defined($params->{$p}));

  return;
}

sub run_sitewise_dNdS
{
  my $self = shift;
  print ("Sitewise_dNdS::run_sitewise_dNdS\n") if ($self->debug);

  return undef unless (defined($self->{protein_tree}));
  $self->{starttime} = time()*1000;

  my $slrexe = $self->{'slr_executable'};
  unless (-e $slrexe) {
    $slrexe = "/software/ensembl/compara/bin/Slr_ensembl";
  }
  throw("can't find an slr executable to run\n") unless(-e $slrexe);

  my $aln = $self->{'cdna'};
  my $tree_string = $self->{'newick'};

  use Bio::EnsEMBL::Compara::TreeAlignUtils;
  my $tree = Bio::EnsEMBL::Compara::TreeAlignUtils::load_tree_from_string($tree_string);
  throw("can't find cds_aln\n") if ( ! $aln );
  throw("can't find tree\n") if ( ! $tree );

  # Reorder the alignment according to the tree
  my $sorted_aln = Bio::EnsEMBL::Compara::TreeAlignUtils::sort_aln_by_tree($aln,$tree);

  # Get the number of leaf nodes.
  my $ct = 1;
  foreach my $node ($tree->get_leaf_nodes) {
      $ct++;
  }
  
  my $tmpdir = $self->worker_temp_directory;
  
  # Output the tree file.
  # We need to add a line with the num of leaves ($ct-1) and the num of trees (1)
  my $treeout = Bio::TreeIO->new('-format' => 'newick',
                                 '-file'     => ">$tmpdir/tree");
  $treeout->_print(sprintf("%d 1\n",($ct-1)));
  $treeout->write_tree($tree);
  $treeout->close();
  
  # For each masking step...
  my $results;
  my @mask_steps = @{$self->{'mask_steps'}};
  foreach my $param (@mask_steps) {
      # Mask the alignment.
      print " --> MASKING: t_coffee, $param\n" if ($self->debug);

      my $masked_aln = $self->mask_t_coffee($param,$sorted_aln);
      my $filename = ">$tmpdir/aln";

      # Count up the total number of non-Gap and non-N residues. Remember, this is a DNA alignment.
      my $nucleotide_count = 0;
      foreach my $seqI ($masked_aln->each_seq) {
	  my $seq = $seqI->seq;
	  $seq =~ s/[-N]//g;
	  $nucleotide_count += length($seq);
      }

      print "Nucleotide count: $nucleotide_count\n";
      
      # Output to a file.
      my $tmpdir = $self->worker_temp_directory;
      my $alnout = Bio::AlignIO->new
	  ('-format'      => 'phylip',
	   '-file'          => $filename,
	   '-interleaved' => 0,
	   '-idlinebreak' => 1,
	   '-idlength'    => $aln->maxdisplayname_length + 1);
      $alnout->write_aln($masked_aln);
      $alnout->close();
      undef $alnout;
      
      # now let's print the ctl file.
      my $slr_ctl = "$tmpdir/slr.ctl";
      open(SLR, ">$slr_ctl") or $self->throw("cannot open $slr_ctl for writing");
      print SLR "seqfile\: aln\n";
      print SLR "treefile\: tree\n";
      my $outfile = "slr.res";
      print SLR "outfile\: $outfile\n";
      $self->{'saturated'} = 0.001;
      print SLR "saturated\: ". $self->{'saturated'} . "\n";
      print SLR "gencode\: ". $self->{'gencode'} . "\n";
      print SLR "aminof\: 1\n"; # aminof
      close(SLR);
      
      eval {
	  my $cwd = cwd();
	  my $exit_status;
	  chdir($tmpdir);
	  my $run;
	  my $quiet = '';
	  $quiet = ' 2>/dev/null' unless ($self->debug);

	  # Run the SLR command and gather the output.
	  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(1);
	  open($run, "$slrexe $quiet |") or $self->throw("Cannot open exe $slrexe");
	  my @output = <$run>;
	  $exit_status = close($run);
	  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(0);

	  foreach my $line (@output) {
	      print $line;
	  }
	  
	  if ( (grep { /is saturated/ } @output)) {
	      print ("  -> Tree is dS saturated!\n") if ($self->debug);
	      # Find the highest branchwise dS value.
	      my @branch_lengths = ();
	      foreach my $line (grep { /is saturated/ } @output) {
		  # SLR will list all the branches saturated above our threshold.
		  $line =~ /length = (\S+)/;
		  push @branch_lengths, sprintf("%2.3f",$1);
	      }
	      
	      # Use the kmedian subroutine (at bottom of this script) to gather the min, median, max.
	      my ($min,$max,$median) = kmedian(\@branch_lengths);
	      # Calculate the mean.
	      my $sum = 0;
	      foreach my $val (@branch_lengths) {
		  $sum += $val;
	      }
	      my $mean = $sum / scalar(@branch_lengths);

	      print "$max  $min $median $mean\n";
	      my $values = {'max' => $max,
			    'min' => $min,
			    'median' => $median,
			    'mean' => $mean,
			    'nucleotide_count' => $nucleotide_count
			};
	      $results->{$param} = $values;
	  }
	  if ( (grep { /\berr(or)?: /io } @output)  || !$exit_status) {
	      warn("There was an error running SLR!!");
	  }      
	  chdir($cwd);
      };
  }

  $self->{'results'} = $results;
}


sub check_job_fail_options
  {
    my $self = shift;

    if ($self->input_job->retry_count >= 2) {
      $self->dataflow_output_id($self->input_id, 2);
      $self->input_job->update_status('FAILED');

      if ($self->{'protein_tree'}) {
        $self->{'protein_tree'}->release_tree;
        $self->{'protein_tree'} = undef;
      }
      throw("Job failed >=3 times: try something else and FAIL it");
    }
  }


sub mask_t_coffee
{
    my $self = shift;
    my $threshold = shift;
    my $cdna_aln = shift;

    # Collect the hash_ref of score strings from the database.
    my $table = $self->{'input_table'};
    my $score_table = $table . '_score';
    my $tree = $self->{'protein_tree'};
    
    my $pta = $self->{'comparaDBA'}->get_ProteinTreeAdaptor;

    my $hash_ref;
    foreach my $leaf (@{$tree->get_all_leaves}) {
	# Grab the score line for each leaf node.
	my $id = $leaf->stable_id; # Must be stable_id to match the aln object.
	my $member_id = $leaf->member_id;
	my $cmd = "SELECT cigar_line FROM $score_table where member_id=$member_id;";
	my $sth = $pta->prepare($cmd);
	$sth->execute();
	my $data = $sth->fetchrow_hashref();
	$sth->finish();
	my $scores = $data->{'cigar_line'};

	#print $id."\t".$scores ."\n";

	# Convert the protein mask into a DNA-level mask by repeating each char 3 times.
	my @arr = split(//,$scores);
	@arr = map { ($_ . '') x 3 } @arr;
	$scores = join("",@arr);
	$hash_ref->{$id} = $scores;
    }
    my $new_aln = $cdna_aln->mask_below_score($threshold,$hash_ref);
    return $new_aln;
}

sub mask_gblocks
  {
    my $self = shift;
    my $gmin = shift;
    my $aln  = shift;

    printf("Sitewise_dNdS::run_gblocks\n") if($self->debug);

    throw("Sitewise_dNdS : error getting Peptide SimpleAlign") unless (defined($aln));

    my $aln_length = $aln->length;
    my $tree_id = $self->{'protein_tree'}->node_id;
    my $tmpdir = $self->worker_temp_directory;
    my $filename = "$tmpdir". "$tree_id.fasta";
    my $tmpfile = Bio::AlignIO->new
      (-file => ">$filename",
       -format => 'fasta');
    $tmpfile->write_aln($aln);
    $tmpfile->close;
    my $min_leaves_gblocks = int(($self->{'protein_tree'}->num_leaves+1) * $gmin + 0.5);
    my $cmd = "echo -e \"o\n$filename\n\nb\n2\n$min_leaves_gblocks\n5\n5\ng\nm\nq\n\" | /software/ensembl/compara/bin/Gblocks 2>/dev/null 1>/dev/null";
    my $ret = system("$cmd");
    open FLANKS, "$filename-gb.htm" or die "$!\n";
    my $segments_string;
    while (<FLANKS>) {
      chomp $_;
      next unless ($_ =~ /Flanks/);
      $segments_string = $_;
      last;
    }
    close FLANKS;
    $segments_string =~ s/Flanks\: //g;
    $segments_string =~ s/\s+$//g;

    $self->{flanks} = $segments_string;
    $self->{'protein_tree'}->store_tag('Gblocks_flanks', $segments_string);
  }


sub store_sitewise_dNdS
{
  my $self = shift;

  my $results = $self->{'results'};
  my $root_id = $self->{'protein_tree'}->node_id;
  
  ####

  my @mask_steps = @{$self->{'mask_steps'}};
  foreach my $param (@mask_steps) {
      my $obj = $results->{$param}; # Object holds min, max, median, mean.

      my $sth;
      my $table = 'sitewise_sweep';
      if (defined($self->{'output_table'})) {
	  $table = $self->{'output_table'};
      }
      
      my $mask_level = $param;
      my $max_dS = $obj->{'max'};
      my $median_dS = $obj->{'median'};
      my $mean_dS = $obj->{'mean'};
      my $nucleotide_count = $obj->{'nucleotide_count'};

      print "$max_dS $median_dS $mean_dS $nucleotide_count\n";
      
      $sth = $self->{'comparaDBA'}->dbc->prepare
	  ("REPLACE INTO $table
                           (node_id,
                            mask_level,
                            max_dS,
                            median_dS,
                            mean_dS,
                            nucleotide_count
                            ) VALUES (?,?,?,?,?,?)");
      $sth->execute($root_id,$mask_level,$max_dS,$median_dS,$mean_dS,$nucleotide_count);
      $sth->finish();
      print "Wrote $root_id $param to table!\n";
  }

  return undef;
}


#---------------------------------------------------------------
#-- kmedian.pl - Function to calculate minimum, maximum, and median of
#-- an array of numbers.  Runs in expected linear time and reasonable
#-- care has been taken in its implementation.
#--
#-- This code is released into the public domain.  Original
#-- author is Joseph B. Kopena (tjkopena@cs.drexel.edu), December 2006.
#---------------------------------------------------------------

#---------------------------------------------------------------
#-------------------------------------------------------
#-- (min, max, median) kmedian(array reference)
#--
#-- Determine the minimum, maximum, and median values for the given
#-- array of numbers.  Duplicate values in the array are acceptable.
#-- The order of the array will not be preserved---a copy is not made,
#-- and neither will it be fully sorted.
sub kmedian {

    my $a = shift;

    # Bail early if there's nothing going on... -tj
    if ($#{@{$a}} == -1) {
        # Should throw an error or something... -tj
	return (0, 0, 0);
    } elsif ($#{@{$a}} == 0) {
	return ($a->[0], $a->[0], $a->[0]);
    } elsif ($#{@{$a}} == 1) {
	if ($a->[0] <= $a->[1]) {
	    return ($a->[0], $a->[1], ($a->[0]+$a->[1])/2);
	} else {
	    return ($a->[1], $a->[0], ($a->[0]+$a->[1])/2);
	}
      # end there is more than 0, 1, or 2 elements
    }

    # mys cost, so push them until they're needed. -tj
    my $left = 0;
    my $right = $#{@{$a}};

    # Pick a random location for the first pivot index, which we'll
    # keep in $k so we can increment and not reassign the first
    # $pivotIndex. -tj
    my $k = int(rand($#{@{$a}}+1)), my $pivotIndex = 0, my $i;

    # Initialize min and max to a value in the table.
    # $a->[$k] will not be compared in the following loop, so
    # double bonus for using it as the initializer. -tj
    my $min = $a->[$k];
    my $max = $a->[$k];

    swap(\$a->[$k], \$a->[$right]);

    # This whole section is unrolled instead of using kpartition()
    # because we only have to do the min/max comparisons on the first
    # partitioning loop. -tj
    for ($i = 0; $i < $right; $i++) {
	if ($a->[$i] <= $a->[$right]) {
	    if ($a->[$i] < $min) {
		$min = $a->[$i];
	    }
	    swap(\$a->[$pivotIndex], \$a->[$i]);
	    $pivotIndex++;
	} elsif ($a->[$i] > $max) {
	    $max = $a->[$i];
	}

      # end looping over for initial partition
    }

    swap(\$a->[$right], \$a->[$pivotIndex]);

    # This algorithm is based off general selection, but we just want
    # the median. -tj
    $k = int($right/2);

    while ($pivotIndex != $k) {
	if ($k < $pivotIndex) {
	    $right = $pivotIndex - 1;
	} else {
	    $left = $pivotIndex + 1;
	}

        # Partition the subarray.
	swap(\$a->[$right], \$a->[$left+int(rand(($right-$left+1)))]);
	$pivotIndex = $left;
	for ($i = $left; $i < $right; $i++) {
	    if ($a->[$i] <= $a->[$right]) {
		swap(\$a->[$pivotIndex], \$a->[$i]);
		$pivotIndex++;
	    }
	}
	swap(\$a->[$right], \$a->[$pivotIndex]);

      # end looping through the array
    }

    if ($#{@{$a}} % 2) {
      # It's even---$#a = length-1.  Find the successor of the value
      # at the current pivotIndex and then return # the mean of the
      # two as the median. -tj

	$k = $max; # Note that max will either be on the right of
                   # the pivot or will be the pivot. -tj
	for ($i = $pivotIndex+1; $i <= $#{@{$a}}; $i++) {
	    if ($a->[$i] < $k) {
		$k = $a->[$i];
	    }
	}

	return ($min, $max, ($a->[$pivotIndex] + $k)/2);
    }

    return ($min, $max, $a->[$pivotIndex]);

# end kmedian
}

#-------------------------------------------------------
sub swap {
    my $a = ${$_[0]};
${$_[0]} = ${$_[1]};
${$_[1]} = $a;
}


1;
