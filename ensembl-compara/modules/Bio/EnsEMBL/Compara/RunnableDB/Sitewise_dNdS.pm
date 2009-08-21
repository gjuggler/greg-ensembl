package Bio::EnsEMBL::Compara::RunnableDB::Sitewise_dNdS;

use strict;
use Time::HiRes qw(time gettimeofday tv_interval sleep);
use Cwd;
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::NestedSet;
use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Hive::Process;

our @ISA = qw(Bio::EnsEMBL::Hive::Process);

#
# Some global-ish variables.
#
my $dba;
my $pta;

# INPUT FILES / OBJECTS.
my $tree;
my $input_cdna;
my $input_aa;
my $params;

# OUTPUT FILES / OBJECTS / STATES.
my $results;
my @sitewise_omegas;
my $cluster_too_small;
my $dont_write_output;
my $aa_map;

sub debug {1;}

sub fetch_input {
  my( $self) = @_;
  
  # Load up the Compara DBAdaptor.
#  if ($self->dba) {
#    $dba = $self->dba;
#  } else {
    $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-DBCONN=>$self->db->dbc);
#  }
  $pta = $dba->get_ProteinTreeAdaptor;
  
  ### DEFAULT PARAMETERS ###
  $params = {
    input_table            => 'protein_tree_member',
    output_table           => 'sitewise_aln',
    parameter_set_id       => 1,

    alignment_quality_filtering => 1,
    sequence_quality_filtering => 1,

    parameter_sets         => '1',

    # SLR Parameters
    gencode                => 'universal',
    aminof                 => 1,                    # Options: 0, 1, 2 (default 1)
    codonf                 => 0,                    # Options: 0, 1, 2 (default 0)
    freqtype               => 0,                    # Options: 0, 1, 2 (default 0)

    # Actions
    action                 => 'sitewise',           # Which action to perform.
                                                    # 'sitewise' - Infer sitewise omega values
                                                    # 'reoptimise' - Reoptimise branch lengths using codon model
    };
  
  # For aminof, codonf, and freqtype, see the SLR readme.txt for more info.

  #########################

  print "TEMP: ".$self->worker_temp_directory."\n";
  print "PARAMS: ".$self->parameters."\n";
  print "INPUT_ID: ".$self->input_id."\n";

  # Fetch parameters from the two possible locations. Input_id takes precedence!

  $params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_string($params,$self->parameters);
  $params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_string($params,$self->input_id);

  #########################

  $self->check_job_fail_options;
  #$self->check_if_exit_cleanly;

  $dba->dbc->disconnect_when_inactive(1);
}

sub run {
  my $self = shift;
  
  my $node_id = $params->{'node_id'};
  my $sth = $dba->dbc->prepare("select * from sitewise_aln where node_id=$node_id limit 1;");
  $sth->execute();
  if ($sth->fetchrow_arrayref) {
    #print "HAS SOMETHING ALREADY!\n";
    #$dont_write_output = 1;
    $sth->finish;
    #return undef;
  }

  my @param_sets = split(",",$params->{'parameter_sets'});
  delete $params->{'parameter_sets'};
  foreach my $param_set (@param_sets) {
    print "PARAM SET: $param_set\n";
    my $param_set_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_param_set($dba->dbc,$param_set);
    my $new_params = Bio::EnsEMBL::Compara::ComparaUtils->replace_params($params,$param_set_params);
    #print "PARAMS: ".Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($new_params)."\n";

    # Load the tree. This logic has been delegated to ComparaUtils.
    $tree = Bio::EnsEMBL::Compara::ComparaUtils->get_tree_for_comparative_analysis($dba,$new_params);
    $tree = $tree->minimize_tree;
    
    # Think of reasons why we want to fail the job.
    my @leaves = $tree->leaves;
    if (scalar(@leaves) <= 3) {
      #$self->fail_job("Tree contains <= 3 leaves -- that's way too small for sitewise analysis!");
      next;
    } elsif (scalar(@leaves) > 300) {
      #$self->fail_job("Tree contains more than 300 leaves -- that's too big!!");
      next;
    }
    foreach my $leaf (@leaves) {
      if ($leaf->distance_to_parent > 100) {
	$self->fail_job("Tree contains a massive branch length -- probably an error in b.l. optimization!");
      }
    }
    
    #$tree->print_tree;
    #eval {
      $self->run_with_params($new_params,$tree);
      $tree->release_tree;
    #};
  }
}

sub run_with_params {
  my $self = shift;
  my $params = shift;
  my $tree = shift;
  $self->check_job_fail_options;

  print "Getting alignments...\n";
  my $aa_aln = $tree->get_SimpleAlign();
  my $cdna_aln = $tree->get_SimpleAlign(-cdna => 1);

  $input_aa = Bio::EnsEMBL::Compara::ComparaUtils->fetch_masked_alignment($aa_aln,$cdna_aln,$tree,$params,0);
  $input_cdna = Bio::EnsEMBL::Compara::ComparaUtils->fetch_masked_alignment($aa_aln,$cdna_aln,$tree,$params,1);

  exit(0);

  #if ($self->input_job->retry_count > 0) {
    #my $arr_ref = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns($input_aa);
    #my @arr = @{$arr_ref};
    #$input_aa = $arr[0];
    #$aa_map = $arr[1];
    #$input_cdna = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns_in_threes($input_cdna);
  #}

  Bio::EnsEMBL::Compara::AlignUtils->pretty_print($input_aa,{length => 150});
  Bio::EnsEMBL::Compara::AlignUtils->pretty_print($input_cdna,{length => 150});

  if ($params->{'action'} eq 'gblocks') {
    my $blocks = $self->run_gblocks($tree,$input_aa,$params);
    $self->store_gblocks($tree,$input_aa,$blocks);
  } else {
    my $results = $self->run_sitewise_dNdS($tree,$input_cdna,$params);
    $self->store_sitewise_dNdS($results,$tree,$params);
  }
}

sub run_sitewise_dNdS
{
  my $self = shift;
  my $tree = shift;
  my $cdna_aln = shift;
  my $params = shift;

  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
  
  # LOAD VARIABLES FROM PARAMS.
  my $slrexe = $params->{'slr_executable'};
  my $gencode = $params->{'gencode'};
#  my $saturated = $params->{'saturate_ds'};
  my $aminof = $params->{'aminof'};
  my $codonf = $params->{'codonf'};
  my $freqtype = $params->{'freqtype'};

#  $slrexe = "/nfs/users/nfs_g/gj1/bin/Slr_ensembl" if (! -e $slrexe);
#  $slrexe = "/nfs/users/nfs_g/gj1/bin/Slr" if (! -e $slrexe);
#  $slrexe = "/nfs/users/nfs_g/gj1/bin/Slr2" if (! -e $slrexe);
#  $slrexe = "/nfs/users/nfs_g/gj1/bin/Slr_Linux_static" if (! -e $slrexe);
  $slrexe = "/nfs/users/nfs_g/gj1/build_scratch/slr/bin/Slr";

  throw("can't find an slr executable to run\n") if (!-e $slrexe);

  # Reorder the alignment according to the tree
  $cdna_aln = Bio::EnsEMBL::Compara::AlignUtils->sort_by_tree($cdna_aln,$treeI);
  
  my $num_leaves = scalar(@{$tree->get_all_leaves});
  my $tmpdir = $self->worker_temp_directory;

  # CLEAN UP OLD RESULTS FILES.
  unlink "$tmpdir/slr.res";
  unlink "$tmpdir/tree";
  unlink "$tmpdir/aln";
  unlink "$tmpdir/slr.ctl";

  # OUTPUT THE ALIGNMENT.
  my $alnout = Bio::AlignIO->new
    ('-format'      => 'phylip',
     '-file'          => ">$tmpdir/aln",
     '-interleaved' => 0,
     '-idlinebreak' => 1,
     '-idlength'    => $cdna_aln->maxdisplayname_length + 1);
  $alnout->write_aln($cdna_aln);
  $alnout->close();
  
  # OUTPUT THE TREE.
  my $tree_newick = Bio::EnsEMBL::Compara::TreeUtils->to_newick($tree);
  open(OUT,">$tmpdir/tree");
  print OUT sprintf("%d 1\n",$num_leaves);
  print OUT $tree_newick . "\n";
  close(OUT);
  
  # OUTPUT THE CTL FILE.
  my $slr_ctl = "$tmpdir/slr.ctl";
  open(SLR, ">$slr_ctl") or $self->throw("cannot open $slr_ctl for writing");
  print SLR "seqfile\: aln\n";
  print SLR "treefile\: tree\n";
  my $outfile = "slr.res";
  print SLR "outfile\: $outfile\n";
#  print SLR "saturated\: $saturated\n";
  print SLR "gencode\: $gencode\n";
  print SLR "aminof\: $aminof\n";
  print SLR "codonf\: $codonf\n";
  print SLR "skipsitewise\: 1\n" if ($params->{'action'} eq 'reoptimise');
  print SLR "freqtype\: $freqtype\n";
  print SLR "seed\: 1\n";
  close(SLR);
  
  my $rc = 1;
  my $results;
  my $error_string;
  {
    my $cwd = cwd();
    chdir($tmpdir);
    my $exit_status = 0;

    #
    # Run SLR and grab the output.
    #
    my $run;
    my $quiet = '';
    $quiet = ' 2>/dev/null' unless ($self->debug);

    #print "RETRY COUNT: ".$self->input_job->retry_count."\n";
    #my $prefix = "";
    #if ($self->input_job->retry_count > 0) {
    my $prefix= "export MALLOC_CHECK_=1;";
    #}

    print "Running!\n";
    open($run, "$prefix $slrexe |") or $self->throw("Cannot open exe $slrexe");
    my @output = <$run>;

    $exit_status = close($run);
    #$dba->dbc->disconnect_when_inactive(0);
    $error_string = (join('',@output));

    #print $error_string."\n";

    # GJ 2009-02-22: Updated and fleshed out the saturated tree handling.
    if ( (grep { /is saturated/ } @output)) {
      print ("  -> Tree is dS saturated!\n") if ($self->debug);
      
      # Find the highest branchwise dS value.
      my $max = $params->{'saturated_ds'};
      foreach my $line (grep { /is saturated/ } @output) {
	$line =~ /length = (\S+)/;
	$max = $1 if ($1 > $max);
      }

      $results->{'tree_is_saturated'} = 1;
      $results->{'max_ds_branch'} = $max;
      $results->{'trees'} = [];
      if (-e "$tmpdir/subtrees.out") {
	my $treeio = Bio::TreeIO->new
	  ('-format' => 'newick',
	   '-file'   => "$tmpdir/subtrees.out");
	while ( my $tree = $treeio->next_tree ) {
	  print "Found saturated tree!\n";
	  push(@{$results->{'trees'}}, $tree);
	}
      }
      chdir($cwd);
      return $results;
    }
    
    foreach my $outline (@output) {
      if ($outline =~ /lnL = (\S+)/) {
	$results->{"lnL"} = $1;
      }
      if ($outline =~ /kappa = (\S+)/) {
	$results->{'kappa'} = $1;
      }
      if ($outline =~ /omega = (\S+)/) {
	$results->{'omega'} = $1;
      }
    }

    if ( (grep { /\berr(or)?: /io } @output)  || !$exit_status) {
      print "@output\n";
      throw("There was an error running SLR!");
      $rc = 0;
    }

    # Reoptimise action.
    if ($params->{'action'} eq 'reoptimise') {
      my $next_line_is_it = 0;
      foreach my $outline (@output) {
	if ($next_line_is_it) {
	  my $new_tree = $outline;
	  my $new_pt = Bio::EnsEMBL::Compara::TreeUtils->from_newick($new_tree);
	  print "Original tree   : ".$tree->newick_format."\n";
	  print "Reoptimised tree: ".$new_pt->newick_format."\n";
	  
	  # Store new branch lengths back in the original protein_tree_node table.
	  foreach my $leaf ($tree->leaves) {
	    my $new_leaf = $new_pt->find_leaf_by_name($leaf->name);
	    $leaf->distance_to_parent($new_leaf->distance_to_parent);
	    $leaf->store;
	  }
	    print "  -> Branch lengths updated!\n";
	    $dont_write_output = 1;
	    return;
	}
	# We know that the new tree will show up just after the LnL line.
	$next_line_is_it = 1 if ($outline =~ /lnL =/);
      }
    }
    
    if (!-e "$tmpdir/$outfile") {
      throw("Error running SLR (no output file!!)\n");
    }

    #eval {
    open RESULTS, "$tmpdir/$outfile" or die "couldnt open results file: $!";
    my $okay = 0;
    my @sites;
    my $type = '';
    my $note = '';
    while (<RESULTS>) {
      $type = '';
      $note = '';
      chomp $_;
      print $_."\n";
      if ( /^\#/ ) {
	next;
      }

      # Parse the notes.
      if ( /All gaps/i ) {
	$note = 'all_gaps';
      } elsif ( /Constant/i ) {
	$note = 'constant';
      } elsif ( /Single char/i ) {
	$note = 'single_char';
      } elsif ( /Synonymous/i ) {
	$note = 'synonymous';
      }	elsif (/\!/) {                # Make sure to parse the random note last, because it's sometimes shown alongside other notes.
	$note = 'random';
      }
      
      # Parse the pos/sel test results.
      if ( /\+\+\+\+\s+/ ) {
	$type = 'positive4';
      } elsif ( /\+\+\+\s+/ ) {
	$type = 'positive3';
      } elsif ( /\+\+\s+/ ) {
	$type = 'positive2';
      } elsif ( /\+\s+/ ) {
	$type = 'positive1';
      } elsif ( /\-\-\-\-\s+/ ) {
	$type = 'negative4';
      } elsif ( /\-\-\-\s+/ ) {
	$type = 'negative3';
      } elsif ( /\-\-\s+/ ) {
	$type = 'negative2';
      } elsif ( /\-\s+/ ) {
	$type = 'negative1';
      }
      if ( /^\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/ ) {
	my @arr = ($1,$2,$3,$4,$5,$6,$7,$8,$9,10,$type,$note);
	push @sites, \@arr;
	#print join(",",@arr)."\n";
      } else {
	warn("error parsing the results: $_\n");
      }
    }
    $results->{'sites'} = \@sites;
    close RESULTS;
    #};
    if ( $@ ) {
      #warn($error_string);
      warn($@);
    }
    chdir($cwd);
  }
  
  $results->{'rc'} = $rc;
  return $results;
}


sub write_output {
  my $self = shift;
}


sub check_job_fail_options
{
  my $self = shift;

  #print "Checking retry count: ".$self->input_job->retry_count."\n";
  if ($self->input_job->retry_count > 3) {
    $self->fail_job("Job failed >3 times!");
  }
}

sub fail_job {
  my $self = shift;
  my $reason = shift;

  $reason = "None" unless (defined $reason);

  my $input_job = $self->input_job;
  if (defined $input_job) {
    $input_job->update_status('FAILED');
    #print "Failing job! Reason given: $reason\n";
    $dont_write_output = 1;
    $self->DESTROY;
    throw("Failing job! Reason given: $reason");
  }
}

sub run_gblocks
{
  my $self = shift;
  my $tree = shift;
  my $aln  = shift;
  my $params = shift;

  my $defaults = {
    t    => 'p',        # Type of sequence (p=protein,c=codon,d=dna)
    b3   => '5000',       # Max # of contiguous nonconserved positions
    b4   => '3',        # Minimum length of a block
    b5   => 'a',        # Allow gap positions (n=none, h=with half,a=all)
  };
  my $use_params = Bio::EnsEMBL::Compara::ComparaUtils->replace_params($defaults,$params);

  printf("Sitewise_dNdS::run_gblocks\n") if($self->debug);

  my $aln_length = $aln->length;
  my $tmpdir = $self->worker_temp_directory;
  my $filename = "$tmpdir". "gblocks_aln.fasta";
  my $tmpfile = Bio::AlignIO->new
    (-file => ">$filename",
     -format => 'fasta');
  $tmpfile->write_aln($aln);
  $tmpfile->close;

  my @leaves = $tree->leaves;
  my $num_leaves = scalar(@leaves);
  my $min_leaves_gblocks = int(($num_leaves+1)/2 + 0.5);
  #Example command: Gblocks 2138.fasta -t=p -b2=20 -b3=50 -b4=3 -b5=a -p=s
  my $cmd = sprintf("Gblocks %s -t=%s -b1=%s -b2=%s -b3=%s -b4=%s -b5=%s -p=s\n",
		    $filename,
		    $use_params->{'t'},
		    $min_leaves_gblocks,
		    $min_leaves_gblocks,
		    $use_params->{'b3'},
		    $use_params->{'b4'},
		    $use_params->{'b5'});
  print "GBLOCKS: $cmd\n";
  my $ret = system("$cmd");
  open FLANKS, "$filename-gb.txts" or die "$!\n";
  my $segments_string;
  while (<FLANKS>) {
    chomp $_;
    next unless ($_ =~ /Flanks:/);
    $segments_string = $_;
    last;
  }
  close FLANKS;
  $segments_string =~ s/Flanks\: //g;
  $segments_string =~ s/\s+$//g;

  print $segments_string . "\n";
  return $segments_string;
}


sub store_gblocks {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $blocks_str = shift;

  $_ = $blocks_str;
  my @bs = /\[(.+?)\]/g;

  my $sth = $tree->adaptor->prepare
    ("REPLACE INTO gblocks_aln
                           (parameter_set_id,
                            aln_start,
                            aln_end,
                            node_id
                            ) VALUES (?,?,?,?)"
     );

  my %sites;
  foreach my $b (@bs) {
    my ($start,$end) = split(" ",$b);
      $sth->execute($params->{'parameter_set_id'},
		    $start,
		    $end,
		    $tree->node_id
		    );
      sleep(0.05);
  }
  $sth->finish();
}

sub store_sitewise_dNdS
{
  my $self = shift;
  my $results = shift;
  my $tree = shift;
  my $params = shift;

  my $tree_node_id = 0;
  $tree_node_id = $tree->subroot->node_id if (defined $tree->subroot);
  my $node_id = $tree->node_id;

  foreach my $site (@{$results->{'sites'}}) {
    sleep(0.08);

    # Site  Neutral  Optimal   Omega    lower    upper LRT_Stat    Pval     Adj.Pval    Q-value Result Note
    # 1     4.77     3.44   0.0000   0.0000   1.4655   2.6626 1.0273e-01 8.6803e-01 1.7835e-02  --     Constant;
    # 0     1        2      3        4        5        6      7          8          9           10     11
    my ($site, $neutral, $optimal, $omega, $lower, $upper, $lrt_stat, $pval, $adj_pval, $q_value, $type, $note) = @{$site};

    my $mapped_site = $site;
    my $original_site = $site;
    if (defined $aa_map) {
      $mapped_site = $aa_map->{$site};
      if (defined $mapped_site) {
	print "Returned column: $site   mapped: $mapped_site\n";
	$original_site = $site;
	$site = $mapped_site;
      }
    }

    # These two cases are pretty useless, so we'll skip 'em.
    next if ($note eq 'all_gaps');
    next if ($note eq 'single_char');

    my $nongaps = Bio::EnsEMBL::Compara::AlignUtils->get_nongaps_at_column($input_aa,$original_site);
    next if ($nongaps == 0);

    printf("Site: %s  nongaps: %d  omegas: %3f (%3f - %3f)  note: %s \n",$site,$nongaps,$omega,$lower,$upper,$note);

    my $table = 'sitewise_aln';

    my $parameter_set_id = 0;
    if (defined($params->{'parameter_set_id'})) {
      $parameter_set_id = $params->{'parameter_set_id'};
    }

    my $sth = $tree->adaptor->prepare
      ("REPLACE INTO $table
                           (parameter_set_id,
                            aln_position,
                            node_id,
                            tree_node_id,
                            omega,
                            omega_lower,
                            omega_upper,
                            lrt_stat,
                            ncod,
                            type,
                            note) VALUES (?,?,?,?,?,?,?,?,?,?,?)");
    $sth->execute($parameter_set_id,
		  $site,
		  $node_id,
		  $tree_node_id,
		  $omega,
		  $lower,
		  $upper,
		  $lrt_stat,
		  $nongaps,
		  $type,
		  $note);
    $sth->finish();
  }
}

sub DESTROY {
    my $self = shift;
    $tree->release_tree if ($tree);
    $tree = undef;
    $self->SUPER::DESTROY if $self->can("SUPER::DESTROY");
}


1;
