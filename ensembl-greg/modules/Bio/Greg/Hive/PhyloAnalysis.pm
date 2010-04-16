package Bio::Greg::Hive::PhyloAnalysis;

use strict;
use Bio::EnsEMBL::Utils::Argument;
use Bio::EnsEMBL::Utils::Exception;

use Time::HiRes qw(sleep);
use Cwd;
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::NestedSet;

use Bio::Greg::Codeml;

use base ('Bio::Greg::Hive::Process');

#
# Some global-ish variables.
#

sub debug {1;}

sub fetch_input {
  my( $self) = @_;
    
  ### DEFAULT PARAMETERS ###
  my $params = {
    alignment_table            => 'protein_tree_member',
    omega_table           => 'sitewise_omega',
    parameter_set_id       => 1,

    sitewise_store_opt_tree         => 1,
    sitewise_store_gaps             => 0,
    sitewise_minimum_leaf_count     => 3,
    sitewise_strip_gaps             => 0,
    sitewise_parameter_sets         => 'all',
    analysis_action                 => 'slr',                # Which action(s) to perform. Space-delimited.
                                                    # 'slr' - SLR sitewise omegas.
                                                    # 'paml_sitewise' - PAML sitewise omegas.
                                                    # 'wobble' - SLR_wobble test.
                                                    # 'slr_reoptimise' - reoptimise with SLR
                                                    # 'paml_reoptimise' - reoptimise with PAML
                                                    # 'paml_lrt' - likelihood ratio test with PAML
                                                    # 'hyphy_dnds' - Calculate dn/ds with HYPHY, along with confidence intervals.

    # SLR Parameters
    slr_gencode                => 'universal',
    slr_aminof                 => 1,                    # Options: 0, 1, 2 (default 1)
    slr_codonf                 => 0,                    # Options: 0, 1, 2 (default 0)
    slr_freqtype               => 0,                    # Options: 0, 1, 2 (default 0)
    slr_skipsitewise           => 0,                    # Options: 0, 1 (default 0)
    slr_malloc_check           => 0,

    # PAML Parameters
    #paml_model                  => 'M3',                 # Used for Bayes Empirical Bayes sitewise analysis.
    #paml_model_b                => 'M7',                 # Used for the likelihood ratio tests.

    alignment_score_mask_character => 'N',
    sequence_quality_mask_character_character => 'N'

    };
  
  $self->load_all_params($params);

  $self->fail_on_single_member();
}

sub run {
  my $self = shift;  

  $self->compara_dba->dbc->disconnect_when_inactive(1);

  my $tree = $self->get_tree;

  print $tree->newick_format."\n";

  # Think of reasons why we want to fail the job.
  my @leaves = $tree->leaves;
  if (scalar(@leaves) < $self->param('sitewise_minimum_leaf_count')) {
    my $value = sprintf "Too small (%d < %d)",scalar(@leaves),$self->param('sitewise_minimum_leaf_count');
    $self->store_tag("slr_skipped",$value);
    next;
  } elsif (scalar(@leaves) > 300) {
    my $value = sprintf "Too big (%d > %d)",scalar(@leaves),300;
    $self->store_tag("slr_skipped",$value);
    next;
  }
  foreach my $leaf (@leaves) {
    if ($leaf->distance_to_parent > 100) {
      $leaf->distance_to_parent(4);
      #$self->fail_job("Tree contains a massive branch length -- probably an error in b.l. optimization!");
    }
  }
  
  #eval {
    $self->run_with_params($self->params,$tree);
  #};
  if ($@) {
    print "ERROR - Trying with removed gaps...\n";
    $self->param('sitewise_strip_gaps',1);
    
    eval {
      $self->run_with_params($self->param,$tree);
    };
    if ($@) {
      print "ERROR - Trying with malloc_check...\n";
      
      $self->param('sitewise_malloc_check',1);
      
      eval {
	$self->run_with_params($self->params,$tree);
      };
      if ($@) {
	my $param_set = $self->param('parameter_set_id');
	print "ERROR - GIVING UP! Parameter set $param_set\n";
	$self->store_tag('slr_error',1);
	sleep(2);
      }
    }
  }
  $tree->release_tree;
}

sub run_with_params {
  my $self = shift;
  my $params = shift;
  my $tree = shift;

  print "Getting alignments...\n";

  foreach my $leaf ($tree->leaves) {
    #print $leaf->member_id."\n";
  }

  my $input_aa = $self->get_aln;
  my $input_cdna = $self->get_cdna_aln;

  my ($slim_cdna,$cdna_new_to_old,$cdna_old_to_new) = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns_in_threes($input_cdna);
  my ($slim_aa,$aa_new_to_old,$aa_old_to_new) = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns($input_aa);

  my $input_aa = $slim_aa;
  my $input_cdna = $slim_cdna;
  $self->param('input_cdna',$slim_cdna);
  $self->param('input_aa',$slim_aa);
  $self->param('aln_map_aa',$aa_new_to_old);
  $self->param('aln_map_cdna',$cdna_new_to_old);

  Bio::EnsEMBL::Compara::AlignUtils->pretty_print($input_aa,{length => 200});
  Bio::EnsEMBL::Compara::AlignUtils->pretty_print($input_cdna,{length => 200});

  my $action = $self->param('analysis_action');

  if ($action =~ m/slr/i) {
    $self->run_sitewise_dNdS($tree,$input_cdna,$self->params);
    sleep(2);
  } elsif ($action =~ m/paml/i) {
    $self->run_paml($tree,$input_cdna,$self->params);
  } elsif ($action =~ m/wobble/i) {
    $self->param('slr_wobble',0);
    my $results_nowobble = $self->run_wobble($tree,$input_cdna,$self->params);
    $self->store_tag("lnl_nowobble",$results_nowobble->{'lnL'});

    sleep(2);

    $self->param('slr_wobble',1);
    my $results_wobble = $self->run_wobble($tree,$input_cdna,$self->params);    
    $self->store_tag("lnl_wobble",$results_wobble->{'lnL'});
  } elsif ($action =~ m/hyphy_dnds/i) {
    $self->run_hyphy($tree,$input_cdna,$self->params);
  } elsif ($action =~ m/xrate_indels/i) {
    $self->run_xrate_indels($tree,$input_cdna,$self->params);
  } elsif ($action =~ m/indelign/i) {
    $self->run_indelign($tree,$input_cdna,$self->params);
  }

}

sub run_hyphy {
  my $self = shift;
  my $tree = shift;
  my $cdna_aln = shift;
  my $params = shift;

  my $cwd = cwd();
  my $tmpdir = $self->worker_temp_directory;
  chdir($tmpdir);

  # OUTPUT THE ALIGNMENT.
  my $aln_f = $tmpdir. "aln.fa";
  open(OUT,">$aln_f");
  foreach my $seq ($cdna_aln->each_seq) {
    my $name = $seq->id;
    print OUT ">$name\n";
    print OUT $seq->seq."\n";
  }
  close(OUT);

  # OUTPUT THE TREE.
  my $tree_f = $tmpdir . "tree.nh";
  my $tree_newick = Bio::EnsEMBL::Compara::TreeUtils->to_newick($tree);
  open(OUT,">$tree_f");
  print OUT $tree_newick . "\n";
  close(OUT);

  my $control_f = $tmpdir . "control.txt";
  my $model_filename = "gy94.bf";
  open(OUT,">$control_f");
  print OUT $aln_f."\n";
  print OUT $tree_f . "\n";
  print OUT $model_filename . "\n";
  close(OUT);

  my $batch_file = "~gj1/src/greg-ensembl/ensembl-greg/scripts/hyphy/fit_codon_model.bf";
  
  my $cmd = "HYPHY $batch_file < $control_f";
  print $cmd."\n";
  my @output = `$cmd`;
  print "@output\n";

  my $omega_lo = 0;
  my $omega_hi = 0;
  my $omega = 0;
  foreach my $line (@output) {
    $omega_lo = $1 if ($line =~ m/omega_lo:(.*)/);
    $omega_hi = $1 if ($line =~ m/omega_hi:(.*)/);
    $omega = $1 if ($line =~ m/omega:(.*)/);
  }
  
  print "$omega_lo $omega $omega_hi\n";

  my $ps = $params->{parameter_set_id};
  $self->store_tag("hyphy_dnds",$omega);
  $self->store_tag("hyphy_dnds_lo",$omega_lo);
  $self->store_tag("hyphy_dnds_hi",$omega_hi);

  chdir($cwd);

}

sub run_indelign {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $params = shift;

  my ($ins,$del,$ins_rate,$del_rate) = Bio::EnsEMBL::Compara::AlignUtils->indelign($aln,$tree,$self->params,$self->worker_temp_directory);

  print "ins: $ins_rate del: $del_rate\n";

  my $ps = $params->{parameter_set_id};
  $self->store_tag("indelign_ins",$ins_rate);
  $self->store_tag("indelign_del",$del_rate);
}

sub run_xrate_indels {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $params = shift;

  my $cwd = cwd();
  my $tmpdir = $self->worker_temp_directory;
  #my $tmpdir = "/tmp/xrate/";
  chdir($tmpdir);

  # OUTPUT THE ALIGNMENT.
  my $aln_f = $tmpdir. "aln";
  my $tree_f = $tmpdir. "tree";
  my $index_f = $tmpdir."index.txt";

  use Bio::Annotation::AnnotationFactory;
  my $factory = Bio::Annotation::AnnotationFactory->new(-type => 'Bio::Annotation::SimpleValue');

  $tree->distance_to_parent(0);
  #$tree = Bio::EnsEMBL::Compara::TreeUtils->scale($tree,1);
  my $tree_newick = Bio::EnsEMBL::Compara::TreeUtils->to_newick($tree);
  my $rfann = $factory->create_object(-value => $tree_newick, 
				      -tagname => 'NH');
  my $coll = $aln->annotation;
  map {$_->description('')} $aln->each_seq;
  $coll->add_Annotation('custom', $rfann);
  my $alnout = Bio::AlignIO->new
    ('-format'      => 'stockholm',
     '-file'          => ">$aln_f");
  $alnout->write_aln($aln);
  $alnout->close;

  # OUTPUT THE TREE.
  open(OUT,">$tree_f");
  my $string = Bio::EnsEMBL::Compara::TreeUtils->to_newick($tree);
  $string =~ s/\:[0-9\.]+;/;/g;
  print "NEWICK : $string\n";
  print OUT $string."\n";
  close(OUT);

  open(INDEX,">$index_f");
  print INDEX "test1 $tree_f $aln_f\n";
  close(INDEX);

  my $exec_folder = '/nfs/users/nfs_g/gj1/src/dart/bin';
  my $cmd = "${exec_folder}/tkfidem $index_f -log 9";
  print "$cmd\n";
  my @results = `$cmd`;

  my $lambda;
  my $mu;
  foreach my $line (@results) {
    $lambda = $1 if ($line =~ m/lambda:\s*([0-9\.]+)/);
    $mu = $1 if ($line =~ m/\s*mu:\s*([0-9\.]+)/);
  }

  print "lambda: $lambda mu: $mu\n";

  my $ps = $params->{parameter_set_id};
  $self->store_tag("xrate_ins",$lambda);
  $self->store_tag("xrate_del",$mu);

  #print "@results\n";
}

sub run_wobble {
  my $self = shift;
  my $tree = shift;
  my $cdna_aln = shift;
  my $params = shift;

  if (scalar $tree->leaves > 30) {
    $self->fail_job("Tree too large for wobbly analysis!!");
  }

  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
  
  # LOAD VARIABLES FROM PARAMS.
  my $slrexe = get_slr('executable');
  my $gencode = get_slr('gencode');
  my $aminof = get_slr('aminof');
  my $codonf = get_slr('codonf');
  my $freqtype = get_slr('freqtype');
  my $wobble = get_slr('wobble');

  $slrexe = "/nfs/users/nfs_g/gj1/bin/Slr_wobble";

  # Reorder the alignment according to the tree
  $cdna_aln = Bio::EnsEMBL::Compara::AlignUtils->sort_by_tree($cdna_aln,$treeI);
  
  my $num_leaves = scalar(@{$tree->get_all_leaves});
  my $tmpdir = $self->worker_temp_directory;

  # CLEAN UP OLD RESULTS FILES.
  unlink "$tmpdir/slr.res";
  unlink "$tmpdir/tree";
  unlink "$tmpdir/aln";
  unlink "$tmpdir/slr.ctl";

  my $tree_newick = Bio::EnsEMBL::Compara::TreeUtils->to_newick($tree);

  # Map the stable_ids into shorter integers, to stay below Slr_wobble's max label length.
  my $tree_map;
  my $i=0;
  my @leaves = $tree->leaves;
  foreach my $seq ($cdna_aln->each_seq) {
    $i++;
    my $id = $seq->id;
    map {$tree_newick =~ s/$id/$i/g if ($_->stable_id eq $id)} @leaves;
  }

  Bio::EnsEMBL::Compara::AlignUtils->pretty_print($cdna_aln,{length => 50});

  #$cdna_aln = $cdna_aln->slice(1,100);

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
  print SLR "gencode\: $gencode\n";
  print SLR "aminof\: $aminof\n";
  print SLR "codonf\: $codonf\n";
  print SLR "freqtype\: $freqtype\n";
  print SLR "seed\: 1\n";
  print SLR "wobble\: $wobble\n";
  close(SLR);

# omega<-c(4.25610800e-03,9.13459526e-02,4.01739002e-01,1.19324070e+00,3.61980510e+00)
# Pomega<-c(2.00000000e-01,2.00000000e-01,2.00000000e-01,2.00000000e-01,2.00000000e-01)
# wobble<-c(1.53800820e-02)
# Pwobble<-c(1.00000000e+00)
# lnL     4.8308754989342551e+03
#   0     3.604307e-01
# #Sitewise log-likelihoods
# Penalty = -2.290851e+00

  my $cwd = cwd();
  chdir($tmpdir);

  my $results;
  my $error_string;
  {
    my $exit_status = 0;

    my $prefix="";
#    my $prefix= "export MALLOC_CHECK_=1;";
    print "Running: $slrexe\n";
    my $run;
    open($run, "$prefix $slrexe |") or $self->throw("Cannot open exe $slrexe");
    my @output;
    while (<$run>) {
      next if ($_ =~ 'Unrecognised');
      next if ($_ =~ 'Odd gapping');
      print $_;
      push @output, $_;
    }

    foreach my $outline (@output) {
      if ($outline =~ /lnL\s+(\S+)/) {
	$results->{'lnL'} = $1;
      }
      if ($outline =~ /^omega<-c\((.*)\)/) {
	my @omega_cats = split(',',$1);
	$results->{'omega_cats'} = \@omega_cats;
      }
      if ($outline =~ /^wobble<-c\((.*)\)/) {
	my @wobble_cats = split(',',$1);
	$results->{'wobble_cats'} = \@wobble_cats;
      }
    }

    print " -> LnL: ".$results->{'lnL'}."\n";
    print " -> OMEGA: ". @{$results->{'omega_cats'}}."\n";
    print " -> WOBBLE: ". @{$results->{'wobble_cats'}}."\n";
  }


  open RESULTS, "$tmpdir/$outfile" or die "couldnt open results file: $!";
  my @sites;
  while (<RESULTS>) {
    print $_;

    if ( /(\S+)\s+(\S+)/ ) {
      my $site = $1;
      my $val = $2;
      my @arr = 0 x 12;
      $arr[3] = $val;
      push @sites,\@arr;
    }
  }
  $results->{'sites'} = \@sites;

  chdir $cwd;
  return $results;
}


sub run_sitewise_dNdS
{
  my $self = shift;
  my $tree = shift;
  my $cdna_aln = shift;
  my $params = shift;

  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
  
  # LOAD VARIABLES FROM PARAMS.
  my $slrexe = $params->{slr_executable};
  my $gencode = $params->{slr_gencode};
  my $aminof = $params->{slr_aminof};
  my $codonf = $params->{slr_codonf};
  my $freqtype = $params->{slr_freqtype};
  my $skipsitewise = $params->{slr_skipsitewise};

  $slrexe = "/nfs/users/nfs_g/gj1/bin/Slr_gj1" if (! -e $slrexe);
  $slrexe = "/homes/greg/bin/Slr" if (! -e $slrexe);

  throw("can't find an slr executable to run\n") if (!-e $slrexe);

  # Reorder the alignment according to the tree
  #$cdna_aln = Bio::EnsEMBL::Compara::AlignUtils->sort_by_tree($cdna_aln,$treeI);
  
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
  print SLR "gencode\: $gencode\n";
  print SLR "aminof\: $aminof\n";
  print SLR "codonf\: $codonf\n";
  print SLR "freqtype\: $freqtype\n";
  print SLR "skipsitewise\: $skipsitewise\n";
  print SLR "seed\: 1\n";
  close(SLR);
  
  my $results;
  
  my $cwd = cwd();
  chdir($tmpdir);
  
  #
  # Run SLR and grab the output.
  #
  my $run;
  my $quiet = '';
  #$quiet = ' 2>/dev/null' unless ($self->debug);
  
  #print "RETRY COUNT: ".$self->input_job->retry_count."\n";
  my $prefix = "";
  if ($params->{sitewise_malloc_check}) {
    $prefix= "export MALLOC_CHECK_=1;";
  } else {
    $prefix= "unset MALLOC_CHECK_;";
  }
  
  print "Running!\n";
  open($run, "$prefix $slrexe |") or $self->throw("Cannot open exe $slrexe");
  my @output = <$run>;
  
  my $exited_well = close($run);
  if (!$exited_well) {
    throw("Slr didn't exit well!");
  }

  foreach my $outline (@output) {
    if ($outline =~ /lnL = (\S+)/) {
      $results->{"lnL"} = $1;
    }
    if ($outline =~ /kappa = (\S+)/i) {
      $results->{'kappa'} = $1;
    }
    if ($outline =~ /omega = (\S+)/i) {
      $results->{'omega'} = $1;
    }
    }
  
  if ( (grep { /\berr(or)?: /io } @output)) {
    throw("There was an error running SLR!");
  }
  
  # Collect SLR's reoptimized tree.
  my $next_line_is_it = 0;
  my $new_pt;
  foreach my $outline (@output) {
    if ($next_line_is_it) {
      my $new_tree = $outline;
      $new_pt = Bio::EnsEMBL::Compara::TreeUtils->from_newick($new_tree);
      print "Original tree   : ".$tree->newick_format."\n";
      print "Reoptimised tree: ".$new_pt->newick_format."\n";
      last;
    }
    # We know that the new tree will show up just after the LnL line.
    $next_line_is_it = 1 if ($outline =~ /lnL =/);
  }
  
  # Reoptimise action.
  if ($self->param('analysis_action') =~ m/reoptimise/i) {
    # Store new branch lengths back in the original protein_tree_node table.
    foreach my $leaf ($tree->leaves) {
      my $new_leaf = $new_pt->find_leaf_by_name($leaf->name);
      $leaf->distance_to_parent($new_leaf->distance_to_parent);
      $leaf->store;
    }
    print "  -> Branch lengths updated!\n";
    return;
  }
  
  # if the setting is set, store the optimized tree as a protein_tree_tag.
  if ($params->{sitewise_store_opt_tree}) {
    my $tree_key = "slr_tree_".$params->{parameter_set_id};
    $self->store_tag($tree_key,$new_pt->newick_format);
  }  
  
  if (!$skipsitewise) {
    if (!-e "$tmpdir/$outfile") {
      throw("Error running SLR (no output file!!)\n");
    }
    open RESULTS, "$tmpdir/$outfile" or die "couldnt open results file: $!";
    my $okay = 0;
    my @sites;
    my $type = '';
    my $note = '';
    while (<RESULTS>) {
      $type = '';
      $note = '';
      chomp $_;
#      print $_."\n";
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
      } elsif (/\!/) {                # Make sure to parse the random note last, because it's sometimes shown alongside other notes.
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
  }
  chdir($cwd);
  
  # Store the results in the database.
  my $kappa_key = "slr_kappa";
  my $omega_key = "slr_omega";
  my $lnl_key = "slr_lnL";
  
  print "$omega_key ".$results->{omega}."\n";
  
  $self->store_tag($kappa_key,$results->{'kappa'});
  $self->store_tag($omega_key,$results->{'omega'});
  $self->store_tag($lnl_key,$results->{'lnL'});
  $self->store_sitewise($results,$tree,$params);
}



sub run_paml {
  my $self = shift;
  my $tree = shift;
  my $input_cdna = shift;
  my $params = shift;

  print $tree->newick_format()."\n";

  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
  $treeI->get_root_node->branch_length(0);

  if ($self->param('analysis_action') =~ m/reoptimize/i) {
    # Not sure if this works. GJ 2009-10-09
    ### Reoptimize branch lengths.
    print $tree->newick_format."\n";
    my $new_treeI = Bio::Greg::Codeml->get_m0_tree($treeI,$input_cdna);
    my $new_pt = Bio::EnsEMBL::Compara::TreeUtils->from_treeI($new_treeI);
    print $new_pt->newick_format."\n";

    # Store new branch lengths back in the original protein_tree_node table.
    foreach my $leaf ($tree->leaves) {
      my $new_leaf = $new_pt->find_leaf_by_name($leaf->name);
      $leaf->distance_to_parent($new_leaf->distance_to_parent);
      $leaf->store;
      print "  -> Branch lengths updated!\n";
    }
  }

  if ($self->param('analysis_action') =~ m/lrt/i) {

    ### Perform a likelihood ratio test between two models.
    my $model_b = $params->{'paml_model'} || $params->{'paml_model_b'};
    my $model_a = $params->{'paml_model_a'};
    
    $model_a =~ s/m//i;
    $model_b =~ s/m//i;

    my ($twice_lnL,$codeml_a,$codeml_b) = Bio::Greg::Codeml->NSsites_ratio_test($treeI,$input_cdna,$model_a,$model_b,$self->worker_temp_directory);
    my $test_label = sprintf("PAML LRT M%s-M%s",$model_a,$model_b);
    $self->store_tag($test_label,$twice_lnL);
  }

  if ($self->param('analysis_action') =~ m/sitewise/i) {
    ### Perform a BEB sitewise analysis of omegas.

    # Should we be scaling the tree? Dunno... GJ 2009-10-09
    #Bio::EnsEMBL::Compara::TreeUtils->scale($treeI,3);
    
    my $model = $params->{'paml_model'} || $params->{'paml_model_b'};
    $model =~ s/m//i;
    my $codeml_params = {
      NSsites => $model,
      ncatG => 10
    };
    my $codeml = Bio::Greg::Codeml->new( -params => $codeml_params,
                                         -alignment => $input_cdna, 
                                         -tree => $treeI,
                                         -tempdir => $self->worker_temp_directory
                                         );
    $codeml->save_tempfiles(1);
    $codeml->run();

    # Use this block for testing the parsing.
    #open(IN,"/tmp/cIIZtxAvlt_codeml/rst");
    #my @supps = <IN>;
    #close(IN);
    #my ($naive,$bayes,$se,$prob,$pos) = $codeml->extract_empirical_bayes(\@supps);

    my ($bayes,$se,$prob,$pos) = $codeml->extract_empirical_bayes();
    my $paml_results = $bayes;
    
    my $results;
    my @slr_style_sites;
    foreach my $site (1 .. $input_cdna->length/3) {
      my $omega = $paml_results->{$site};
      my $se = $se->{$site};
      my $prob = $prob->{$site};

      next if (!defined $omega);
      my $lower = $omega - $se;
      $lower = 0 if ($lower < 0);
      my $upper = $omega + $se;

      my $result = '';
      $result = 'positive1' if ($prob > 0.5);
      $result = 'positive4' if ($pos->{$site});
      
      # Just for reference: slr style output.
      # Site  Neutral  Optimal   Omega    lower    upper LRT_Stat    Pval     Adj.Pval    Q-value Result Note
      # 1     4.77     3.44   0.0000   0.0000   1.4655   2.6626 1.0273e-01 8.6803e-01 1.7835e-02  --     Constant;
      # 0     1        2      3        4        5        6      7          8          9           10     11

      my @slr_style = ('') x 11;
      $slr_style[0] = $site;
      $slr_style[3] = $omega;
      $slr_style[4] = $lower;
      $slr_style[5] = $upper;
      $slr_style[6] = $prob || 0;
      $slr_style[10] = $result;
      
      push @slr_style_sites, \@slr_style;
    }
    $results->{'sites'} = \@slr_style_sites;

    $self->store_sitewise($results,$tree,$params);
  }
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
    $self->DESTROY;
    throw("Failing job! Reason given: $reason");
  }
}

sub store_sitewise {
  my $self = shift;
  my $results = shift;
  my $tree = shift;
  my $params = shift;

  my $node_id = $self->param('data_id');

  my $table = 'sitewise_aln';
  $table = $params->{'omega_table'} if ($params->{'omega_table'});
  
  my $parameter_set_id = 0;
  $parameter_set_id = $params->{'parameter_set_id'} if (defined($params->{'parameter_set_id'}));

  my $aln_map_aa = $self->param('aln_map_aa');
  my $input_aa = $self->param('input_aa');
  
  foreach my $site (@{$results->{'sites'}}) {
    #sleep(0.08);

    # Site  Neutral  Optimal   Omega    lower    upper LRT_Stat    Pval     Adj.Pval    Q-value Result Note
    # 1     4.77     3.44   0.0000   0.0000   1.4655   2.6626 1.0273e-01 8.6803e-01 1.7835e-02  --     Constant;
    # 0     1        2      3        4        5        6      7          8          9           10     11
    my ($site, $neutral, $optimal, $omega, $lower, $upper, $lrt_stat, $pval, $adj_pval, $q_value, $type, $note) = @{$site};

    my $mapped_site = $site;
    my $original_site = $site;
    if (defined $aln_map_aa) {
      $mapped_site = $aln_map_aa->{$site};
      if (defined $mapped_site) {
	$original_site = $site;
	$site = $mapped_site;
      } else {
	die "Something went wrong! NO mapped site for $site!\n";
      }
    }

    # These two cases are pretty useless, so we'll skip 'em.
    if (!$self->param('sitewise_store_gaps')) {
      next if ($note eq 'all_gaps');
      next if ($note eq 'single_char');
    }

    my $nongaps = Bio::EnsEMBL::Compara::AlignUtils->get_nongaps_at_column($input_aa,$original_site);
    next if ($nongaps == 0);

    printf("Site: %s  nongaps: %d  omegas: %3f (%3f - %3f) type: %s note: %s \n",$site,$nongaps,$omega,$lower,$upper,$type,$note);

    my $sth = $tree->adaptor->prepare
      ("REPLACE INTO $table
                           (parameter_set_id,
                            aln_position,
                            node_id,
                            omega,
                            omega_lower,
                            omega_upper,
                            lrt_stat,
                            ncod,
                            type,
                            note) VALUES (?,?,?,?,?,?,?,?,?,?)");
    $sth->execute($parameter_set_id,
		  $site,
		  $node_id,
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
    $self->SUPER::DESTROY if $self->can("SUPER::DESTROY");
}


1;