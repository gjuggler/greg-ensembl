package Bio::EnsEMBL::Compara::RunnableDB::Sitewise_dNdS;

use strict;
#use Time::HiRes qw(time gettimeofday tv_interval sleep);
use Cwd;
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::NestedSet;
use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Hive::Process;

use Bio::Greg::Codeml;
use Bio::Greg::ProcessUtils;

our @ISA = qw(Bio::EnsEMBL::Hive::Process Bio::Greg::ProcessUtils);

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
  if ($self->{dba}) {
    $dba = $self->dba;
  } else {
    $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-DBCONN=>$self->db->dbc);
  }
  $pta = $dba->get_ProteinTreeAdaptor;
  
  ### DEFAULT PARAMETERS ###
  $params = {
    alignment_table            => 'protein_tree_member',
    omega_table           => 'sitewise_omega',
    parameter_set_id       => 1,

    sitewise_store_gaps             => 1,
    sitewise_parameter_sets         => 'all',
    sitewise_action                 => 'slr',                # Which action(s) to perform. Space-delimited.
                                                    # 'slr' - SLR sitewise omegas.
                                                    # 'paml_sitewise' - PAML sitewise omegas.
                                                    # 'wobble' - SLR_wobble test.
                                                    # 'slr_reoptimise' - reoptimise with SLR
                                                    # 'paml_reoptimise' - reoptimise with PAML
                                                    # 'paml_lrt' - likelihood ratio test with PAML

    # SLR Parameters
    slr_gencode                => 'universal',
    slr_aminof                 => 1,                    # Options: 0, 1, 2 (default 1)
    slr_codonf                 => 0,                    # Options: 0, 1, 2 (default 0)
    slr_freqtype               => 0,                    # Options: 0, 1, 2 (default 0)

    # PAML Parameters
    #paml_model                  => 'M3',                 # Used for Bayes Empirical Bayes sitewise analysis.
    #paml_model_b                => 'M7',                 # Used for the likelihood ratio tests.

    alignment_score_mask_character => 'N',
    sequence_quality_mask_character_character => 'N'

    };
  
  # For aminof, codonf, and freqtype, see the SLR readme.txt for more info.

  #########################
  print "TEMP: ".$self->worker_temp_directory."\n";

  my $p_params = $self->get_params($self->parameters);
  my $i_params = $self->get_params($self->input_id);
  my $node_id = $i_params->{'protein_tree_id'} || $i_params->{'node_id'};
  my $t_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_tree_tags($dba,$node_id);
  
  $params = $self->replace_params($params,$p_params,$i_params,$t_params);
  #########################

  $dba->dbc->disconnect_when_inactive(1);
}

sub run {
  my $self = shift;  

  my $node_id = $params->{'node_id'};
  my $param_set_string = get('parameter_sets');
  print "Param sets: $param_set_string\n";

  my @param_sets;
  if ($param_set_string eq 'all') {
    my $query = qq^select distinct(parameter_set_id) FROM parameter_set order by parameter_set_id;^;
    @param_sets = @{$dba->dbc->db_handle->selectcol_arrayref($query)};
  } else {
    @param_sets = split(",",$param_set_string);
  }

  delete $params->{'parameter_sets'};
  foreach my $param_set (@param_sets) {
    my $param_set_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_param_set($dba->dbc,$param_set);
    my $new_params = Bio::EnsEMBL::Compara::ComparaUtils->replace_params($params,$param_set_params);

    Bio::EnsEMBL::Compara::ComparaUtils->hash_print($new_params);

    # Load the tree. This logic has been delegated to ComparaUtils.
    $tree = Bio::EnsEMBL::Compara::ComparaUtils->get_tree_for_comparative_analysis($dba,$new_params);
    $tree = $tree->minimize_tree;
    
    # Think of reasons why we want to fail the job.
    my @leaves = $tree->leaves;
    if (scalar(@leaves) <= 3) {
      next;
    } elsif (scalar(@leaves) > 300) {
      next;
    }
    foreach my $leaf (@leaves) {
      if ($leaf->distance_to_parent > 100) {
	$leaf->distance_to_parent(4);
	#$self->fail_job("Tree contains a massive branch length -- probably an error in b.l. optimization!");
      }
    }
    
    eval {
    $self->run_with_params($new_params,$tree);
      $tree->release_tree;
    };
  }
}

sub run_with_params {
  my $self = shift;
  my $params = shift;
  my $tree = shift;

  #$self->check_job_fail_options;

  print "Getting alignments...\n";
  my $aa_aln = $tree->get_SimpleAlign();
  my $cdna_aln = $tree->get_SimpleAlign(-cdna => 1);

  foreach my $leaf ($tree->leaves) {
    #print $leaf->member_id."\n";
  }

  $input_aa = Bio::EnsEMBL::Compara::ComparaUtils->fetch_masked_alignment($aa_aln,$cdna_aln,$tree,$params,0);
  $input_cdna = Bio::EnsEMBL::Compara::ComparaUtils->fetch_masked_alignment($aa_aln,$cdna_aln,$tree,$params,1);

  Bio::EnsEMBL::Compara::AlignUtils->pretty_print($input_aa,{length => 100});
  Bio::EnsEMBL::Compara::AlignUtils->pretty_print($input_cdna,{length => 100});

  if (get('action') =~ m/slr/i) {
    $self->run_sitewise_dNdS($tree,$input_cdna,$params);
  } elsif (get('action') =~ m/paml/i) {
    $self->run_paml($tree,$input_cdna,$params);
  } elsif (get('action') =~ m/wobble/i) {
    $params->{'slr_wobble'} = 0;
    my $results_nowobble = $self->run_wobble($tree,$input_cdna,$params);
    $tree->store_tag("lnl_nowobble",$results_nowobble->{'lnL'});

    sleep(2);

    $params->{'slr_wobble'} = 1;
    my $results_wobble = $self->run_wobble($tree,$input_cdna,$params);    
    $tree->store_tag("lnl_wobble",$results_wobble->{'lnL'});
  }

}

sub get {
  my $key = shift;
  return $params->{'sitewise_'.$key};
}

sub get_slr {
  my $key = shift;
  return $params->{'slr_'.$key};
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

  my $rc = 1;
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

    $exit_status = close($run);
    
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
  my $slrexe = get_slr('executable');
  my $gencode = get_slr('gencode');
  my $aminof = get_slr('aminof');
  my $codonf = get_slr('codonf');
  my $freqtype = get_slr('freqtype');

#  $slrexe = "/nfs/users/nfs_g/gj1/bin/Slr_ensembl" if (! -e $slrexe);
#  $slrexe = "/nfs/users/nfs_g/gj1/bin/Slr" if (! -e $slrexe);
#  $slrexe = "/nfs/users/nfs_g/gj1/bin/Slr2" if (! -e $slrexe);
#  $slrexe = "/nfs/users/nfs_g/gj1/bin/Slr_Linux_static" if (! -e $slrexe);
  $slrexe = "/nfs/users/nfs_g/gj1/build_scratch/slr/bin/Slr";
  $slrexe = "/homes/greg/bin/Slr" if (! -e $slrexe);

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
  print SLR "gencode\: $gencode\n";
  print SLR "aminof\: $aminof\n";
  print SLR "codonf\: $codonf\n";
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

    print $error_string."\n";
    
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

    if ( (grep { /\berr(or)?: /io } @output)  || !$exit_status) {
      print "@output\n";
      throw("There was an error running SLR!");
      $rc = 0;
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
    if (get('action') =~ m/reoptimise/i) {
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

  # Store the results in the database.
  my $kappa_key = "slr_kappa_".$params->{parameter_set_id};
  my $omega_key = "slr_omega_".$params->{parameter_set_id};
  my $lnl_key = "slr_lnL_".$params->{parameter_set_id};
    
  $tree->store_tag($kappa_key,$results->{'kappa'});
  $tree->store_tag($omega_key,$results->{'omega'});
  $tree->store_tag($lnl_key,$results->{'lnL'});
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

  if (get('action') =~ m/reoptimize/i) {
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

  if (get('action') =~ m/lrt/i) {

    ### Perform a likelihood ratio test between two models.
    my $model_b = $params->{'paml_model'} || $params->{'paml_model_b'};
    my $model_a = $params->{'paml_model_a'};
    
    $model_a =~ s/m//i;
    $model_b =~ s/m//i;

    my ($twice_lnL,$codeml_a,$codeml_b) = Bio::Greg::Codeml->NSsites_ratio_test($treeI,$input_cdna,$model_a,$model_b,$self->worker_temp_directory);
    my $test_label = sprintf("PAML LRT M%s-M%s",$model_a,$model_b);
    my $tags;
    $tags->{$test_label} = $twice_lnL;
    Bio::EnsEMBL::Compara::ComparaUtils->store_tags($tree,$tags);
  }

  if (get('action') =~ m/sitewise/i) {
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
    $dont_write_output = 1;
    $self->DESTROY;
    throw("Failing job! Reason given: $reason");
  }
}

sub store_sitewise {
  my $self = shift;
  my $results = shift;
  my $tree = shift;
  my $params = shift;

  my $tree_node_id = 0;
  $tree_node_id = $tree->subroot->node_id if (defined $tree->subroot);
  my $node_id = $tree->node_id;

  my $table = 'sitewise_aln';
  $table = $params->{'omega_table'} if ($params->{'omega_table'});
  
  my $parameter_set_id = 0;
  $parameter_set_id = $params->{'parameter_set_id'} if (defined($params->{'parameter_set_id'}));
  
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
	$original_site = $site;
	$site = $mapped_site;
      }
    }

    # These two cases are pretty useless, so we'll skip 'em.
    if (!get('store_gaps')) {
      next if ($note eq 'all_gaps');
      next if ($note eq 'single_char');
    }

    my $nongaps = Bio::EnsEMBL::Compara::AlignUtils->get_nongaps_at_column($input_aa,$original_site);
    next if ($nongaps == 0);

    #printf("Site: %s  nongaps: %d  omegas: %3f (%3f - %3f) type: %s note: %s \n",$site,$nongaps,$omega,$lower,$upper,$type,$note);

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
    $tree->release_tree if ($tree);
    $tree = undef;
    $self->SUPER::DESTROY if $self->can("SUPER::DESTROY");
}


1;
