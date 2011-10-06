package Bio::Greg::Hive::PhyloAnalysis;

use strict;
use Bio::EnsEMBL::Utils::Argument;
use Bio::EnsEMBL::Utils::Exception;

use Time::HiRes qw(sleep);
use Cwd;
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::NestedSet;

use File::Path;

use Bio::Greg::Codeml;

use base ('Bio::Greg::Hive::Process');

#
# Some global-ish variables.
#

sub debug { 1; }

sub param_defaults {
  my $self = shift;
  my $params = {
    alignment_table  => 'protein_tree_member',
    omega_table      => 'sitewise_omega',
    parameter_set_id => 0,

    sitewise_filter_order => 'before',

    genewise_table_output => 1,
    sitewise_table_output => 1,
    genewise_table => 'dnds_genes',

    sitewise_store_opt_tree     => 1,
    sitewise_skip_gaps         => 0,
    sitewise_skip_single_chars         => 0,
    sitewise_minimum_leaf_count => 2,
    sitewise_strip_gaps         => 0,
    sitewise_parameter_sets     => 'all',
    analysis_action             => 'hyphy_dnds',    # Which action(s) to perform. Space-delimited.
                                                    # 'slr' - SLR sitewise omegas.
                                                    # 'paml_sitewise' - PAML sitewise omegas.
                                                    # 'wobble' - SLR_wobble test.
                                                    # 'slr_reoptimise' - reoptimise with SLR
                                                    # 'paml_reoptimise' - reoptimise with PAML
                                                    # 'paml_lrt' - likelihood ratio test with PAML
         # 'hyphy_dnds' - Calculate dn/ds with HYPHY, along with confidence intervals.

    # SLR Parameters
    slr_gencode      => 'universal',
    slr_aminof       => 1,             # Options: 0, 1, 2 (default 1)
    slr_codonf       => 0,             # Options: 0, 1, 2 (default 0)
    slr_freqtype     => 0,             # Options: 0, 1, 2 (default 0)
    slr_skipsitewise => 0,             # Options: 0, 1 (default 0)
    slr_malloc_check => 0,
    
    # PAML Parameters
    paml_model                  => 'M8',                 # Used for Bayes Empirical Bayes sitewise analysis.


    paml_model_b                => 'M7',                 # Used for the likelihood ratio tests.
    paml_model_a                => 'M8',
    paml_branch_model => 0,
    paml_fix_kappa => 0,
    paml_kappa => 2,
    paml_fix_omega => 0,
    paml_omega => 0.8,
    paml_fix_alpha => 1,
    paml_alpha => 0,

    # HYPHY Parameters
    hyphy_codon_model => 'mg94',
  };
  return $params;
}

sub fetch_input {
  my ($self) = @_;

  $self->load_all_params;
}

use base ('Bio::Greg::Hive::CollectSitewiseStats', 'Bio::Greg::StatsCollectionUtils');

sub get_genewise_structure {
  my $self = shift;

  my $structure = {
    data_id   => 'int',
    node_id   => 'int',
    parameter_set_id => 'int',

    method => 'char32',
    dnds => 'float',
    lnL => 'float',
    kappa => 'float',
    tree => 'string',

    n_leaves => 'int',
    aln_length => 'int',
    aln_site_count => 'int',

    unique_keys => 'data_id,node_id,parameter_set_id'
  };

  return $structure;
}


sub run {
  my $self = shift;

  my $tree = $self->get_tree;

  print $tree->newick_format . "\n";

  # Think of reasons why we want to fail the job.
  my @leaves = $tree->leaves;
  if ( scalar(@leaves) < $self->param('sitewise_minimum_leaf_count') ) {
    my $value = sprintf "Too small (%d < %d)", scalar(@leaves),
      $self->param('sitewise_minimum_leaf_count');
#    $self->store_tag( "slr_skipped", $value );
    return;
  } elsif ( scalar(@leaves) > 300 ) {
    my $value = sprintf "Too big (%d > %d)", scalar(@leaves), 300;
#    $self->store_tag( "slr_skipped", $value );
    return;
  } elsif ( $self->check_tree_aln < 0 ) {
#    $self->store_tag( "slr_skipped", "Tree or align doesn't look good!" );
    print "Skipping: tree or align doesn't look good!\n";
    return;
  }
  foreach my $leaf (@leaves) {
    if ( $leaf->distance_to_parent > 100 ) {
      $leaf->distance_to_parent(4);
    }
  }

  #  eval {
  $self->run_with_params( $self->params, $tree );

  #  };
  if ($@) {
    print "ERROR - Trying with removed gaps...\n";
    $self->param( 'sitewise_strip_gaps', 1 );

    eval { $self->run_with_params( $self->param, $tree ); };
    if ($@) {
      print "ERROR - Trying with malloc_check...\n";

      $self->param( 'sitewise_malloc_check', 1 );

      eval { $self->run_with_params( $self->params, $tree ); };
      if ($@) {
        print "ERROR - GIVING up! \n";
#        $self->store_tag( 'slr_error', 1 );
        $self->fail_and_die;
        sleep(2);
      }
    }
  }
  $tree->release_tree;

}

sub run_with_params {
  my $self   = shift;
  my $params = shift;
  my $tree   = shift;

  print "Getting alignments...\n";

  foreach my $leaf ( $tree->leaves ) {

    #print $leaf->member_id."\n";
  }

  my $input_aa   = $self->get_aln;
  my $input_cdna = $self->get_cdna_aln;


  if ( $self->param('sitewise_filter_order') eq 'after' ) {
    $self->param( 'aa_filtered',   $input_aa );
    $self->param( 'cdna_filtered', $input_cdna );

    # Disable alignment score filtering and re-fetch the input alignment.
    $self->param( 'alignment_score_filtering', 0 );
    $input_aa   = $self->get_aln;
    $input_cdna = $self->get_cdna_aln;
  }

#  my ($slim_cdna,$cdna_new_to_old,$cdna_old_to_new) = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns_in_threes($input_cdna);
#  my ($slim_aa,$aa_new_to_old,$aa_old_to_new) = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns($input_aa);
#  my $input_aa = $slim_aa;
#  my $input_cdna = $slim_cdna;
#  $self->param('input_cdna',$slim_cdna);
#  $self->param('input_aa',$slim_aa);
#  $self->param('aln_map_aa',$aa_new_to_old);
#  $self->param('aln_map_cdna',$cdna_new_to_old);

  $self->param( 'input_aa',   $input_aa );
  $self->param( 'input_cdna', $input_aa );

  Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $input_aa,   { length => 100,full=>1 } );
  Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $input_cdna, { length => 100,full=>1 } );
  my $action = $self->param('analysis_action');

  $self->compara_dba->dbc->disconnect_when_inactive(1);

  if ( $action =~ m/slr/i ) {
    $self->run_sitewise_dNdS( $tree, $input_cdna, $self->params );
    sleep(2);
  } elsif ( $action =~ m/paml/i ) {
    $self->run_paml( $tree, $input_cdna, $self->params );
  } elsif ( $action =~ m/wobble/i ) {
    $self->param( 'slr_wobble', 0 );
    my $results_nowobble = $self->run_wobble( $tree, $input_cdna, $self->params );
#    $self->store_tag( "lnl_nowobble", $results_nowobble->{'lnL'} );
    sleep(2);

    $self->param( 'slr_wobble', 1 );
    my $results_wobble = $self->run_wobble( $tree, $input_cdna, $self->params );
#    $self->store_tag( "lnl_wobble", $results_wobble->{'lnL'} );
  } elsif ( $action =~ m/hyphy_dnds/i ) {
    $self->run_hyphy( $tree, $input_cdna, $self->params );
  } elsif ( $action =~ m/xrate_indels/i ) {
    $self->run_xrate_indels( $tree, $input_cdna, $self->params );
  } elsif ( $action =~ m/indelign/i ) {
    $self->run_indelign( $tree, $input_cdna, $self->params );
  }

  $self->compara_dba->dbc->disconnect_when_inactive(0);
}

sub run_hyphy {
  my $self     = shift;
  my $tree     = shift;
  my $cdna_aln = shift;
  my $params   = shift;

  my $cwd    = cwd();
  my $tmpdir = $self->worker_temp_directory;
  chdir($tmpdir);

  # OUTPUT THE ALIGNMENT.
  my $aln_f = $tmpdir . "aln.fa";
  open( OUT, ">$aln_f" );
  foreach my $seq ( $cdna_aln->each_seq ) {
    my $name = $seq->id;
    print OUT ">$name\n";
    print OUT $seq->seq . "\n";
  }
  close(OUT);

  # OUTPUT THE TREE.
  my $tree_f      = $tmpdir . "tree.nh";
  my $tree_newick = Bio::EnsEMBL::Compara::TreeUtils->to_newick($tree);
  open( OUT, ">$tree_f" );
  print OUT $tree_newick . "\n";
  close(OUT);

  my $control_f = $tmpdir . "control.txt";
  open( OUT, ">$control_f" );
  print OUT $aln_f . "\n";
  print OUT $tree_f . "\n";
  close(OUT);

  my $hyphy_base = "/nfs/users/nfs_g/gj1/src/greg-ensembl/ensembl-greg/scripts/hyphy";

  my $batch_file = "${hyphy_base}/fit_codon_model.bf";
  open( IN, $batch_file );
  my @lines = <IN>;
  close(IN);
  my $batch_string = join( "", @lines );

  sub replace_string {
    my $string = shift;
    my $search = shift;
    my $file   = shift;

    open( IN, $file );
    my @lines = <IN>;
    close(IN);
    my $file_str = join( "", @lines );

    my $i = index( $string, $search );
    substr( $string, $i, length($search), $file_str );
    return $string;
  }

  my $codon_model = $self->param('hyphy_codon_model');
  $codon_model = 'mg94' unless (defined $codon_model);
  my $model_bf    = "${hyphy_base}/${codon_model}.bf";
  $batch_string = replace_string( $batch_string, '[MODEL_BF]', $model_bf );

  # Write the new batch file.
  my $batch_f = $tmpdir . "batch.bf";
  open( OUT, ">$batch_f" );
  print OUT $batch_string . "\n";
  close(OUT);

  # Copy all the other HYPHY crap to the temp directory.
  system("cp ${hyphy_base}/*.* $tmpdir");

  my $cmd = "HYPHY $batch_f < $control_f";
  print $cmd. "\n";
  my @output = ();
  my @output = `$cmd`;
  print "@output\n";

  my $omega_lo = 0;
  my $omega_hi = 0;
  my $omega    = 0;
  foreach my $line (@output) {
    $omega_lo = $1 if ( $line =~ m/omega_lo:(.*)/ );
    $omega_hi = $1 if ( $line =~ m/omega_hi:(.*)/ );
    $omega    = $1 if ( $line =~ m/omega:(.*)/ );
  }

  print "$omega_lo $omega $omega_hi\n";

#  $self->store_tag( "hyphy_dnds",    $omega );
#  $self->store_tag( "hyphy_dnds_lo", $omega_lo );
#  $self->store_tag( "hyphy_dnds_hi", $omega_hi );

  $self->param('hyphy_dnds',$omega);
  $self->param('hyphy_dnds_lo',$omega_lo);
  $self->param('hyphy_dnds_hi',$omega_hi);

  chdir($cwd);

}

sub run_indelign {
  my $self   = shift;
  my $tree   = shift;
  my $aln    = shift;
  my $params = shift;

  my ( $ins, $del, $ins_rate, $del_rate ) =
    Bio::EnsEMBL::Compara::AlignUtils->indelign( $aln, $tree, $self->params,
    $self->worker_temp_directory );

  print "ins: $ins_rate del: $del_rate\n";

#  $self->store_tag( "indelign_ins", $ins_rate );
#  $self->store_tag( "indelign_del", $del_rate );
}

sub run_xrate_indels {
  my $self   = shift;
  my $tree   = shift;
  my $aln    = shift;
  my $params = shift;

  my $cwd    = cwd();
  my $tmpdir = $self->worker_temp_directory;

  #my $tmpdir = "/tmp/xrate/";
  chdir($tmpdir);

  # OUTPUT THE ALIGNMENT.
  my $aln_f   = $tmpdir . "aln";
  my $tree_f  = $tmpdir . "tree";
  my $index_f = $tmpdir . "index.txt";

  use Bio::Annotation::AnnotationFactory;
  my $factory = Bio::Annotation::AnnotationFactory->new( -type => 'Bio::Annotation::SimpleValue' );

  $tree->distance_to_parent(0);

  #$tree = Bio::EnsEMBL::Compara::TreeUtils->scale($tree,1);
  my $tree_newick = Bio::EnsEMBL::Compara::TreeUtils->to_newick($tree);
  my $rfann       = $factory->create_object(
    -value   => $tree_newick,
    -tagname => 'NH'
  );
  my $coll = $aln->annotation;
  map { $_->description('') } $aln->each_seq;
  $coll->add_Annotation( 'custom', $rfann );
  my $alnout = Bio::AlignIO->new(
    '-format' => 'stockholm',
    '-file'   => ">$aln_f"
  );
  $alnout->write_aln($aln);
  $alnout->close;

  # OUTPUT THE TREE.
  open( OUT, ">$tree_f" );
  my $string = Bio::EnsEMBL::Compara::TreeUtils->to_newick($tree);
  $string =~ s/\:[0-9\.]+;/;/g;
  print "NEWICK : $string\n";
  print OUT $string . "\n";
  close(OUT);

  open( INDEX, ">$index_f" );
  print INDEX "test1 $tree_f $aln_f\n";
  close(INDEX);

  my $exec_folder = '/nfs/users/nfs_g/gj1/src/dart/bin';
  my $cmd         = "${exec_folder}/tkfidem $index_f -log 9";
  print "$cmd\n";
  my @results = `$cmd`;

  my $lambda;
  my $mu;
  foreach my $line (@results) {
    $lambda = $1 if ( $line =~ m/lambda:\s*([0-9\.]+)/ );
    $mu     = $1 if ( $line =~ m/\s*mu:\s*([0-9\.]+)/ );
  }

  print "lambda: $lambda mu: $mu\n";

#  $self->store_tag( "xrate_ins", $lambda );
#  $self->store_tag( "xrate_del", $mu );

  #print "@results\n";
}

sub run_wobble {
  my $self     = shift;
  my $tree     = shift;
  my $cdna_aln = shift;
  my $params   = shift;

  if ( scalar $tree->leaves > 30 ) {
    $self->fail_job("Tree too large for wobbly analysis!!");
  }

  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);

  # LOAD VARIABLES FROM PARAMS.
  my $slrexe   = get_slr('executable');
  my $gencode  = get_slr('gencode');
  my $aminof = 1;
  $aminof   = get_slr('aminof') if (defined get_slr('aminof'));
  my $codonf   = get_slr('codonf');
  my $freqtype = get_slr('freqtype');
  my $wobble   = get_slr('wobble');

  $slrexe = "/nfs/users/nfs_g/gj1/bin/Slr_wobble";

  # Reorder the alignment according to the tree
  $cdna_aln = Bio::EnsEMBL::Compara::AlignUtils->sort_by_tree( $cdna_aln, $treeI );

  my $num_leaves = scalar( @{ $tree->get_all_leaves } );
  my $tmpdir     = $self->worker_temp_directory;

  # CLEAN UP OLD RESULTS FILES.
  #unlink "$tmpdir/slr.res";
  #unlink "$tmpdir/tree";
  #unlink "$tmpdir/aln";
  #unlink "$tmpdir/slr.ctl";

  my $tree_newick = Bio::EnsEMBL::Compara::TreeUtils->to_newick($tree);

  # Map the stable_ids into shorter integers, to stay below Slr_wobble's max label length.
  my $tree_map;
  my $i      = 0;
  my @leaves = $tree->leaves;
  foreach my $seq ( $cdna_aln->each_seq ) {
    $i++;
    my $id = $seq->id;
    map { $tree_newick =~ s/$id/$i/g if ( $_->stable_id eq $id ) } @leaves;
  }

  Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $cdna_aln, { length => 50 } );

  # OUTPUT THE ALIGNMENT.
  my $alnout = Bio::AlignIO->new(
    '-format'      => 'phylip',
    '-file'        => ">$tmpdir/aln",
    '-interleaved' => 0,
    '-idlinebreak' => 1,
    '-idlength'    => $cdna_aln->maxdisplayname_length + 1
  );
  $alnout->write_aln($cdna_aln);
  $alnout->close();

  # OUTPUT THE TREE.
  open( OUT, ">$tmpdir/tree" );
  print OUT sprintf( "%d 1\n", $num_leaves );
  print OUT $tree_newick . "\n";
  close(OUT);

  # OUTPUT THE CTL FILE.
  my $slr_ctl = "$tmpdir/slr.ctl";
  open( SLR, ">$slr_ctl" ) or $self->throw("cannot open $slr_ctl for writing");
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

    my $prefix = "";

    #    my $prefix= "export MALLOC_CHECK_=1;";
    print "Running: $slrexe\n";
    my $run;
    open( $run, "$prefix $slrexe |" ) or $self->throw("Cannot open exe $slrexe");
    my @output;
    while (<$run>) {
      next if ( $_ =~ 'Unrecognised' );
      next if ( $_ =~ 'Odd gapping' );
      print $_;
      push @output, $_;
    }

    foreach my $outline (@output) {
      if ( $outline =~ /lnL\s+(\S+)/ ) {
        $results->{'lnL'} = $1;
      }
      if ( $outline =~ /^omega<-c\((.*)\)/ ) {
        my @omega_cats = split( ',', $1 );
        $results->{'omega_cats'} = \@omega_cats;
      }
      if ( $outline =~ /^wobble<-c\((.*)\)/ ) {
        my @wobble_cats = split( ',', $1 );
        $results->{'wobble_cats'} = \@wobble_cats;
      }
    }

    print " -> LnL: " . $results->{'lnL'} . "\n";
    print " -> OMEGA: " . @{ $results->{'omega_cats'} } . "\n";
    print " -> WOBBLE: " . @{ $results->{'wobble_cats'} } . "\n";
  }

  open RESULTS, "$tmpdir/$outfile" or die "couldnt open results file: $!";
  my @sites;
  while (<RESULTS>) {
    print $_;

    if (/(\S+)\s+(\S+)/) {
      my $site = $1;
      my $val  = $2;
      my @arr  = 0 x 12;
      $arr[3] = $val;
      push @sites, \@arr;
    }
  }
  $results->{'sites'} = \@sites;

  chdir $cwd;
  return $results;
}

sub run_sitewise_analysis {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pep_aln = shift;

  my $result_arrayref;

  my $params = $self->replace($self->param_defaults, $self->params);
  #$self->set_params($params);

  my $action = $self->param('analysis_action');
  if ($action eq 'slr') {
    $result_arrayref = $self->run_sitewise_dNdS($tree, $aln, $params);
  } elsif ($action eq 'paml' || $action eq 'paml_m2' || $action eq 'paml_m8') {
    $result_arrayref = $self->run_paml_sitewise($tree, $aln, $params);
  }
  return $result_arrayref;
}

sub run_sitewise_dNdS {
  my $self     = shift;
  my $tree     = shift;
  my $cdna_aln = shift;
  my $params   = shift;

  $params = $self->replace($self->params, $params);

  print "  running SLR...\n" if ($self->debug);

  if ($tree->isa("Bio::Tree::TreeI")) {
    $tree = Bio::EnsEMBL::Compara::TreeUtils->from_treeI($tree);
  }
  $tree = Bio::EnsEMBL::Compara::TreeUtils->unroot($tree);
  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);

  # LOAD VARIABLES FROM PARAMS.
  my $slrexe       = $params->{slr_executable};
  my $gencode      = $params->{slr_gencode};
  my $aminof       = $params->{slr_aminof};
  my $codonf       = $params->{slr_codonf};
  my $freqtype     = $params->{slr_freqtype};
  my $skipsitewise = $params->{slr_skipsitewise};
  my $icode = $params->{icode};

  my $gencode = 'universal';
  $gencode = 'mammalian' if ($icode == 1);

  $slrexe = "Slr" if (!-e $slrexe);
  $slrexe = "/nfs/users/nfs_g/gj1/bin/Slr_Linux_shared" if ( !-e $slrexe );
  $slrexe = "/homes/greg/bin/Slr"              if ( !-e $slrexe );
  $slrexe = "Slr"                              if ( !-e $slrexe );

  # Reorder the alignment according to the tree
  #$cdna_aln = Bio::EnsEMBL::Compara::AlignUtils->sort_by_tree($cdna_aln,$treeI);

  my $num_leaves = scalar( @{ $tree->get_all_leaves } );
  my $tmpdir     = $self->worker_temp_directory;

  # CLEAN UP OLD RESULTS FILES.
  unlink "$tmpdir/slr.res" if ( -e "$tmpdir/slr.res" );
  unlink "$tmpdir/tree"    if ( -e "$tmpdir/tree" );
  unlink "$tmpdir/aln"     if ( -e "$tmpdir/aln" );
  unlink "$tmpdir/slr.ctl" if ( -e "$tmpdir/slr.ctl" );

  # OUTPUT THE ALIGNMENT.
  my $alnout = Bio::AlignIO->new(
    '-format'      => 'phylip',
    '-file'        => ">$tmpdir/aln",
    '-interleaved' => 0,
    '-idlinebreak' => 1,
    '-idlength'    => $cdna_aln->maxdisplayname_length + 1
  );
  $alnout->write_aln($cdna_aln);
  $alnout->close();

  # OUTPUT THE TREE.
  my $tree_copy = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree);
  map {$_->name('') if (!$_->is_leaf) } $tree_copy->nodes;

  my $tree_newick = Bio::EnsEMBL::Compara::TreeUtils->to_newick($tree_copy);
  #print "Newick: [$tree_newick]\n";
  open( OUT, ">$tmpdir/tree" );
  print OUT sprintf( "%d 1\n", $num_leaves );
  print OUT $tree_newick . "\n";
  close(OUT);

  # OUTPUT THE CTL FILE.
  my $slr_ctl = "$tmpdir/slr.ctl";
  open( SLR, ">$slr_ctl" ) or $self->throw("cannot open $slr_ctl for writing");
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

  #my $asdf = `cat $tmpdir/slr.ctl`;
  #print $asdf."\n";
  #$asdf = `ls $tmpdir`;
  #print $asdf."\n";

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
  if ( $params->{sitewise_malloc_check} ) {
    $prefix = "export MALLOC_CHECK_=1;";
  } else {
    $prefix = "unset MALLOC_CHECK_;";
  }

  #print "Running SLR!\n";
  #print "$prefix $slrexe\n";
  open( $run, "$prefix $slrexe |" ) or $self->throw("Cannot open exe $slrexe");
  my @output = <$run>;

  my $exited_well = close($run);
  if ( !$exited_well ) {
    throw("Slr didn't exit well!");
  }

  if ( ( grep { /\berr(or)?: /io } @output ) ) {
    throw("There was an error running SLR ('error' in output)!");
  }

  # Concatenate the SLR output to the end of our @output array.
  open(IN,"$outfile");
  while (<IN>) {
    push @output,$_;
  }
  close(IN);

  chdir($cwd);
  return \@output;
}


sub parse_sitewise_output_lines {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pep_aln = shift;
  my $output_arrayref = shift;

  # The results hash contains objects keyed by the amino acid site position
  my $results_hash;
  
  my $action = $self->param('analysis_action');
  if ($action eq 'slr') {
    $results_hash = $self->parse_slr_output($tree, $aln, $pep_aln, $output_arrayref);
  } elsif ($action eq 'paml' || $action eq 'paml_m2' || $action eq 'paml_m8') {
    $results_hash = $self->parse_paml_output($tree, $aln, $pep_aln, $output_arrayref);
  }

  return $results_hash;
}

sub parse_sitewise_file {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pep_aln = shift;
  my $output_file = shift;

  open(IN,"$output_file");
  my @output = <IN>;
  close(IN);

  return $self->parse_sitewise_output_lines($tree, $aln, $pep_aln, \@output);
}

sub parse_slr_output {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pep_aln = shift;
  my $output_lines_arrayref = shift;

  my $params = $self->params;
  my @output = @$output_lines_arrayref;

  my $results;

  foreach my $outline (@output) {
    if ( $outline =~ /lnL = (\S+)/ ) {
      $results->{"lnL"} = $1;
    }
    if ( $outline =~ /kappa = (\S+)/i ) {
      $results->{'kappa'} = $1;
    }
    if ( $outline =~ /omega = (\S+)/i ) {
      $results->{'omega'} = $1;
    }
  }

  # Collect SLR's reoptimized tree.
  my $next_line_is_it = 0;
  my $new_pt;
  foreach my $outline (@output) {
    if ($next_line_is_it) {
      my $new_tree = $outline;
      $new_pt = Bio::EnsEMBL::Compara::TreeUtils->from_newick($new_tree);
      #print "Reoptimised tree: " . $new_pt->newick_format . "\n";
      last;
    }
    # We know that the new tree will show up just after the LnL line.
    $next_line_is_it = 1 if ( $outline =~ /lnL =/ );
  }
  $results->{'slr_tree'} = $new_pt;
  $results->{'slr_tree_string'} = $new_pt->newick_format;

  my $skipsitewise = $params->{slr_skipsitewise};
  my $tmpdir     = $self->worker_temp_directory;

  my $okay = 0;
  my @sites;
  my $type = '';
  my $note = '';
  my $random = '';
  foreach my $line (@output) {
    $_ = $line;
    $type = '';
    $note = '';
    $random = '';
    chomp $_;
    
    #      print $_."\n";
    if (/^\#/) {
      next;
    }
    
    # Parse the notes.
    if (/All gaps/i) {
      $note = 'all_gaps';
    } elsif (/Constant/i) {
      $note = 'constant';
    } elsif (/Single char/i) {
      $note = 'single_char';
    } elsif (/Synonymous/i) {
      $note = 'synonymous';
    }
    
    if (/\!/)
    { # Make sure to parse the random note last, because it's sometimes shown alongside other notes.
      $random = 'random';
    }
    
    # Parse the pos/sel test results.
    if (/\+\+\+\+\s+/) {
      $type = 'positive4';
    } elsif (/\+\+\+\s+/) {
      $type = 'positive3';
    } elsif (/\+\+\s+/) {
      $type = 'positive2';
    } elsif (/\+\s+/) {
      $type = 'positive1';
    } elsif (/\-\-\-\-\s+/) {
      $type = 'negative4';
    } elsif (/\-\-\-\s+/) {
      $type = 'negative3';
    } elsif (/\-\-\s+/) {
      $type = 'negative2';
    } elsif (/\-\s+/) {
      $type = 'negative1';
    }
    if (/^\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
      my @arr = ( $1, $2, $3, $4, $5, $6, $7, $8, $9, 10, $type, $note, $random );
      
      my (
        $site,     $neutral, $optimal,  $omega,   $lower, $upper,
        $lrt_stat, $pval,    $adj_pval, $q_value, $type,  $note, $random
        ) = @arr;
      my $nongaps = Bio::EnsEMBL::Compara::AlignUtils->get_nongaps_at_column( $pep_aln, $site );
      
      my $signed_lrt = $lrt_stat;
      $signed_lrt = -$signed_lrt if ($omega < 1);

      my $obj = {
        aln_position => $site,
        ncod => $nongaps,
        nongap_count => $nongaps,
        omega => $omega,
        omega_lower => $lower,
        omega_upper => $upper,
        lrt_stat => $signed_lrt,
        signed_lrt => $signed_lrt,
        type => $type,
        note => $note,
        random => $random
      };
      $results->{$site} = $obj;

    } else {
      #warn("Doesn't look like a result line: [$_]\n");
    }
  }

  $self->param('slr_kappa',$results->{kappa});
  $self->param('slr_dnds',$results->{omega});
  $self->param('slr_lnL',$results->{lnL});

  return $results;
}

sub run_paml_sitewise {
  my $self = shift;
  my $tree = shift;
  my $input_cdna = shift;
  my $params = shift;

  print "  running PAML...\n" if ($self->debug);

  my $tmp = $self->worker_temp_directory;
  #mkpath([$tmp]);

  foreach my $file (qw(mlc rst rub lnf)) {
  #  if (-e "${tmp}/${file}") {
  #    unlink "${tmp}/${file}";
  #  }
  }

  # PAML has a problem with sequences that start with a number... we need to map
  # the IDs to avoid this.
  my $needs_translating = 0;
  foreach my $leaf ($tree->leaves) {
    my $id = $leaf->name;
    if ($id =~ m/^[0-9]+$/g) {
      $needs_translating = 1;
    }
  }
  if ($needs_translating) {
    die("Needs translating!");
    my $leaf_map;
    my $seq_map;
    map {$leaf_map->{$_->node_id} = 's_'.$_->name;} $tree->leaves;
    map {$seq_map->{$_->id} = 's_'.$_->id} $input_cdna->each_seq;
    $tree = Bio::EnsEMBL::Compara::TreeUtils->translate_ids($tree, $leaf_map);
    $input_cdna = Bio::EnsEMBL::Compara::AlignUtils->translate_ids($input_cdna, $seq_map);
  }

  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
  #$treeI->get_root_node->branch_length(0);

  #$self->pretty_print($input_cdna);
  #print $treeI->ascii;

  # Output the tree and alignment files
  my $aln_file = ">${tmp}/aln.phy";
  my $aln_out = Bio::AlignIO->new(
      '-format'      => 'phylip',
      '-file'          => $aln_file,
      '-interleaved' => 0,
      '-idlength'    => $input_cdna->maxdisplayname_length() + 1
  );
  $aln_out->write_aln($input_cdna);

  my $tree_file = "${tmp}/tree.nh";
  Bio::EnsEMBL::Compara::TreeUtils->to_file($tree, $tree_file);

  # Create a codeml.ctl file here.

  my $sites_model = $params->{paml_model};
  $sites_model =~ s/M//gi;
  my $branch_model = $params->{paml_branch_model};
  my $fix_kappa = $params->{paml_fix_kappa};
  my $kappa = $params->{paml_kappa};
  my $fix_omega = $params->{paml_fix_omega};
  my $omega = $params->{paml_omega};
  my $fix_alpha = $params->{paml_fix_alpha};
  my $alpha = $params->{paml_alpha};

  if ($params->{analysis_action} eq 'paml_m2') {
    $sites_model = '2';
  } elsif ($params->{analysis_action} eq 'paml_m8') {
    $sites_model = '8';
  }

  my $ctl_string = qq^
seqfile = aln.phy
treefile = tree.nh

model = ${branch_model}
NSsites = ${sites_model}
fix_kappa = ${fix_kappa}
kappa = ${kappa}
fix_omega = ${fix_omega}
omega = ${omega}

fix_alpha = 1
alpha = 0
outfile = mlc
noisy = 0
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
clock = 0
icode = 0
Malpha = 0
getSE = 0
RateAncestor = 1
fix_blength = 1
method = 0
^;  
  print $ctl_string."\n";
  my $cwd = cwd();
  chdir $tmp;
  open(OUT, ">codeml.ctl");
  print OUT $ctl_string."\n";
  close(OUT);

  system("codeml");

  open(IN, "mlc");
  my @main_lines = <IN>;
  close(IN);
  open(IN, "rst");
  my @supp_lines = <IN>;
  close(IN);

  my @results_lines;
  #push @results_lines, split("\n", $ctl_string);
  push @results_lines, @main_lines;
  push @results_lines, @supp_lines;

  #print "@results_lines\n";

  chdir $cwd;
  return \@results_lines;
}

sub parse_paml_output {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pep_aln = shift;
  my $output_lines_arrayref = shift;

  print "    parsing paml output\n";

  my $codeml = Bio::Greg::Codeml->new();
  my $codeml_results = $codeml->parse_results($aln, $pep_aln, $output_lines_arrayref);

  my $normalized_results;
  foreach my $site ( 1 .. $pep_aln->length ) {
    my $obj = $codeml_results->{$site};
    #next unless (defined $obj);

    my $omega = $obj->{omega};
    my $se = $obj->{se};
    my $prob = $obj->{prob_gt_one};
    my $pos = $obj->{is_positive_site};

    my $lower = $omega - $se;
    my $upper = $omega + $se;
      
    my $type = '';
    $type = 'positive1' if ( $omega > 1 && $prob > 0.5 );
    $type = 'positive2' if ( $omega > 1 && $prob > 0.9 );
    $type = 'positive4' if ( $pos );
    $type = 'negative1' if ( $omega < 1 && $prob < 0.5 );
    $type = 'negative2' if ( $omega < 1 && $prob < 0.1 );

    my $nongaps = Bio::EnsEMBL::Compara::AlignUtils->get_nongaps_at_column( $pep_aln, $site );

    my $obj = {
      aln_position => $site,
      ncod => $nongaps,
      nongap_count => $nongaps,
      omega => $omega,
      omega_lower => $lower,
      omega_upper => $upper,
      lrt_stat => $prob,
      type => $type
    };
    
    $normalized_results->{$site} = $obj;
  }

  return $normalized_results;
}

sub run_paml {
  my $self       = shift;
  my $tree       = shift;
  my $input_cdna = shift;
  my $params     = shift;

  print $tree->newick_format() . "\n";

  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
  $treeI->get_root_node->branch_length(0);

  if ( $self->param('analysis_action') =~ m/lrt/i ) {

    ### Perform a likelihood ratio test between two models.
    my $model_b = $params->{'paml_model'} || $params->{'paml_model_b'};
    my $model_a = $params->{'paml_model_a'};

    $model_a =~ s/m//i;
    $model_b =~ s/m//i;

    my ( $twice_lnL, $codeml_a, $codeml_b ) =
      Bio::Greg::Codeml->NSsites_ratio_test( $treeI, $input_cdna, $model_a, $model_b,
      $self->worker_temp_directory );
    my $test_label = sprintf( "PAML LRT M%s-M%s", $model_a, $model_b );
    $self->param($test_label,$twice_lnL);
#    $self->store_tag( $test_label, $twice_lnL );
  }

  if ( $self->param('analysis_action') =~ m/sitewise/i ) {
    ### Perform a BEB sitewise analysis of omegas.

    # Should we be scaling the tree? Dunno... GJ 2009-10-09
    #Bio::EnsEMBL::Compara::TreeUtils->scale($treeI,3);

    my $model = $params->{'paml_model'} || $params->{'paml_model_b'};
    $model =~ s/m//i;
    my $codeml_params = {
      NSsites => $model,
      ncatG   => 10
    };
    my $codeml = Bio::Greg::Codeml->new(
      -params    => $codeml_params,
      -alignment => $input_cdna,
      -tree      => $treeI,
      -tempdir   => $self->worker_temp_directory
    );
    $codeml->save_tempfiles(1);
    $codeml->run();

  }
}


sub store_sitewise {
  my $self    = shift;
  my $tree    = shift;
  my $pep_aln = shift;
  my $results = shift;
  my $params  = shift;

  my $node_id = $self->node_id;
  my $data_id = $self->data_id;

  my $table = 'sites';
  $table = $params->{'omega_table'} if ( $params->{'omega_table'} );
  $self->create_table_from_params( $self->dbc, $table, $self->sitewise_table_structure);

  foreach my $i (1 .. $pep_aln->length) {
    my $site = $results->{$i};

    my $note = $site->{note};
    # These two cases are pretty useless, so we'll skip 'em.
    if ( $self->param('sitewise_skip_gaps') ) {
      #next if ( $note eq 'all_gaps' );
      next if ($site->{nongap_count} == 0);
    }
    if ( $self->param('sitewise_skip_single_chars') ) {
      #next if ( $note eq 'single_char' );
      next if ($site->{nongap_count} == 1);
    }
    
    my $p = $self->replace($site, {
      data_id => $data_id,
      aln_position => $i
                           });
    $self->store_params_in_table($self->dbc, $table, $p);
  }
}


sub sitewise_table_structure {
  my $self = shift;

  return {
    data_id => 'int',
    aln_position => 'int',
    omega => 'float',
    omega_lower => 'float',
    omega_upper => 'float',
    lrt_stat => 'float',
    nongap_count => 'int',
    type => 'char16',
    note => 'char16',
    unique_keys => 'data_id,aln_position'
  };
}

1;
