#!/usr/bin/env perl
use strict;

use Getopt::Long;
use File::Path;
use Bio::AlignIO;
use Bio::TreeIO;

use Bio::Greg::Hive::Align;
use Bio::Greg::Hive::PhyloAnalysis;
use Bio::Greg::PrimateHIV::WindowAnalysis;
use File::Copy;

my $arg = $ARGV[0];

my @prefixes = ('IFITM1','IFITM2','IFITM3','IFITM5');

foreach my $p (@prefixes) {
  print "$p\n";
#  align_file($p.'.fasta',$p.'_aligned.fasta',$p.'.nh');
  slr_file($p.'_aligned.fasta',$p.'.nh',$p.'.res',$p.'.txt');  
}

sub align_file {
  my $in_f = shift;
  my $out_f = shift;
  my $out_tree = shift;

  my $in = new Bio::AlignIO(-file => $in_f);
  my $aln = $in->next_aln;

  my $params = {alignment_prank_codon_model => 1};
  my $prank = new Bio::Greg::Hive::Align;
  my $prank_aln = $prank->align_with_prank($aln,undef,$params);

  my $out = new Bio::AlignIO(-file => '>'.$out_f, -format => 'fasta');
  $out->write_aln($prank_aln);
  
  copy('/tmp/process_tmp/output.2.dnd',$out_tree);
}


sub slr_file {
  my $in_f = shift;
  my $in_tree = shift;
  my $out_f = shift;
  my $out_table = shift;

  my $in = new Bio::AlignIO(-file => $in_f);
  my $aln = $in->next_aln;
  my $pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($aln);  

  $aln = Bio::EnsEMBL::Compara::AlignUtils->filter_stop_codons($aln);

  $in = new Bio::TreeIO(-file => $in_tree);
  my $tree = $in->next_tree;
  
  my $slr = new Bio::Greg::PrimateHIV::WindowAnalysis;

  my $params = {
    output_to_file => $out_f,
    window_output_file => $out_table
  };
  print "Running Slr...\n";

  my $slr_results;
  if (!-e $out_f) {
    $slr_results = $slr->run_sitewise_dNdS($tree,$aln,$params);
  } else {
    open(IN,"$out_f");
    my @output = <IN>;
    close(IN);
    $slr_results = $slr->parse_slr_output(\@output,$params);
  }
  my $hash = $slr->results_to_psc_hash($slr_results,$pep_aln);

  my ($ref_seq) = grep {$_->id =~ m/enst0/gi} $aln->each_seq;

  # Store windowed p-values.
  my @row_output;
  foreach my $size (10, 30, 50, 100, 9999) {
    my @rows = $slr->run_with_windows($size,$size/2,$aln,$tree,$ref_seq,$hash);
    push @row_output, @rows;
  }
  
  $slr->output_rows_to_file(\@row_output,$out_table);
}
