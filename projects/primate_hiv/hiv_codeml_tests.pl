#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use File::Path;
use File::Basename;

use Bio::Greg::Codeml;
use FreezeThaw qw(freeze thaw cmpStr safeFreeze cmpStrHard);

use base (
  'Bio::Greg::StatsCollectionUtils',
  'Bio::Greg::Hive::Process',
  'Bio::Greg::Hive::Align',
  'Bio::Greg::Hive::AlignmentScores',
  'Bio::Greg::Hive::CountSubstitutions',
  'Bio::Greg::Hive::PhyloAnalysis',
);

my ( $aln_f ) = undef;
GetOptions(
  'aln=s' => \$aln_f
);
die("No input file given!") unless (defined $aln_f);

my ($f,$d) = fileparse($aln_f);
my $id = $f;
$id =~ s/\..*//;
print "$id\n";

my $aln = Bio::EnsEMBL::Compara::AlignUtils->from_file($aln_f);
$aln = Bio::EnsEMBL::Compara::AlignUtils->translate_ensembl($aln);

my $process = new Bio::Greg::Hive::Process;
$process->pretty_print($aln);

my $tree_str = qq^
(((((((Human, Chimpanzee), Gorilla), Orangutan), Macaque), Marmoset), Tarsier), (MouseLemur, Bushbaby));
^;
my $tree = Bio::EnsEMBL::Compara::TreeUtils->from_newick($tree_str);
my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
#print $treeI->ascii;

foreach my $leaf ($treeI->leaves) {
  foreach my $i ('', '_1', '_2', '_3') {
    my $new_child = new $leaf;
    $new_child->name($leaf->name.$i);

    $leaf->add_child($new_child);
  }
  $leaf->name('');
}
#print $treeI->ascii;

$tree = Bio::EnsEMBL::Compara::TreeUtils->from_treeI($treeI);
$tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_aln($tree, $aln);
$treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);

print $treeI->to_newick."\n";

my $tmp = $process->worker_temp_directory;
print "$tmp\n";

my $m0_f = _save_f('m0', 'txt');
if (!-e $m0_f) {
  ### Get M0 results.
  my $params = {
    model => 0,
    fix_blength => 0,
    method => 0,
    getSE => 1,
    Small_Diff => 1e-7,
    verbose => 1,
    debug => 1
  };

  print $m0_f."\n";
  print "  calculating M0...\n";
  my $m0 = Bio::Greg::Codeml->branch_model_likelihood( $treeI, $aln, $tmp, $params );
  my $lines = $m0->{lines};

  _out($m0_f, $lines);
}

print "  loading M0...\n";
my $m0 = _in($m0_f);
my $m0_tree = Bio::Greg::Codeml->parse_codeml_results($m0);
Bio::EnsEMBL::Compara::TreeUtils->transfer_branchlengths($m0_tree, $treeI);
print $treeI->to_newick."\n";

my $m7_f = _save_f('m7', 'txt');
if (!-e $m7_f) {
  ### Get M0 results.
  my $params = {
    model => 0,
    NSsites => 7,
    fix_blength => 0,
    method => 0,
    getSE => 1,
    Small_Diff => 1e-7
  };

  print $m7_f."\n";
  print "  calculating M7...\n";
  my $m7 = Bio::Greg::Codeml->branch_model_likelihood( $treeI, $aln, $tmp, $params );
  my $lines = $m7->{lines};
  _out($m7_f, $lines);
}

my $m8_f = _save_f('m8', 'txt');
if (!-e $m8_f) {
  ### Get M0 results.
  my $params = {
    model => 0,
    NSsites => 8,
    fix_blength => 0,
    method => 0,
    getSE => 1,
    Small_Diff => 1e-7
  };

  print $m8_f."\n";
  print "  calculating M8...\n";
  my $m8 = Bio::Greg::Codeml->branch_model_likelihood( $treeI, $aln, $tmp, $params );
  my $lines = $m8->{lines};
  _out($m8_f, $lines);
}
print "  loading M8...\n";
my $m8 = _in($m8_f);
my $m8_tree = Bio::Greg::Codeml->parse_codeml_results($m8);
Bio::EnsEMBL::Compara::TreeUtils->transfer_branchlengths($m8_tree, $treeI);
print $treeI->to_newick."\n";

rmtree($tmp);

sub _out {
  my $file = shift;
  my $lines = shift;

  open(OUT, ">$file");
  print OUT join("", @$lines) . "\n";
  close(OUT);
}

sub _in {
  my $file = shift;
  open(IN, $file);
  my @lines = <IN>;
  close(IN);
  return \@lines;
}

sub _save_f {
  my $name = shift;
  my $ext = shift;

  my $fld = '/nfs/users/nfs_g/gj1/src/greg-ensembl/projects/primate_hiv/manuscript_codeml';
  return $fld . '/' . $id . '_' . $name . '.' . $ext;
}
