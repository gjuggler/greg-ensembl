#!/usr/bin/env perl

use warnings;
use strict;
use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::TreeUtils;
use File::Path;
use File::Basename;

my $C = "Bio::EnsEMBL::Compara::ComparaUtils";
my $T = "Bio::EnsEMBL::Compara::TreeUtils";

my $str = qq^(((((((((((((((
((Human:0.006591,Chimp:0.006639):0.002184,Gorilla:0.009411):0.009942,
Orangutan:0.018342):0.014256,Rhesus:0.036199):0.021496,
Marmoset:0.066389):0.056911,Tarsier:0.135169):0.011307,
(Mouse_lemur:0.091452,Bushbaby:0.128984):0.035463):0.015304,
TreeShrew:0.183583):0.004688,(((((Mouse:0.083220,Rat:0.090564):0.196605,
Kangaroo_rat:0.209532):0.022555,Guinea_Pig:0.223415):0.009828,
Squirrel:0.146894):0.025042,
(Rabbit:0.116009,Pika:0.198295):0.100037):0.015355):0.020666,
(((Alpaca:0.105252,(Dolphin:0.064182,Cow:0.121911):0.025111):0.039691,
((Horse:0.107726,(Cat:0.097971,Dog:0.100888):0.049486):0.006252,
(Microbat:0.141155,Megabat:0.111787):0.033187):0.004179):0.011699,
(Hedgehog:0.220580,Shrew:0.266859):0.056117):0.021065):0.023276,
(((Elephant:0.083775,Rock_hyrax:0.152633):0.026190,Tenrec:0.240221):0.049905,
(Armadillo:0.115179,Sloth:0.096272):0.052373):0.006713):0.232748,
Opossum:0.325899):0.072430,Platypus:0.453916):0.109903,
((Chicken:0.166386,Zebra_finch:0.170717):0.199763,
Lizard:0.509545):0.108130):0.166150,X_tropicalis:0.852482):0.300396,
(((Tetraodon:0.224774,Fugu:0.205294):0.191836,
(Stickleback:0.313967,Medaka:0.478451):0.058404):0.322824,
Zebrafish:0.731166):0.155214):0.511293,Lamprey:0.511293);
^;
my $v_out = tree_subset($str,['Human','Chimp','Rhesus','Mouse','Rat','Guinea_Pig','Horse','Dog','Cow']);
print $v_out."\n";

exit(0);

# From http://hgdownload.cse.ucsc.edu/goldenPath/hg18/phastCons44way/vertebrate.mod
my $v_tree_str = "(((((((((((((((((hg18:0.006591,panTro2:0.006639):0.002184,gorGor1:0.009411):0.009942,ponAbe2:0.018342):0.014256,rheMac2:0.036199):0.021496,calJac1:0.066389):0.056911,tarSyr1:0.135169):0.011307,(micMur1:0.091452,otoGar1:0.128984):0.035463):0.015304,tupBel1:0.183583):0.004688,(((((mm9:0.083220,rn4:0.090564):0.196605,dipOrd1:0.209532):0.022555,cavPor3:0.223415):0.009828,speTri1:0.146894):0.025042,(oryCun1:0.116009,ochPri2:0.198295):0.100037):0.015355):0.020666,(((vicPac1:0.105252,(turTru1:0.064182,bosTau4:0.121911):0.025111):0.039691,((equCab2:0.107726,(felCat3:0.097971,canFam2:0.100888):0.049486):0.006252,(myoLuc1:0.141155,pteVam1:0.111787):0.033187):0.004179):0.011699,(eriEur1:0.220580,sorAra1:0.266859):0.056117):0.021065):0.023276,(((loxAfr2:0.083775,proCap1:0.152633):0.026190,echTel1:0.240221):0.049905,(dasNov2:0.115179,choHof1:0.096272):0.052373):0.006713):0.232748,monDom4:0.325899):0.072430,ornAna1:0.453916):0.109903,((galGal3:0.166386,taeGut1:0.170717):0.199763,anoCar1:0.509545):0.108130):0.166150,xenTro2:0.852482):0.300396,(((tetNig1:0.224774,fr2:0.205294):0.191836,(gasAcu1:0.313967,oryLat2:0.478451):0.058404):0.322824,danRer5:0.731166):0.155214):0.511293,petMar1:0.511293);";
my $pl_tree_str = "(((((((((((hg18:0.006591,panTro2:0.006639):0.002184,gorGor1:0.009411):0.009942,ponAbe2:0.018342):0.014256,rheMac2:0.036199):0.021496,calJac1:0.066389):0.056911,tarSyr1:0.135169):0.011307,(micMur1:0.091452,otoGar1:0.128984):0.035463):0.015304,tupBel1:0.183583):0.004688,(((((mm9:0.083220,rn4:0.090564):0.196605,dipOrd1:0.209532):0.022555,cavPor3:0.223415):0.009828,speTri1:0.146894):0.025042,(oryCun1:0.116009,ochPri2:0.198295):0.100037):0.015355):0.020666,(((vicPac1:0.105252,(turTru1:0.064182,bosTau4:0.121911):0.025111):0.039691,((equCab2:0.107726,(felCat3:0.097971,canFam2:0.100888):0.049486):0.006252,(myoLuc1:0.141155,pteVam1:0.111787):0.033187):0.004179):0.011699,(eriEur1:0.220580,sorAra1:0.266859):0.056117):0.021065):0.023276,(((loxAfr2:0.083775,proCap1:0.152633):0.026190,echTel1:0.240221):0.049905,(dasNov2:0.115179,choHof1:0.096272):0.052373):0.006713);";
my $p_tree_str = "(((((((hg18:0.006591,panTro2:0.006639):0.002184,gorGor1:0.009411):0.009942,ponAbe2:0.018342):0.014256,rheMac2:0.036199):0.021496,calJac1:0.066389):0.056911,tarSyr1:0.135169):0.011307,(micMur1:0.091452,otoGar1:0.128984):0.035463);";
my $a_bglobin = "((xenlaev:0.51040247, xentrop:0.78663390):0.001, ((duck:0.17345744, chicken:0.05771191):0.57801375, (marsupial:0.89258965, ((hamster:0.36642328, (rat:0.13974134, mouse:0.19356138):0.10343856):0.28360696, ((elephseal:0.21571060, (pig:0.32152498, (cow:0.06157821, sheep:0.09541233):0.13588222):0.06570764):0.05974793, ((bushbaby:0.14974872, (hare:0.02840006, rabbit:0.06031739):0.15985391):0.01637450, (human:0.17351861, tarsier:0.22127351):0.06466335):0.07090812):0.05240999):0.33978458):0.19069468):1.86481009);";
my $a_art = "((3:0.2,(Human:0.1,2:0.1):0.1):0.05,(4:0.2,(5:0.1,6:0.1):0.1):0.05);";

my $pfam = read_file("trees/pfam.nhx");


tree_stats($v_tree_str,'ucsc vertebrate');
tree_stats($pl_tree_str,'ucsc placental');
tree_stats($p_tree_str,'ucsc primate');
tree_stats($a_bglobin,'anisimova bglobin');
tree_stats($a_art,'anisimova artificial');
tree_stats($pfam,'pfam family');

sub read_file {
  my $file = shift;
  open(IN,$file);
  my @lines = <IN>;
  close(IN);
  map {chomp} @lines;
  return join("",@lines);
}

sub tree_subset {
  my $tree_str = shift;
  my $keep_species = shift;

  my @keepers = @{$keep_species};

  my $tree = $T->from_newick($tree_str);
  $tree = $T->keep_members_by_method_call($tree,\@keepers,'name');

  use Bio::EnsEMBL::Compara::ComparaUtils;

  return $tree->newick_format;
}

sub tree_stats {
  my $tree_str = shift;
  my $label = shift;
  my $tree = $T->from_newick($tree_str);

  print "$label: {\n";
  print " mean_path:".$T->mean_path($tree)."\n";
  print " total_dist:".$T->total_distance($tree)."\n";
  print " max_dist:".$T->max_distance($tree)."\n";
  print "}\n";
}
