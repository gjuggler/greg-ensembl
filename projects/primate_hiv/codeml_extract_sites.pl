#!/usr/bin/env perl

use warnings;
use strict;

my $filename = $ARGV[0];
open(IN, $filename);
my @lines = <IN>;
close(IN);

print "pos\tdnds\tprob_gt_one\n";

my $within_bayes = 0;
my $number_of_bayes_classes = 0;
foreach my $line (@lines) {
  chomp $line;
  $line = strip($line);

  next if ($line eq '');
  next if ($line =~ m/refer to /i);

  $within_bayes = 1 if ($line =~ m/Bayes Empirical Bayes \(BEB\) probabilities/i);
  $within_bayes = 0 if ($within_bayes && $line =~ m/Positively selected sites/i);
  $within_bayes = 0 if ($within_bayes && $line =~ m/Ancestral reconstruction by CODONML/i);

  if ($within_bayes && $line =~ m/probabilities for (\d+) classes/i) {
    $number_of_bayes_classes = $1;
    print STDERR "  NUMBER OF CLASSES: $number_of_bayes_classes\n";
  }

  if ($within_bayes) {
    my @bits = split( /\s+/, $line );
    my $pos = shift @bits;
    my $residue = shift @bits;

    my $se              = pop @bits;
    my $plus_minus_sign = pop @bits;
    my $omega           = pop @bits;

    my $prob_gt_one = 0;
    my $last_class_posterior = $bits[ $number_of_bayes_classes - 1 ];
    my $prob_gt_one = $last_class_posterior;
    
    printf "%d\t%.4f\t%.4f\n", $pos, $omega, $prob_gt_one;
  }
  

}

sub strip {
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}
