#!/usr/bin/env perl

use warnings;
use strict;

my $filename = $ARGV[0];
open(IN, $filename);
my @lines = <IN>;
close(IN);

print join("\t", qw^pos id aa_from aa_to codon_from codon_to confidence^) . "\n";

my $within_subs = 0;
my $node_id;
my $parent_id;
my $branch_line;
foreach my $line (@lines) {
  chomp $line;
  $line = strip($line);

  next if ($line eq '');
  next if ($line =~ m/refer to /i);

  $within_subs = 1 if ($line =~ m/Summary of changes along branches/i);
  $within_subs = 0 if ($within_subs && $line =~ m/List of extant and reconstructed sequences/i);

  if ($within_subs) {
    # Branch 5:    9..5  (ENSGGOP00000012863)  (n= 2.0 s=10.0)
    if ($line =~ m/Branch (\d+):\s+(\d+)\.\.(\d+)\s*(\S*)?\s*(\S*)?/gi) {
      $branch_line = $line;
      $parent_id = $2;
      $node_id = $3;
      if ($4 ne '' && !($4 =~ m/n=/i)) {
        $node_id = substr($4, 1, length($4)-2);
      }
    }
    if ($line =~ m/\s*(\d+) (\S+) \((\S)\) (.*) -> (\S+) \((\S)\)/g) {
      my $obj = {
        id => $node_id,
        parent_id => $parent_id,
        pos => $1,
        codon_a => $2,
        codon_b => $5,
        aa_a => $3,
        aa_b => $6,
        confidence => $4,
        line => $line,
        branch_line => $branch_line
      };
      next if ($obj->{codon_b} eq '---');
      next if ($obj->{codon_a} eq '---');

      printf "%d\t%s\t%s\t%s\t%s\t%s\t%.3f\n", 
                  $obj->{pos}, $obj->{id}, $obj->{aa_a}, $obj->{aa_b},
                  $obj->{codon_a}, $obj->{codon_b}, $obj->{confidence};
    }
  }  
}

sub strip {
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}

