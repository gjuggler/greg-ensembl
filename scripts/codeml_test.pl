#!/usr/bin/perl -w 
use strict;
use Bio::Greg::Codeml;

my $file = $ARGV[0];
my $file2 = $ARGV[1];

my @lines = Bio::Greg::Codeml::file_array($file);
my @lines2 = Bio::Greg::Codeml::file_array($file2);
my @all = (@lines,@lines2);

#Bio::Greg::Codeml->extract_params(\@all);

my $tree = Bio::Greg::Codeml->parse_codeml_results(\@all);

foreach my $node ($tree->get_nodes) {
#  next if (!$node->is_Leaf);

  my $subs = $node->get_tag_values('substitutions');
#  print $node->id."\n";
  foreach my $key (sort {$a <=> $b} keys %$subs) {
    my $subst = $subs->{$key};
#    printf "  %s\n", $subst->{line};
  }
  printf "%s %s %s\n", $node->id, $node->get_tag_values('dN/dS'), $node->get_tag_values('dN/dS_se');
}
