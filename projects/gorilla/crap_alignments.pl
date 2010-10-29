#!/usr/bin/env perl
use warnings;
use strict;

use Bio::Greg::EslrUtils;
use Bio::EnsEMBL::Registry;

my $url = 'mysql://ensadmin:ensembl@ens-research:3306/gj1_gor_58';
my $clean = 1;

my $args = Bio::Greg::EslrUtils->url_to_hashref($url);
Bio::EnsEMBL::Registry->load_registry_from_multiple_dbs(
  {
    -host => 'ens-livemirror',
    -user => 'ensro'
  }
  );
my $compara = Bio::EnsEMBL::Registry->get_DBAdaptor('multi','compara');
my $mba = $compara->get_MemberAdaptor;
my $ha = $compara->get_HomologyAdaptor;
my $gdba = $compara->get_GenomeDBAdaptor;

my $ref = 9606;
my $ref_gdb = $gdba->fetch_by_taxon_id($ref);

my $mutation_run_length_counts;

foreach my $other (9593,9598,9600,9544) {
  my $other_gdb = $gdba->fetch_by_taxon_id($other);

  $mutation_run_length_counts->{$other} = [];

  my @members = @{$mba->fetch_all_by_source_taxon('ENSEMBLGENE',$ref)};
  #@members = @members[1..500];

  my $i=0;
  foreach my $ref_member (@members) {
    $i++;
    print STDERR "$i / ".scalar(@members)."\n";
    #print $ref_member->stable_id."\n";
    #print $other_gdb->name."\n";
    my @homologies = @{$ha->fetch_all_by_Member_paired_species($ref_member,$other_gdb->name)};
    foreach my $homology (@homologies) {
      my $desc = $homology->description;

      if ($desc eq 'ortholog_one2one') {
	my $aln = $homology->get_SimpleAlign(-cdna => 1);
        $aln = $aln->remove_gaps(undef,1);

	my @substitution_run_lengths = get_substitution_runs($aln);

	foreach my $run_length (@substitution_run_lengths) {
          print join(",",$other,$run_length)."\n";

	  # Count up from 1 to $run_length, adding a value to the hash for each length.
	  foreach my $i (1 .. $run_length) {
	    if (!defined $mutation_run_length_counts->{$other}->[$i]) {
	      $mutation_run_length_counts->{$other}->[$i] = 0;
	    } else {
	      $mutation_run_length_counts->{$other}->[$i]++;
	    }
 	  }
        }
      }
    }
  }
}

#my @lengths = (1..5,10,20,30,50,100);
#print join("\t",'species',@lengths)."\n";
#foreach my $species (keys %$mutation_run_length_counts) {
#  my @arr = @{$mutation_run_length_counts->{$species}};
#  my @clean_subs;
#  foreach my $i (@lengths) {
#    my $count = $arr[$i] || 0;
#    push @clean_subs, $count;
#  }
#  print join("\t",$species,@clean_subs) . "\n";
#}


sub get_substitution_runs {
  my $aln = shift;

  my ($seq1,$seq2) = $aln->each_seq;
  my $str1 = $seq1->seq;
  my $str2 = $seq2->seq;

  my $sub_run_length = 0;
  my @sub_runs;
  for (my $i=0; $i <= $aln->length; $i++) {
    my $c1 = substr($str1,$i-1,1);
    my $c2 = substr($str2,$i-1,1);

    if ($c1 ne $c2 && $c1 ne '-' && $c2 ne '-') {
      $sub_run_length++;
    } elsif ($sub_run_length > 0) {
      push @sub_runs, $sub_run_length;
      $sub_run_length = 0;
    }
  }
  return @sub_runs;
}
