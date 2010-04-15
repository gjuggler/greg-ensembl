#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use Bio::EnsEMBL::Registry;

#my $chr;
#my $start;
#my $end;
#GetOptions('chr=s' => \$chr,,
#           'start=s' => \$start,
#           'end=s' => \$end,
#           'strand=s' => \$strand
#	   );

Bio::EnsEMBL::Registry->load_registry_from_db(-host => 'ensdb-archive',-user => 'ensro',-port => 5304);
my $dbea = Bio::EnsEMBL::Registry->get_adaptor("human","core","dbentry");
my $txa = Bio::EnsEMBL::Registry->get_adaptor("human","core","transcript");
my $ga = Bio::EnsEMBL::Registry->get_adaptor("human","core","gene");
my $slice_a = Bio::EnsEMBL::Registry->get_adaptor("human","core","slice");


my $i=0;
while (<>) {
  chomp;
  print STDERR $i++." ".$_."\n";
  my ($chr,$start,$end,$strand) = split("\t",$_);
  my $stable_id = gene_id_from_region($chr,$start,$end,$strand);
  print STDERR "  -> $stable_id\n" if ($stable_id ne 'NA');
  print "$stable_id\n";
}

sub gene_id_from_region {
  my $chr = shift;
  my $start = shift;
  my $end = shift;
  my $strand_input = shift;

  $chr =~ s/chr//g; # Get rid of "chr"

  my $strand = 1;
  $strand = -1 if ($strand_input eq '-');

  # Get the forward strand slice, so rev-strand genes show up as strand==-1
  my $slice = $slice_a->fetch_by_region('chromosome',$chr,$start,$end,1);

  my @genes = @{$slice->get_all_Genes};
  foreach my $gene (@genes) {
    #printf STDERR "g:%s s:%d input_s:%d\n",$gene->stable_id,$gene->strand,$strand;
    if ($gene->strand == $strand) {
      return $gene->stable_id;
    }
  }
  return "NA";
}
