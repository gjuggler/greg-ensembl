#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use Bio::EnsEMBL::Registry;

Bio::EnsEMBL::Registry->load_registry_from_db(-host => 'ens-livemirror',-user=>'ensro');
my $dbea = Bio::EnsEMBL::Registry->get_adaptor("human","core","dbentry");
my $txa = Bio::EnsEMBL::Registry->get_adaptor("human","core","transcript");
my $ga = Bio::EnsEMBL::Registry->get_adaptor("human","core","gene");
my $sa = Bio::EnsEMBL::Registry->get_adaptor("human","core","slice");

my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                             -dbname => 'homo_sapiens_core_54_36p',
                                             -user => 'anonymous',
                                             -pass => '',
                                             -host => 'ensembldb.ensembl.org',
                                             -port => 5306,
                                            );
my $old_slice_adaptor = $dba->get_SliceAdaptor;

my $file = "variant_coordinates.txt";

open(IN,"$file");
my @lines = <IN>;
close(IN);

print join("\t",'var_id','dna_36','position_37','dna_37','gene_id','transcript_id','protein_id','position_peptide')."\n";

open(OUT,">failures.txt");

foreach my $line (@lines) {
  chomp $line;
  #print "line: [$line]\n";
  next if ($line =~ m/var_id/gi);
  
  my $var_id = $line;
  my ($chr,$pos,$nuc) = split("_",$line);

  #next unless ($chr eq '11' && $pos eq '5642842');
  
  my $slice = $sa->fetch_by_region( "chromosome", $chr, $pos, $pos, 1,
                                    "NCBI36");
  my $old_slice = $old_slice_adaptor->fetch_by_region("chromosome",$chr,$pos,$pos,1);
  my $old_slice_seq = $old_slice->seq;

  print STDERR $slice->name, "\n";
  my $proj_segments = $slice->project("chromosome", "GRCh37");
  if (scalar(@$proj_segments) > 1) {
    print OUT join("\t",$var_id,"Multiple segments found when projecting to GRCh37")."\n";
    next;
  }
  foreach my $proj_segment (@$proj_segments) {
    my $proj_slice = $proj_segment->to_Slice;
    print STDERR $proj_slice->name."\n";
    my $new_pos = $proj_slice->start;
    #my @txs = @{$txa->fetch_all_by_Slice($proj_slice)};
    my @txs = @{$proj_slice->get_all_Transcripts};
    if (scalar(@txs) == 0) {
      print OUT join("\t",$var_id,"No transcripts found within region")."\n";
      next;
    }
    foreach my $tx (@txs) {
      if ($tx->translateable_seq eq '') {
        #print STDERR " Pseudogene: $var_id\n";
        next;
      } else {
        #print STDERR " Gene: ".$tx->stable_id."\n";
      }

      # Re-fetch the tx so it's based on the genome coordinates.
      $tx = $txa->fetch_by_stable_id($tx->stable_id);

      my $tx_id = $tx->stable_id;
      my $gene = $ga->fetch_by_transcript_stable_id($tx_id);
      my $gene_id = $gene->stable_id;
      my $tl_id = $tx->translation->stable_id;

      my $tm = $tx->get_TranscriptMapper;
      my @pep_coords = $tm->genomic2pep($new_pos,$new_pos,$tx->strand);
      @pep_coords = grep {$_->isa('Bio::EnsEMBL::Mapper::Coordinate')} @pep_coords;
      if (scalar(@pep_coords) == 0) {
        print OUT join("\t",$var_id,"No transcripts found within region")."\n";
        next;
      }
      foreach my $coord (@pep_coords) {
        #print STDERR $coord->start."\n";
        my $pep_pos = $coord->start;
        my $new_dna = $proj_slice->seq;
        print join("\t",$var_id,$old_slice_seq,$new_pos,$new_dna,$gene_id,$tx_id,$tl_id,$pep_pos)."\n";
        print STDERR "  ".join("\t",$var_id,$old_slice_seq,$new_pos,$new_dna,$gene_id,$tx_id,$tl_id,$pep_pos)."\n";
      }
    }
  }
}

close(OUT);
