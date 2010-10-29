#!/usr/bin/perl -w
 
use strict;
 
use Bio::SeqIO;
use Bio::Seq::Quality;
 
use Getopt::Long;
 
die "pass a fasta and a fasta-quality file\n"
  unless @ARGV;
 
 
my ($seq_infile,$qual_infile)
  = (scalar @ARGV == 1) ?($ARGV[0], "$ARGV[0].qual") : @ARGV;
 
## Create input objects for both a seq (fasta) and qual file

my $in_seq_obj;
my $in_qual_obj;

if ($seq_infile =~ m/gz/) {
  $in_seq_obj =
    Bio::SeqIO->new( -file   => "gunzip -c $seq_infile |",
                     -format => 'fasta',
    );
} else {
  $in_seq_obj =
    Bio::SeqIO->new( -file   => $seq_infile,
                     -format => 'fasta',
    );
}

if ($qual_infile =~ m/gz/) {
  $in_qual_obj =
    Bio::SeqIO->new( -file   => "gunzip -c $qual_infile |",
                     -format => 'qual',
    );
} else {
  $in_qual_obj =
    Bio::SeqIO->new( -file   => $qual_infile,
                     -format => 'qual',
    );
}

my $out_fastq_obj =
  Bio::SeqIO->new( -format => 'fastq'
  );
 

while (1){
  ## create objects for both a seq and its associated qual
  my $seq_obj  = $in_seq_obj->next_seq || last;
  my $qual_obj = $in_qual_obj->next_seq;
 
  die "foo!\n"
    unless
      $seq_obj->id eq
      $qual_obj->id;
 

  #print $seq_obj->seq."\n";
  #print $qual_obj->qual."\n";

  ## Here we use seq and qual object methods feed info for new BSQ
  ## object.
  my $bsq_obj =
    Bio::Seq::Quality->
    new( -id   => $seq_obj->id,
              -seq  => $seq_obj->seq,
              -qual => $qual_obj->qual,
    );
 
  ## and print it out.
  $out_fastq_obj->write_fastq($bsq_obj);
}
