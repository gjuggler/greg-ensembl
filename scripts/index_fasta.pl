use Bio::DB::Qual;
use Bio::DB::Fasta;
use Bio::Greg::IndexedFasta;

my ($file) = $ARGV[0];
print $file."\n";

my $is_qual = ($file =~ m/qual/);

my $ind = new Bio::Greg::IndexedFasta;
my $sep = '';
$sep = ' ' if ($is_qual);

$ind->load_indexed_fasta($file,1,{sep => $sep});

#my $chr = 'a';
#2Bmy $s = 1;
#my $e = 5;
#my $seq = $ind->get_sequence_region($chr,$s,$e,{use_subindex => 0,sep => $sep});
#print "$seq\n";
#my $seq = $ind->get_sequence_region($chr,$s,$e,{sep => $sep});
#print "$seq\n";

#exit(0);


my @keys = $ind->get_all_ids;
foreach my $key (@keys) {
  print $key."\n";
  if ($file =~ m/qual/) {
    my $seq = $ind->get_qual_region($key,1,10);
    print "  $seq\n";
  } else {
    my $seq = $ind->get_sequence_region($key,1,10);
    print "  $seq\n";
  }
}

exit(0);

my $db;
if ($file =~ m/qual/) {
  $db = Bio::DB::Qual->new($file,(-debug => 1));
} elsif ($file =~ m/base/) {
  $db = Bio::DB::Fasta->new($file);
}
 
my @ids = $db->get_all_ids;
foreach my $id (@ids) {
  $id = 'chr22';
  print "  $id\n";
  if ($file =~ m/qual/) {
    $obj = $db->get_Qual_by_id($id);
    print "[".join("-",@{$obj->subqual(1 => 10)})."]\n";
  } elsif ($file =~ m/base/) {
    $obj = $db->get_Seq_by_id($id);
    print $obj->subseq(1 => 10)."\n";
  }
}
