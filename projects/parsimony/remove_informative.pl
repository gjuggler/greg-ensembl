#!/usr/bin/env perl
#use strict;
use Getopt::Long;
use File::Path;
use File::Basename;
use Bio::SimpleAlign;
use Bio::AlignIO;

my $file = "";
my $file_out = "uninformative.phy";
GetOptions('file=s' => \$file,
	   'out=s' => \$file_out
	   );

my $alnio = Bio::AlignIO->new(-file => $file,
			      -format => 'phylip');

my $out = Bio::AlignIO->new(-file => ">$file_out",
			    -format => 'phylip'
			    );

print is_col_informative("TAATGGGGGG")."\n";
print is_col_informative("TATGGGGGGG")."\n";
print is_col_informative("TACGGGGGGG")."\n";

while (my $aln = $alnio->next_aln()) {
  $aln = remove_uninformative_sites($aln);
  $out->write_aln($aln);
}

exit(0);

sub remove_uninformative_sites {
  my $aln = shift;

  my @keeper_columns;  
  foreach my $col (1 .. $aln->length) {
    my $column_string = column_string($aln,$col);
    if (!is_col_informative($column_string)) {
      push @keeper_columns,$col;
      #print $column_string."  ".$col."\n";
    } else {
      #print $column_string."\n";
    }
  }

  my $new_aln = new $aln;
  foreach my $seq ($aln->each_seq) {
    my $new_seq = new Bio::LocatableSeq(-id=>$seq->id);
    my $sequence = extract_cols($seq->seq,\@keeper_columns);
    $new_seq->seq($sequence);
    $new_aln->add_seq($new_seq);
  }
  return $new_aln;
}

sub extract_cols {
  my $seq = shift;
  my $cols = shift;

  my $str = "";
  foreach my $col (@{$cols}) {
    $str .= substr($seq,$col,1);
  }
  return $str;
}

sub column_string {
  my $aln = shift;
  my $col = shift;

  my $str = "";
  foreach my $seq ($aln->each_seq) {
    my $char = substr($seq->seq,$col,1);
    $str .= $char unless ($char eq '-');
  }
  return $str;
}

sub is_col_informative {
  my $str = shift;

  # Returns 1 if this column is informative.
  # From Yang's CME book: For a site to be parsimony-informative, at least two characters have to
  # be observed, each at least twice.
  my $char_obs;

  foreach my $char (split("",$str)) {
    if (!defined $char_obs->{$char}) {
      $char_obs->{$char} = 1;      
    } else {
      $char_obs->{$char} = $char_obs->{$char} + 1;
    }
  }

  my $num_chars_above_two = 0;
  foreach my $char_val (values %{$char_obs}) {
    $num_chars_above_two++ if ($char_val >= 2);
    return 1 if ($num_chars_above_two == 2);
  }
  return 0;
}
