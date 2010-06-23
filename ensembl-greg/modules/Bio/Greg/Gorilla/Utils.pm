package Bio::Greg::Gorilla::Utils;

use strict;

sub taxon_letter {
  my $class = shift;
  my $taxon = shift;
  
  my $taxon_to_letter = {
    9606 => 'H',
    9598 => 'C',
    9593 => 'G',
    9600 => 'O'
  };
  my $letter_to_taxon = {
    9606 => 'H',
    9598 => 'C',
    9593 => 'G',
    9600 => 'O'
  };

  if (defined $taxon_to_letter->{$taxon}) {
    return $taxon_to_letter->{$taxon};
  } elsif (defined $letter_to_taxon->{$taxon}) {
    return $letter_to_taxon->{$taxon};
  }
}

1;
