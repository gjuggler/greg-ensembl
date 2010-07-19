package Bio::Greg::Gorilla::Utils;

use strict;

sub taxon_short {
  my $class = shift;
  my $taxon = shift;
  
  my $taxon_to_short = {
    9606 => 'Hsap',
    9598 => 'Ptro',
    9593 => 'Ggor',
    9600 => 'Ppyg',
    9544 => 'Mmul'
  };
  my $short_to_taxon = {
    'Hsap' => 9606,
    'Ptro' => 9598,
    'Ggor' => 9593,
    'Ppyg' => 9600,
    'Mmul' => 9544
  };
  if (defined $taxon_to_short->{$taxon}) {
    return $taxon_to_short->{$taxon};
  } elsif (defined $short_to_taxon->{$taxon}) {
    return $short_to_taxon->{$taxon};
  }
}

sub taxon_letter {
  my $class = shift;
  my $taxon = shift;
  
  my $taxon_to_letter = {
    9606 => 'H',
    9598 => 'C',
    9593 => 'G',
    9600 => 'O',
    9544 => 'M'
  };
  my $letter_to_taxon = {
    'H' => 9606,
    'C' => 9598,
    'G' => 9593,
    'O' => 9600,
    'Hsap' => 9606,
    'Ptro' => 9598,
    'Ggor' => 9593,
    'Ppyg' => 9600,
    'Mmul' => 9544
  };

  if (defined $taxon_to_letter->{$taxon}) {
    return $taxon_to_letter->{$taxon};
  } elsif (defined $letter_to_taxon->{$taxon}) {
    return $letter_to_taxon->{$taxon};
  }
}

1;
