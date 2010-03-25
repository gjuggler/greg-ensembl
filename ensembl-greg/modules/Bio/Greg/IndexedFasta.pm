package Bio::Greg::IndexedFasta;
use warnings;
use strict;

sub new {
  my ( $class, @args ) = @_;
  ## Allows to create a new object from an existing one with $object->new
  $class = ref($class) if ( ref($class) );
  my $self = $class->alloc(@args);
  $self->init;
  return $self;
}

sub alloc {
  my ( $class, @args ) = @_;
  my $self = {};
  bless $self, $class;
  return $self;

}

sub init {
  my $self = shift;
  $self->{'_position_hash'} = undef;
  $self->{'_filehandle'}    = undef;
}

sub position_hash {
  my ( $self, $hash ) = @_;

  if ( defined $hash ) {
    $self->{'_position_hash'} = $hash;
  }

  return $self->{'_position_hash'};
}

sub filehandle {
  my ( $self, $fh ) = @_;

  if ( defined $fh ) {
    $self->{'_filehandle'} = $fh;
  }
  return $self->{'_filehandle'};
}

sub get_sequence {
  my $self = shift;
  my $key  = shift;

  my $filehandle = $self->filehandle;

  if ( !defined $filehandle ) {
    die("No fasta file was loaded! Dying...");
  }

  my $position_hash = $self->position_hash;

  my $position = $position_hash->{$key};
  if ( !defined $position ) {

    #die ("No sequence found for $key!");
    return undef;
  }

  seek( $filehandle, $position, 0 );
  my $line;
  my $seq = "";
  while ( $line = <$filehandle> ) {
    if ( $line =~ m/\d/ ) {

      # Quality scores -- keep a single space between each one.
      $line =~ s/\s+/ /g;
    } else {

      # DNA seq -- no spaces.
      $line =~ s/\s+//g;
    }
    if ( $line =~ /^>/ ) {
      if ( $seq eq "" ) {

        #$line =~ s/^>\S+//;
        $line = "";
      } else {
        last;
      }
    }
    $seq .= $line;

    #print $seq."\n";
  }
  return $seq;
}

sub get_sequence_region {
  my $self  = shift;
  my $key   = shift;
  my $start = shift;
  my $end   = shift;

  my $filehandle = $self->filehandle;
  if ( !defined $filehandle ) {
    die("No fasta file was loaded! Dying...");
  }

  my $position_hash = $self->position_hash;
  my $position      = $position_hash->{$key};
  if ( !defined $position ) {

    #$key =~ s/genescaffold/contig/i;
    #$key =~ s/scaffold/contig/i;
    $position = $position_hash->{$key};

    #die ("No sequence found for $key!") if (!defined $position);
    return undef if ( !defined $position );
  }

  my $current_position = 0;
  my $numbers          = 0;

  seek( $filehandle, $position, 0 );
  my $line;
  my $seq = "";
  while ( $line = <$filehandle> ) {
    if ( $line =~ /^>/ ) {
      next if ( $seq eq "" );
      return $seq;
    }
    $numbers = 1 if ( $line =~ m/\d/ );

    my @toks = split( /\s/, $line );
    for ( my $i = 0 ; $i < scalar(@toks) ; $i++ ) {
      if ( $current_position == $end ) {

        #print $seq."\n";
        return $seq;
      }
      if ( $current_position >= $start ) {
        $seq .= $toks[$i];
        $seq .= " " if ($numbers);
      }
      $current_position++;
    }
  }
  return undef;
}

sub load_indexed_fasta {
  my $self      = shift;
  my $orig_file = shift;

  my $index_file = $orig_file . ".ind";

  # Create the index if it doesn't exist yet.
  if ( !-e $index_file ) {
    $self->index_fasta($orig_file);
  }

  # Load the index into memory.
  open( IN, $index_file );
  my $line;
  my $position_hash;
  while ( $line = <IN> ) {
    chomp $line;
    my @toks = split( "\t", $line );
    $position_hash->{ $toks[0] } = $toks[1];
  }
  close(IN);

  # Open the original file and keep the filehandle hanging around.
  my $filehandle;
  open( $filehandle, $orig_file );

  #use IO::Uncompress::Gunzip;
  #$filehandle = new IO::Uncompress::Gunzip $orig_file;

  $self->filehandle($filehandle);
  $self->position_hash($position_hash);
}

sub DESTROY {
  my $self = shift;

  if ( defined $self->filehandle ) {
    close( $self->filehandle );
    $self->{'_filehandle'} = undef;
  }
  if ( defined $self->position_hash ) {
    $self->{'_position_hash'} = undef;
  }
}

sub index_fasta {
  my $class   = shift;
  my $file_in = shift;

  printf "Indexing fasta file %s...\n", $file_in;

  #open(OUT,">${file_in}.ind");
  #print OUT "Indexing in progress...\n";
  #close(OUT);

  open( IN, "${file_in}" );

  #use IO::Uncompress::Gunzip;
  #my $in = new IO::Uncompress::Gunzip $file_in;

  my $line;
  my $prev_tell = 0;
  my $tell      = 0;
  my $qual_hash = {};
  while ( $line = <IN> ) {
    chomp $line;
    $tell = tell(IN);
    if ( substr( $line, 0, 1 ) eq ">" ) {
      $line =~ s/\s//g;
      $qual_hash->{ substr( $line, 1 ) } = $prev_tell;
    }
    $prev_tell = $tell;
  }
  close(IN);

  open( OUT, ">${file_in}.ind" );
  my @keys = keys( %{$qual_hash} );
  map { print OUT sprintf( "%s\t%d\n", $_, $qual_hash->{$_} ) } @keys;
  close(OUT);

  printf "Done!\n";
  undef $qual_hash;

}

1;
