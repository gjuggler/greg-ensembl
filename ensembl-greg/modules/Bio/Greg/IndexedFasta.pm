package Bio::Greg::IndexedFasta;
use warnings;
use strict;

sub new {
  my ( $class, @args ) = @_;
  ## Allows to create a new object from an existing one with $object->new
  $class = ref($class) if ( ref($class) );
  my $self = $class->alloc(@args);
  $self->init(@args);
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
  my $file = shift;
  $self->{'_position_hash'} = undef;
  $self->{'_filehandle'}    = undef;

  #print STDERR "Loading index for [$file]\n";
  $self->load_indexed_fasta($file) if (defined $file && -e $file);
}

sub position_hash {
  my ( $self, $hash ) = @_;

  if ( defined $hash ) {
    $self->{'_position_hash'} = $hash;
  }

  return $self->{'_position_hash'};
}

sub get_all_ids {
  my $self = shift;
  return keys %{ $self->position_hash };
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

sub has_key {
  my $self = shift;
  my $key = shift;

  my $hash = $self->position_hash;
  #print "($key)\n";
  return 1 if ($hash->{$key});
  return 0;
}

sub get_qual_region {
  my $self  = shift;
  my $key   = shift;
  my $start = shift;
  my $end   = shift;

  return $self->get_sequence_region( $key, $start, $end, { sep => ' ' } );
}

sub get_sequence_region {
  my $self   = shift;
  my $key    = shift;
  my $start  = shift;
  my $end    = shift;
  my $params = shift || {};

  my $sep = $params->{sep};
  $sep = '' unless ( defined $sep );
  my $use_subindex = $params->{use_subindex};
  $use_subindex = 1 unless ( defined $use_subindex );

  my $filehandle = $self->filehandle;
  if ( !defined $filehandle ) {
    die("No fasta file was loaded! Dying...");
  }

  my $position_hash = $self->position_hash;

  my $index_freq = 100_000;
  my $position;
  my $orig_key = $key;
  if ($start > $index_freq && $use_subindex) {
    # We'll make use of the sub-sequence hashing to go within our sequence.
    my $least_index = sprintf("%d",($start/$index_freq));
    my $start_from = ($least_index*$index_freq);
    my $new_id = $key."_".$start_from;
    #print "New id:[$new_id]\n";
    $position = $position_hash->{$new_id};
    if (defined $position) {
      #print "Used sub-index!\n";
      #print "New start: $start_from $position\n";
      $start = $start - $start_from;
      $end = $end - $start_from;
    }
  } else {
    # Use the boring old sequence identifier.
    $position = $position_hash->{$key};
  }

  if ( !defined $position ) {
    return undef;
  }

  my $current_position = 0;

  seek( $filehandle, $position, 0 );
  my $line;
  my $seq = "";
  my @toks;
  my $skip_count = 0;
  while ( $line = <$filehandle> ) {
    chomp $line;
    if ( $line =~ /^>/ ) {
      next if ( $seq eq "" );
      #print $skip_count."\n";
      return $seq;
    }

    #    $line =~ s/^\s+//;
    #    $line =~ s/\s+$//;

    if ( $sep eq ' ' ) {
      @toks = split( ' ', $line );
    } else {
      @toks = split( '', $line );
    }

    for ( my $i = 0 ; $i < scalar(@toks) ; $i++ ) {
      next if ( $toks[$i] eq ' ' || $toks[$i] eq '' );
      $current_position++;
      $skip_count++;
      if ( $current_position == $end ) {
        #print "$skip_count\n";
        return $seq . $toks[$i];
      }
      if ( $current_position >= $start ) {
        $seq .= $toks[$i] . $sep;
      }
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
    $self->index_fasta( $orig_file, @_ );
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
  my $class         = shift;
  my $file_in       = shift;
  my $index_subseqs = shift;
  my $params        = shift || {};

  $index_subseqs = 0 unless ( defined $index_subseqs );
  my $sep = $params->{sep};
  $sep = '' unless ( defined $sep );
  my $index_spacing = $params->{index_spacing};
  $index_spacing = 100_000 unless ( defined $index_spacing );

  printf "Indexing fasta file %s...\n", $file_in;
  print "  (with subseqs each $index_spacing)\n" if ($index_subseqs);

  open( IN, "${file_in}" );

  my $line;
  my $cur_seq_id;
  my $prev_tell  = 0;
  my $tell       = 0;
  my $char_count = 0;
  my $char_str   = '';
  my @toks;
  my $qual_hash = {};
  while ( $line = <IN> ) {
    chomp $line;
    $tell = tell(IN);
    if ( substr( $line, 0, 1 ) eq ">" ) {
      $line =~ s/\s//g;
      $line =~ s/\[.*\]//g;
      $cur_seq_id = substr( $line, 1 );
      print "  $cur_seq_id $tell\n";
      $qual_hash->{$cur_seq_id} = $tell;
      $char_count = 0;
    } else {
      if ($index_subseqs) {
        @toks = split( $sep, $line );
        for ( my $i = 0 ; $i < scalar(@toks) ; $i++ ) {
          next if ( $toks[$i] eq ' ' || $toks[$i] eq '' );
          my $tok = $toks[$i];
          $char_count++;
          if ( $char_count % $index_spacing == 0 ) {
            {
              # Index a specific region in the file.
              # The current tell() call brought us to the end of this line,
              # so subtract what's left of this line from the current tell() to get
              # the right position to start from.
              use bytes;
              my $so_far = join( $sep, @toks[ 0 .. $i ] );
              my $total_length = length($line);
              my $so_far_length = length($so_far);
              my $cur_location = $tell - $total_length + $so_far_length;
              die("Format error: no sequence ID encountered before sequence data!") unless (defined $cur_seq_id);
              my $subindex_key = $cur_seq_id.'_'.$char_count;
              $qual_hash->{$subindex_key} = $cur_location;
              print "  $subindex_key $cur_location\n";
            }
          }
        }
      }
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
