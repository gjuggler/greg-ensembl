=head1 NAME - Bio::EnsEMBL::Utils::SeqDumper

=head1 SYNOPSIS

  $seq_dumper = Bio::EnsEMBL::Utils::SeqDumper;

  #don't dump snps or repeats
  $seq_dumper->disable_feature_type('repeat');
  $seq_dumper->disable_feature_type('variation');

  #dump EMBL format to STDOUT
  $seq_dumper->dump($slice, 'EMBL');

  #dump GENBANK format to a file
  $seq_dumper->dump($slice, 'GENBANK', 'out.genbank');

  #dump FASTA format to a file
  $seq_dumper->dump($slice, 'FASTA', 'out.fasta');

=head1 DESCRIPTION

  A relatively simple and lite-weight flat file dumper for Ensembl slices.
  The memory efficiency could be improved and this is currently
  not very good for dumping very large sequences such as whole chromosomes.

=head1 CONTACT

  Contact the Ensembl development list with questions: <ensembl-dev@ebi.ac.uk>

=cut

use strict;

package Bio::EnsEMBL::Utils::SeqDumper;

use IO::File;
use vars qw(@ISA);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

#keys must be uppercase
my $DUMP_HANDLERS = 
  { 'FASTA'     => \&dump_fasta,
    'EMBL'      => \&dump_embl,
    'GENBANK'   => \&dump_genbank };

my @COMMENTS = 
  ('This sequence was annotated by the Ensembl system. Please visit ' .
   'the Ensembl web site, http://www.ensembl.org/ for more information.',

   'All feature locations are relative to the first (5\') base ' .
   'of the sequence in this file.  The sequence presented is '.
   'always the forward strand of the assembly. Features ' .
   'that lie outside of the sequence contained in this file ' .
   'have clonal location coordinates in the format: ' .
   '<clone accession>.<version>:<start>..<end>',

   'The /gene indicates a unique id for a gene, /note="transcript_id=..."' . 
   ' a unique id for a transcript, /protein_id a unique id for a peptide ' .
   'and note="exon_id=..." a unique id for an exon. These ids are ' .
   'maintained wherever possible between versions.',

   'All the exons and transcripts in Ensembl are confirmed by ' .
   'similarity to either protein or cDNA sequences.');


=head2 new

  Arg [1]    : none
  Example    : $seq_dumper = Bio::EnsEMBL::Utils::SeqDumper->new;
  Description: Creates a new SeqDumper 
  Returntype : Bio::EnsEMBL::Utils::SeqDumper
  Exceptions : none
  Caller     : general

=cut

sub new {
  my ($caller, $slice) = @_;

  my $class = ref($caller) || $caller;

  my $feature_types = {'gene'        => 1,
		       'genscan'     => 1,
		       'repeat'      => 1,
		       'similarity'  => 1,
		       'variation'   => 1,
		       'contig'      => 1,
		       'marker'      => 1,
		       'estgene'     => 0,
                       'vegagene'    => 0};

  my $self = bless {'feature_types' => $feature_types}, $class;

  return $self;
}



=head2 enable_feature_type

  Arg [1]    : string $type
  Example    : $seq_dumper->enable_feature_type('similarity');
  Description: Enables the dumping of a specific type of feature
  Returntype : none
  Exceptions : warn if invalid feature type is passed,
               thrown if no feature type is passed
  Caller     : general

=cut

sub enable_feature_type {
  my ($self, $type) = @_;

  $type || throw("type arg is required");

  if(exists($self->{'feature_types'}->{$type})) {
    $self->{'feature_types'}->{$type} = 1;
  } else {
    warning("unknown feature type '$type'\n" .
	  "valid types are: " . join(',', keys %{$self->{'feature_types'}})); 
  }
}



=head2 attach_database

  Arg [1]    : string name
  Arg [2]    : Bio::EnsEMBL::DBSQL::DBAdaptor
  Example    : $seq_dumper->attach_database('estgene', $estgene_db);
  Description: Attaches a database to the seqdumper that can be used to 
               dump data which is external to the ensembl core database.
               Currently this is necessary to dump est genes and vega genes
  Returntype : none
  Exceptions : thrown if incorrect argument is supplied
  Caller     : general

=cut

sub attach_database {
  my ($self, $name, $db) = @_;

  $name || throw("name arg is required");
  unless($db && ref($db) && $db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    throw("db arg must be a Bio::EnsEMBL::DBSQL::DBAdaptor not a [$db]");
  }

  $self->{'attached_dbs'}->{$name} = $db;
}



=head2 get_database

  Arg [1]    : string $name
  Example    : $db = $seq_dumper->get_database('vega');
  Description: Retrieves a database that has been attached to the 
               seqdumper via the attach database call.
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : thrown if incorrect argument is supplied
  Caller     : dump_feature_table

=cut

sub get_database {
  my ($self, $name) = @_;

  $name || throw("name arg is required");
  
  return $self->{'attached_dbs'}->{$name};
}



=head2 remove_database

  Arg [1]    : string $name 
  Example    : $db = $seq_dumper->remove_database('estgene');
  Description: Removes a database that has been attached to the seqdumper
               via the attach database call.  The database that is removed
               is returned (or undef if it did not exist).
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : thrown if incorrect argument is supplied
  Caller     : general

=cut

sub remove_database {
  my ($self, $name) = @_;

  $name || throw("name arg is required");

  if(exists $self->{'attached_dbs'}->{$name}) {
    return delete $self->{'attached_dbs'}->{$name};
  }

  return undef;
}


=head2 disable_feature_type

  Arg [1]    : string $type
  Example    : $seq_dumper->disable_feature_type('genes');
  Description: Disables the dumping of a specific type of feature
  Returntype : none
  Exceptions : warn if an invalid feature type is passed,
               thrown if no feature type is passed
  Caller     : general

=cut

sub disable_feature_type {
  my ($self, $type) = @_;
  
  $type || throw("type arg is required");

  if(exists($self->{'feature_types'}->{$type})) {
    $self->{'feature_types'}->{$type} = 0;
  } else {
    warning("unknown feature type '$type'\n" .
	    "valid types are: " . join(',', keys %{$self->{'feature_types'}}));
  }
}



=head2 is_enabled

  Arg [1]    : string $type 
  Example    : do_something() if($seq_dumper->is_enabled('gene'));
  Description: checks if a specific feature type is enabled
  Returntype : none
  Exceptions : warning if invalid type is passed, 
               thrown if no type is passed 
  Caller     : general

=cut

sub is_enabled {
  my ($self, $type) = @_;

  $type || throw("type arg is required");

  if(exists($self->{'feature_types'}->{$type})) {
    return $self->{'feature_types'}->{$type};
  } else {
    warning("unknown feature type '$type'\n" .
	   "valid types are: " . join(',', keys %{$self->{'feature_types'}}));
  }
}


=head2 dump

  Arg [1]    : Bio::EnsEMBL::Slice slice
               The slice to dump
  Arg [1]    : string $format
               The name of the format to dump
  Arg [2]    : (optional) $outfile
               The name of the file to dump to. If no file is specified STDOUT
               is used
  Example    : $seq_dumper->dump($slice, 'EMBL');
  Description: Dumps a region of a genome specified by the slice argument into
               an outfile of the format $format
  Returntype : none
  Exceptions : thrown if slice or format args are not supplied
  Caller     : general

=cut


sub dump {
  my ($self, $slice, $format, $outfile, $seq) = @_;

  $format || throw("format arg is required");
  $slice  || throw("slice arg is required");

  my $dump_handler = $DUMP_HANDLERS->{uc($format)};

  unless($dump_handler) {
    throw("No dump handler is defined for format $format\n");
  }


  my $FH = IO::File->new;;
  if($outfile) {
    $FH->open(">>$outfile") or throw("Could not open file $outfile");
  } else {
    $FH = \*STDOUT;
    #mod_perl did not like the following
    #$FH->fdopen(fileno(STDOUT), "w") 
    #  or throw("Could not open currently selected output filehandle " .
    #		      "for writing");
  }
  
  &$dump_handler($self, $slice, $FH, $seq);

  $FH->close if ($outfile); #close if we were writing to a file
}



=head2 dump_embl

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : IO::File $FH
  Arg [3]    : optional sequence string
  Example    : $seq_dumper->dump_embl($slice, $FH);
  Description: Dumps an EMBL flat file to an open file handle
  Returntype : none
  Exceptions : none
  Caller     : dump

=cut

sub dump_embl {
  my $self = shift;
  my $slice = shift;
  my $FH   = shift;
  my $SEQ = shift;

  my $len = $slice->length;

  my $version;
  my $acc;

  my $cs = $slice->coord_system();
  my $name_str = $cs->name() . ' ' . $slice->seq_region_name();
  $name_str .= ' ' . $cs->version if($cs->version);

  my $start = $slice->start;
  my $end   = $slice->end;

  #determine if this slice is the entire seq region
  #if it is then we just use the name as the id
  my $slice_adaptor = $slice->adaptor();
  my $full_slice =
    $slice->adaptor->fetch_by_region($cs->name,
                                    $slice->seq_region_name,
                                    undef,undef,undef,
                                    $cs->version);


  my $entry_name = $slice->seq_region_name();



  if($full_slice->name eq $slice->name) {
    $name_str .= ' full sequence';
    $acc = $slice->seq_region_name();
    my @acc_ver = split(/\./, $acc);
    if(@acc_ver == 2) {
      $acc = $acc_ver[0];
      $version = $acc_ver[0] . '.'. $acc_ver[1];
    } elsif(@acc_ver == 1 && $cs->version()) {
      $version = $acc . '.'. $cs->version();
    } else {
      $version = $acc;
    }
  } else {
    $name_str .= ' partial sequence';
    $acc = $slice->name();
    $version = $acc;
  }

  $acc = $slice->name();



  #line breaks are allowed near the end of the line on ' ', "\t", "\n", ',' 
  $: = (" \t\n-,");

  #############
  # dump header
  #############

  my $EMBL_HEADER = 
'@<   ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<~
';

  #ID and moltype
  # HTG = High Throughput Genome division, probably most suitable
  #       and it would be hard to come up with another appropriate division
  #       that worked for all organisms (e.g. plants are in PLN but human is
  #       in HUM).
  my $VALUE = "$entry_name    standard; DNA; HTG; $len BP.";
  $self->write($FH, $EMBL_HEADER, 'ID', $VALUE);  
  $self->print( $FH, "XX\n" );

  #Accession
  $self->write($FH, $EMBL_HEADER, 'AC', $acc);
  $self->print( $FH, "XX\n" );

  #Version
  $self->write($FH, $EMBL_HEADER, 'SV', $version);
  $self->print( $FH, "XX\n" );

  #Date
  $self->write($FH, $EMBL_HEADER, 'DT', $self->_date_string);
  $self->print( $FH, "XX\n" );

  my $species   = $slice->adaptor->db->get_MetaContainer->get_Species();

  #Description
  $self->write($FH, $EMBL_HEADER, 'DE', $species->binomial .
               " $name_str $start..$end annotated by Ensembl");
  $self->print( $FH, "XX\n" );

  #key words
  $self->write($FH, $EMBL_HEADER, 'KW', '.');
  $self->print( $FH, "XX\n" );

  #Species
  my $species_name = $species->binomial();
  if(my $cn = $species->common_name()) {
    $species_name .= " ($cn)";
  }

  $self->write($FH, $EMBL_HEADER, 'OS', $species_name);

  #Classification
  my @cls = $species->classification;
  shift @cls; #shift off species name
  $self->write($FH, $EMBL_HEADER, 'OC', join('; ', reverse(@cls)) . '.');
  $self->print( $FH, "XX\n" );
  
  #References (we are not dumping refereneces)

  #Database References (we are not dumping these)

  #comments
  foreach my $comment (@COMMENTS) {
    $self->write($FH, $EMBL_HEADER, 'CC', $comment);
    $self->print( $FH, "XX\n" );
  }

  ####################
  #DUMP FEATURE TABLE
  ####################
  $self->print( $FH, "FH   Key             Location/Qualifiers\n" );

  my $FEATURE_TABLE = 
'FT   ^<<<<<<<<<<<<<<<^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<~
';
  $self->_dump_feature_table($slice, $FH, $FEATURE_TABLE);  

  #write an XX after the feature tables
  $self->print( $FH, "XX\n" );

  ###################
  #DUMP SEQUENCE
  ###################

  if(!defined($SEQ)){
    $SEQ = $slice->seq();
  }
#  my $SEQ     = $slice->seq();
  my $length  = length($SEQ);
  my $a_count = $SEQ =~ tr/aA/aA/;
  my $c_count = $SEQ =~ tr/cC/cC/;
  my $t_count = $SEQ =~ tr/tT/tT/;
  my $g_count = $SEQ =~ tr/gG/gG/;
  my $other_count = $length - $a_count - $c_count - $t_count - $g_count;

  my $value = "Sequence $length BP; $a_count A; $c_count C; " .
    "$g_count G; $t_count T; $other_count other;";
  $self->write($FH, $EMBL_HEADER, 'SQ', $value);

  $self->write_embl_seq($FH, \$SEQ);


  $self->print( $FH, "//\n" );

  # Set formatting back to normal
  $: = " \n-";
}




=head2 dump_genbank

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : IO::File $FH
  Example    : $seq_dumper->dump_genbank($slice, $FH);
  Description: Dumps a GENBANK flat file to an open file handle
  Returntype : none
  Exceptions : none
  Caller     : dump

=cut

sub dump_genbank {
  my ($self, $slice, $FH, $SEQ) = @_;

  #line breaks are allowed near the end of the line on ' ', "\t", "\n", ',' 
  $: = " \t\n-,";

  my $GENBANK_HEADER = 
'^<<<<<<<<<  ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<~~
';

  my $GENBANK_SUBHEADER =
'  ^<<<<<<<  ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<~~
';

  my $GENBANK_FT =
'     ^<<<<<<<<<<<<<< ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<~~
';

  my $version;
  my $acc;

  my $cs = $slice->coord_system();

  my $name_str = $cs->name() . ' ' . $slice->seq_region_name();

  $name_str .= ' ' . $cs->version if($cs->version);

  #determine if this slice is the entire seq region
  #if it is then we just use the name as the id
  my $slice_adaptor = $slice->adaptor();
  my $full_slice =
    $slice->adaptor->fetch_by_region($cs->name,
                                    $slice->seq_region_name,
                                    undef,undef,undef,
                                    $cs->version);


  my $entry_name = $slice->seq_region_name();

  if($full_slice->name eq $slice->name) {
    $name_str .= ' full sequence';
    $acc = $slice->seq_region_name();
    my @acc_ver = split(/\./, $acc);
    if(@acc_ver == 2) {
      $acc = $acc_ver[0];
      $version = $acc_ver[0] . $acc_ver[1];
    } elsif(@acc_ver == 1 && $cs->version()) {
      $version = $acc . $cs->version();
    } else {
      $version = $acc;
    }
  } else {
    $name_str .= ' partial sequence';
    $acc = $slice->name();
    $version = $acc;
  }

  $acc = $slice->name();     # to keep format consistent for all

  my $length = $slice->length;
  my $start = $slice->start();
  my $end   = $slice->end();

  my $date = $self->_date_string;

  my $species = $slice->adaptor->db->get_MetaContainer->get_Species;


  #LOCUS
  my $tag   = 'LOCUS';
  my $value = "$entry_name $length bp DNA HTG $date";
  $self->write($FH, $GENBANK_HEADER, $tag, $value);

  #DEFINITION
  $tag   = "DEFINITION";
  $value = $species->binomial . 
    " $name_str $start..$end reannotated via EnsEMBL";
  $self->write($FH, $GENBANK_HEADER, $tag, $value);

  #ACCESSION
  $self->write($FH, $GENBANK_HEADER, 'ACCESSION', $acc);

  #VERSION
  $self->write($FH, $GENBANK_HEADER, 'VERSION', $version);

  # KEYWORDS
  $self->write($FH, $GENBANK_HEADER, 'KEYWORDS', '.');

  # SOURCE
  $self->write($FH, $GENBANK_HEADER, 'SOURCE', $species->common_name());

  #organism
  my @cls = $species->classification();
  shift @cls;
  $self->write($FH, $GENBANK_SUBHEADER, 'ORGANISM', $species->binomial);
  $self->write($FH, $GENBANK_SUBHEADER, '', join('; ', reverse @cls) . ".");

  #refereneces

  #comments
  foreach my $comment (@COMMENTS) {
    $self->write($FH, $GENBANK_HEADER, 'COMMENT', $comment);
  }

  ####################
  # DUMP FEATURE TABLE
  ####################
  $self->print( $FH, "FEATURES             Location/Qualifiers\n" );
  $self->_dump_feature_table($slice, $FH, $GENBANK_FT);

  ####################
  # DUMP SEQUENCE
  ####################

  if(!defined($SEQ)){
    $SEQ = $slice->seq();
  }
#  my $SEQ       = $slice->seq();
  my $a_count = $SEQ =~ tr/aA/aA/;
  my $c_count = $SEQ =~ tr/cC/cC/;
  my $t_count = $SEQ =~ tr/tT/tT/;
  my $g_count = $SEQ =~ tr/gG/gG/;
  my $bp_length = length($SEQ);
  my $other_count = $bp_length - $a_count - $c_count - $t_count - $g_count;

  $tag   = 'BASE COUNT';
  $value = "$a_count a $c_count c $g_count g $t_count t";
  $value .= " $other_count n" if($other_count);
  $self->write($FH, $GENBANK_HEADER, $tag, $value);
  $self->print( $FH, "ORIGIN\n" );

  $self->write_genbank_seq($FH, \$SEQ);

  $self->print( $FH, "//\n" );

  # Set formatting back to normal
  $: = " \n-";
}



=head2 _dump_feature_table

  Arg [1]    : Bio::EnsEMBL::Slice slice
  Example    : none
  Description: Helper method used to dump feature tables used in EMBL, FASTA,
               GENBANK.  Assumes formating of file handle has been setup
               already to use $FEAT and $VALUE values.
  Returntype : none
  Exceptions : none
  Caller     : internal

=cut

sub _dump_feature_table {
  my $self   = shift;
  my $slice  = shift;
  my $FH     = shift;
  my $FORMAT = shift;

  #use only the core database to dump features (except for bloody snps)
  my $lite = $slice->adaptor->db->remove_db_adaptor('lite');

  my $meta = $slice->adaptor->db->get_MetaContainer;
  my $species = $meta->get_Species;

  #lump file handle and format string together for simpler method calls
  my @ff = ($FH, $FORMAT);
  my $value;

  #source
  my $classification = join(', ', $species->classification);
  $self->write(@ff,'source', "1.." . $slice->length());
  $self->write(@ff,''      , '/organism="'.$species->binomial . '"');
  $self->write(@ff,''      , '/db_xref="taxon:'.$meta->get_taxonomy_id().'"');

  #
  # Transcripts & Genes
  #
  my @gene_slices;
  if($self->is_enabled('gene')) {
    push @gene_slices, $slice;
  }

  # Retrieve slices of other database where we need to pull genes from

  my $gene_dbs = {'vegagene' => 'vega',
                  'estgene'  => 'estgene'};

  foreach my $gene_type (keys %$gene_dbs) {
    if($self->is_enabled($gene_type)) {
      my $db = $self->get_database($gene_dbs->{$gene_type});
      if($db) {
        my $sa = $db->get_SliceAdaptor();
        push @gene_slices, $sa->fetch_by_name($slice->name());
      } else {
        warning("A [". $gene_dbs->{$gene_type} ."] database must be " .
                "attached to this SeqDumper\n(via a call to " .
                "attach_database) to retrieve genes of type [$gene_type]");
      }
    }
  }

  foreach my $gene_slice (@gene_slices) {
    foreach my $gene (@{$gene_slice->get_all_Genes(undef,undef, 1)}) {
      $value = $self->features2location( [$gene] );
      $self->write( @ff, 'gene', $value );
      $self->write( @ff, "", '/gene='.$gene->stable_id() );


      if(defined($gene->display_xref)){
	$self->write( @ff, "",'/locus_tag="'.$gene->display_xref->display_id.'"');
      }
      my $desc = $gene->description;
      if(defined($desc) and $desc ne ""){
	$desc =~ s/\"//; 
	$self->write( @ff, "", '/note="'.$gene->description.'"');
      }



      foreach my $transcript (@{$gene->get_all_Transcripts}) {
        my $translation = $transcript->translation;

        # normal transcripts get dumped differently than pseudogenes
        if($translation) {
          #normal transcript
          $value = $self->features2location($transcript->get_all_Exons);
          $self->write(@ff, 'mRNA', $value);
          $self->write(@ff,''   , '/gene="'.$gene->stable_id().'"');
          $self->write(@ff,''
                       ,'/note="transcript_id='.$transcript->stable_id().'"');

          # ...and a CDS section
          $value = 
            $self->features2location($transcript->get_all_translateable_Exons);
          $self->write(@ff,'CDS', $value);
          $self->write(@ff,''   , '/gene="'.$gene->stable_id().'"'); 
          $self->write(@ff,'', '/protein_id="'.$translation->stable_id().'"');
          $self->write(@ff,''
                       ,'/note="transcript_id='.$transcript->stable_id().'"');

          foreach my $dbl (@{$transcript->get_all_DBLinks}) {
            $value = '/db_xref="'.$dbl->dbname().':'.$dbl->display_id().'"';
            $self->write(@ff, '', $value);
          }

          $value = '/translation="'.$transcript->translate()->seq().'"';
          $self->write(@ff, '', $value);
        } else {
          #pseudogene
          $value = $self->features2location($transcript->get_all_Exons);
          $self->write(@ff, 'misc_RNA', $value);
          $self->write(@ff,''   , '/gene="'.$gene->stable_id().'"');
          foreach my $dbl (@{$transcript->get_all_DBLinks}) {
            $value = '/db_xref="'.$dbl->dbname().':'.$dbl->primary_id().'"';
            $self->write(@ff, '', $value);
          }
          $self->write(@ff,''   , '/note="'.$transcript->biotype().'"');
          $self->write(@ff,''
                       ,'/note="transcript_id='.$transcript->stable_id().'"');
        }
      }
    }

    # exons
    foreach my $gene (@{$gene_slice->get_all_Genes(undef,undef,1)}) {
      foreach my $exon (@{$gene->get_all_Exons}) {
        $self->write(@ff,'exon', $self->features2location([$exon]));
        $self->write(@ff,''    , '/note="exon_id='.$exon->stable_id().'"');
      }
    }
  }

  #
  # genscans
  #
  if($self->is_enabled('genscan')) {
    my @genscan_exons;
    foreach my $transcript(@{$slice->get_all_PredictionTranscripts(undef,1)}) {
      my $exons = $transcript->get_all_Exons();
      push @genscan_exons, @$exons;
      $self->write(@ff, 'mRNA', $self->features2location($exons));
      $self->write(@ff, '', '/product="'.$transcript->translate()->seq().'"');
      $self->write(@ff, '', '/note="identifier='.$transcript->stable_id.'"');
      $self->write(@ff, '', '/note="Derived by automated computational' .
		   ' analysis using gene prediction method:' .
		   $transcript->analysis->logic_name . '"');
    }
  }

  #
  # snps
  #
  if($self->is_enabled('variation') && $slice->can('get_all_VariationFeatures')) {
#    $slice->adaptor->db->add_db_adaptor('lite', $lite) if $lite;

    foreach my $snp (@{$slice->get_all_VariationFeatures}) {
      my $ss = $snp->start;
      my $se = $snp->end;
      #skip snps that hang off edge of slice
      next if($ss < 1 || $se > $slice->length); 

      $self->write(@ff, 'variation', "$ss..$se");
      $self->write(@ff, ''         , '/replace="'.$snp->allele_string.'"'); 
      #$self->write(@ff, ''         , '/evidence="'.$snp->status.'"'); 
      my $rs_id = $snp->variation_name();
      my $db = $snp->source();
#      foreach my $link ($snp->each_DBLink) {
#        my $id = $link->primary_id;
#        my $db = $link->database;
        $self->write(@ff, '', "/db_xref=\"$db:$rs_id\""); 
#      }
    }

#    $slice->adaptor->db->remove_db_adaptor('lite') if $lite;
  }

  #
  # similarity features
  #
  if($self->is_enabled('similarity')) {
    foreach my $sim (@{$slice->get_all_SimilarityFeatures}) {
      $self->write(@ff, 'misc_feature', $self->features2location([$sim]));
      $self->write(@ff, ''       , '/note="match: '.$sim->hseqname.
		  ' : '.$sim->hstart.'..'.$sim->hend.'('.$sim->hstrand.')"');
    }
  }

  #
  # repeats
  #
  if($self->is_enabled('repeat')) {
    my $rfs = $slice->get_all_RepeatFeatures();

    foreach my $repeat (@$rfs) {
      $self->write(@ff, 'repeat_region', $self->features2location([$repeat]));
      $self->write(@ff, ''    , '/note="' . $repeat->repeat_consensus->name.
		   ' repeat: matches ' . $repeat->hstart.'..'.$repeat->hend .
		   '('.$repeat->hstrand.') of consensus"');
    }

  }

  #
  # markers
  #
  if($self->is_enabled('marker') && $slice->can('get_all_MarkerFeatures')) {
    foreach my $mf (@{$slice->get_all_MarkerFeatures}) {
      $self->write(@ff, 'STS', $self->features2location([$mf]));
      if($mf->marker->display_MarkerSynonym) {
        $self->write(@ff, ''   , '/standard_name="' .
                     $mf->marker->display_MarkerSynonym->name . '"');
      }


      #grep out synonyms without a source
      my @synonyms = @{$mf->marker->get_all_MarkerSynonyms};
      @synonyms = grep {$_->source } @synonyms;
      foreach my $synonym (@synonyms) {
        $self->write(@ff, '', '/db_xref="'.$synonym->source.
                     ':'.$synonym->name.'"');
      }
      $self->write(@ff, '', '/note="map_weight='.$mf->map_weight.'"');
    }
  }

  #
  # contigs
  #
  if($self->is_enabled('contig')) {
    foreach my $segment (@{$slice->project('seqlevel')}) {
      my ($start, $end, $slice) = @$segment;
      $self->write(@ff, 'misc_feature',
                   $start .'..'. $end);
      $self->write(@ff, '', '/note="contig '.$slice->seq_region_name .
		   ' ' . $slice->start . '..' . $slice->end .
		    '(' . $slice->strand . ')"');
    }
  }

  $slice->adaptor->db->add_db_adaptor('lite', $lite) if $lite;

}



=head2 dump_fasta

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : IO::File $FH
  Example    : $seq_dumper->dump_fasta($slice, $FH);
  Description: Dumps an FASTA flat file to an open file handle
  Returntype : none
  Exceptions : none
  Caller     : dump

=cut

sub dump_fasta {
  my $self = shift;
  my $slice = shift;
  my $FH   = shift;

  my $id       = $slice->seq_region_name;
  my $seqtype  = 'dna';
  my $idtype   = $slice->coord_system->name;
  my $location = $slice->name;
  my $start = 1;
  my $end = $slice->length();

  my $header = ">$id $seqtype:$idtype $location\n";
  $self->print( $FH, $header );

  #set the formatting to FASTA
  my $FORMAT = '^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
';

  #chunk the sequence in 60kb chunks to use less memory
  my $cur = $start;
  while($cur < $end) {
    my $to = $cur + 59_999;
    $to = $end if($to > $end); 
    my $seq = $slice->subseq($cur, $to);
    $cur = $to + 1;
    $self->write($FH, $FORMAT, $seq);
  }
}



=head2 features2location

  Arg [1]    : listref of Bio::EnsEMBL::SeqFeatures
  Example    : $location = $self->features2location(\@features);
  Description: Constructs an EMBL location string from a list of features
  Returntype : string
  Exceptions : none
  Caller     : internal

=cut

sub features2location {
  my $self = shift;
  my $features = shift;

  my @join = ();

  foreach my $f (@$features) {
    my $slice = $f->slice;
    my $start = $f->start();
    my $end   = $f->end();
    my $strand = $f->strand();

    if($start >= 1 && $end <= $slice->length) {
      #this feature in on a slice and doesn't lie outside the boundary
	
      if($strand == 1) {
        push @join, "$start..$end";
      } else {
        push @join, "complement($start..$end)";
      }
    } else {
      my @fs = ();
      #this feature is outside the boundary of the dump,
      # yet implemented and 'seqlevel' is guaranteed to be 1step
      my $projection = $f->project('seqlevel');
      foreach my $segment (@$projection) {
        my $slice = $segment->[2];
        my $slc_start = $slice->start();
        my $slc_end   = $slice->end();
        my $seq_reg   = $slice->seq_region_name();
        if($slice->strand == 1) {
          push @join, "$seq_reg:$slc_start..$slc_end";
        } else {
          push @join, "complement($seq_reg:$slc_start..$slc_end)";
        }
      }
    }
  }

  my $out = join ',', @join;

  if(scalar @join > 1) {
    $out = "join($out)";
  }

  return $out;
}


sub _date_string {
  my $self = shift;

  my ($sec, $min, $hour, $mday,$mon, $year) = localtime(time());

  my $month = ('JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL',
	       'AUG', 'SEP', 'OCT', 'NOV', 'DEC')[$mon];
  $year += 1900;
  
  return "$mday-$month-$year";
}


sub write {
  my ($self, $FH, $FORMAT, @values) = @_;
  
  #while the last value still contains something
  while($values[-1] ne '') {
    formline($FORMAT, @values);
    $self->print( $FH, $^A );
    $^A = '';
  }
}

sub write_genbank_seq {
  my $self = shift;
  my $FH  = shift;
  my $seq = shift;
  my $base_total = shift;

  $base_total ||= 0;

  my $GENBANK_SEQ = 
'@>>>>>>>> ^<<<<<<<<< ^<<<<<<<<< ^<<<<<<<<< ^<<<<<<<<< ^<<<<<<<<< ^<<<<<<<<<~
';

  my $total = -59 + $base_total;
  #keep track of total and print lines of 60 bases with spaces every 10bp
  while($$seq) {
    $total += 60; 
    formline($GENBANK_SEQ,$total, $$seq, $$seq, $$seq, $$seq, $$seq, $$seq);
    $self->print( $FH, $^A );
    $^A = '';
  }
}

sub write_embl_seq {
  my $self = shift;
  my $FH   = shift;
  my $seq  = shift;
  my $base_total = shift;

  $base_total ||= 0;

  my $EMBL_SEQ = 
'     ^<<<<<<<<< ^<<<<<<<<< ^<<<<<<<<< ^<<<<<<<<< ^<<<<<<<<< ^<<<<<<<<<@>>>>>>>>>~
';
  #keep track of total and print lines of 60 bases with spaces every 10bp
  my $length = length($$seq);
  my $total = $length - $base_total;
  while($$seq) {
    $total -= 60;
    $total = 0 if($total < 0);
    formline($EMBL_SEQ, 
	     $$seq, $$seq, $$seq, $$seq, $$seq, $$seq, 
	     $length - $total);
    $self->print( $FH, $^A );
    $^A = '';
  }
}

sub print {
  my( $self, $FH, $string ) = @_;
  print $FH $string;
}

1;   
