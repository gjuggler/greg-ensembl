package Bio::Greg::Gorilla::MapGorillaSubs;

use strict;
use Time::HiRes qw(sleep);

use base (
  'Bio::Greg::Hive::Process'
);

my $TU = 'Bio::EnsEMBL::Compara::TreeUtils';
my $CU = 'Bio::EnsEMBL::Compara::ComparaUtils';
my $AU = 'Bio::EnsEMBL::Compara::AlignUtils';

sub param_defaults {
  return {
    job_role => 'fan_jobs'
  };
}

sub fetch_input {
  my ($self) = @_;

  $self->create_table_from_params( $self->compara_dba, 'species_subs',
    $self->get_gene_stats_def );

}

sub fan_job {
  my $self = shift;
  my $line_lo = shift;
  my $line_hi = shift;

  my $output_params = {
    job_role => 'map_subs',
    line_lo => $line_lo,
    line_hi => $line_hi
  };
 
  my $output_id = $self->string_to_hash( $self->input_id ); 
  $output_id = $self->replace($output_id,$output_params);
  my ($output_job_id) = @{ $self->dataflow_output_id($output_id, 99)};
  print " --> Flowed out to job $output_job_id\n";

  $self->hash_print($output_id);
  sleep(0.1);

}

sub run {
  my $self = shift;

  my $scratch           = Bio::Greg::EslrUtils->scratchDirectory;
  my $subs_file = "$scratch/gorilla/Gorilla.NS_St.csv";

  my $i = 0;
  my $bin = 0;
  my $bin_size = 5;

  open(IN,$subs_file) or die $!;
  <IN>;
  while(<IN>){
    chomp;
    my $line = $_;
    print $line."\n";
    my @tokens = split(/\t/,$line);

    my $chr = $tokens[0];
    my $start = $tokens[1];
    my $ts_id = $tokens[2];
    my $effect = $tokens[3];
    my $aa = $tokens[4];

    if ($self->param('job_role') eq 'map_subs') {
      if ($i >= $self->param('line_lo') && $i < $self->param('line_hi')) {
        $self->map_gorilla_subs($chr,$start,$ts_id,$effect,$aa);
      }
    } elsif ($bin >= $bin_size) {
      $self->fan_job($i-$bin_size,$i);
      $bin = 0;
    }

    $i++;
    $bin++;
  }
  close(IN);

  if ($bin > 0) {
    $self->fan_job($i-$bin,$i);
  }
}

sub map_gorilla_subs {
  my $self = shift;
  my ($sub_chr,$sub_start,$ts_id,$sub_effect,$sub_aa) = @_;

  $self->load_registry;


  my $ts_adaptor = Bio::EnsEMBL::Registry->get_adaptor( "Gorilla", "Core", "Transcript" );
  my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor( "Gorilla", "Core", "Gene" );
  my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor( "Gorilla", "Core", "Slice" );

  my $compara_dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
    -url => 'mysql://ensadmin:ensembl@ensdb-archive:5304/ensembl_compara_58' );
#  my $compara_dba = $self->compara_dba;
  my $as_a       = $compara_dba->get_AlignSliceAdaptor;
  my $mlss_a       = $compara_dba->get_MethodLinkSpeciesSetAdaptor;
  my $mlss = $mlss_a->fetch_by_method_link_type_species_set_name( 'EPO', 'primates' );

  my $slice = $slice_adaptor->fetch_by_region(undef,$sub_chr,$sub_start-3,$sub_start+3);

  print "$ts_id\n";
  print $slice->seq."\n";
  print "  ".$sub_aa."\n";

  my $ts_index = 0;
  my $sub_ts_index = 0;
  my $ts_seq = '';

  my $ts = $ts_adaptor->fetch_by_stable_id($ts_id);
  my $gene = $gene_adaptor->fetch_by_transcript_stable_id($ts_id);

  my $pep = $ts->translate;
  my $pep_str = $pep->seq;
  my $pep_length = length($pep_str);

  print "Tx length: $pep_length aa's\n";

  $self->param('name',$gene->external_name);
  $self->param('description',$gene->description);

  print "Strand: ".$ts->strand."\n";
  my @exons = @{ $ts->get_all_translateable_Exons };
  foreach my $exon (@exons) {
    my $start     = $exon->coding_region_start($ts);
    my $end       = $exon->coding_region_end($ts);
    my $strand    = $exon->strand;
    my $frame     = $exon->frame;
    my $phase     = $exon->phase;
    my $end_phase = $exon->end_phase;
    
    my $off_phase_start;
    if ($ts_index == 0 && $exon->phase > 0) {
      print "OFF_PHASE_START!!!\n";
      $off_phase_start = 1;
      print $exon->stable_id." ".$exon->phase."\n";
      if ($exon->strand == -1) {
        $end = $end + $exon->phase;
      } else {
        $start = $start - $exon->phase;
      }
    }
    
    my $slice = $exon->slice->sub_Slice( $start, $end, $strand );
    #print $slice->seq."\n";
    
    for (my $i=1; $i <= $slice->length; $i++) {
      $ts_index++;
      
      my $sub_slice = $slice->sub_Slice($i,$i);
      my $nuc = $sub_slice->seq;
      $ts_seq .= $nuc;
      
      # Most important bit: Flag the transcript CDNA index of the mutation.
      if ($sub_start == $sub_slice->start) {
        $sub_ts_index = $ts_index;
        print "TS index: $sub_ts_index\n";
      }
      
      if ($ts_index % 3 == 0) {
        my $codon = substr($ts_seq,$ts_index-3,3);
        
        my $aa = new Bio::LocatableSeq(-seq=>$codon)->translate->seq;
        my $aa_index = $ts_index / 3;
        
        my $d = $sub_start - $sub_slice->start;
        if ($d < 2 && $d > -2) {
          print "$codon $aa\n";
          my $mutation_slice = $exon->slice->sub_Slice($sub_start,$sub_start,$strand);
          my $align_slice = $as_a->fetch_by_Slice_MethodLinkSpeciesSet($mutation_slice,$mlss,0);          
          my $sa = $align_slice->get_SimpleAlign;
          $sa = $self->remove_ancestrals($sa);
          $self->pretty_print($sa);          
        } else {
#          print "$codon $aa\n";
        }
      }
    }
  }

  print "Getting aln...\n";
  my ($cdna_aln,$aa_aln,$extra_info) = Bio::EnsEMBL::Compara::ComparaUtils->genomic_aln_for_transcript($compara_dba,$ts,{quality_threshold => 0,species_set => 'primates'});
#  $self->pretty_print($cdna_aln,{full => 1,width=>200});

  # Find the codon which this index is within.
  my $codon_index = $sub_ts_index - (($sub_ts_index-1) % 3);
  my $sub_codon_position = ($sub_ts_index - 1) % 3;

  my $sub_aln = $cdna_aln->slice($sub_ts_index,$sub_ts_index,1);
  my $codon_aln = $cdna_aln->slice($codon_index,$codon_index+2,1);
  my $aa_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($codon_aln);

  $codon_aln = $self->remove_ancestrals($codon_aln);
  $aa_aln = $self->remove_ancestrals($aa_aln);

  my $dna_string = Bio::EnsEMBL::Compara::AlignUtils->get_column_string($sub_aln,1);
  my $aa_string = Bio::EnsEMBL::Compara::AlignUtils->get_column_string($aa_aln,1);
  print "$aa_string $dna_string\n";

  # Get the majority amino acid allele from the ancestral alignments.
  my $ancestral_aa = $self->get_majority_allele($aa_string,$sub_aa);

  my @a_b = split("/",$sub_aa);
  my $ancestral_allele = 0;
  $ancestral_allele = 1 if ($a_b[0] eq $ancestral_aa);
  $ancestral_allele = 2 if ($a_b[1] eq $ancestral_aa);

  my $params = $self->params;
  my $data = {
    data_id => 12345,
    chr => $sub_chr,
    start => $sub_start,
    ts_id => $ts_id,
    aa_length => $pep_length,

    effect => $sub_effect,
    aa => $sub_aa,
    
    sub_pos => $sub_ts_index,
    sub_codon_pos => $sub_codon_position,
    aln_aa_string => $aa_string,
    aln_dna_string => $dna_string,
    
    ancestral_aa => $ancestral_aa,
    ancestral_allele => $ancestral_allele
  };
  $params = $self->replace($params,$data);
  $self->store_params_in_table( $self->dbc, 'species_subs', $params );

  my $aln_id = "${ts_id}.${sub_chr}.${sub_start}";

  my $file = $self->save_aln(
    $codon_aln, {
      id        => $aln_id,
      filename  => $aln_id,
      subfolder => 'species_subs'
    }
    );
}

sub get_majority_allele {
  my $self = shift;
  my $column_string = shift;
  my $orig_allele_string = shift;
   
  my @orig_alleles = split("/",$orig_allele_string);

  $column_string =~ s/-//g; # No gaps.

  return '' if (length($column_string) == 0);

  my @chars = split('',$column_string);

  my %histogram;
  $histogram{$_}++ for @chars;

  my @sorted_chars = (sort { $histogram{$b} <=> $histogram{$a} } keys %histogram);

  # We have a list of amino acids sorted by popularity in the alignment column.
  # Go through and return the most popular a.a. that matches either one of the
  # original alleles.
  foreach my $char (@sorted_chars) {
    if ($char ne '*' && $orig_allele_string =~ m/$char/) {
      return $char;
    }
  }
  return $sorted_chars[0];
}

sub remove_ancestrals {
  my $self = shift;
  my $aln = shift;

  my @seqs;

  foreach my $seq ($aln->each_seq) {
    next if ($seq->id =~ m/ancestral/i);
    push @seqs, $seq;
  }

  my $new_aln = new $aln;
  foreach my $seq (@seqs) {
    $new_aln->add_seq($seq);
  }
  return $new_aln;
}

sub get_gene_stats_def {
  my $self = shift;

  my $params = {
    ts_id => 'char32',
    chr => 'string',
    start => 'int',
    aa_length => 'int',

    name => 'string',
    description => 'string',

    effect => 'string',
    aa => 'char16',

    sub_pos => 'int',
    sub_codon_pos => 'int',
    aln_aa_string => 'char32',
    aln_dna_string => 'char32',

    ancestral_aa => 'char8',
    ancestral_allele => 'int'
  };

}


sub write_output {
  my $self = shift;
}
