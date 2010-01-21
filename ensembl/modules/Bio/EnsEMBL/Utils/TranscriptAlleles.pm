=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

TranscriptAlleles - A utility class used to obtain information about the
relationships between a transcript and Alleles

=head1 SYNOPSIS

  use Bio::EnsEMBL::Utils::TranscriptAlleles;

  # get the peptide variations caused by a set of Alleles

  %variations = %{
    Bio::EnsEMBL::Utils::TranscriptAlleles::get_all_peptide_variations(
      $transcript, $alleles ) };

=head1 DESCRIPTION

This is a utility class which can be used to find consequence type of an
AlleleFeature in a transcript, and to determine the amino acid changes
caused by the AlleleFeature in the Transcript


=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::TranscriptAlleles;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::ConsequenceType;
use vars qw(@ISA @EXPORT_OK);

use Data::Dumper;

@ISA = qw(Exporter);

@EXPORT_OK = qw(&get_all_ConsequenceType &type_variation);


=head2 get_all_ConsequenceType

  Arg [1]    : $transcript the transcript to obtain the peptide variations for
  Arg [2]    : $alleles listref of AlleleFeatures
  Example    : $consequence_types = get_all_ConsequenceType($transcript, \@alleles);
               foreach my $ct (@{$consequence_types}){
                  print "Allele : ", $ct->allele_string, " has a consequence type of :",$ct->type;
                  print " and is affecting the transcript with ",@{$ct->aa_alleles}, "in position ", 
		              $ct->aa_start,"-", $ct->aa_end if (defined $ct->aa_alleles);
		  print "\n";
	      }
  Description: Takes a list of AlleleFeatures and a Transcritpt, and return a list
                     of ConsequenceType of the alleles in the given Transcript
  Returntype : listref of Bio::EnsEMBL::Variation::ConsequenceType
  Exceptions : none
  Caller     : general 

=cut

sub get_all_ConsequenceType {
  my $transcript = shift;
  my $alleles = shift;

  if(!ref($transcript) || !$transcript->isa('Bio::EnsEMBL::Transcript')) {
    throw('Bio::EnsEMBL::Transcript argument is required.');
  }

  if(!ref($alleles) || (ref($alleles) ne 'ARRAY')) {
    throw('Reference to a list of Bio::EnsEMBL::Variation::AlleleFeature objects is required');
  }


  my @alleles_ordered = sort { $a->start <=> $b->start} @$alleles; #sort the alleles by the genomic position
  my @same_codon; #contains up to 3 allele features, that are in the same codon, but each position can contain more than 1 allele
  my @out; #array containing the consequence types of the alleles in the transcript
  foreach my $allele (@alleles_ordered) {
#  foreach my $allele (@{$alleles}) {
    #get consequence type of the AlleleFeature
    # my $new_allele = $allele->transform('chromosome');
    #my $consequence_type = Bio::EnsEMBL::Variation::ConsequenceType->new($transcript->dbID(),'',$allele->start,$allele->end,$allele->strand,[$allele->allele_string]);
    ### REAL HACK BY js5 because something is borked in TranscriptMapper
    ### This relies on the Allele being of the form i.e. a SNP! [ACGT-](/[ACGT-])+
    ### The rest don't work anyway until we have a AlignStrainSlice
    ### MUST BE SORTED....

    #we have to consider het alleles
    my $allele_string;
    if ($allele->allele_string =~ /[\|\\\/]/){
      my @alleles = split /[\|\\\/]/,$allele->allele_string;
      if ($alleles[0] ne $allele->ref_allele_string){
	$allele_string = $alleles[0];
	  }
      else{
	$allele_string = $alleles[1];
      }
    }
    else{
      $allele_string = $allele->allele_string;  
    }    
    my $opposite_strand = 0; #to indicate wether transcript and allele and in different strands
    my $transcript_allele = $allele_string;
    if( $transcript->strand != $allele->strand ) {
      $transcript_allele =~tr/ACGT/TGCA/;
      $opposite_strand = 1;
    }

    my $consequence_type = Bio::EnsEMBL::Variation::ConsequenceType->new($transcript->dbID(),'',$allele->start, $allele->end, $transcript->strand, [$transcript_allele]);
    #calculate the consequence type of the Allele if different from the reference Allele
    #if (($opposite_strand && $allele->ref_allele_string eq $allele_string) || (!$opposite_strand && $allele->ref_allele_string eq $allele_string)){      	#same allele as reference, there is no consequence, called SARA
    if ($allele->ref_allele_string eq $allele_string) {      	#same allele as reference, there is no consequence, called SARA
      	#same allele as reference, there is no consequence, called SARA
	#we have to calculate if there are more than 2 in the same codon
	empty_codon(\@out,\@same_codon);	
	$consequence_type->type('SARA');
	push @out, $consequence_type;
	next;
    }
	
    my $ref_consequences = type_variation($transcript,"",$consequence_type);
    if ($allele->start != $allele->end){
	empty_codon(\@out,\@same_codon);
	#do not calculate for indels effects of 2 or more in same codon
	push @out, @{$ref_consequences};
	next;
    }

    my $new_consequence = shift @{$ref_consequences};
    if (! defined $new_consequence ) {
	empty_codon(\@out,\@same_codon);
      push @out, $consequence_type; # should be empty
      next;
    }

    if ( !defined $new_consequence->aa_start){
	empty_codon(\@out,\@same_codon);
	push @out, $new_consequence;
	next;
    }
    #first element of the codon
    if (!defined $same_codon[0]){
	push @{$same_codon[0]}, $new_consequence; #goes to the first position
	next;
    }
    #for alleles with aa effect, find out if they are in the same codon
    if ($same_codon[-1]->[0]->aa_start == $new_consequence->aa_start){
	#they are in the same codon, find out if it is the same position
	if ($same_codon[-1]->[0]->start == $new_consequence->start){
	    #it is the same position
	    push @{$same_codon[-1]},$new_consequence; #push in the last 
	}
	else{	    
	    push @{$same_codon[$#same_codon + 1]},$new_consequence; #this is a new element in the codon
	}

    }
    else{
	#if there is more than one element in the same_codon array, calculate the effect of the codon
	if (@same_codon > 1){
	    calculate_same_codon(\@same_codon);
	}
	map {push @out, @{$_}} @same_codon;
	@same_codon = ();
	push @{$same_codon[0]}, $new_consequence; #push the element not in the same codon
    }
  }
  #add last consequence_type
  empty_codon(\@out,\@same_codon);

  return \@out;
}

sub empty_codon{
    my $out = shift;
    my $same_codon = shift;

    if (@{$same_codon} == 1){
	map {push @{$out}, @{$_}} @{$same_codon};
    }
    elsif (@{$same_codon} > 1){
	calculate_same_codon($same_codon);
	map {push @{$out}, @{$_}} @{$same_codon};    
    }
    @{$same_codon} = ();
}

# recalculates the effect of 2 or 3 SNPs in the same codon
sub calculate_same_codon{
    my $same_codon = shift;
    my $new_codon;
    my $old_aa;
    my $codon_table = Bio::Tools::CodonTable->new;
    if (@{$same_codon} == 3){
	#if there are 3 alleles in the same codon
	map {$new_codon .= @{$_->[0]->alleles};$old_aa = $_->[0]->aa_alleles()->[0]} @{$same_codon};
    }
    else{
	#if there are 2 alleles affecting the same codon
	my $first_pos = ($same_codon->[0]->[0]->cdna_start -1) % 3; #position of the first allele in the codon
	my $second_pos = ($same_codon->[1]->[0]->cdna_start -1)% 3; #position of the second allele in the codon
	if ($first_pos == 0){
	    #codon starts with first allele
	    $new_codon = $same_codon->[0]->[0]->alleles->[0]; #first base in the codon
	    if ($second_pos == 1){
		$new_codon .= $same_codon->[1]->[0]->alleles->[0]; #second base in the codon
		$new_codon .= substr($same_codon->[1]->[0]->codon,2,1); #third base in the codon
	    }
	    else{
		$new_codon .= substr($same_codon->[1]->[0]->codon,1,1); #second base in the codon
		$new_codon .= $same_codon->[1]->[0]->alleles->[0]; #third base in the codon
	    }
	}
	else{
	    #alleles are in position 1 and 2 in the codon
	    $new_codon = substr($same_codon->[1]->[0]->codon,0,1); #first base in the codon
	    $new_codon .= $same_codon->[0]->[0]->alleles->[0]; #second base in the codon
	    $new_codon .= $same_codon->[1]->[0]->alleles->[0]; #third base in the codon
	}
	$old_aa = $same_codon->[0]->[0]->aa_alleles->[0];	
    }
    #calculate the new_aa
    my $new_aa = $codon_table->translate($new_codon);
    #and update the aa_alleles field in all the codons
    foreach my $codon (@{$same_codon}){
	map {$_->aa_alleles([$old_aa,$new_aa])} @{$codon};
    }

}
#
# Classifies a variation which is in the vicinity of a transcript
#
sub type_variation {
  my $tr  = shift;
  my $g   = shift;
  my $var = shift;

  my $UPSTREAM = 5000;
  my $DOWNSTREAM = 5000;

  #empty type first in the case of recursive call
  $var->empty_type if defined $var->type;
  
  if (!$var->isa('Bio::EnsEMBL::Variation::ConsequenceType')){
      throw("Not possible to calculate the consequence type for ",ref($var)," : Bio::EnsEMBL::Variation::ConsequenceType object expected");
  }

  if (($var->start < $tr->start - $UPSTREAM) || ($var->end  > $tr->end + $DOWNSTREAM)){
    #since the variation is more than UPSTREAM and DOWNSTREAM of the transcript, there is no effect in the transcript
    return [];
  }

  if (!$tr->translation()) {#for other biotype rather than coding/IG genes
    # check if the variation is completely outside the transcript:

    if($var->end < $tr->start()) {
      $var->type( ($tr->strand() == 1) ? 'UPSTREAM' : 'DOWNSTREAM' );
      return [$var];
    }
    if($var->start > $tr->end()) {
      $var->type( ($tr->strand() == 1) ? 'DOWNSTREAM' : 'UPSTREAM' );
      return [$var];
    }
    if ($var->start >= $tr->start() and $var->end <= $tr->end()) {#within the transcript
      if ($tr->biotype() eq "miRNA") {
	my ($attribute) = @{$tr->get_all_Attributes('miRNA')};
	#the value is the mature miRNA coordinate within miRNA transcript
	if ( $attribute->value =~ /(\d+)-(\d+)/ ) {
	  my @mapper_objs = $tr->cdna2genomic($1, $2, $tr->strand);#transfer cdna value to genomic coordinates
	  foreach my $obj ( @mapper_objs ){#Note you can get more than one mature seq per miRNA
	    if( $obj->isa("Bio::EnsEMBL::Mapper::Coordinate")){
	      if ($var->start >= $obj->start() and $var->end <= $obj->end()) {
		$var->type("WITHIN_MATURE_miRNA");
		return [$var];
	      }
	    }
	  }
	}
      }
      $var->type("WITHIN_NON_CODING_GENE");
      return [$var];
    }
  }

  my $tm = $tr->get_TranscriptMapper();
  my @coords = $tm->genomic2cdna($var->start,
                                 $var->end,
                                 $var->strand);

  # Handle simple cases where the variation is not split into parts.
  # Call method recursively with component parts in complicated case.
  # E.g. a single multi-base variation may be both intronic and coding


  if(@coords > 1) {
    my @out;
    #this will be a new type, complex_indel
    $var->type('COMPLEX_INDEL');
    return [$var];
#     foreach my $c (@coords) {
#       my %new_var = %{$var};
#       $new_var{'end'} = $var->start + $c->length() - 1;
#       $var->start( $new_var{'end'} + 1);
#       #empty the type before re-run
#       $var->empty_type ;  
#       push @out, @{type_variation($tr, $g, bless \%new_var, ref($var))};
#     }
#    return \@out;


  }

  my @coords_splice_2 = $tm->genomic2cdna($var->start -2, $var->end +2, $var->strand);
  my @coords_splice_3 = $tm->genomic2cdna($var->start -3, $var->end +3, $var->strand);
  my @coords_splice_8 = $tm->genomic2cdna($var->start -8, $var->end +8, $var->strand);

  my ($splice_site_2, $splice_site_3, $splice_site_8);

  if (scalar @coords_splice_2 >1) {
    $splice_site_2=1;
  }
  elsif (scalar @coords_splice_3 >1) {
    $splice_site_3=1;
  }
  elsif (scalar @coords_splice_8 >1) {
    $splice_site_8=1;
  }
  

  my $c = $coords[0];
  if($c->isa('Bio::EnsEMBL::Mapper::Gap')) {

    # check if the variation is completely outside the transcript:

    if($var->end < $tr->start()) {
      $var->type( ($tr->strand() == 1) ? 'UPSTREAM' : 'DOWNSTREAM' );
      return [$var];
    }
    if($var->start > $tr->end()) {
      $var->type( ($tr->strand() == 1) ? 'DOWNSTREAM' : 'UPSTREAM' );
      return [$var];
    }

    # variation must be intronic since mapped to cdna gap, but is within
    # transcript, note that ESSENTIAL_SPLICE_SITE only consider first (AG) and last (GT) 2 bases inside the intron.
    # if variation is in intron, we need to check the lenth of intron, if it's shoter than 6, we call it SYNONYMOUS_CODING rather then INTRONIC

    foreach my $intron (@{$tr->get_all_Introns()}) {
      if ($intron->length <=5) {#the length of frameshift intron could be 1,2,4,5 bases
	if ($var->start>=$intron->start and $var->end<=$intron->end) {
	  #this is a type of SYNONYMOUS_CODING since changes happen in frameshift intron, which don't change exon structure
	  $var->type('SYNONYMOUS_CODING');
	  return [$var];
	}
      }
    }
    #if it's not in frameshift intron, then it's in normal intron
    $var->type('INTRONIC');

    if ($splice_site_2) {
      $var->type('ESSENTIAL_SPLICE_SITE');
    }
    elsif ($splice_site_3 or $splice_site_8) {
      $var->type('SPLICE_SITE');
    }
    return [$var];
  }

  #now variation must be in exons, the first 3 bs into exon could be splice_site

  if ($splice_site_2 or $splice_site_3) {
    $var->type('SPLICE_SITE');
  }
  
  $var->cdna_start( $c->start() );
  $var->cdna_end( $c->end() );

  @coords = $tm->genomic2cds($var->start, $var->end,$var->strand);

  if(@coords > 1) {
#    my @out;
      #this is a new type, complex_indel
      $var->type('COMPLEX_INDEL');
      return [$var];
#     foreach my $c (@coords) {
#       my %new_var = %{$var};
#       $new_var{'end'} = $var->start + $c->length() - 1;
#       $var->start( $new_var{'end'} + 1);
#       #empty the type before re-run       
#       $var->empty_type ;
#       push @out, @{type_variation($tr, $g, bless \%new_var, ref($var))};
#     }
#     return \@out;
  }

  $c = $coords[0];

  if($c->isa('Bio::EnsEMBL::Mapper::Gap')) {
    # mapped successfully to CDNA but not to CDS, must be UTR

    if($var->end < $tr->coding_region_start()) {
      $var->type( ($tr->strand() == 1) ? '5PRIME_UTR' : '3PRIME_UTR' );
    }
    elsif($var->start > $tr->coding_region_end()) {
      $var->type( ($tr->strand() == 1) ? '3PRIME_UTR' : '5PRIME_UTR');
    }
    else {
      throw('Unexpected: CDNA variation which is not in CDS is not in UTR');
    }
    return [$var];
  }
  
  $var->cds_start( $c->start() );
  $var->cds_end( $c->end() );

  @coords = $tm->genomic2pep($var->start, $var->end, $var->strand);


  if(@coords != 1 || $coords[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
    throw("Unexpected: Could map to CDS but not to peptide coordinates.");
  }

  $c = $coords[0];

  $var->aa_start( $c->start());
  $var->aa_end( $c->end());

  apply_aa_change($tr, $var);

  return [$var];
}

#
# Determines the effect of a coding variation on the peptide sequence
#

sub apply_aa_change {
  my $tr = shift;
  my $var = shift;

  #my $peptide = $tr->translate->seq();
  #to consider stop codon as well
  my $mrna = $tr->translateable_seq();
  my $mrna_seqobj = Bio::Seq->new( -seq        => $mrna,
				   -moltype    => "dna",
				   -alphabet   => "dna");

  my ($attrib) = @{$tr->slice()->get_all_Attributes('codon_table')}; #for mithocondrial dna it is necessary to change the table

  my $codon_table;
  $codon_table = $attrib->value() if($attrib);
  $codon_table ||= 1; # default vertebrate codon table 

  my $peptide = $mrna_seqobj->translate(undef,undef,undef,$codon_table)->seq;
  my $len = $var->aa_end - $var->aa_start + 1;
  my $old_aa = substr($peptide, $var->aa_start -1 , $len);

  my $codon_cds_start = $var->aa_start * 3 - 2;
  my $codon_cds_end   = $var->aa_end   * 3;
  my $codon_len = $codon_cds_end - $codon_cds_start + 1;

  my @alleles = @{$var->alleles};

  my $var_len = $var->cds_end - $var->cds_start + 1;

  my @aa_alleles = ($old_aa);

  #here could generate multi type if have multi-allele change: "ACTAGT/-/T"
  foreach my $a (@alleles) {
    $a =~ s/\-//;
    my $cds = $tr->translateable_seq();
    if($var_len != length($a)) {
      if(abs(length($a) - $var_len) % 3) {
        # frameshifting variation, do not set peptide_allele string
        # since too complicated and could be very long
	
	$var->type('FRAMESHIFT_CODING');
        return [$var];
      }

      if($codon_len == 0) { # insertion
        $aa_alleles[0] = '-';
        $old_aa    = '-';
      }
    }

    my $new_aa;

    substr($cds, $var->cds_start-1, $var_len) = $a;
	
    my $codon_str = substr($cds, $codon_cds_start-1, $codon_len + length($a)-$var_len);
    
    $var->codon($codon_str); #add the codon to the ConsequenceType object
    my $codon_seq = Bio::Seq->new(-seq      => $codon_str,
				  -moltype  => 'dna',
				  -alphabet => 'dna');

    $new_aa = $codon_seq->translate(undef,undef,undef,$codon_table)->seq();

    if(length($new_aa)<1){
      $new_aa='-';
    }
    
    if(uc($new_aa) ne uc($old_aa)) {
      push @aa_alleles, $new_aa;
      if ($new_aa =~ /\*/) {
	$var->type('STOP_GAINED');
      }
      elsif ($old_aa =~ /\*/) {
	$var->type('STOP_LOST');
      }
    }
  }
  
  #note if type is already defined as SOTP_GAINED OR STOP_LOST, then even @aa_alleles > 1, we are not given type
  # of 'NON_SYNONYMOUS_CODING'
  if(@aa_alleles > 1) {
    if (!$var->type or (join ' ',@{$var->type}) !~ /STOP/) {
      $var->type('NON_SYNONYMOUS_CODING');
    }
  }
  else {
    $var->type('SYNONYMOUS_CODING');
  }

  $var->aa_alleles(\@aa_alleles);
}


1;
