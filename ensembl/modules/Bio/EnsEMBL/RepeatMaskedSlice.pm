#
# Ensembl module for Bio::EnsEMBL::RepeatMaskedSlice
#
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::RepeatMaskedSlice - Arbitary Slice of a genome

=head1 SYNOPSIS

   $sa = $db->get_SliceAdaptor;

   $slice = $sa->fetch_by_region('chromosome', 'X', 1_000_000, 2_000_000);

   $repeat_masked_slice = $slice->get_repeatmasked_seq();

   #get repeat masked sequence:
   my $dna = $repeat_masked_slice->seq();
   $dna    = $repeat_masked_slice->subseq(1, 1000);


=head1 DESCRIPTION

This is a specialised B<Bio::EnsEMBL::Slice> class that is used to 
retrieve repeat masked genomic sequence rather than normal genomic 
sequence.

=head1 AUTHOR - Graham McVicker

=head1 CONTACT

This modules is part of the Ensembl project http://www.ensembl.org

Questions can be posted to the ensembl-dev mailing list:
ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut


package Bio::EnsEMBL::RepeatMaskedSlice;

use strict;
use warnings;

use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Utils::Exception;

use vars qw(@ISA);

@ISA = ('Bio::EnsEMBL::Slice');

# The BLOCK_PWR is the lob_bin of the chunksize where you want your repeat features
# to be retreived. This will create repeat feature retrieval calls that are likely 
# to be on the same slice and hopefully create cache hits and less database traffic
my $BLOCK_PWR = 18;



=head2 new

  Arg [-REPEAT_MASK] : The logic name of the repeats to be used for masking.
                      If not provided, all repeats in the database are used.
  Arg [...]  : Named superclass arguments. See B<Bio::EnsEMBL::Slice>.
  Example    : my $slice = Bio::EnsEMBL::RepeatMaskedSlice->new
                  (-START  => $start,
                   -END    => $end,
                   -STRAND => $strand,
                   -SEQ_REGION_NAME => $seq_region,
                   -SEQ_REGION_LENGTH => $seq_region_length,
                   -COORD_SYSTEM  => $cs,
                   -ADAPTOR => $adaptor,
                   -REPEAT_MASK => ['repeat_masker'],
                   -SOFT_MASK => 1,
                   -NOT_DEFAULT_MASKING_CASES => {"repeat_class_SINE/MIR" => 1,
                                                  "repeat_name_AluSp" => 0});
  Description: Creates a Slice which behaves exactly as a normal slice but
               that returns repeat masked sequence from the seq method.
  Returntype : Bio::EnsEMBL::RepeatMaskedSlice
  Exceptions : none
  Caller     : RawComputes (PredictionTranscript creation code).
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($logic_names, $soft_mask, $not_default_masking_cases) = rearrange(['REPEAT_MASK',
                                                                         'SOFT_MASK',
                                                                         'NOT_DEFAULT_MASKING_CASES'], @_);

  my $self = $class->SUPER::new(@_);


  $logic_names ||= [''];
  if(ref($logic_names) ne 'ARRAY') {
    throw("Reference to list of logic names argument expected.");
  }

  $self->{'repeat_mask_logic_names'} = $logic_names;
  $self->{'soft_mask'} = $soft_mask;
  $self->{'not_default_masking_cases'} = $not_default_masking_cases;
  $self->{'not_default_masking_cases'} ||= {};
  
  return $self;
}


=head2 repeat_mask_logic_names

  Arg [1]    : reference to list of strings $logic_names (optional)
  Example    : $rm_slice->repeat_mask_logic_name(['repeat_masker']);
  Description: Getter/Setter for the logic_names of the repeats that are used
               to mask this slices sequence.
  Returntype : reference to list of strings
  Exceptions : none
  Caller     : seq() method
  Status     : Stable

=cut

sub repeat_mask_logic_names {
  my $self = shift;

  if(@_) {
    my $array = shift;
    if(ref($array) ne 'ARRAY') {
      throw('Reference to list of logic names argument expected.');
    }
  }
  
  return $self->{'repeat_mask_logic_names'};
}


=head2 soft_mask

  Arg [1]    : boolean $soft_mask (optional)
  Example    : $rm_slice->soft_mask(0);
  Description: Getter/Setter which is used to turn on/off softmasking of the
               sequence returned by seq.
  Returntype : boolean
  Exceptions : none
  Caller     : seq() method
  Status     : Stable

=cut

sub soft_mask {
  my $self = shift;
  $self->{'soft_mask'} = shift if(@_);
  return $self->{'soft_mask'} || 0;
}

=head2 not_default_masking_cases

  Arg [1]    : hash reference $not_default_masking_cases (optional, default is {})
               The values are 0 or 1 for hard and soft masking respectively
               The keys of the hash should be of 2 forms
               "repeat_class_" . $repeat_consensus->repeat_class,
                e.g. "repeat_class_SINE/MIR"
               "repeat_name_" . $repeat_consensus->name
                e.g. "repeat_name_MIR"
               depending on which base you want to apply the not default masking either 
               the repeat_class or repeat_name. Both can be specified in the same hash
               at the same time, but in that case, repeat_name setting has priority over 
               repeat_class. For example, you may have hard masking as default, and 
               you may want soft masking of all repeat_class SINE/MIR,
               but repeat_name AluSp (which are also from repeat_class SINE/MIR)
  Example    : $rm_slice->not_default_masking_cases({"repeat_class_SINE/MIR" => 1,
                                                     "repeat_name_AluSp" => 0});
  Description: Getter/Setter which is used to escape some repeat class or name from the default 
               masking in place. 
  Returntype : hash reference
  Exceptions : none
  Caller     : seq() and subseq() methods
  Status     : Stable

=cut

sub not_default_masking_cases {
  my $self = shift;
  $self->{'not_default_masking_cases'} = shift if (@_);
  return $self->{'not_default_masking_cases'};
}

=head2 seq

  Arg [1]    : none
  Example    : print $rmslice->seq(), "\n";
  Description: Retrieves the entire repeat masked sequence for this slice.
               See also the B<Bio::EnsEMBL::Slice> implementation of this 
               method.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seq {
  my $self = shift;
  #
  # get all the features
  #
  my $logic_names = $self->repeat_mask_logic_names();
  my $soft_mask   = $self->soft_mask();
  my $not_default_masking_cases = $self->not_default_masking_cases;

  my $repeats = [];

  foreach my $l (@$logic_names) {
    push @{$repeats}, @{$self->get_all_RepeatFeatures($l)};
  }

  #
  # get the dna
  #
  my $dna = $self->SUPER::seq(@_);

  #
  # mask the dna
  #
  $self->_mask_features(\$dna,$repeats,$soft_mask,$not_default_masking_cases);
  return $dna;
}



=head2 subseq

  Arg [1]    : none
  Example    : print $rmslice->subseq(1, 1000);
  Description: Retrieves a repeat masked sequence from a specified subregion
               of this slice.  See also the B<Bio::EnsEMBL::Slice> 
               implementation of this method.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub subseq {
  my $self   = shift;
  my $start  = shift;
  my $end    = shift;
  my $strand = shift;

  #
  # get all the features
  #
  my $logic_names = $self->repeat_mask_logic_names();
  my $soft_mask   = $self->soft_mask();
  my $not_default_masking_cases = $self->not_default_masking_cases;


  # If frequent subseqs happen on repeatMasked sequence this results in
  # a lot of feature retrieval from the database. To avoid this, features
  # are only retrieved from subslices with fixed space boundaries. 
  # The access happens in block to make cache hits more likely

  # The blocksize can be defined on the top of this module.

  my $seq_region_slice = $self->seq_region_Slice();
  my $block_min = ($self->start()-1) >> $BLOCK_PWR;
  my $block_max = ($self->end()-1) >> $BLOCK_PWR;

  my $repeats = [];

  my $sub_start = ($block_min << $BLOCK_PWR)+1;
  my $sub_end = ($block_max+1)<<$BLOCK_PWR;
  if ($sub_end > $seq_region_slice->length) {
    $sub_end =  $seq_region_slice->length ;
  }

  my $subslice = $seq_region_slice->sub_Slice( $sub_start, $sub_end);

  foreach my $l (@$logic_names) {
    push @{$repeats}, @{$subslice->get_all_RepeatFeatures($l)};
  }

  #
  # get the dna
  #
  my $subsequence_slice = $self->sub_Slice( $start, $end, $strand );
  my $dna = $subsequence_slice->seq();
  #
  # mask the dna
  #
  $subsequence_slice->_mask_features(\$dna,$repeats,$soft_mask,$not_default_masking_cases);

  return $dna;
}


1;
