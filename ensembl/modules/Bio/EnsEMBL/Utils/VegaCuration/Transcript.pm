=head1 LICENSE

  Copyright (c) 1999-2010 The European Bioinformatics Institute and
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

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::VegaCuration::Transcript;

use strict;
use warnings;
no warnings 'uninitialized';
use vars qw(@ISA);

use Bio::EnsEMBL::Utils::VegaCuration::Gene;
use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Utils::VegaCuration::Gene);


=head2 find_non_overlaps

   Args       : arrayref of B::E::Transcripts
   Example    : find_non_overlaps($all_transcripts)
   Description: identifies any non-overlapping transcripts
   Returntype : array refs of stable IDs
   Exceptions : none

=cut

sub find_non_overlaps {
  my $self = shift;
  my ($all_transcripts) = @_;
  my $non_overlaps = [];
  foreach my $transcript1 (@{$all_transcripts}) {
    foreach my $transcript2 (@{$all_transcripts}) {
      if ($transcript1->end < $transcript2->start) {
	push @{$non_overlaps}, $transcript1->stable_id;
	push @{$non_overlaps}, $transcript2->stable_id;
      }
    }
  }
  return $non_overlaps;
}

=head2 check_remarks_and_update names

   Arg[1]     : B::E::Gene (with potentially duplicated transcript names)
   Arg[2]     : counter 1 (no. of patched genes)
   Arg[3]     : counter 2 (no. of patched transcripts)
   Example    : $support->update_names($gene,\$c1,\$c2)
   Description: - checks remarks and patches transcripts with identical names according to
                CDS and length
                - adds remark to gene if there arefragmented gene/transcript_remarks
   Returntype : true | false (depending on whether patched or not), counter1, counter2

=cut

sub check_remarks_and_update_names {
  my $self = shift;
  my ($gene,$gene_c,$trans_c) = @_;
  my $action = ($self->param('dry_run')) ? 'Would add' : 'Added';
  my $aa  = $gene->adaptor->db->get_AttributeAdaptor;
  my $dbh = $gene->adaptor->db->dbc->db_handle;

  #get list of IDs that have previously been sent to annotators
  my $seen_genes = $self->get_havana_fragmented_loci_comments;

  my $gsi    = $gene->stable_id;
  my $gid    = $gene->dbID;
  my $g_name;
  my $study_more = 1;
  eval {
    $g_name = $gene->display_xref->display_id;
  };	
  if ($@) {
    $g_name = $gene->get_all_Attributes('name')->[0]->value;
  }
  my $gene_remark = 'This locus has been annotated as fragmented because either there is not enough evidence covering the whole locus to identify the exact exon structure of the transcript, or because the transcript spans a gap in the assembly';
  my $attrib = [
    Bio::EnsEMBL::Attribute->new(
      -CODE => 'remark',
      -NAME => 'Remark',
      -DESCRIPTION => 'Annotation remark',
      -VALUE => $gene_remark,
    ) ];

  #get existing gene and transcript remarks
  my %remarks;
  foreach my $type ('remark','hidden_remark') {
    $remarks{$type}->{'gene'} = [ map {$_->value} @{$gene->get_all_Attributes($type)} ];
    foreach my $trans (@{$gene->get_all_Transcripts()}) {
      my $tsi = $trans->stable_id;
      push @{$remarks{$type}->{'transcripts'}}, map {$_->value} @{$trans->get_all_Attributes('remark')};
    }
  }

  #if any of the remarks identify this gene as being known by Havana as being fragmented...
  if ( (grep {$_ =~ /fragmen/i }
        @{$remarks{'hidden_remark'}->{'gene'}},
	@{$remarks{'remark'}->{'gene'}},
	@{$remarks{'remark'}->{'transcripts'}}, 
	@{$remarks{'hidden_remark'}->{'transcripts'}} ) ) {
    if (grep { $_ eq $gene_remark} @{$remarks{'remark'}->{'gene'}}) {
      $self->log_verbose("Fragmented loci annotation remark for gene $gsi already exists\n");
    }
    #add gene_attrib
    else {
      if (! $self->param('dry_run') ) {
	$aa->store_on_Gene($gid,$attrib);
      }			
      $self->log("$action correctly formatted fragmented loci annotation remark for gene $gsi\n");
    }
    $study_more = 0;
  }
  #log if it's been reported before since the gene should have a remark.
  elsif ($seen_genes->{$gsi} eq 'fragmented') {
    $self->log_warning("PREVIOUS: $action correctly formatted fragmented loci annotation remark for gene $gsi (has previously been OKeyed by Havana as being fragmented but has no Annotation remark, please add one!)\n");
    #add gene_attrib anyway.
    if (! $self->param('dry_run') ) {
      $aa->store_on_Gene($gid,$attrib);
    }
  }

  ##patch transcript names according to length and CDS
  $gene_c++;

  #separate coding and non_coding transcripts
  my $coding_trans = [];
  my $noncoding_trans = [];
  foreach my $trans ( @{$gene->get_all_Transcripts()} ) {
    if ($trans->translate) {
      push @$coding_trans, $trans;
    }
    else {
      push @$noncoding_trans, $trans;
    }
  }

  #sort transcripts coding > non-coding, then on length
  my $c = 0;
  $self->log("\nPatching names according to CDS and length:\n",1);
  foreach my $array_ref ($coding_trans,$noncoding_trans) {
    foreach my $trans ( sort { $b->length <=> $a->length } @$array_ref ) {
      $trans_c++;
      my $tsi = $trans->stable_id;
      my $t_name;
      eval {
	$t_name = $trans->display_xref->display_id;
      };	
      if ($@) {
	$t_name = $trans->get_all_Attributes('name')->[0]->value;
      }
      $c++;
      my $ext = sprintf("%03d", $c);
      my $new_name = $g_name.'-'.$ext;
      $self->log(sprintf("%-20s%-3s%-20s", "$t_name ", "-->", "$new_name")."\n",1);
      if (! $self->param('dry_run')) {
	
	# update transcript display xref
	$dbh->do(qq(UPDATE  xref x, external_db edb
                                SET     x.display_label  = "$new_name"
                                WHERE   x.external_db_id = edb.external_db_id
                                AND     x.dbprimary_acc  = "$tsi"
                                AND     edb.db_name      = "Vega_transcript"));
      }
    }
  }
  return ($study_more,$gene_c,$trans_c);
}

=head2 check_names_and_overlap

   Arg[1]     : arayref of arrayrefs of duplicated names
   Arg[2]     : B::E::Gene (with potentially duplicated transcript names)
   Arg[3]     : FH (to log new duplicates)
   Example    : $support->check_names_and_overlap($transcripts,$gene,$fh)
   Description: checks pairs of transcripts identified as having duplicate Vega names:
                - to see if they have identical names in loutre (shouldn't have)
                - distinguish between overlapping and non overlapping transcripts
   Returntype : none

=cut

sub check_names_and_overlap {
  my $self = shift;
  my ($transcript_info,$gene,$n_flist_fh) = @_;
  my $ta  = $gene->adaptor->db->get_TranscriptAdaptor;
  my $gsi = $gene->stable_id;
  my $g_name = $gene->get_all_Attributes('name')->[0]->value;
  foreach my $set (values %{$transcript_info} ) {
    next if (scalar @{$set} == 1);
    my $transcripts = [];
    my $all_t_names;
    my %ids_to_names;
    foreach my $id1 (@{$set}) {
      my ($name1,$tsi1) = split /\|/, $id1;
      $ids_to_names{$tsi1} = $name1;
      $all_t_names .= "$tsi1 [$name1] ";
      my $t = $ta->fetch_by_stable_id($tsi1);
      push @{$transcripts}, $t;
    }

    my $non_overlaps;
    eval {
      $non_overlaps = $self->find_non_overlaps($transcripts);
    };
    if ($@) {
      $self->log_warning("Problem looking for overlapping transcripts for gene $gsi (is_current = 0 ?). Skipping this bit\n");
    }

    #if the transcripts don't overlap
    elsif (@{$non_overlaps}) {
      my $tsi_string;
      foreach my $id (@{$non_overlaps}) {
	my $string = " $id [ $ids_to_names{$id} ] ";
	$tsi_string .= $string;
      }

      $self->log_warning("NEW: Non-overlapping: $gsi ($g_name) has non-overlapping transcripts ($tsi_string) with duplicated Vega names, and it has no \'Annotation_remark- fragmented_loci\' on the gene or \'\%fragmen\%\' remark on any transcripts. Neither has it been OKeyed by Havana before. Transcript names are being patched but this needs checking by Havana.\n");
      #log gsi (to be sent to Havana)
      print $n_flist_fh "$gsi\n";
    }
    #...otherwise if the transcripts do overlap
    else {
      $self->log_warning("NEW: Overlapping: $gsi ($g_name) has overlapping transcripts ($all_t_names) with Vega duplicated names and it has no \'Annotation_remark- fragmented_loci\' on the gene or \'\%fragmen\%\' remark on any transcripts. Neither has it been OKeyed by Havana before. Transcript names are being patched but this could be checked by Havana if they were feeling keen.\n");
      print $n_flist_fh "$gsi\n";
    }
  }
}		

=head2 get_havana_fragmented_loci_comments

   Args       : none
   Example    : my $results = $support->get_havana_fragmented_loci_comments
   Description: parses the HEREDOC containing Havana comments in this module
   Returntype : hashref

=cut

sub get_havana_fragmented_loci_comments {
	my $seen_genes;
	while (<DATA>) {
		next if /^\s+$/ or /#+/;
		my ($obj,$comment) = split /=/;
		$obj =~ s/^\s+|\s+$//g;
		$comment =~ s/^\s+|\s+$//g;
		$seen_genes->{$obj} = $comment;
	}
	return $seen_genes;
}



#details of genes with duplicated transcript names that have already been reported to Havana
#identified as either fragmented or as being OK to patch
__DATA__

OTTMUSG00000005478 = fragmented
OTTMUSG00000001936 = fragmented
OTTMUSG00000017081 = fragmented
OTTMUSG00000011441 = fragmented
OTTMUSG00000013335 = fragmented
OTTMUSG00000011654 = fragmented
OTTMUSG00000001835 = fragmented
OTTHUMG00000035221 = fragmented
OTTHUMG00000037378 = fragmented
OTTHUMG00000060732 = fragmented
OTTHUMG00000132441 = fragmented
OTTHUMG00000031383 = fragmented
OTTHUMG00000012716 = fragmented
OTTHUMG00000031102 = fragmented
OTTHUMG00000148816 = fragmented
OTTHUMG00000149059 = fragmented
OTTHUMG00000149221 = fragmented
OTTHUMG00000149326 = fragmented
OTTHUMG00000149644 = fragmented
OTTHUMG00000149574 = fragmented
OTTHUMG00000058101 = fragmented

OTTHUMG00000150119 = OK
OTTHUMG00000149850 = OK
OTTHUMG00000058101 = OK
OTTHUMG00000058907 = OK

OTTMUSG00000011654 = fragmented
OTTMUSG00000019369 = fragmented
OTTMUSG00000017081 = fragmented
OTTMUSG00000001835 = fragmented
OTTMUSG00000011499 = fragmented
OTTMUSG00000013335 = fragmented
OTTMUSG00000008023 = fragmented
OTTMUSG00000019369 = fragmented


OTTMUSG00000022266
OTTMUSG00000006697





OTTMUSG00000012302 =
OTTMUSG00000013368 =
OTTMUSG00000015766 =
OTTMUSG00000016025 =
OTTMUSG00000001066 =
OTTMUSG00000016331 =
OTTMUSG00000006935 =
OTTMUSG00000007263 =
OTTMUSG00000000304 =
OTTMUSG00000009150 =
OTTMUSG00000008023 =
OTTMUSG00000017077 =
OTTMUSG00000003440 =
OTTMUSG00000016310 =
OTTMUSG00000026199 =
OTTMUSG00000028423 =
OTTMUSG00000007427 =
