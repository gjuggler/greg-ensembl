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

Bio::EnsEMBL::IdMapping::ScoreBuilder - score builder base class

=head1 SYNOPSIS

This class is not instantiated. Please see subclasses for usage examples
(e.g.  GeneScoreBuilder).

=head1 DESCRIPTION

This is the base class for the score builders used in the stable Id
mapping application. It contains methods which are used by more than one
ScoreBuilder.

=head1 METHODS

  create_shrinked_matrix
  internal_id_rescore
  log_matrix_stats

=cut

package Bio::EnsEMBL::IdMapping::ScoreBuilder;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::BaseObject;
our @ISA = qw(Bio::EnsEMBL::IdMapping::BaseObject);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);
use Bio::EnsEMBL::IdMapping::ScoredMappingMatrix;


=head2 create_shrinked_matrix

  Arg[1]      : Bio::EnsEMBL::Idmapping::ScoredMappingMatrix $matrix - a scoring
                matrix
  Arg[2]      : Bio::EnsEMBL::Idmapping::MappingList $mappings - mappings
  Arg[3]      : String $cache_file - base name of a cache file (extension '.ser'
                will be added automatically) for the returned matrix
  Example     : my $new_scores = $score_builder->create_shrinked_matrix(
                  $gene_scores, $mappings, "gene_matrix1");
  Description : Create a shrinked scoring matrix which doesn't contain entries
                which were already mapped. It also logs how many new mappings
                were added in this process.
  Return type : Bio::EnsEMBL::IdMapping::ScoredMappingMatrix
  Exceptions  : thrown on wrong or missing arguments
  Caller      : InternalIdMapper plugin
  Status      : At Risk
              : under development

=cut

#
sub create_shrinked_matrix {
  my $self = shift;
  my $matrix = shift;
  my $mappings = shift;
  my $cache_file = shift; # base name, extension '.ser' will be added

  # argument checks
  unless ($matrix and
          $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }
  
  unless ($mappings and
          $mappings->isa('Bio::EnsEMBL::IdMapping::MappingList')) {
    throw('Need a gene Bio::EnsEMBL::IdMapping::MappingList.');
  }

  throw('Need a cache file name.') unless ($cache_file);

  my $dump_path = path_append($self->conf->param('basedir'), 'matrix');
  $cache_file .= '.ser';

  my $shrinked_matrix = Bio::EnsEMBL::IdMapping::ScoredMappingMatrix->new(
    -DUMP_PATH   => $dump_path,
    -CACHE_FILE  => $cache_file,
    -AUTO_LOAD   => 1,
  );

  # if we already found a saved matrix, just return it
  if ($shrinked_matrix->loaded) {
  
    $self->logger->info("Read existing scoring matrix from $cache_file.\n");
  
  } else {
    
    # create lookup hashes for sources and targets in the MappingList
    my %sources = ();
    my %targets = ();

    foreach my $entry (@{ $mappings->get_all_Entries }) {
      $sources{$entry->source} = 1;
      $targets{$entry->target} = 1;
    }

    # add all entries to shrinked matrix which are not in the MappingList
    foreach my $entry (@{ $matrix->get_all_Entries }) {
      unless ($sources{$entry->source} or $targets{$entry->target}) {
        $shrinked_matrix->add_Entry($entry);
      }
    }

  }

  # log shrinking stats
  $self->logger->info('Sources '.$matrix->get_source_count.' --> '.
    $shrinked_matrix->get_source_count."\n");
  $self->logger->info('Targets '.$matrix->get_target_count.' --> '.
    $shrinked_matrix->get_target_count."\n");
  $self->logger->info('Entries '.$matrix->get_entry_count.' --> '.
    $shrinked_matrix->get_entry_count."\n");
  $self->logger->info('New mappings: '.$mappings->get_entry_count."\n\n");

  return $shrinked_matrix;
}


=head2 internal_id_rescore

  Arg[1]      : Bio::EnsEMBL::Idmapping::ScoredMappingMatrix $matrix - a scoring
                matrix
  Example     : $score_builder->internal_id_rescore($gene_scores);
  Description : Rescore ambiguous mappings based on internal Ids. This is the
                last disambiguation step and is only useful if objects with the
                same internal Id were used in source and target dbs (e.g. in
                patch builds or if objects were copied from source to target).

                If a source and target gene have the same internal Id and there
                are mappings to other target genes then these *other* mappings
                are penalised.
  Return type : none
  Exceptions  : thrown on wrong or missing argument
  Caller      : InternalIdMapper plugins
  Status      : At Risk
              : under development

=cut

sub internal_id_rescore {
  my $self = shift;
  my $matrix = shift;

  unless ($matrix and
          $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  my $i = 0;

  foreach my $source (@{ $matrix->get_all_sources }) {

    my @entries = sort { $b <=> $a }
      @{ $matrix->get_Entries_for_source($source) };

    # nothing to do if we only have one mapping
    next unless (scalar(@entries) > 1);

    # only penalise if mappings are ambiguous
    next unless ($entries[0]->score == $entries[1]->score);

    # only penalise if one source id == target id where score == best score
    my $ambiguous = 0;
    
    foreach my $e (@entries) {
      if ($e->target == $source and $e->score == $entries[0]) {
        $ambiguous = 1;
      }
    }

    next unless ($ambiguous);

    # now penalise those where source id != target id and score == best score
    foreach my $e (@entries) {
      if ($e->target != $source and $e->score == $entries[0]) {
        $matrix->set_score($source, $e->target, ($e->score * 0.8));
        $i++;
      }
    }

  }
  
  $self->logger->debug("Scored entries with internal ID mismatch: $i\n", 1);
}


=head2 log_matrix_stats

  Arg[1]      : Bio::EnsEMBL::Idmapping::ScoredMappingMatrix $matrix - a scoring
                matrix
  Example     : $score_builder->log_matrix_stats;
  Description : Logs scoring matrix statistics (number of entries, min/max/avg
                scores).
  Return type : none
  Exceptions  : thrown on wrong or missing argument
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub log_matrix_stats {
  my $self = shift;
  my $matrix = shift;

  unless ($matrix and
          $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('You must provide a ScoredMappingMatrix.');
  }

  my $fmt1 = "%-40s%10.0f\n";
  my $fmt2 = "%-40s%10.5f\n";
  
  $self->logger->info(sprintf($fmt1, "Scoring matrix entries:",
    $matrix->get_entry_count), 1);
  
  $self->logger->info(sprintf($fmt1, "Scoring matrix sources:",
    $matrix->get_source_count), 1);
  
  $self->logger->info(sprintf($fmt1, "Scoring matrix targets:",
    $matrix->get_target_count), 1);
  
  $self->logger->info(sprintf($fmt2, "Average score:",
    $matrix->get_average_score), 1);
  
  my ($min, $max) = @{ $matrix->get_min_max_scores };
  $self->logger->info(sprintf($fmt2, "Min. score:", $min), 1);
  $self->logger->info(sprintf($fmt2, "Max. score:", $max), 1);
}


1;

