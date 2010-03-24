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

Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblGeneGeneric - default Ensembl
InternalIdMapper implementation for genes

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblGeneGeneric;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::InternalIdMapper::BaseMapper;
our @ISA = qw(Bio::EnsEMBL::IdMapping::InternalIdMapper::BaseMapper);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);


#
# basic mapping
#
sub init_basic {
  my $self = shift;
  my $num = shift;
  my $gsb = shift;
  my $mappings = shift;
  my $gene_scores = shift;

  $self->logger->info("Basic gene mapping...\n", 0, 'stamped');

  $mappings = $self->basic_mapping($gene_scores, "gene_mappings$num");
  $num++;
  my $new_scores = $gsb->create_shrinked_matrix($gene_scores, $mappings,
    "gene_matrix$num");

  return ($new_scores, $mappings);
}


#
# build the synteny from unambiguous mappings
#
sub synteny {
  my $self = shift;
  my $num = shift;
  my $gsb = shift;
  my $mappings = shift;
  my $gene_scores = shift;

  unless ($gene_scores->loaded) {
    $self->logger->info("Synteny Framework building...\n", 0, 'stamped');
    my $dump_path = path_append($self->conf->param('basedir'), 'mapping');
    my $sf = Bio::EnsEMBL::IdMapping::SyntenyFramework->new(
      -DUMP_PATH    => $dump_path,
      -CACHE_FILE   => 'synteny_framework.ser',
      -LOGGER       => $self->logger,
      -CONF         => $self->conf,
      -CACHE        => $self->cache,
    );
    $sf->build_synteny($mappings);

    # use it to rescore the genes
    $self->logger->info("\nSynteny assisted mapping...\n", 0, 'stamped');
    $gene_scores = $sf->rescore_gene_matrix_lsf($gene_scores);

    # checkpoint
    $gene_scores->write_to_file;
  }

  my $new_mappings = $self->basic_mapping($gene_scores, "gene_mappings$num");
  $num++;
  my $new_scores = $gsb->create_shrinked_matrix($gene_scores, $new_mappings,
    "gene_matrix$num");

  return ($new_scores, $new_mappings); 
}


#
# rescore with simple scoring function and try again
#
sub best_transcript {
  my $self = shift;
  my $num = shift;
  my $gsb = shift;
  my $mappings = shift;
  my $gene_scores = shift;
  my $transcript_scores = shift;

  $self->logger->info("Retry with simple best transcript score...\n", 0, 'stamped');
  
  unless ($gene_scores->loaded) {
    $gsb->simple_gene_rescore($gene_scores, $transcript_scores);
    $gene_scores->write_to_file;
  }
  
  my $new_mappings = $self->basic_mapping($gene_scores, "gene_mappings$num");
  $num++;
  my $new_scores = $gsb->create_shrinked_matrix($gene_scores, $new_mappings,
    "gene_matrix$num");

  return ($new_scores, $new_mappings); 
}


#
# rescore by penalising scores between genes with different biotypes  
#
sub biotype {
  my $self = shift;
  my $num = shift;
  my $gsb = shift;
  my $mappings = shift;
  my $gene_scores = shift;

  $self->logger->info("Retry with biotype disambiguation...\n", 0, 'stamped');
  
  unless ($gene_scores->loaded) {
    $gsb->biotype_gene_rescore($gene_scores);
    $gene_scores->write_to_file;
  }

  my $new_mappings = $self->basic_mapping($gene_scores, "gene_mappings$num");
  $num++;
  my $new_scores = $gsb->create_shrinked_matrix($gene_scores, $new_mappings,
    "gene_matrix$num");

  return ($new_scores, $new_mappings); 
}


#
# selectively rescore by penalising scores between genes with different
# internalIDs  
#
sub internal_id {
  my $self = shift;
  my $num = shift;
  my $gsb = shift;
  my $mappings = shift;
  my $gene_scores = shift;

  $self->logger->info("Retry with internalID disambiguation...\n", 0, 'stamped');
  
  unless ($gene_scores->loaded) {
    $gsb->internal_id_rescore($gene_scores);
    $gene_scores->write_to_file;
  }

  my $new_mappings = $self->basic_mapping($gene_scores, "gene_mappings$num");
  $num++;
  my $new_scores = $gsb->create_shrinked_matrix($gene_scores, $new_mappings,
    "gene_matrix$num");

  return ($new_scores, $new_mappings); 
}


1;

