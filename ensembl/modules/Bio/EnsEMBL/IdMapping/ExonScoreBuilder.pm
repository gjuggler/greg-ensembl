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

Combines ExonScoreBuilder, ExonDirectMapper and ExonerateRunner from
Java application.

=head1 METHODS

=cut

package Bio::EnsEMBL::IdMapping::ExonScoreBuilder;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::ScoreBuilder;
our @ISA = qw(Bio::EnsEMBL::IdMapping::ScoreBuilder);

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(parse_bytes path_append);
use Bio::EnsEMBL::IdMapping::ScoredMappingMatrix;


#
# exon scoring is done in two steps:
# 1. map exons by overlap (if a common coord_system exists)
# 2. map remaining and poorly scoring exons using exonerate
#
sub score_exons {
  my $self = shift;

  $self->logger->info("-- Scoring exons...\n\n", 0, 'stamped');

  # score using overlaps, then exonerate
  my $matrix = $self->overlap_score;
  my $exonerate_matrix = $self->exonerate_score($matrix);

  # log stats before matrix merging
  $self->logger->info("\nOverlap scoring matrix:\n");
  $self->log_matrix_stats($matrix);
  $self->logger->info("\nExonerate scoring matrix:\n");
  $self->log_matrix_stats($exonerate_matrix);
  
  # merge matrices
  $self->logger->info("\nMerging scoring matrices...\n", 0, 'stamped');
  $matrix->merge($exonerate_matrix);
  $self->logger->info("Done.\n\n", 0, 'stamped');

  # debug logging
  if ($self->logger->loglevel eq 'debug') {
    $matrix->log('exon', $self->conf->param('basedir'));
  }

  # log stats of combined matrix
  $self->logger->info("Combined scoring matrix:\n");
  $self->log_matrix_stats($matrix);
  
  $self->logger->info("\nDone with exon scoring.\n\n", 0, 'stamped');

  return $matrix;
}


#
# direct mapping by overlap (if common coord_system exists)
#
sub overlap_score {
  my $self = shift;

  my $dump_path = path_append($self->conf->param('basedir'), 'matrix');
  
  my $matrix = Bio::EnsEMBL::IdMapping::ScoredMappingMatrix->new(
    -DUMP_PATH   => $dump_path,
    -CACHE_FILE  => 'exon_overlap_matrix.ser',
  );

  my $overlap_cache = $matrix->cache_file;

  if (-s $overlap_cache) {
    
    # read from file
    $self->logger->info("Reading exon overlap scoring matrix from file...\n", 0, 'stamped');
    $self->logger->debug("Cache file $overlap_cache.\n", 1);
    $matrix->read_from_file;
    $self->logger->info("Done.\n", 0, 'stamped');
    
  } else {
    
    # build scoring matrix
    $self->logger->info("No exon overlap scoring matrix found. Will build new one.\n");

    if ($self->cache->highest_common_cs) {
      $self->logger->info("Overlap scoring...\n", 0, 'stamped');
      $matrix = $self->build_overlap_scores($matrix);
      $self->logger->info("Done.\n", 0, 'stamped');
    }

    # write scoring matrix to file
    $matrix->write_to_file;

  }

  return $matrix;
}


#
# map the remaining exons using exonerate
#
sub exonerate_score {
  my $self = shift;
  my $matrix = shift;

  unless ($matrix and
          $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  my $dump_path = path_append($self->conf->param('basedir'), 'matrix');

  my $exonerate_matrix = Bio::EnsEMBL::IdMapping::ScoredMappingMatrix->new(
    -DUMP_PATH   => $dump_path,
    -CACHE_FILE  => 'exon_exonerate_matrix.ser',
  );

  my $exonerate_cache = $exonerate_matrix->cache_file;

  if (-s $exonerate_cache) {

    # read from file
    $self->logger->info("Reading exonerate matrix from file...\n", 0, 'stamped');
    $self->logger->debug("Cache file $exonerate_cache.\n", 1);
    $exonerate_matrix->read_from_file;
    $self->logger->info("Done.\n", 0, 'stamped');

  } else {

    # build scoring matrix
    $self->logger->info("No exonerate matrix found. Will build new one.\n");

    # dump exons to fasta files
    my $dump_count = $self->dump_filtered_exons($matrix);

    if ($dump_count) {
      # run exonerate
      $self->run_exonerate;
      
      # parse results
      $self->parse_exonerate_results($exonerate_matrix);

    } else {

      $self->logger->info("No source and/or target exons dumped, so don't need to run exonerate.\n");

    }
    
    # write scoring matrix to file
    $exonerate_matrix->write_to_file;

  }

  return $exonerate_matrix;
}

#
# Algorithm:
# Get a lists of exon containers for source and target. Walk along both lists,
# set a flag when you first encounter an exon (i.e. it starts). Record all
# alternative exons until you encounter the exon again (exon ends), then score
# against all alternative exons you've recorded.
#
sub build_overlap_scores {
  my $self = shift;
  my $matrix = shift;

  unless ($matrix and
          $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  # get sorted list of exon containers
  $self->logger->info("Reading sorted exons from cache...\n", 1, 'stamped');

  my @source_exons = $self->sort_exons(
    [values %{ $self->cache->get_by_name('exons_by_id', 'source') }]
  );
  my @target_exons = $self->sort_exons(
    [values %{ $self->cache->get_by_name('exons_by_id', 'target') }]
  );

  $self->logger->info("Done.\n", 1, 'stamped');

  # get first source and target exon container
  my $source_ec = shift(@source_exons);
  my $target_ec = shift(@target_exons);

  my %source_overlap = ();
  my %target_overlap = ();
  
  $self->logger->info("Scoring...\n", 1, 'stamped');

  while ($source_ec or $target_ec) {

    my $add_source = 0;
    my $add_target = 0;

    # compare exon containers
    if ($source_ec and $target_ec) {
      my $cmp = $self->compare_exon_containers($source_ec, $target_ec);
      $add_source = 1 if ($cmp <= 0);
      $add_target = 1 if ($cmp >= 0);
    } elsif ($source_ec) {
      $add_source = 1;
    } else {
      $add_target = 1;
    }

    if ($add_source) {
      if ($source_overlap{$source_ec->[0]}) {
        # remove exon from list of overlapping source exons to score target
        # against
        delete $source_overlap{$source_ec->[0]};
      } else {
        # add exon to list of overlapping source exons to score target against
        $source_overlap{$source_ec->[0]} = $source_ec->[0];

        # score source exon against all target exons in current overlap list
        foreach my $target_exon (values %target_overlap) {
          next if (defined($matrix->get_score(
            $source_ec->[0]->id, $target_exon->id)));

          $self->calc_overlap_score($source_ec->[0], $target_exon,
            $matrix);
        }
      }

      # get next source exon container
      $source_ec = shift(@source_exons);
    }

    if ($add_target) {
      if ($target_overlap{$target_ec->[0]}) {
        # remove exon from list of overlapping target exons to score source
        # against
        delete $target_overlap{$target_ec->[0]};
      } else {
        # add exon to list of overlapping target exons to score source against
        $target_overlap{$target_ec->[0]} = $target_ec->[0];

        # score target exon against all source exons in current overlap list
        foreach my $source_exon (values %source_overlap) {
          next if (defined($matrix->get_score(
            $source_exon->id, $target_ec->[0]->id)));

          $self->calc_overlap_score($source_exon, $target_ec->[0], $matrix);
        }
      }

      # get next target exon container
      $target_ec = shift(@target_exons);
    }
  }

  $self->logger->info("Done.\n", 1, 'stamped');

  return $matrix;
}


#
# Return a list of exon containers, sorted by seq_region_name, then location
# (where location is either start-1 or end, so each exon is in the list twice).
# An exon container is a listrefs of a TinyExon object and its location. This
# implements the ExonSortContainer in the java application.
#
sub sort_exons {
  my $self = shift;
  my $exons = shift;

  return
    sort { ($a->[0]->common_sr_name cmp $b->[0]->common_sr_name)
           || ($a->[1] <=> $b->[1]) }
      (map { [$_, $_->common_start - 1] } @$exons),
      (map { [$_, $_->common_end] } @$exons);
}


sub compare_exon_containers {
  my $self = shift;
  my $e1 = shift;
  my $e2 = shift;

  return ( ($e1->[0]->common_sr_name cmp $e2->[0]->common_sr_name) ||
           ($e1->[1] <=> $e2->[1]) );
}

#
# Calculates overlap score between two exons. Its done by dividing overlap
# region by exons sizes. 1.0 is full overlap on both exons. Score of at least
# 0,5 are added to the exon scoring matrix.
#
sub calc_overlap_score {
  my $self = shift;
  my $source_exon = shift;
  my $target_exon = shift;
  my $matrix = shift;

  my ($start, $end);

  # don't score if exons on different strand
  return unless ($source_exon->strand == $target_exon->strand);
  
  # determine overlap start
  if ($source_exon->start > $target_exon->start) {
    $start = $source_exon->start;
  } else {
    $start = $target_exon->start;
  }
  
  # determine overlap end
  if ($source_exon->end < $target_exon->end) {
    $end = $source_exon->end;
  } else {
    $end = $target_exon->end;
  }

  #
  # calculate score, which is defined as average overlap / exon length ratio
  #
  my $overlap = $end - $start + 1;
  my $source_length = $source_exon->end - $source_exon->start + 1;
  my $target_length = $target_exon->end - $target_exon->start + 1;
  
  my $score = ($overlap/$source_length + $overlap/$target_length)/2;

  # penalise by 10% if phase if different
  $score *= 0.9 if ($source_exon->phase != $target_exon->phase);

  # add score to scoring matrix if it's at least 0.5
  if ($score >= 0.5) {
    $matrix->add_score($source_exon->id, $target_exon->id, $score);
  }
}


sub run_exonerate {
  my $self = shift;

  my $source_file = $self->exon_fasta_file('source');
  my $target_file = $self->exon_fasta_file('target');
  my $source_size = -s $source_file;
  my $target_size = -s $target_file;

  # check if fasta files exist and are not empty
  unless ($source_size and $target_size) {
    throw("Can't find exon fasta files.");
  }

  # create an empty lsf log directory
  my $logpath = path_append($self->logger->logpath, 'exonerate');
  system("rm -rf $logpath") == 0 or
    $self->logger->error("Unable to delete lsf log dir $logpath: $!\n");
  system("mkdir -p $logpath") == 0 or
    $self->logger->error("Can't create lsf log dir $logpath: $!\n");

  # delete exonerate output from previous runs
  my $dump_path = $self->cache->dump_path;

  opendir(DUMPDIR, $dump_path) or
    $self->logger->error("Can't open $dump_path for reading: $!");

  while (defined(my $file = readdir(DUMPDIR))) {
    next unless /exonerate_map\.\d+/;

    unlink("$dump_path/$file") or
      $self->logger->error("Can't delete $dump_path/$file: $!");
  }
  
  closedir(DUMPDIR);

  # determine number of jobs to split task into
  my $bytes_per_job = $self->conf->param('exonerate_bytes_per_job') || 250000;
  my $num_jobs = $self->conf->param('exonerate_jobs');
  $num_jobs ||= int($source_size/$bytes_per_job + 1);
  
  my $percent = ($self->conf->param('exonerate_threshold') || 0.5) * 100;
  my $lsf_name = 'idmapping_exonerate_'.time;
  my $exonerate_path = $self->conf->param('exonerate_path');
  my $exonerate_extra_params = $self->conf->param('exonerate_extra_params');

  #
  # run exonerate jobs using lsf
  #
  my $exonerate_job = qq{$exonerate_path } .
    qq{--query $source_file --target $target_file } .
    q{--querychunkid $LSB_JOBINDEX } .
    qq{--querychunktotal $num_jobs } .
    q{--model affine:local -M 900 --showalignment FALSE --subopt no } .
    qq{--percent $percent } .
    $self->conf->param('exonerate_extra_params') . " " .
    q{--ryo 'myinfo: %qi %ti %et %ql %tl\n' } .
    qq{| grep '^myinfo:' > $dump_path/exonerate_map.\$LSB_JOBINDEX} . "\n";
  
  $self->logger->info("Submitting $num_jobs exonerate jobs to lsf.\n");
  $self->logger->debug("$exonerate_job\n\n");

  local *BSUB;
  open BSUB, "|bsub -J$lsf_name\[1-$num_jobs\] -o $logpath/exonerate.\%I.out"
    or $self->logger->error("Could not open open pipe to bsub: $!\n");

  print BSUB $exonerate_job;
  $self->logger->error("Error submitting exonerate jobs: $!\n")
    unless ($? == 0); 
  close BSUB;

  # submit dependent job to monitor finishing of exonerate jobs
  $self->logger->info("Waiting for exonerate jobs to finish...\n", 0, 'stamped');

  my $dependent_job = qq{bsub -K -w "ended($lsf_name)" -q small } .
    qq{-o $logpath/exonerate_depend.out /bin/true};

  system($dependent_job) == 0 or
    $self->logger->error("Error submitting dependent job: $!\n");

  $self->logger->info("All exonerate jobs finished.\n", 0, 'stamped');

  #
  # check results
  #
  my @missing;
  my @error;
  
  for (my $i = 1; $i <= $num_jobs; $i++) {
  
    # check that output file exists
    my $outfile = "$dump_path/exonerate_map.$i";
    push @missing, $outfile unless (-f "$outfile");

    # check no errors occurred
    my $errfile = "$logpath/exonerate.$i.err";
    push @error, $errfile if (-s "$errfile");
  }

  if (@missing) {
    $self->logger->info("Couldn't find all exonerate output files. These are missing:\n");
    foreach (@missing) {
      $self->logger->info("$_\n", 1);
    }

    exit(1);
  }

  if (@error) {
    $self->logger->info("One or more exonerate jobs failed. Check these error files:\n");
    foreach (@error) {
      $self->logger->info("$_\n", 1);
    }

    exit(1);
  }

}


sub exon_fasta_file {
  my $self = shift;
  my $type = shift;

  throw("You must provide a type.") unless $type;

  return $self->cache->dump_path."/$type.exons.fasta";
}


sub dump_filtered_exons {
  my $self = shift;
  my $matrix = shift;

  unless ($matrix and
          $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('You must provide a ScoredMappingMatrix.');
  }

  # write exons to fasta files
  my $source_count = $self->write_filtered_exons('source', $matrix);
  my $target_count = $self->write_filtered_exons('target', $matrix);

  # return true if both source and target exons were written; otherwise we
  # don't need to run exonerate
  return (($source_count > 0) and ($target_count > 0));
}


sub write_filtered_exons {
  my $self = shift;
  my $type = shift;
  my $matrix = shift;

  throw("You must provide a type.") unless $type;
  unless ($matrix and
          $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('You must provide a ScoredMappingMatrix.');
  }

  $self->logger->info("\nDumping $type exons to fasta file...\n", 0, 'stamped');

  # don't dump exons shorter than this
  my $min_exon_length = $self->conf->param('min_exon_length') || 15;

  # counters
  my $total_exons = 0;
  my $dumped_exons = 0;

  # filehandle for fasta files
  my $fh;
  my $file = $self->exon_fasta_file($type);
  open($fh, '>', $file) or throw("Unable to open $file for writing: $!");

  # loop over exons, dump sequence to fasta file if longer than threshold and
  # score < 1
  EXON:
  foreach my $eid (sort { $b <=> $a }
                   keys %{ $self->cache->get_by_name('exons_by_id', $type) }) {

    my $exon = $self->cache->get_by_key('exons_by_id', $type, $eid);

    $total_exons++;

    # skip if exon shorter than threshold
    next EXON if ($exon->length < $min_exon_length);

    # skip if overlap score with any other exon is 1
    if ($type eq 'source') {
      foreach my $target (@{ $matrix->get_targets_for_source($eid) }) {
        next EXON if ($matrix->get_score($eid, $target) > 0.9999);
      }
    } else {
      foreach my $source (@{ $matrix->get_sources_for_target($eid) }) {
        next EXON if ($matrix->get_score($source, $eid) > 0.9999);
      }
    }

    # write exon to fasta file
    print $fh '>', $eid, "\n", $exon->seq, "\n";

    $dumped_exons++;

  }

  close($fh);

  # log
  my $fmt = "%-30s%10s\n";
  my $size = -s $file;
  $self->logger->info(sprintf($fmt, 'Total exons:', $total_exons), 1);
  $self->logger->info(sprintf($fmt, 'Dumped exons:', $dumped_exons), 1);
  $self->logger->info(sprintf($fmt, 'Dump file size:', parse_bytes($size)), 1);
  $self->logger->info("Done.\n\n", 0, 'stamped');

  return $dumped_exons;
}


sub parse_exonerate_results {
  my $self = shift;
  my $exonerate_matrix = shift;

  unless ($exonerate_matrix and
      $exonerate_matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('You must provide a ScoredMappingMatrix.');
  }

  $self->logger->info("Parsing exonerate results...\n", 0, 'stamped');

  # loop over all result files
  my $dump_path = $self->cache->dump_path;
  my $num_files = 0;
  my $num_lines = 0;

  opendir(DUMPDIR, $dump_path) or
    $self->logger->error("Can't open $dump_path for reading: $!");

  while (defined(my $file = readdir(DUMPDIR))) {
    next unless $file =~ /exonerate_map\.\d+/;

    $num_files++;

    open(F, '<', "$dump_path/$file");

    while (<F>) {
      $num_lines++;
      chomp;

      # line format:
      # myinfo: source_id target_id match_length source_length target_length
      my (undef, $source_id, $target_id, $match_length, $source_length,
          $target_length) = split;

      my $score = 0;

      if ($source_length == 0 or $target_length == 0) {
        $self->logger->warning("Alignment length is 0 for $source_id/$target_id.\n");
      } else {
        $score = 2 * $match_length / ($source_length + $target_length);
      }

      $exonerate_matrix->add_score($source_id, $target_id, $score);

    }

    close(F);
  }
  
  closedir(DUMPDIR);

  $self->logger->info("Done parsing $num_lines lines from $num_files result files.\n", 0, 'stamped');

  return $exonerate_matrix;
}


sub non_mapped_transcript_rescore {
  my $self = shift;
  my $matrix = shift;
  my $transcript_mappings = shift;

  # argument checks
  unless ($matrix and
      $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix of exons.');
  }

  unless ($transcript_mappings and
          $transcript_mappings->isa('Bio::EnsEMBL::IdMapping::MappingList')) {
    throw('Need a Bio::EnsEMBL::IdMapping::MappingList of transcripts.');
  }

  # create of lookup hash of mapped source transcripts to target transcripts
  my %transcript_lookup = map { $_->source => $_->target }
    @{ $transcript_mappings->get_all_Entries };

  my $i = 0;

  foreach my $entry (@{ $matrix->get_all_Entries }) {

    my @source_transcripts = @{ $self->cache->get_by_key(
      'transcripts_by_exon_id', 'source', $entry->source) };
    my @target_transcripts = @{ $self->cache->get_by_key(
      'transcripts_by_exon_id', 'target', $entry->target) };

    my $found_mapped = 0;

    TR:
    foreach my $source_tr (@source_transcripts) {
      foreach my $target_tr (@target_transcripts) {
        my $mapped_target = $transcript_lookup{$source_tr->id};

        if ($mapped_target and ($mapped_target == $target_tr->id)) {
          $found_mapped = 1;
          last TR;
        }
      }
    }
    
    unless ($found_mapped) {
      $matrix->set_score($entry->source, $entry->target, ($entry->score * 0.8));
      $i++;
    }
  }

  $self->logger->debug("Scored exons in non-mapped transcripts: $i\n", 1);
}


1;

