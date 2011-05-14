package Bio::Greg::Gorilla::ILSCounts;

use strict;
use File::Path;
use Bio::Greg::Codeml;
use FreezeThaw qw(freeze thaw cmpStr safeFreeze cmpStrHard);
use Data::Types qw(:all);

use base (
  'Bio::Greg::StatsCollectionUtils',
  'Bio::Greg::Hive::Process',
);

sub param_defaults {
  return {
  };
}

sub fetch_input {
  my ($self) = @_;

  # Fetch parameters from all possible locations.
  $self->load_all_params();

  # Create tables if necessary.
  $self->create_table_from_params( $self->compara_dba, 'sites',
    $self->sites_table );
  foreach my $chunk_width ($self->chunk_widths) {
    $self->create_table_from_params( $self->compara_dba, 'chunks_'.$chunk_width,
                                     $self->chunks_table );
  }
}

sub chunk_widths {
  my $self = shift;

  return (5, 10, 20, 50, 100);
}

sub run {
  my $self = shift;

  $self->param('data_id', $self->param('chr_name').$self->param('chr_start').$self->param('chr_end'));

  my $slice_a = Bio::EnsEMBL::Registry->get_adaptor( 'Human', 'core', 'slice' );
  my $chr_name = $self->param('chr_name');
  my $chr_start = $self->param('chr_start');
  my $chr_end = $self->param('chr_end');
  my $slice = $slice_a->fetch_by_region(undef, $chr_name, $chr_start, $chr_end);  
  
  my $aln = $self->_load_aln($slice);

  foreach my $width ($self->chunk_widths) {
    $self->_collect_chunks($slice, $aln, $width);
  }
  #$self->_collect_patterns($slice, $aln);
}

sub _collect_chunks {
  my $self = shift;
  my $slice = shift;
  my $aln = shift;
  my $chunk_width = shift;

  my $real_chunk_width = $chunk_width * 1_000;

  # Collect all genes in the vicinity.
  my @genes = @{$slice->get_all_Genes};
  @genes = grep {$_->biotype eq 'protein_coding'} @genes;
  @genes = map {$_->transform('chromosome')} @genes;
  @genes = sort {$a->start <=> $b->start} @genes;

  my @exons;
  my $exon_to_gene;
  foreach my $gene (@genes) {
    my @cur_exons = @{$gene->get_all_Exons};
    foreach my $e (@cur_exons) {
      $exon_to_gene->{$e} = $gene;
      push @exons, $e;
    }
  }
  @exons = sort {$a->start <=> $b->start} @exons;
  print "Genes: " . scalar(@genes)."\n";
  print "Exons: " . scalar(@exons)."\n";
  
  my $chunk_start = 1;
  while ($chunk_start + $slice->start < $slice->end) {
    my $chunk_end = $chunk_start + $real_chunk_width;
    $chunk_end = $aln->length if ($chunk_end > $aln->length);
        
    my $i = int($chunk_start + $real_chunk_width / 2);
    my $chr_mid = $slice->start + $i - 1;

    my $nearest_gene_distance = $slice->length;;
    my $nearest_gene;
    foreach my $gene (@genes) {
      #print "  ". $gene->external_name."  ".$gene->start."\n";
      my $s = $gene->start;
      my $e = $gene->end;
      my $distance = $chr_mid - $s;
      if (abs($chr_mid - $e) < abs($distance)) {
        # Closer to the end.
        $distance = $chr_mid - $e;
      }
      if ($chr_mid > $s && $chr_mid < $e) {
        # within the feature.
        $distance = 0;
      }
      if (!defined $nearest_gene_distance || abs($distance) < abs($nearest_gene_distance)) {
        $nearest_gene = $gene;
        $nearest_gene_distance = $distance;
      }
    }

    my $nearest_exon_distance;
    my $nearest_exon;
    foreach my $exon (@exons) {
      #print "  ". $gene->external_name."  ".$gene->start."\n";
      my $s = $exon->start;
      my $e = $exon->end;
      my $distance = $chr_mid - $s;
      if (abs($chr_mid - $e) < abs($distance)) {
        # Closer to the end coordinate.
        $distance = $chr_mid - $e;
      }
      if ($chr_mid > $s && $chr_mid < $e) {
        # within the feature.
        $distance = 0;
      }
      if (!defined $nearest_exon_distance || abs($distance) < abs($nearest_exon_distance)) {
        $nearest_exon = $exon;
        $nearest_exon_distance = $distance;
      }
    }


    my $gene_dnds = -1;
    my $gene_slr = -1;
    my $gene_name = '';
    my $gene_id = '';
    if (defined $nearest_gene) {
      $gene_name = $nearest_gene->external_name;
      $gene_id = $nearest_gene->stable_id;
    }

    my $exon_id = '';
    if (defined $nearest_exon) {
      $exon_id = $nearest_exon->stable_id;
    }

    my $sth = $self->dbc->prepare("SELECT * from gj1_gorilla.genes where gene_name='$gene_name' limit 1;");
    $sth->execute;
    while ( my $obj = $sth->fetchrow_hashref ) {
      $gene_dnds = $obj->{m0_dnds};
      $gene_slr = $obj->{slr_dnds};
    }

    my $patterns = $self->collect_patterns_from_aln($aln->slice($chunk_start, $chunk_end));

    if (defined $patterns) {
      my $params = $self->replace($patterns, {
        data_id => $self->param('data_id'),
        chr_name => $slice->seq_region_name,
        chr_start => $slice->start + $chunk_start,
        chr_end => $slice->start + $chunk_end,
        gene_dist => $nearest_gene_distance,      
        gene_id => $gene_id,
        gene_name => $gene_name,
        gene_dnds => $gene_dnds,
        gene_slr => $gene_slr,
        exon_dist => $nearest_exon_distance,
        exon_id => $exon_id
                                  });

      my $table_name = 'chunks_'.$chunk_width;
      $self->store_params_in_table($self->dbc, $table_name, $params);
      print "$nearest_gene_distance $nearest_exon_distance $gene_name $gene_dnds " . $patterns->{n_101}."\n";
    }

    $chunk_start += int($real_chunk_width/3);
  }
}

sub collect_patterns_from_aln {
  my $self = shift;
  my $aln = shift;

  my $pattern_hash;
  foreach my $i (0,1) {
    foreach my $j (0,1) {
      foreach my $k (0,1) {
        $pattern_hash->{'n_'.$i.$j.$k} = 0;
      }
    }
  }

  my @seqs = $aln->each_seq;

  my ($h_seq) = grep {$_->id == 9606} @seqs;
  my ($c_seq) = grep {$_->id == 9598} @seqs;
  my ($g_seq) = grep {$_->id == 9593} @seqs;
  my ($o_seq) = grep {$_->id == 9600} @seqs;
  return undef if (!defined $h_seq || !defined $c_seq || !defined $g_seq || !defined $o_seq);

  return undef if (scalar(@seqs) > 6);

  my $h_str = $h_seq->seq;
  my $c_str = $c_seq->seq;
  my $g_str = $g_seq->seq;
  my $o_str = $o_seq->seq;

  my $h = '';
  my $c = '';
  my $g = '';
  my $o = '';
  foreach my $i (1 .. $aln->length) {
    $h = substr($h_str, $i-1, 1);
    $c = substr($c_str, $i-1, 1);
    $g = substr($g_str, $i-1, 1);
    $o = substr($o_str, $i-1, 1);
    my $char = '-';
    next if ($h eq $char || $g eq $char || $c eq $char || $o eq $char);
    $char = 'N';
    next if ($h eq $char || $g eq $char || $c eq $char || $o eq $char);
    $char = '.';
    next if ($h eq $char || $g eq $char || $c eq $char || $o eq $char);
    $char = '';
    next if ($h eq $char || $g eq $char || $c eq $char || $o eq $char);
    
    my $pattern = join('', map {if ($_ eq $o) {0;} else {1;}} ($h, $c, $g));

    next if ($pattern eq '000' || $pattern eq '111');
    
    $pattern_hash->{'n_'.$pattern}++;
  }
  return $pattern_hash;
}

sub _collect_patterns {
  my $self = shift;
  my $slice = shift;
  my $aln = shift;

  # Collect all genes in the vicinity.
  my @genes = @{$slice->get_all_Genes};
  @genes = grep {$_->biotype eq 'protein_coding'} @genes;
  @genes = map {$_->transform('chromosome')} @genes;
  @genes = sort {$a->start <=> $b->start} @genes;

  my @seqs = $aln->each_seq;

  my ($h_seq) = grep {$_->id == 9606} @seqs;
  my ($c_seq) = grep {$_->id == 9598} @seqs;
  my ($g_seq) = grep {$_->id == 9593} @seqs;
  my ($o_seq) = grep {$_->id == 9600} @seqs;
  return if (!defined $h_seq || !defined $c_seq || !defined $g_seq || !defined $o_seq);
  my $h_str = $h_seq->seq;
  my $c_str = $c_seq->seq;
  my $g_str = $g_seq->seq;
  my $o_str = $o_seq->seq;

  my $h = '';
  my $c = '';
  my $g = '';
  my $o = '';
  foreach my $i (1 .. $aln->length) {
    $h = substr($h_str, $i-1, 1);
    $c = substr($c_str, $i-1, 1);
    $g = substr($g_str, $i-1, 1);
    $o = substr($o_str, $i-1, 1);
    my $char = '-';
    next if ($h eq $char || $g eq $char || $c eq $char || $o eq $char);
    $char = 'N';
    next if ($h eq $char || $g eq $char || $c eq $char || $o eq $char);
    $char = '.';
    next if ($h eq $char || $g eq $char || $c eq $char || $o eq $char);
    $char = '';
    next if ($h eq $char || $g eq $char || $c eq $char || $o eq $char);
    
    my $pattern = join('', map {if ($_ eq $o) {0;} else {1;}} ($h, $c, $g));

    next if ($pattern eq '000' || $pattern eq '111');

    my $chr_location = $slice->start + $i - 1;
    my $nearest_distance;
    my $nearest_gene;
    foreach my $gene (@genes) {
      #print "  ". $gene->external_name."  ".$gene->start."\n";
      my $s = $gene->start;
      my $e = $gene->end;
      my $distance = $chr_location - $s;
      if ($distance > 0) {
        if ($chr_location > $e) {
          $distance = $chr_location - $e;
        } else {
          # Site inside the gene.
          $distance = 0;
        }
      }
      if (!defined $nearest_distance || abs($distance) < abs($nearest_distance)) {
        $nearest_gene = $gene;
        $nearest_distance = $distance;
      }
    }

    printf "%-10s %-10s %-10s  %s\n", $chr_location, $nearest_distance, $nearest_gene->external_name, $pattern;

    my $params = {
      data_id => $self->param('data_id'),
      pattern => $pattern,
      chr_name => $slice->seq_region_name,
      chr_start => $chr_location,
      gene_dist => $nearest_distance,      
      gene_id => $nearest_gene->stable_id,
      gene_name => $nearest_gene->external_name
    };
    $self->store_params_in_table($self->dbc, 'sites', $params);
  }  

}

sub _load_aln {
  my $self = shift;
  my $slice = shift;
    
  my $cu = 'Bio::EnsEMBL::Compara::ComparaUtils';

  my $f = $self->_save_file('aln', 'fasta');
  if (!-e $f->{full_file} || $self->param('force_recalc')) {
    print "  fetching EPO alignment...\n";
    my $mirror_dba = $self->compara_dba;
    my $as_a       = $mirror_dba->get_AlignSliceAdaptor;
    my $mlss_a     = $mirror_dba->get_MethodLinkSpeciesSetAdaptor;    
    my $species_set = 'primates';
    my $mlss;
    if ($mlss_a->can('fetch_by_method_link_type_species_set_name')) {
      $mlss = $mlss_a->fetch_by_method_link_type_species_set_name( 'EPO', $species_set );
    } else {      
      my $mlss_list = $mlss_a->fetch_all_by_method_link_type('EPO');
      foreach my $cur_mlss (@$mlss_list) {
        my $name = $cur_mlss->name;
        $mlss = $cur_mlss if ($name =~ m/$species_set/gi);
      }
    }
    
    my $align_slice = $as_a->fetch_by_Slice_MethodLinkSpeciesSet( $slice, $mlss, 0 );
    my $as_aln = $align_slice->get_SimpleAlign;
    my $aln = new Bio::SimpleAlign;
    foreach my $a_s_slice (@{$align_slice->get_all_Slices}) {
      my $gdb       = $a_s_slice->genome_db;
      next if ( $gdb->name =~ m/ancestral/gi );
      my $bioseq = Bio::LocatableSeq->new( 
        -seq => $a_s_slice->seq,
        id => $gdb->taxon->taxon_id 
        );
      $aln->add_seq($bioseq);
    }
    
    if (0) {
      $aln = new Bio::SimpleAlign;
      my @a_s_slices = @{ $align_slice->get_all_Slices };
      
      foreach my $a_s_slice (@a_s_slices) {
        my $aln_start = $a_s_slice->start;
        my $aln_end   = $a_s_slice->end;
        my $gdb       = $a_s_slice->genome_db;
        
        next if ( $gdb->name =~ m/ancestral/gi );
        print $gdb->name."\n";
        
        my $cur_slice_adaptor =
          Bio::EnsEMBL::Registry->get_adaptor( $gdb->taxon->ensembl_alias, 'core', 'slice' );
        
        my @orig_positions;
        for ( my $i = 1 ; $i <= $a_s_slice->length ; $i++ ) {
          my ( $slice, $orig_position ) = $a_s_slice->get_original_seq_region_position($i);
          push @orig_positions, [ $slice, $orig_position ];
        }
        my $slice_ranges = 
          $cu->combine_slice_positions_into_ranges( \@orig_positions );
        
        my @slices = ();
        foreach my $slice_range (@$slice_ranges) {
          my $start = $slice_range->{start};
          my $end   = $slice_range->{end};
          my $slice = $slice_range->{slice};
          #print "$start $end ".$slice->name."\n";
          if ( $end < $start ) {
            my $tmp = $end;
            $end   = $start;
            $start = $tmp;
          }
          if ( $slice->seq_region_name =~ m/gap/i ) {
            push @slices, $slice;
          } else {
            my $sub_slice =
              $cur_slice_adaptor->fetch_by_region( undef, $slice->seq_region_name, $start, $end,
                                                   $slice->strand );
            push @slices, $sub_slice;
          }
        }
        #print "Slice ranges done\n";
        
        my $cache_object;
        my @quals = ();
        my $i=0;
        foreach my $p_slice (@slices) {
          $i++;
          my @cur_quals;
          if ( $p_slice->seq_region_name =~ m/gap/i ) {
            @cur_quals = (99) x $p_slice->length;
          } else {
            my ($dna,$qual) = 
              Bio::EnsEMBL::Compara::ComparaUtils->get_bases_for_slice( $p_slice, $gdb, $cache_object);
            @cur_quals = split( ' ', $qual );
          }
          push @quals, @cur_quals;
        }      
        
        my $quality_threshold = 30;
        my $seq_str = $a_s_slice->seq;
        my $bioseq = Bio::LocatableSeq->new( -seq => $seq_str, id => $gdb->taxon->taxon_id );
        
        if (scalar(@quals) > 0) {
          my $filtered =
            $cu->filter_alignment_seq_by_qual_array( $bioseq, \@quals, $quality_threshold );
          my $n_filtered = $cu->count_filtered_sites( $seq_str, $filtered );
          $self->store_param($gdb->name."_filtered", $n_filtered);
          $bioseq->seq($filtered);
        }
        $aln->add_seq($bioseq);
      }
    }
    
    $self->pretty_print($aln);
    Bio::EnsEMBL::Compara::AlignUtils->to_file($aln, $f->{full_file});  
  }
  print "  loading aln from file\n";
  my $aln = Bio::EnsEMBL::Compara::AlignUtils->from_file($f->{full_file});
  
  $self->pretty_print($aln);  
  return $aln;
}

sub sites_table {
  return {
    data_id => 'int',
    chr_name => 'char8',
    chr_start => 'int',
    pattern => 'char4',
    gene_dist => 'int',
    gene_id => 'char32',
    gene_name => 'char32',
#    gene_dnds => 'float',
#    recomb_rate => 'float',
    unique_keys => 'chr_name,chr_start',
    extra_keys => 'pattern'
  };
}

sub chunks_table {
  return {
    data_id => 'int',
    chr_name => 'char8',
    chr_start => 'int',
    chr_end => 'int',
    n_100 => 'int',
    n_010 => 'int',
    n_001 => 'int',
    n_110 => 'int',
    n_101 => 'int',
    n_011 => 'int',
    gene_dist => 'int',
    gene_id => 'char32',
    gene_name => 'char32',
    gene_dnds => 'float',
    gene_slr => 'float',
    exon_dist => 'int',
    exon_id => 'char32',
    unique_keys => 'chr_name,chr_start',
    extra_keys => 'gene_dnds,gene_dist'
  };
}

sub _save_file {
  my $self = shift;
  my $filename_base = shift;
  my $ext = shift;

  my $id = $self->param('chr_name').$self->param('chr_start').$self->param('chr_end');

  my $filename = "${id}_${filename_base}";

  my $file_params = {
    id => $id,
    filename => $filename,
    extension => $ext,
    subfolder => 'data',
  };

  my $file_obj = $self->save_file($file_params);

  if ($self->param('data_tarball') && !-e $file_obj->{full_file}) {
    # Load the file from the tarball into our temp directory.
    #my $tarball = $self->param('data_tarball');
    #my $tmp = $self->worker_temp_directory;
    #my $rel_f = $file_obj->{rel_file};

    #my $cmd = qq^tar -zxvf $tarball $rel_f^;
    #print $cmd."\n";
    #$file_obj->{full_file} = $tmp . $filename;
  }

  if (!defined $self->param('data_prefix')) {
    $self->store_param('data_prefix', $file_obj->{hash_folder});
  }

  return $file_obj;
}

1;
