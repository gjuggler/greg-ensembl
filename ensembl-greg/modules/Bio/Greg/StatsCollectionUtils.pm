package Bio::Greg::StatsCollectionUtils;

use strict;

use Bio::EnsEMBL::Compara::NestedSet;

our @ISA = qw();

sub max_path {
  my $class = shift;
  my $tree  = shift;
  my $max   = $tree->max_distance;
  return sprintf "%.3f", $max;
}

sub tree_length {
  my $class = shift;
  my $tree  = shift;
  my $dist  = 0;
  map { $dist += $_->distance_to_parent } $tree->nodes;
  return sprintf "%.3f", $dist;
}

sub mean_path {
  my $class = shift;
  my $tree  = shift;
  my $dist  = 0;
  map { $dist += $class->dist_to_root($_); } $tree->leaves;
  $dist = $dist / scalar( $tree->leaves );
  return sprintf "%.3f", $dist;
}

sub dist_to_root {
  my $class = shift;
  my $leaf  = shift;

  my $d = $leaf->distance_to_parent;
  my $p = $leaf->parent;
  while ($p) {
    $d += $p->distance_to_parent;
    $p = $p->parent;
  }
  return $d;
}

sub max_branch {
  my $class = shift;
  my $tree  = shift;

  my $max = 0;
  foreach my $node ( $tree->nodes ) {
    $max = $node->distance_to_parent if ( $node->distance_to_parent > $max );
  }
  return $max;
}

sub mean_branch {
  my $class = shift;
  my $tree  = shift;

  my $total = 0;
  foreach my $node ( $tree->nodes ) {
    $total += $node->distance_to_parent;
  }
  return $total / scalar( $tree->nodes );
}

# GJ 2010-03-16
# Calculate the per-species member count within this gene family. Counted for all
# species that have at least one member in the family.
sub mean_copy_count {
  my $class = shift;
  my $tree  = shift;

  my $taxon_hash;
  foreach my $leaf ( $tree->leaves ) {
    my $tx_id = $leaf->taxon_id;
    if ( !defined $taxon_hash->{$tx_id} ) {
      $taxon_hash->{$tx_id} = 1;
    } else {
      $taxon_hash->{$tx_id} += 1;
    }
  }

  my $sum = 0;
  foreach my $taxon ( keys %$taxon_hash ) {
    $sum += $taxon_hash->{$taxon};
  }
  my $num_taxa = scalar keys %$taxon_hash;

  return 0 if ( $num_taxa == 0 );
  my $mean = $sum / $num_taxa;

  return $mean;
}

sub duplication_count {
  my $class = shift;
  my $tree = shift;

  my $dup_sum = 0;
  my $node_sum = 0;
  foreach my $node ($tree->nodes) {
    my $is_duplication = $node->get_tagvalue("Duplication") || 0;
    my $is_dubious = $node->get_tagvalue("dubious_dup") || 0;
    if ($is_duplication && !$is_dubious) {
      $dup_sum++;
    }
    $node_sum++;
  }
  return $dup_sum;
}

sub duplication_fraction {
  my $class = shift;
  my $tree = shift;

  my $dup_sum = 0;
  my $node_sum = 0;
  foreach my $node ($tree->nodes) {
    my $is_duplication = $node->get_tagvalue("Duplication") || 0;
    my $is_dubious = $node->get_tagvalue("dubious_dup") || 0;
    if ($is_duplication && !$is_dubious) {
      $dup_sum++;
    }
    $node_sum++;
  }
  return $dup_sum / $node_sum;
}


sub seq_length_mean {
  my $class = shift;
  my $tree  = shift;

  my $seq_len = 0;
  map { $seq_len += $_->seq_length } $tree->leaves;
  return $seq_len / scalar( $tree->leaves );
}

sub gc_content_mean {
  my $class = shift;
  my $tr    = shift;

  my @seqs;
  if ( $tr =~ /ProteinTree/i ) {
    foreach my $leaf ( $tr->leaves ) {
      my $tx  = $leaf->transcript;
      my $seq = $tx->seq->seq;
      push @seqs, $seq;
    }
  } elsif ( $tr =~ /Align/i ) {
    foreach my $seq ( $tr->each_seq ) {
      push @seqs, $seq->seq;
    }
  }

  my $sum_gc = 0;
  foreach my $seq (@seqs) {
    $seq =~ s/[nx]//gi;

    my $total_len = length($seq);
    $seq =~ s/[at]//gi;
    my $gc_content = length($seq) / $total_len;
    $sum_gc += $gc_content;
  }
  my $avg_gc = $sum_gc / scalar(@seqs);
}

sub get_tag_hash {
  my $class  = shift;
  my $dbc    = shift;
  my $params = shift;

  # Index tags by the aln_position and tag.
  my $table   = "sitewise_tag";
  my $pset    = $params->{'parameter_set_id'};
  my $node_id = $params->{'node_id'};

  my $cmd = qq^SELECT aln_position,tag,value
    FROM $table WHERE ( parameter_set_id=$pset or parameter_set_id=0 ) and node_id=$node_id
    ^;

  my $tag_hash = {};

  my $sth = $dbc->prepare($cmd);
  $sth->execute;
  while ( my $obj = $sth->fetchrow_hashref ) {
    my $aln_position = $obj->{'aln_position'};
    $tag_hash->{$aln_position} = {} if (!defined $tag_hash->{$aln_position});
    $tag_hash->{$aln_position}->{$obj->{'tag'}} = $obj->{'value'};
  }

  return $tag_hash;
}

sub get_psc_hash {
  my $class  = shift;
  my $dbc    = shift;
  my $params = shift;

  my $table   = $params->{'omega_table'};
  my $pset    = $params->{'parameter_set_id'};
  my $node_id = $params->{'node_id'};
  
  my $CLEAN_WHERE = qq^
    AND o.note != 'random' AND o.omega_upper < 99
    ^;

  my $cmd     = qq^SELECT aln_position,omega,omega_lower,omega_upper,lrt_stat,ncod,type,note 
    FROM $table o WHERE parameter_set_id=$pset and node_id=$node_id $CLEAN_WHERE
    ^;

  if ( $params->{genome} ) {
    $cmd = qq^SELECT * from $table o, sitewise_genome g WHERE o.parameter_set_id=$pset AND 
      o.node_id=$node_id
      AND o.node_id=g.node_id
      AND o.aln_position=g.aln_position $CLEAN_WHERE^;
  }

  if ($params->{filtered} ) {
    my $filter_value = $params->{mammals_alignment_filtering_value} || 1;
    print "Filtering value: $filter_value\n";

    # Filter on alignment columns that pass Pollard et al's filtering criteria.
    $cmd = qq^SELECT * from $table o, sitewise_tag t  WHERE
      o.parameter_set_id=$pset AND 
      o.node_id=$node_id
      AND o.node_id=t.node_id
      AND o.aln_position=t.aln_position
      AND t.tag="FILTER" AND t.value >= $filter_value $CLEAN_WHERE;
      ^;
    # TODO: Join with the sitewise_tag table where tag = "FILTER".
  }

  print $cmd."\n";

  my $sth = $dbc->prepare($cmd);
  $sth->execute;
  my $id_field = 'aln_position';
  my $obj = $sth->fetchall_hashref($id_field);
  $sth->finish;
  return $obj;
}

sub max_lrt {
  my $class             = shift;
  my $psc_hash          = shift;

  my @obj_array = map { $psc_hash->{$_} } keys %$psc_hash;
  my @signed_lrt = map { if ($_->{omega} > 0) {$_->{lrt_stat}} else {-$_->{lrt_stat}}} @obj_array;
  
  my $max_lrt_stat = -10000;
  foreach my $lrt (@signed_lrt) {
    $max_lrt_stat = $lrt if ($lrt > $max_lrt_stat);
  }

  # Correct by the number of viable sites.
  #$max_lrt_stat = $max_lrt_stat / scalar(@signed_lrt);

  return $max_lrt_stat;
}

sub psc_count {
  my $class             = shift;
  my $psc_hash          = shift;
  my $include_weak_pscs = shift;

  my @obj_array = map { $psc_hash->{$_} } keys %$psc_hash;
  my @psc_objs;
  if ($include_weak_pscs) {
    @psc_objs = grep { $_->{type} =~ /positive[1234]/ } @obj_array;
  } else {
    @psc_objs = grep { $_->{type} =~ /positive[34]/ } @obj_array;
  }

  return scalar(@psc_objs);
}

sub sitewise_count {
  my $class             = shift;
  my $psc_hash          = shift;

  my @obj_array = map { $psc_hash->{$_} } keys %$psc_hash;
  return scalar(@obj_array);
}

sub omega_median {
  my $class = shift;
  my $hash  = shift;

  my @obj_array = map { $hash->{$_} } keys %$hash;

  my @omega_values = map {$_->{omega}} @obj_array;
  
  @omega_values = sort {$a <=> $b} @omega_values;

  my $median = $omega_values[scalar(@omega_values)/2];

  return 'NA' if ( scalar @obj_array == 0 );
  return sprintf "%.3f", $median;
}

sub omega_average {
  my $class = shift;
  my $hash  = shift;

  my @obj_array = map { $hash->{$_} } keys %$hash;

  my $omega_total = 0;
  map { $omega_total += $_->{omega} } @obj_array;

  return 'NA' if ( scalar @obj_array == 0 );
  return sprintf "%.3f", $omega_total / scalar(@obj_array);
}

sub omega_average_exclude_pscs {
  my $class = shift;
  my $hash  = shift;

  my @obj_array = map { $hash->{$_} } keys %$hash;
  @obj_array = grep { $_->{omega_lower} < 1 } @obj_array;

  my $omega_total = 0;
  map { $omega_total += $_->{omega} } @obj_array;

  return 'NA' if ( scalar @obj_array == 0 );
  return sprintf "%.3f", $omega_total / scalar(@obj_array);
}

sub site_count {
  my $class = shift;
  my $aln   = shift;

  my $site_count = 0;
  foreach my $seq ( $aln->each_seq ) {
    my $str = $seq->seq;
    $str =~ s/-//g;
    $site_count += length($str);
  }
  return $site_count;
}

sub unfiltered_site_count {
  my $class = shift;
  my $aln   = shift;

  my $site_count = 0;
  foreach my $seq ( $aln->each_seq ) {
    my $str = $seq->seq;
    $str =~ s/[-X]//g;
    $site_count += length($str);
  }
  return $site_count;
}

sub cpg_obs_exp_mean {
  my $class = shift;
  my $aln   = shift;

  my $seq;
  my $running_total;
  foreach my $seq_obj ( $aln->each_seq ) {
    $seq .= $seq_obj->seq;
  }

  $seq =~ s/[-X]//gi;    # Remove gaps and filtered sites.

  my $len = length($seq);
  my $gs  = $seq;
  $gs =~ s/[atc]//gi;
  my $cs = $seq;
  $cs =~ s/[atg]//gi;
  my $g_rate = length($gs) / length($seq);
  my $c_rate = length($cs) / length($seq);

  my $exp_count = $g_rate * $c_rate * $len;

  my $cpgs      = $seq;
  my @matches   = $cpgs =~ /(cg)/gi;
  my $obs_count = scalar(@matches);

  return $obs_count / $exp_count;
}

sub mysql_getval {
  my $class = shift;
  my $tree  = shift;
  my $cmd   = shift;

  my $dbc = $tree->adaptor->dbc;
  my $sth = $dbc->prepare($cmd);
  $sth->execute();
  my $val = @{ $sth->fetchrow_arrayref() }[0];
  $val = 'NA' unless ( defined $val );
  $sth->finish;
  return $val;
}

1;
