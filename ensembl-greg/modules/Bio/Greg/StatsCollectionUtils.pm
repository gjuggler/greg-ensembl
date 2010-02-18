package Bio::Greg::StatsCollectionUtils;

use strict;

use Bio::EnsEMBL::Compara::NestedSet;

our @ISA = qw();

sub max_path {
  my $class = shift;
  my $tree = shift;
  my $max = $tree->max_distance;
  return sprintf "%.3f", $max;
}

sub tree_length {
  my $class = shift;
  my $tree = shift;
  my $dist = 0;
  map {$dist += $_->distance_to_parent} $tree->nodes;
  return sprintf "%.3f", $dist;
}

sub mean_path {
  my $class = shift;
  my $tree = shift;
  my $dist = 0;
  map {$dist += $class->dist_to_root($_);} $tree->leaves;
  $dist = $dist / scalar($tree->leaves);
  return sprintf "%.3f", $dist;
}

sub dist_to_root {
  my $class = shift;
  my $leaf = shift;

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
  my $tree = shift;

  my $max = 0;
  foreach my $node ($tree->nodes) {
    $max = $node->distance_to_parent if ($node->distance_to_parent > $max);
  }
  return $max;
}

sub mean_branch {
  my $class = shift;
  my $tree = shift;

  my $total = 0;
  foreach my $node ($tree->nodes) {
    $total += $node->distance_to_parent;
  }
  return $total / scalar($tree->nodes);
}

sub seq_length_mean {
  my $class = shift;
  my $tree = shift;
  
  my $seq_len = 0;
  map {$seq_len += $_->seq_length} $tree->leaves;
  return $seq_len / scalar($tree->leaves);
}

sub gc_content_mean {
  my $class = shift;
  my $tr = shift;
  
  my @seqs;
  if ($tr =~ /ProteinTree/i) {
    foreach my $leaf ($tr->leaves) {
      my $tx = $leaf->transcript;
      my $seq = $tx->seq->seq;
      push @seqs, $seq;
    }
  } elsif ($tr =~ /Align/i) {
    foreach my $seq ($tr->each_seq) {
      push @seqs,$seq->seq;
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

sub get_psc_hash {
  my $class = shift;
  my $dbc = shift;
  my $params = shift;

  my $table = $params->{'omega_table'};
  my $pset = $params->{'parameter_set_id'};
  my $node_id = $params->{'node_id'};
  my $cmd = qq^SELECT aln_position,omega,omega_lower,omega_upper,lrt_stat,ncod,type,note 
    FROM $table WHERE parameter_set_id=$pset and node_id=$node_id
    AND ncod > 4 AND note != 'random' AND omega_upper > omega^;

  if ($params->{genome}) {
    $cmd = qq^SELECT * from $table o, sitewise_genome g WHERE o.parameter_set_id=$pset AND 
      o.node_id=$node_id AND ncod > 4 AND note != 'random' AND omega_upper > omega AND o.node_id=g.node_id
      AND o.aln_position=g.aln_position^;
  }

  print "$cmd\n";
  my $sth = $dbc->prepare($cmd);
  $sth->execute;
  my $obj = $sth->fetchall_hashref('aln_position');
  $sth->finish;
  return $obj;
}

sub psc_count {
  my $class = shift;
  my $psc_hash = shift;
  my $include_weak_pscs = shift;

  my @obj_array = map {$psc_hash->{$_}} keys %$psc_hash;
  my @psc_objs;
  if ($include_weak_pscs) {
    @psc_objs = grep {$_->{type} =~ /positive[1234]/} @obj_array;
  } else {  
    @psc_objs = grep {$_->{type} =~ /positive[34]/} @obj_array;
  }

  return scalar(@psc_objs);
}


sub omega_average {
  my $class = shift;
  my $hash = shift;

  my @obj_array = map {$hash->{$_}} keys %$hash;
  
  my $omega_total = 0;
  map {$omega_total += $_->{omega}} @obj_array;

  return 'NA' if (scalar @obj_array == 0);
  return sprintf "%.3f", $omega_total/scalar(@obj_array);
}

sub omega_average_exclude_pscs {
  my $class = shift;
  my $hash = shift;

  my @obj_array = map {$hash->{$_}} keys %$hash;
  @obj_array = grep {$_->{omega_lower} < 1} @obj_array;

  my $omega_total = 0;
  map {$omega_total += $_->{omega}} @obj_array;

  return 'NA' if (scalar @obj_array == 0);
  return sprintf "%.3f", $omega_total/scalar(@obj_array);
}

sub site_count {
  my $class = shift;
  my $aln = shift;

  my $site_count = 0;
  foreach my $seq ($aln->each_seq) {
    my $str = $seq->seq;
    $str =~ s/-//g;
    $site_count += length($str);
  }
  return $site_count;
}

sub unfiltered_site_count {
  my $class = shift;
  my $aln = shift;

  my $site_count = 0;
  foreach my $seq ($aln->each_seq) {
    my $str = $seq->seq;
    $str =~ s/[-X]//g;
    $site_count += length($str);
  }
  return $site_count;
}

sub mysql_getval {
  my $class = shift;
  my $tree = shift;
  my $cmd = shift;

  my $dbc = $tree->adaptor->dbc;
  my $sth = $dbc->prepare($cmd);
  $sth->execute();
  my $val = @{$sth->fetchrow_arrayref()}[0];
  $val = 'NA' unless (defined $val);
  $sth->finish;
  return $val;
}

1;
