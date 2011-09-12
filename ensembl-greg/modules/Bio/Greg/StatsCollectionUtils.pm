package Bio::Greg::StatsCollectionUtils;

use strict;

use Bio::EnsEMBL::Compara::NestedSet;

use Bio::Greg::EslrUtils;

our @ISA = qw();

sub dawg_lambda {
  my $self = shift;
  my $tree = shift;

  my $params = $self->params;
  $params->{tree} = $tree;
  my $aln = $self->get_cdna_aln($params);

  my $lambda =
    Bio::EnsEMBL::Compara::AlignUtils->dawg_lambda( $aln, $tree, {}, $self->worker_temp_directory );

  return $lambda;
}

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

sub species_with_dups {
  my $class = shift;
  my $tree  = shift;

  my @leaves = $tree->leaves;
  my $taxid_hash;
  map {$taxid_hash->{$_->taxon_id} = 1} @leaves;

  my @dup_list;
  foreach my $taxid (keys %$taxid_hash) {
    my $count = scalar(grep {$_->taxon_id == $taxid} @leaves);
    if ($count > 1) {
      push @dup_list, $taxid;
    }
  }
  return join(' ', @dup_list);
}

sub duplication_count {
  my $class = shift;
  my $tree  = shift;

  my @leaves = $tree->leaves;
  my $taxid_hash;
  map {$taxid_hash->{$_->taxon_id} = 1} @leaves;

  my $dup_count = 0;
  foreach my $taxid (keys %$taxid_hash) {
    my $count = scalar(grep {$_->taxon_id == $taxid} @leaves);
    if ($count > 1) {
      $dup_count++;
    }
  }
  return $dup_count;
}

sub duplication_fraction {
  my $class = shift;
  my $tree  = shift;

  my $dup_sum  = 0;
  my $node_sum = 0;
  foreach my $node ( $tree->nodes ) {
    my $is_duplication = $node->get_tagvalue("Duplication") || 0;
    my $is_dubious     = $node->get_tagvalue("dubious_dup") || 0;
    if ( $is_duplication && !$is_dubious ) {
      $dup_sum++;
    }
    $node_sum++;
  }
  return $dup_sum / $node_sum;
}

sub root_node_gene_count {
  my $class = shift;
  my $tree  = shift;

  #my $root_node_id = $class->root_node_id($tree);
  my $root_node_id = $class->param('orig_node_id') || $class->root_node_id($tree);

  my $pta       = $tree->adaptor;
  my $root_tree = $pta->fetch_node_by_node_id($root_node_id);
  return scalar( $root_tree->leaves );
}

sub root_node_id {
  my $class = shift;
  my $tree  = shift;

  my $root_id = $tree->node_id;
  my $cmd     = qq^SELECT parent.node_id FROM protein_tree_node parent, protein_tree_node child
               WHERE child.left_index BETWEEN parent.left_index AND parent.right_index
               AND child.node_id=$root_id^;
  return $class->mysql_getval( $tree, $cmd );
}

sub seq_length_mean {
  my $class = shift;
  my $tree  = shift;

  my $seq_len = 0;
  map { $seq_len += $_->seq_length } $tree->leaves;
  return $seq_len / scalar( $tree->leaves );
}

sub genomic_gc_content {
  my $self = shift;
  my $member = shift;

  my $gene = $member->get_Gene;
  my $slice = $gene->slice;
  my $full_chromosome;

  # Small spelling difference between Ensembl versions.
  if ($slice->can('seq_region_slice')) {
    $full_chromosome = $slice->seq_region_slice;
  } else {
    $full_chromosome = $slice->seq_region_Slice;    
  }

  my $start = $gene->seq_region_start;
  my $end = $gene->seq_region_end;
  my $strand = $gene->seq_region_strand;

  my $segment = $full_chromosome->sub_Slice($start,$end,$strand);

  my $seq = $segment->seq;

  $seq =~ s/[-nx]//gi;
  my $total_len = length($seq);
  my $at = $seq;
  my $gc = $seq;

  $at =~ s/[gc]//gi;
  $gc =~ s/[at]//gi;

  #die("Error!") if ($seq =~ m/[^gc]/i);
  my $gc_content = length($gc) / (length($gc) + length($at));
  
  return $gc_content;
}

sub gc3_content {
  my $self = shift;
  my $member = shift;

  my $seq = $member->sequence_cds;

  # Get 3rd position nucleotides.
  $seq =~ s/..(.)/$1/gi;

  # Remove nonexistent seqs & gaps.
  $seq =~ s/[nx-]//gi;

  my $total_len = length($seq);
  $seq =~ s/[at]//gi;
  my $gc_content = length($seq) / $total_len;

  return $gc_content;
}

sub gc_content {
  my $self = shift;
  my $member = shift;

  my $seq = $member->sequence_cds;

  $seq =~ s/[nx-]//gi;

  my $total_len = length($seq);
  $seq =~ s/[at]//gi;
  my $gc_content = length($seq) / $total_len;
 
  return $gc_content;
}


sub gc_content_mean {
  my $class = shift;
  my $tr    = shift;

  my @seqs;
  if ( $tr =~ /ProteinTree/i ) {
    foreach my $leaf ( $tr->leaves ) {
      my $seq = $leaf->sequence_cds;

      #my $seq = $tx->seq->seq;
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
  my $pset    = $params->{'parameter_set_id'} || 0;
  my $node_id = $params->{'node_id'};

  my $cmd = qq^SELECT aln_position,tag,value
    FROM $table WHERE ( parameter_set_id=$pset or parameter_set_id=0 ) and node_id=$node_id
    ^;
  print $cmd. "\n";
  my $tag_hash = {};

  my $sth = $dbc->prepare($cmd);
  $sth->execute;
  while ( my $obj = $sth->fetchrow_hashref ) {
    my $aln_position = $obj->{'aln_position'};

    $tag_hash->{$aln_position} = {} if ( !defined $tag_hash->{$aln_position} );
    $tag_hash->{$aln_position}->{ $obj->{'tag'} } = $obj->{'value'};
  }

  return $tag_hash;
}

sub get_gene_name {
  my $class    = shift;
  my $tree     = shift;
  my $taxon_id = shift;

  foreach my $member ( $tree->leaves ) {
    if ( $member->taxon_id == $taxon_id ) {

    }
  }
}

sub get_psc_hash {
  my $class  = shift;
  my $dbc    = shift;
  my $params = shift;

  my $table = $params->{'omega_table'}      || 'sitewise_omega';
  my $pset  = $params->{'parameter_set_id'} || '0';
  my $data_id              = $class->data_id;
  my $node_id              = $class->node_id;
  my $include_crappy_sites = $params->{'get_all_sites'};

  my $CLEAN_WHERE = qq^
    AND o.note != 'random' AND o.omega_upper < 99
    ^;
  if ($include_crappy_sites) {
    $CLEAN_WHERE = '';
  }

  my $cmd;

  if ( $params->{filtered} && $params->{genome} ) {
    print "Both!\n";
    my $filter_value = $params->{alignment_filtering_value} || 1;
    print "Filtering value: $filter_value\n";

    $cmd = qq^SELECT * from
$table o LEFT OUTER JOIN sitewise_genome g ON
  g.node_id=$node_id AND g.aln_position=o.aln_position AND g.parameter_set_id=0
LEFT OUTER JOIN sitewise_tag t ON
  t.node_id=$node_id AND t.aln_position=o.aln_position AND t.parameter_set_id=0
WHERE
  o.parameter_set_id=$pset AND o.data_id=$data_id
  AND t.tag="FILTER" AND t.value >= $filter_value
$CLEAN_WHERE;
      ^;
  } elsif ( $params->{filtered} ) {
    print "Filtered!\n";
    my $filter_value = $params->{alignment_filtering_value} || 1;

    # Filter on alignment columns that pass Pollard et al's filtering criteria.
    $cmd = qq^SELECT * from $table o, sitewise_tag t  WHERE
      o.parameter_set_id=$pset AND 
      (o.data_id=$data_id)
      AND t.node_id=$node_id
      AND o.aln_position=t.aln_position
      AND t.tag="FILTER" AND t.value >= $filter_value $CLEAN_WHERE;
      ^;
  } elsif ( $params->{genome} ) {
    print "Genome!\n";
    $cmd = qq^SELECT * from $table o, sitewise_genome g WHERE o.parameter_set_id=$pset AND 
      (o.data_id=$data_id)
      AND g.node_id=$node_id
      AND o.aln_position=g.aln_position $CLEAN_WHERE^;
  } else {
    $cmd = qq^SELECT aln_position,omega,omega_lower,omega_upper,lrt_stat,ncod,type,note 
    FROM $table o WHERE parameter_set_id=$pset and (data_id=$node_id) $CLEAN_WHERE
    ^;
  }

  print $cmd. "\n";

  my $sth = $dbc->prepare($cmd);
  $sth->execute;
  my $id_field = 'aln_position';
  my $obj      = $sth->fetchall_hashref($id_field);
  $sth->finish;

  printf "SITE COUNT: %d\n", scalar( keys(%$obj) );
  return $obj;
}

sub add_ungapped_branch_lengths {
  my $self = shift;
  my $tree = shift;
  my $pep_aln = shift;
  my $psc_hash = shift;

  print "  calculating ungapped branch length...\n";

  my $new_hash;
  foreach my $i (1 .. $pep_aln->length) {
    $new_hash->{$i} = $psc_hash->{$i} if (defined $psc_hash->{$i});
  }

  my @values;
#  my @values = Bio::EnsEMBL::Compara::AlignUtils->nongap_branch_lengths($tree, $pep_aln);

  my $use_r = 1;
  if ($use_r) {
    my $tmp = $self->worker_temp_directory;
    my $tree_f = "${tmp}/tree.nh";
    my $aln_f = "${tmp}/aln.fasta";

    Bio::EnsEMBL::Compara::TreeUtils->to_file($tree, $tree_f);
    Bio::EnsEMBL::Compara::AlignUtils->to_file($pep_aln, $aln_f);

    my $csr = Bio::Greg::EslrUtils->baseDirectory . "/scripts/collect_sitewise.R";
    my $cmd = qq^
source("${csr}")
library(phylosim)

sim <- PhyloSim()
readTree(sim, "${tree_f}")
readAlignment(sim, "${aln_f}")

sites <- add.nongap.bl(sites, sim\$.phylo, sim\$.alignment)
sites <- sites[order(sites[,'aln_position']),]

for (i in 1:nrow(sites)) {
  print(sites[i,'aln_position'])
  print(sites[i,'nongap.bl'])
}
^;
    @values = $self->_run_r_with_sites($new_hash,$cmd);
  }

  for (my $i=0; $i < scalar(@values); $i++) {
    my $pos = $values[$i];
    $i++;
    my $nongap_bl = $values[$i];

    #print "$pos $nongap_bl\n";
    $psc_hash->{$pos}->{nongap_bl} = $nongap_bl;
  }
  print "Done!\n";
}

sub calculate_fractions {
  my $self = shift;
  my $psc_hash = shift;

  return if ( scalar(keys %$psc_hash) == 0 );

  my $cmd = qq^
  sites <- subset(sites, note != 'random')

  sites[,'pval'] <- 1 - pchisq(sites[,'lrt_stat'],1)
  pos.sites <- subset(sites,omega > 1 & pval < 0.01)
  neg.sites <- subset(sites,omega < 1 & pval < 0.01)
  neutral.sites <- subset(sites,pval > 0.01)

  n <- nrow(sites)
  if (n == 0) { n <- 1 }

  print(c(n,nrow(pos.sites),nrow(neg.sites),nrow(neutral.sites)))
^;
  my @values = $self->_run_r_with_sites($psc_hash,$cmd);
  
  my $n = $values[0];

  $self->param('n_pos',$values[1]);
  $self->param('n_neg',$values[2]);
  $self->param('n_neut',$values[3]);
  $self->param('f_pos',$values[1]/$n);
  $self->param('f_neg',$values[2]/$n);
  $self->param('f_neut',$values[3]/$n);
}


sub combined_pval {
  my $self    = shift;
  my $psc_hash = shift;
  my $method   = shift || 'stouffer';

  return -2 if ( scalar(keys %$psc_hash) == 0 );

  my $combine_p = Bio::Greg::EslrUtils->baseDirectory . "/scripts/combine.p.R";
  my $cmd      = qq^
source("$combine_p")
sites <- subset(sites, note != 'random')

p.values = 1 - pchisq(sites[,'lrt_stat'],1)
p.values[p.values <= 0] = 1e-10
p.values[p.values >= 1] = 1
sites[,'pval'] = p.values

pos.sites = subset(sites,omega>1)
p.value = 999
if (nrow(pos.sites) > 0) {
  comb.p = combine.p(pos.sites[,'pval'],method='$method')
  p.value = as.numeric(comb.p[['p.value']])
} else if (nrow(sites) > 0) {
  p.value = -1
} else {
  p.value = -2
}
print(p.value)
^;
  my @values = $self->_run_r_with_sites($psc_hash,$cmd);
  my $pval = $values[0];
  return $pval;
}

sub _run_r_with_sites {
  my $self = shift;
  my $psc_hash = shift;
  my $cmd = shift;

  my @obj_array = map { $psc_hash->{$_} } keys %$psc_hash;

  my $pval;
  my $header = '';
  my $body = '';
  if ( scalar @obj_array > 0) {

    my $first_obj = $obj_array[0];
    my @keys      = sort keys %$first_obj;
    
    $header = join( "\t", @keys );
    foreach my $obj (@obj_array) {
      my $line = join( "\t", map { $obj->{$_} } @keys );
      $body .= $line . "\n";
    }
  }
  
  my $ug = new Data::UUID;
  my $filename = $ug->create_str;
  my $temp_f = $self->worker_temp_directory . $filename . ".txt";
  open( OUT, ">$temp_f" );
  print OUT $header . "\n";
  
  #print $header."\n";
  print OUT $body . "\n";
  
  #print $body."\n";
  close(OUT);
  
  $cmd = qq^
sites <- read.table(file="$temp_f",sep="\t",header=T)
^ . $cmd;
  
  my @values = Bio::Greg::EslrUtils->get_r_values( $cmd, $self->worker_temp_directory );
  return @values;
}

sub max_lrt {
  my $class    = shift;
  my $psc_hash = shift;

  my @obj_array = map { $psc_hash->{$_} } keys %$psc_hash;
  my @signed_lrt = map {
    if   ( $_->{omega} > 0 ) { $_->{lrt_stat} }
    else                     { -$_->{lrt_stat} }
  } @obj_array;

  my $max_lrt_stat = -10000;
  foreach my $lrt (@signed_lrt) {
    $max_lrt_stat = $lrt if ( $lrt > $max_lrt_stat );
  }

  # Correct by the number of viable sites.
  #$max_lrt_stat = $max_lrt_stat / scalar(@signed_lrt);

  return $max_lrt_stat;
}

sub psc_count {
  my $class           = shift;
  my $psc_hash        = shift;
  my $include_psc_num = shift;

  my @obj_array = map { $psc_hash->{$_} } keys %$psc_hash;
  my @psc_objs;

  if ( $include_psc_num == 1 ) {
    @psc_objs = grep { $_->{type} =~ /positive[1234]/ } @obj_array;
  } elsif ( $include_psc_num == 2 ) {
    @psc_objs = grep { $_->{type} =~ /positive[234]/ } @obj_array;
  } elsif ( $include_psc_num == 3 ) {
    @psc_objs = grep { $_->{type} =~ /positive[34]/ } @obj_array;
  } elsif ( $include_psc_num == 4 ) {
    @psc_objs = grep { $_->{type} =~ /positive[4]/ } @obj_array;
  } elsif ( $include_psc_num == -1 ) {
    @psc_objs = grep { $_->{type} =~ /negative[1234]/ } @obj_array;
  } elsif ( $include_psc_num == -2 ) {
    @psc_objs = grep { $_->{type} =~ /negative[234]/ } @obj_array;
  } elsif ( $include_psc_num == -3 ) {
    @psc_objs = grep { $_->{type} =~ /negative[34]/ } @obj_array;
  } elsif ( $include_psc_num == -4 ) {
    @psc_objs = grep { $_->{type} =~ /negative[4]/ } @obj_array;
  }

  return scalar(@psc_objs);
}

sub sitewise_count {
  my $class    = shift;
  my $psc_hash = shift;

  my @obj_array = map { $psc_hash->{$_} } keys %$psc_hash;
  return scalar(@obj_array);
}

sub omega_median {
  my $class = shift;
  my $hash  = shift;

  my @obj_array = map { $hash->{$_} } keys %$hash;

  my @omega_values = map { $_->{omega} } @obj_array;

  @omega_values = sort { $a <=> $b } @omega_values;

  my $median = $omega_values[ scalar(@omega_values) / 2 ];

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

sub get_hg18_from_hg19 {
  my $self   = shift;
  my $chr    = shift;
  my $pos    = shift;
  my $strand = shift;

  # Get the dba for the reference species.
  my $alias = 'Human';
  my $dba   = Bio::EnsEMBL::Registry->get_DBAdaptor( $alias, 'core' );
  my $asma  = $dba->get_AssemblyMapperAdaptor;
  my $csa   = $dba->get_CoordSystemAdaptor;

  my $new_cs = $csa->fetch_by_name('chromosome');
  my $old_cs = $csa->fetch_by_name( 'chromosome', 'NCBI36' );

  my $asm_mapper = $asma->fetch_by_CoordSystems( $new_cs, $old_cs );

  my @coords = $asm_mapper->map( $chr, $pos, $pos, $strand, $new_cs );
  my ($coord) = @coords;
  if ($coord) {
    return $coord->start;
  }  
}

sub get_coords_from_pep_position {
  my $self   = shift;
  my $member = shift;
  my $pos    = shift;

  my $ref_tx = $member->get_Transcript;

  # Get the dba for the reference species.
  my $alias = $member->taxon->ensembl_alias;
  my $dba   = Bio::EnsEMBL::Registry->get_DBAdaptor( $alias, 'core' );
  my $asma  = $dba->get_AssemblyMapperAdaptor;
  my $csa   = $dba->get_CoordSystemAdaptor;

  if ($member->taxon_id != 9606) {
    warn("Non-human ref, not collecting coordinates!\n");
    return {
    };
  }

  my $new_cs = $csa->fetch_by_name('chromosome');
  my $old_cs = $csa->fetch_by_name( 'chromosome', 'NCBI36' );

  my $asm_mapper = $asma->fetch_by_CoordSystems( $new_cs, $old_cs );

  my @genomic1 = $ref_tx->pep2genomic( $pos, $pos );
  foreach my $coord (@genomic1) {
    next unless ($coord->isa('Bio::EnsEMBL::Mapper::Coordinate'));
    my $chr    = $ref_tx->seq_region_name;
    my $start  = $coord->start;
    my $end    = $coord->end;
    my $strand = $coord->strand;
    
    my $obj = {
      chr_name => $chr,
      chr_start => $start,
      chr_end => $end,
      chr_strand => $strand,
      chr => $chr,
      strand => $strand,
      'hg19_name' => $chr,
      'hg19_pos' => $start,
    };

    # Map to old coordinates.
    eval {
      my @old_coords = $asm_mapper->map( $chr, $start, $end, $strand, $new_cs );
      my ($old_coord) = @old_coords;
      if ($old_coord) {
        $obj = $self->replace($obj,{
          'hg19_name' => $chr,
          'hg18_name' => $chr,
          'hg19_pos'  => $start,
          'hg18_pos'  => $coord->start,
                              });
      }
    };
    return $obj;
  }

  warn("No coordinates found!\n");
  return {};
}

sub get_ref_member {
  my $self = shift;
  my $tree = shift;

  my $ref_member;
  my @members = $tree->leaves;

  # Try for human first.
  ($ref_member) = grep { $_->taxon_id == 9606 } @members;
  
  # Then mouse.
  ($ref_member) = grep { $_->taxon_id == 10090 } @members if (!defined $ref_member);

  # Then dog.
  ($ref_member) = grep { $_->taxon_id == 9615 } @members if (!defined $ref_member);

  # Then, take whatever we can get.
  ($ref_member) = $members[0] if (!defined $ref_member);
  
  $self->param('ref_member_id',$ref_member->stable_id);
  return $ref_member;
}

sub standard_deviation {
  my $self = shift;
  my $numbers_hashref = shift;

  my @numbers = @{$numbers_hashref};

  #Prevent division by 0 error in case you get junk data
  return -1 unless(scalar(@numbers));

  # Step 1, find the mean of the numbers
  my $total1 = 0;
  foreach my $num (@numbers) {
    $total1 += $num;
  }
  my $mean1 = $total1 / (scalar @numbers);

  # Step 2, find the mean of the squares of the differences
  # between each number and the mean
  my $total2 = 0;
  foreach my $num (@numbers) {
    $total2 += ($mean1-$num)**2;
  }
  my $mean2 = $total2 / (scalar @numbers);

  # Step 3, standard deviation is the square root of the
  # above mean
  my $std_dev = sqrt($mean2);
  return $std_dev;
}


1;
