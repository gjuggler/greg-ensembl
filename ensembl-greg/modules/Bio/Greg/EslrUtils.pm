package Bio::Greg::EslrUtils;

use strict;

sub defaultMysqlURL {
  my $class = shift;

  if ($ENV{'USER'} =~ /gj1/) {
    my $url = 'mysql://ensadmin:ensembl@ens-research:3306/gj1_2xmammals_gold';
    print $url."\n";
    return $url;
    # return 'mysql://ensadmin:ensembl@compara2:5316/avilella_compara_homology_53';
  } else {
    my $url = 'mysql://ensadmin:ensembl@127.0.0.1:5425/gj1_2xmammals_gold';
    print $url."\n";
    return $url;
  }
}

sub urlFromConnection {
  my $class = shift;
  my $dbc = shift;

  my $url = sprintf("mysql://%s:%s@%s:%s/%s",
                    $dbc->username,
                    $dbc->password,
                    $dbc->host,
                    $dbc->port,
                    $dbc->dbname);
  return $url;
}

sub loadConfiguration {
  my $class = shift;
  my $params = shift;
  my $file = shift;

  my $obj = {};

  foreach my $key (keys %{$params}) {
    $obj->{$key} = $params->{$key};
  }

  open(IN,$file);
  my @lines = <IN>;
  foreach my $line (@lines) {
    next if $line =~ /^#/;
    
    $line =~ /^(.*?):(.*)/;
    my $key = $1;
    my $val = $2;
    $key =~ s/^\s*//g;
    $val =~ s/^\s*//g;
#    $val =~ s/'/"/g;
    $obj->{$key} = $val;
  }
  close(IN,$file);

  return $obj;
}


sub mapSitewiseToGenome {
  my $class = shift;
  my $tree = shift;
  my $taxon_id = shift;
  my $pos_value = shift; # Hashref where key=position, value=value.

  my $sa = $tree->get_SimpleAlign;
  my $aln_len = $sa->length;

  my @leaves = $tree->leaves;
  my @genomic_coords = ();
  foreach my $leaf (@leaves) {
    next unless ($leaf->taxon_id == $taxon_id);
    eval {
      my ($seq) = $sa->each_seq_with_id($leaf->stable_id);
      my $seq_str = $seq->seq;

      my $tscr = $leaf->get_Transcript;
      print STDERR $tscr->stable_id."\n";
      $tscr = $tscr->transform("chromosome");
      next unless (defined $tscr);
      my $chr = "chr".$tscr->slice->seq_region_name;

      foreach my $i (1..$aln_len) {
        my $char = substr($seq_str,$i-1,1);
        
        my $loc = $seq->location_from_column($i);
        next if (!defined $loc || $loc->location_type() eq 'IN-BETWEEN');
        
        #next if ($char eq '-');
        my @gc_arr = $tscr->pep2genomic($loc->start,$loc->end);
        my $gc = $gc_arr[0];
        next unless ($gc && $gc->isa("Bio::EnsEMBL::Mapper::Coordinate"));

        my $strand = "+";
        $strand = "-" if ($gc->strand == -1);

        my $value = -1;
        $value = $pos_value->{$i} if (exists $pos_value->{$i});
        next if ($value == -1);

        my $start = $gc->start;
        my $end = $gc->end;

        push (@genomic_coords, {
          chr => $chr,
          start => $start,
          end => $end,
          value => $value,
          aln_position => $i,
          char => $char,
          stable_id => $leaf->stable_id,
          node_id => $tree->node_id,
          member_id => $leaf->dbID,
          strand => $strand
        });
      }
    };
    warn() if $@;
  }
  return \@genomic_coords;
}

sub collectDuplicationTags {
  my $class = shift;
  my $tree = shift;
  my $params = shift;

  my $taxon_name = $tree->get_tagvalue('taxon_name');
  my $taxon_id = $tree->get_tagvalue('taxon_id');

  my $sth = $tree->adaptor->prepare("SELECT name FROM ncbi_taxa_name WHERE taxon_id=$taxon_id AND name_class='ensembl timetree mya';");
  $sth->execute;
  my $taxon_mya = $sth->fetchrow_array();

  my $total_bl = sprintf "%.3f", total_distance($tree);
  my $max_bl = sprintf "%.3f", $tree->max_distance;

  # Collect chromosome on human gene.
  my $hum_chr = '';
  my @hum_genes = grep {$_->taxon_id == 9606} $tree->leaves;
  if (scalar(@hum_genes) > 0) {
    my $hum_gene = $hum_genes[0];
    my $gene = $hum_gene->get_Gene;
    $hum_chr = $gene->slice->seq_region_name;
  }

  my $tags = {
    taxon_name => $taxon_name,
    taxon_id => $taxon_id,
    taxon_mya => $taxon_mya,
    bl_total => $total_bl,
    bl_max => $max_bl,
    human_chr => $hum_chr
  };
  my $prefix = 'dupldiv';
  my $mapped_tags;
  map {$mapped_tags->{$prefix.'_'.$_} = $tags->{$_}} keys %{$tags};
  return $mapped_tags;
}

sub collectGeneTags {
  my $class = shift;
  our $tree = shift;
  our $params = shift;

  $tree->re_root;
  my $hash;

  print "Node ID: ".$tree->node_id."\n";

  sub mysql_getval {
    my $cmd = shift;
    my $dbc = $tree->adaptor->dbc;
    my $sth = $dbc->prepare($cmd);
    $sth->execute();
    my $val = @{$sth->fetchrow_arrayref()}[0];
    $val = 'NA' unless (defined $val);
    return $val;
  }

  sub mysql_array {
    my $cmd = shift;
    my $dbc = shift;
    my $sth = $dbc->prepare($cmd);
    $sth->execute;
    my $array_ref = $sth->fetchall_arrayref([0]);
    my @vals = @{$array_ref};
    @vals = map {@{$_}[0]} @vals;  # Some weird mappings to unpack the numbers from the arrayrefs.
    $sth->finish;
    return @vals;
  }

  sub psc_hash {
    my $tree = shift;
    my $table = shift;
    my $pset = shift;
    my $weak = shift;

    my $psc_str = qq^("positive3","positive4")^;
    $psc_str = qq^("positive1","positive2","positive3","positive4")^ if ($weak);
    

    my $node_id = $tree->node_id;
    my $cmd = qq^SELECT aln_position FROM $table sa WHERE node_id=$node_id AND parameter_set_id=$pset
      AND omega_upper > omega AND type != 'random'
      AND ncod >= 4 AND type IN $psc_str;^;
    my @vals = mysql_array($cmd,$tree->adaptor->dbc);
    my $return_hash;
    map {$return_hash->{$_} = 1} @vals;
    return $return_hash;
  }

  sub psc_clusters {
    my $tree = shift;
    my $sa = shift;
    my $table = shift;
    my $pset = shift;
    my $weak = shift;
    my $dbl = shift;

    my $hash = psc_hash($tree,$table,$pset,$weak);
    return 0 unless ($hash);
    my @pscs = keys %{$hash};
    @pscs = sort {$a <=> $b} @pscs;
 #   print "@pscs\n";
    my $len = $sa->length;
    my @sites = (0) x $len;
    
    my $width = $len / 10;
    $width = 5 if ($width < 5);
    $width = 50 if ($width > 50);

    for (my $i=0; $i < scalar(@pscs); $i++) {
      my $psc = $pscs[$i];
      my $lo = $psc - $width;
      my $hi = $psc + $width;
      $lo = 1 if ($lo < 1);
      $hi = $len if ($hi > $len);

      if (!$dbl) {
	my $len = ($hi-$lo);
	splice(@sites,$psc,$len,(1)x$len);
      } else {
	if ($i > 0 && $pscs[$i-1] >= $lo) {
	  my $len = $psc-$pscs[$i-1];
	  splice(@sites,$pscs[$i-1],$len,(1)x$len);
	}
	if ($i < scalar(@pscs)-1 && $pscs[$i+1] <= $hi) {
	  my $len = $pscs[$i+1] - $psc;
	  splice(@sites,$psc,$len,(1)x$len);
	}
      }
    }

    my $str = join("",@sites);
#    print $str."\n";

    my @toks = split(/0+/,$str);
    my $num_clusters = 0;
    foreach my $tok (@toks) {
      $num_clusters++ if (length($tok) > 0);
    }
#    print "COUNT: $num_clusters\n";
    return $num_clusters;
  }

  sub avg_sitewise {
    my $value = shift;
    my $node_id = shift;
    my $table = shift;
    my $parameter_set_id = shift;

    my $cmd = qq^SELECT avg($value) FROM $table 
      WHERE node_id=$node_id AND ncod >= 4 AND type != 'random' AND omega_upper > omega
      AND parameter_set_id=$parameter_set_id
      ^;
    return sprintf("%.4f",mysql_getval($cmd));
  }

  sub num_pscs {
    my ($node_id,$table,$pset) = @_;
    return mysql_getval(qq^ SELECT count(*) FROM $table sa
			WHERE sa.node_id=$node_id AND sa.parameter_set_id=$pset
			AND sa.omega_upper > sa.omega AND sa.type != 'random' AND
			sa.type IN ("positive4","positive3") AND
			sa.ncod >= 4;
			^);
  }

  sub weak_pscs {
    my ($node_id,$table,$pset) = @_;
    return mysql_getval(qq^ SELECT count(*) FROM $table sa
			WHERE sa.node_id=$node_id AND sa.parameter_set_id=$pset
			AND sa.omega_upper > sa.omega AND sa.type != 'random' AND
			sa.type IN ("positive1","positive2","positive3","positive4") AND
			sa.ncod >= 4;
			^);
  }
  
  
  sub gc_content {
    my $tr = shift;

    my $sum_gc = 0;
    foreach my $leaf ($tr->leaves) {
      my $tx = $leaf->transcript;
      my $seq = $tx->seq->seq;
      $seq =~ s/[nx]//gi;

      my $total_len = length($seq);
      $seq =~ s/[at]//gi;
      my $gc_content = length($seq) / $total_len;
      $sum_gc += $gc_content;
    }
    my $avg_gc = $sum_gc / scalar($tr->leaves);
  }
  
  my $node_id = $tree->node_id;
  my $sw = "sitewise_omega";

  $hash->{'leaf_count'} = scalar($tree->leaves);
  $hash->{'node_count'} = scalar($tree->nodes);
  $hash->{'gappiness'} = mysql_getval("SELECT gappiness($node_id)");
  $hash->{'duplications'} = mysql_getval("SELECT num_dups_under_node($node_id)");
  $hash->{'duplication_fraction'} = mysql_getval("SELECT num_dups_under_node($node_id)/node_count($node_id)");
  my @hum_gen = grep {$_->taxon_id==9606} $tree->leaves;
  $hash->{'human_genes'} = scalar(@hum_gen);
  $hash->{'ensp'} = $hum_gen[0]->stable_id if ($hum_gen[0]);
  $hash->{'tree_length_total'} = sprintf "%.3f", total_distance($tree);
  $hash->{'tree_length_max'} = sprintf "%.3f", max_distance($tree);
  $hash->{'tree_length_avg'} = sprintf "%.3f", avg_distance($tree);
  $hash->{'avg_gc'} = gc_content($tree);

  # Alignment stats.
  my $sa = $tree->get_SimpleAlign;
  $hash->{'aln_length'} = sprintf "%d", $sa->length;
  $hash->{'aln_percent_identity'} = sprintf "%.3f", $sa->percentage_identity;

  # Avg seq length.
  my $seq_len=0;
  map {$seq_len += $_->seq_length} $tree->leaves;
  $hash->{'avg_seq_length'} = sprintf "%.3f", $seq_len / scalar($tree->leaves);

  my @param_names = mysql_getval("SELECT parameter_value FROM parameter_set where parameter_name='name' ORDER BY parameter_set_id;");
  my @param_ids = mysql_getval("SELECT parameter_set_id FROM parameter_set where parameter_name='name' ORDER BY parameter_set_id;");
  my %param_hash = zip(@param_names,@param_ids);
  
  foreach my $name (keys %param_hash) {
    my $cl = substr($name,0,2);
    my $num = $param_hash{$name};

    $hash->{$cl.'_lrt'} = avg_sitewise("lrt_stat",$node_id,$sw,$num);
    $hash->{$cl.'_omega'} = avg_sitewise("omega",$node_id,$sw,$num);
    $hash->{$cl.'_pscs'} = num_pscs($node_id,$sw,$num);
    $hash->{$cl.'_weak_pscs'} = weak_pscs($node_id,$sw,$num);

    $hash->{$cl.'_num_clusters'} = psc_clusters($tree,$sa,$sw,$num,0,0);
    $hash->{$cl.'_num_clusters_dbl'} = psc_clusters($tree,$sa,$sw,$num,0,1);
    $hash->{$cl.'_num_weak_clusters'} = psc_clusters($tree,$sa,$sw,$num,1);

    my $param_set = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_param_set($tree->adaptor->dbc,$num);
    $hash->{$cl.'_bl_total'} = '';
    eval {
      $hash->{$cl.'_bl_total'} = subtree_total($tree,$param_set);
    };
    $hash->{$cl.'_leaves'} = subtree_leaves($tree,$param_set);
  }

  my $prefix = "eslr";
  $prefix = $params->{'tag_prefix'} if ($params->{'tag_prefix'});

  my $hash2;
  map {$hash2->{$prefix.'_'.$_} = $hash->{$_}} keys %{$hash};

  my @sorted_keys = sort keys %{$hash2};
  my @sorted_vals = map {"  ".$_ . "=>".$hash2->{$_}} @sorted_keys;
  my $val_string = join("\n",@sorted_vals);
  print "Tags: {\n". $val_string."\n}\n";

  return $hash2;
}


sub avg_distance {
  my $tree = shift;
  my $dist = 0;
  map {$dist += dist_to_root($_);} $tree->leaves;
  #map {$dist += dist_to_root($_);print dist_to_root($_)."\n";} $tree->leaves;
  $dist = $dist / scalar($tree->leaves);
  #print $tree->newick_format."\n";
  #print " -> ".$dist."\n";
  return sprintf "%.4f", $dist;
}

sub dist_to_root {
  my $leaf = shift;

  my $d = $leaf->distance_to_parent;
  my $p = $leaf->parent;
  while ($p) {
    $d += $p->distance_to_parent;
#    print $p . "  " . $p->node_id."  ".$p->distance_to_parent."\n";
    $p = $p->parent;
  }
  return $d;
}

sub total_distance {
  my $tree = shift;
  my $dist = 0;
  map {$dist += $_->distance_to_parent} $tree->nodes;
  return sprintf "%.4f", $dist;
}

sub max_distance {
  my $tree = shift;
  my $max = $tree->max_distance;
  return sprintf "%.4f", $max;
}

sub subtree_avg {
  my $tree = shift;
  my $params = shift;

  $params->{'node_id'} = $tree->node_id;
  my $new_tree = Bio::EnsEMBL::Compara::ComparaUtils->get_tree_for_comparative_analysis($tree->adaptor->db,$params);
  $new_tree->print_tree;
  return avg_distance($new_tree);
}

sub subtree_total {
  my $tree = shift;
  my $params = shift;

  $params->{'node_id'} = $tree->node_id;
  my $new_tree = Bio::EnsEMBL::Compara::ComparaUtils->get_tree_for_comparative_analysis($tree->adaptor->db,$params);
  return total_distance($new_tree);
}

sub subtree_leaves {
  my $tree = shift;
  my $params = shift;

  $params->{'node_id'} = $tree->node_id;
  my $new_tree = Bio::EnsEMBL::Compara::ComparaUtils->get_tree_for_comparative_analysis($tree->adaptor->db,$params);
  return scalar($new_tree->leaves);
}

sub subtree_max {
  my $tree = shift;
  my $params = shift;

  $params->{'node_id'} = $tree->node_id;
  my $new_tree = Bio::EnsEMBL::Compara::ComparaUtils->get_tree_for_comparative_analysis($tree->adaptor->db,$params);
  return max_distance($new_tree);
}


sub run_r {
  my $class = shift;
  my $rcmd = shift;
  my $params = shift;

  my $temp_dir = "/tmp/eslr_rcmd";
  use File::Path;
  mkpath($temp_dir);

  my $temp_in = $temp_dir . "/temp_in.txt";
  #my $temp_out = $temp_dir . "/temp_out.txt";

  open(OUT,">$temp_in");
  print OUT $rcmd."\n";
  close(OUT);
  
  my $vanilla = "--vanilla";
  $vanilla = "--slave" if ($params->{'silent'});

  my $r_cmd = "R";
  if ($params->{'farm'}) {
    $r_cmd = "/software/bin/R-2.9.0";
  } elsif ($params->{'bigmen'}) {
    $r_cmd = "bsub -Is -R'select[mem>10000] rusage[mem=10000]' -M10000000 /software/R-2.9.0/bin/R ";
  } elsif ($ENV{'USER'} =~ /gj1/) {
    $r_cmd = "/software/R-2.9.0/bin/R";
  } else {

  }

  my $rc = system("$r_cmd $vanilla < $temp_in");
  die "R returned an error!" if ($rc);

  unlink($temp_in);
  unlink($temp_dir);
}

sub get_r_values {
  my $class = shift;
  my $rcmd = shift;
  my $temp_dir = shift;

  if (!$temp_dir) {
    $temp_dir = "/tmp/eslr_rcmd";
    use File::Path;
    mkpath($temp_dir);
  }

  use Digest::MD5 qw(md5_hex);
  my $digest = md5_hex($rcmd . rand());
  $digest = substr($digest,0,10);

  my $temp_in = $temp_dir . "/temp_in_".$digest.".txt";
  my $temp_out = $temp_dir . "/temp_out".$digest.".txt";

#  print "TEMP IN: $temp_in\n";
  open(OUT,">$temp_in");
  print OUT $rcmd."\n";
  close(OUT);
  
# cmd to run: /software/R-2.7.1/bin/R CMD BATCH $filename
  my $rc = system("/ebi/research/software/Linux_x86_64/bin/R-2.7.0 --slave --vanilla < $temp_in > $temp_out");
  #my $rc = system("/software/R-2.7.1/bin/R --vanilla --slave < $temp_in ");
  
  open(IN,"$temp_out");
  my @lines = <IN>;
  map {chomp} @lines; 
  close(IN);

  if ($rc) {
    print join("\n",@lines);
    die "R returned an error!";
  }
  
  unlink($temp_in);
  unlink($temp_out);
  #unlink($temp_dir);

  return @lines;
}

sub mysql_array {
  my $class = shift;
  my $dbc = shift;
  my $cmd = shift;

  my $sth = $dbc->prepare($cmd);
  $sth->execute();

  my $array_ref = $sth->fetchall_arrayref([0]);
  $sth->finish;
  my @vals = @{$array_ref};
  @vals = map {@{$_}[0]} @vals;  # Some weird mappings to unpack the numbers from the arrayrefs.
  return @vals;
}

sub get_short_name_for_parameter_set {
  my $class = shift;
  my $parameter_set_id = shift;

  my $id = $parameter_set_id;
  my $str = "";
  $str = "vertebrates" if ($id == 1);
  $str = "no-2x" if ($id == 2);
  $str = "only-2x" if ($id == 3);
  $str = "no-seq-f" if ($id == 4);
  $str = "no-aln-f" if ($id == 5);
  $str = "no-f" if ($id == 6);
  $str = "primates" if ($id == 7);
  $str = "laurasiatheria" if ($id == 9);
  $str = "glires" if ($id == 8);
  return ucfirst($str);
}

1;
