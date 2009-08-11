package Bio::Greg::EslrUtils;

use strict;

sub defaultMysqlURL {
  my $class = shift;

  if ($ENV{'USER'} =~ /gj1/) {
    my $url = 'mysql://ensadmin:ensembl@ens-research:3306/gj1_compara_54';
    print $url."\n";
    return $url;
    # return 'mysql://ensadmin:ensembl@compara2:5316/avilella_compara_homology_53';
  } else {
    my $url = 'mysql://ensadmin:ensembl@127.0.0.1:5425/gj1_compara_54';
    print $url."\n";
    return $url;
  }
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


sub collectGeneTags {
  my $class = shift;
  my $tree = shift;

  $tree->re_root;

  my $hash;

  sub mysql_getval {
    my $cmd = shift;
    my $dbc = $tree->adaptor->dbc;
    my $sth = $dbc->prepare($cmd);
    $sth->execute();
    my $val = @{$sth->fetchrow_arrayref()}[0];
    $val = 'NA' unless (defined $val);
    return $val;
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

  $hash->{'leaf_count'} = scalar($tree->leaves);
  $hash->{'node_count'} = scalar($tree->nodes);
  $hash->{'gappiness'} = mysql_getval("SELECT gappiness($node_id)");
  $hash->{'duplications'} = mysql_getval("SELECT num_dups_under_node($node_id)");
  $hash->{'duplication_fraction'} = mysql_getval("SELECT num_dups_under_node($node_id)/node_count($node_id)");
  my @hum_gen = grep {$_->taxon_id==9606} $tree->leaves;
  $hash->{'human_genes'} = scalar(@hum_gen);
  $hash->{'ensp'} = $hum_gen[0]->stable_id if ($hum_gen[0]);

  my $sw = "sitewise_aln";
  $hash->{'v_lrt'} = avg_sitewise("lrt_stat",$node_id,$sw,1);
  $hash->{'v_omega'} = avg_sitewise("omega",$node_id,$sw,1);
  $hash->{'no2x_omega'} = avg_sitewise("omega",$node_id,$sw,2);
  $hash->{'only2x_omega'} = avg_sitewise("omega",$node_id,$sw,3);
  $hash->{'no_seq_f_omega'} = avg_sitewise("omega",$node_id,$sw,4);
  $hash->{'no_aln_f_omega'} = avg_sitewise("omega",$node_id,$sw,5);
  $hash->{'no_f_omega'} = avg_sitewise("omega",$node_id,$sw,6);
  $hash->{'p_omega'} = avg_sitewise("omega",$node_id,$sw,7);
  $hash->{'g_omega'} = avg_sitewise("omega",$node_id,$sw,8);
  $hash->{'l_omega'} = avg_sitewise("omega",$node_id,$sw,9);

  $hash->{'v_pscs'} = num_pscs($node_id,$sw,1);
  $hash->{'p_pscs'} = num_pscs($node_id,$sw,7);
  $hash->{'g_pscs'} = num_pscs($node_id,$sw,8);
  $hash->{'l_pscs'} = num_pscs($node_id,$sw,9);
  $hash->{'no2x_pscs'} = num_pscs($node_id,$sw,2);
  $hash->{'only2x_pscs'} = num_pscs($node_id,$sw,3);

  

  # Collect total branch lengths of taxonomically-defined subtrees:
  my $v = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_param_set($tree->adaptor->dbc,1);
  print Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($v)."\n";
  my $no2 = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_param_set($tree->adaptor->dbc,2);
  my $only2 = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_param_set($tree->adaptor->dbc,3);
  my $p = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_param_set($tree->adaptor->dbc,7);
  my $g = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_param_set($tree->adaptor->dbc,8);
  my $l = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_param_set($tree->adaptor->dbc,9);

  $hash->{'v_bl_total'} = subtree_total($tree,$v);
  $hash->{'p_bl_total'} = subtree_total($tree,$p);
  $hash->{'l_bl_total'} = subtree_total($tree,$l);
  $hash->{'g_bl_total'} = subtree_total($tree,$g);
  $hash->{'v_bl_max'} = subtree_max($tree,$v);
  $hash->{'p_bl_max'} = subtree_max($tree,$p);
  $hash->{'l_bl_max'} = subtree_max($tree,$l);
  $hash->{'g_bl_max'} = subtree_max($tree,$g);
  $hash->{'v_bl_avg'} = subtree_avg($tree,$v);
  $hash->{'p_bl_avg'} = subtree_avg($tree,$p);
  $hash->{'l_bl_avg'} = subtree_avg($tree,$l);
  $hash->{'g_bl_avg'} = subtree_avg($tree,$g);
  $hash->{'v_leaves'} = subtree_leaves($tree,$v);
  $hash->{'g_leaves'} = subtree_leaves($tree,$p);
  $hash->{'l_leaves'} = subtree_leaves($tree,$l);
  $hash->{'g_leaves'} = subtree_leaves($tree,$g);
  eval {
    $hash->{'no2x_leaves'} = subtree_leaves($tree,$no2);
    $hash->{'only2x_leaves'} = subtree_leaves($tree,$only2);
    $hash->{'no2x_bl_total'} = subtree_total($tree,$no2);
    $hash->{'no2x_bl_avg'} = subtree_avg($tree,$no2);
    $hash->{'only2x_bl_total'} = subtree_total($tree,$only2);
    $hash->{'only2x_bl_avg'} = subtree_avg($tree,$only2);
  };
  if (@$) {
    $hash->{'no2x_bl_total'} = '';
    $hash->{'no2x_bl_avg'} = '';
    $hash->{'only2x_bl_total'} = '';
    $hash->{'only2x_bl_avg'} = '';
  }
  $hash->{'tree_length_total'} = sprintf "%.3f", total_distance($tree);
  $hash->{'tree_length_max'} = sprintf "%.3f", max_distance($tree);
  $hash->{'tree_length_avg'} = sprintf "%.3f", avg_distance($tree);
  
  $hash->{'avg_gc'} = gc_content($tree);

  my $seq_len=0;
  map {$seq_len += $_->seq_length} $tree->leaves;
  $hash->{'avg_seq_length'} = sprintf "%.3f", $seq_len / scalar($tree->leaves);

  my $sa = $tree->get_SimpleAlign;
  $hash->{'aln_length'} = sprintf "%d", $sa->length;
  $hash->{'aln_percent_identity'} = sprintf "%.3f", $sa->percentage_identity;

  my $hash2;
  map {$hash2->{'eslr_'.$_} = $hash->{$_}} keys %{$hash};

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
  if ($ENV{'USER'} =~ /gj1/) {
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

  my $temp_dir = "/tmp/eslr_rcmd";
  use File::Path;
  mkpath($temp_dir);

  my $temp_in = $temp_dir . "/temp_in.txt";
  my $temp_out = $temp_dir . "/temp_out.txt";

  open(OUT,">$temp_in");
  print OUT $rcmd."\n";
  close(OUT);
  
# cmd to run: /software/R-2.7.1/bin/R CMD BATCH $filename
  my $rc = system("/ebi/research/software/Linux_x86_64/bin/R-2.7.0 --slave --vanilla < $temp_in > $temp_out");
  #my $rc = system("/software/R-2.7.1/bin/R --vanilla --slave < $temp_in ");
  die "R returned an error!" if ($rc);
  
  open(IN,"$temp_out");
  my @lines = <IN>;
  map {chomp} @lines; 
  close(IN);
  
  unlink($temp_in);
  unlink($temp_out);
  unlink($temp_dir);

  return @lines;
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
