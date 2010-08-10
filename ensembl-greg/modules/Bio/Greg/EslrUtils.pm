package Bio::Greg::EslrUtils;

use strict;

sub baseDirectory {
  my $class = shift;
  if ( $ENV{'USER'} =~ /gj1/ ) {
    return $ENV{'HOME'} . '/src/greg-ensembl';
  } else {
    return $ENV{'HOME'} . '/lib/greg-ensembl';
  }
}

sub defaultMysqlURL {
  my $class = shift;

  if ( $ENV{'USER'} =~ /gj1/ ) {
    my $url = 'mysql://ensadmin:ensembl@ens-research:3306/gj1_2xmammals_gold';
    print $url. "\n";
    return $url;

    # return 'mysql://ensadmin:ensembl@compara2:5316/avilella_compara_homology_53';
  } else {
    my $url = 'mysql://ensadmin:ensembl@127.0.0.1:5425/gj1_2xmammals_gold';
    print $url. "\n";
    return $url;
  }
}

sub urlFromConnection {
  my $class = shift;
  my $dbc   = shift;

  my $url = sprintf( "mysql://%s:%s@%s:%s/%s",
    $dbc->username, $dbc->password, $dbc->host, $dbc->port, $dbc->dbname );
  return $url;
}

sub mysqlArgsFromConnection {
  my $class = shift;
  my $dbc   = shift;

  my $args = sprintf( "-h%s -u%s -p%s -P%s %s",
    $dbc->host, $dbc->username, $dbc->password, $dbc->port, $dbc->dbname );
  return $args;
}

sub url_to_mysql_args {
  my $class = shift;
  my $url   = shift;

  $url =~ m!mysql://(.*?):(.*?)@(.*?)/(.*)!;

  # Something like mysql -uensadmin -pensembl -hensdb-2-12 -P5106 -A gj1_compara_53
  # $1 = user, $2 = pw, $3 = server+port, $4 = db
  my ( $u, $p, $h, $P, $db ) = undef;
  $u  = $1;
  $p  = $2;
  $h  = $3;
  $P  = 3306;
  $db = $4;
  if ( $h =~ m/(.*?):(.*)/ ) {
    $h = $1;
    $P = $2;
  }

  my $mysqlargs = sprintf( "-u%s -p%s -h%s -P%s -A %s", $u, $p, $h, $P, $db );
  return $mysqlargs;
}

sub url_to_hashref {
  my $class = shift;
  my $url   = shift;

  my $args = $class->url_to_mysql_args($url);

  $args =~ m/-u(.*) -p(.*) -h(.*) -P(.*) -A (.*)/;
  my $data = {
    user     => $1,
    password => $2,
    host     => $3,
    port     => $4,
    database => $5
  };
  return $data;
}

sub loadConfiguration {
  my $class  = shift;
  my $params = shift;
  my $file   = shift;

  my $obj = {};

  foreach my $key ( keys %{$params} ) {
    $obj->{$key} = $params->{$key};
  }

  open( IN, $file );
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
  close( IN, $file );

  return $obj;
}

sub find_member_by_external_id {
  my $class       = shift;
  my $compara_dba = shift;
  my $id          = shift;

  Bio::EnsEMBL::Registry->load_registry_from_multiple_dbs( {
      -host => 'ensembldb.ensembl.org',
      -user => 'anonymous'
    }
  );
  my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "gene" );
  my $member_adaptor = $compara_dba->get_MemberAdaptor;

  my @genes = @{ $gene_adaptor->fetch_all_by_external_name($id) };
  push @genes, $gene_adaptor->fetch_by_stable_id($id) if ( scalar @genes == 0 );

  foreach my $gene (@genes) {
    my $stable_id = $gene->stable_id;

    my $member = $member_adaptor->fetch_by_source_stable_id( undef, $stable_id );
    return $member if ( defined $member );
  }
  return undef;
}

sub mapSitewiseToGenome {
  my $class    = shift;
  my $tree     = shift;
  my $taxon_id = shift;

  my $sa      = $tree->get_SimpleAlign;
  my $aln_len = $sa->length;
  Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $sa, { length => 200 } );

  my @leaves         = $tree->leaves;
  my @genomic_coords = ();
  foreach my $leaf (@leaves) {
    next unless ( $leaf->taxon_id == $taxon_id );
    eval {
      my ($seq) = $sa->each_seq_with_id( $leaf->stable_id );
      my $seq_str = $seq->seq;

      my $tscr = $leaf->get_Transcript;
      print STDERR $tscr->stable_id . "\n";
      $tscr = $tscr->transform("chromosome");
      next unless ( defined $tscr );
      my $chr = "chr" . $tscr->slice->seq_region_name;

      foreach my $i ( 1 .. $aln_len ) {
        my $char = substr( $seq_str, $i - 1, 1 );

        my $loc = $seq->location_from_column($i);

        #print "$i $char $loc\n";
        next if ( !defined $loc || $loc->location_type() eq 'IN-BETWEEN' );

        #next if ($char eq '-');
        my @gc_arr = $tscr->pep2genomic( $loc->start, $loc->end );
        my $gc = $gc_arr[0];
        next unless ( $gc && $gc->isa("Bio::EnsEMBL::Mapper::Coordinate") );

        my $strand = "+";
        $strand = "-" if ( $gc->strand == -1 );

        my $start = $gc->start;
        my $end   = $gc->end;

        push(
          @genomic_coords, {
            chr          => $chr,
            start        => $start,
            end          => $end,
            aln_position => $i,
            char         => $char,
            stable_id    => $leaf->stable_id,
            node_id      => $tree->node_id,
            member_id    => $leaf->dbID,
            strand       => $strand
          }
        );
      }
    };
    warn() if $@;
  }
  return \@genomic_coords;
}

sub collectDuplicationTags {
  my $class  = shift;
  my $tree   = shift;
  my $params = shift;

  my $taxon_name = $tree->get_tagvalue('taxon_name');
  my $taxon_id   = $tree->get_tagvalue('taxon_id');

  my $sth = $tree->adaptor->prepare(
    "SELECT name FROM ncbi_taxa_name WHERE taxon_id=$taxon_id AND name_class='ensembl timetree mya';"
  );
  $sth->execute;
  my $taxon_mya = $sth->fetchrow_array();

  my $total_bl = sprintf "%.3f", total_distance($tree);
  my $max_bl   = sprintf "%.3f", $tree->max_distance;

  # Collect chromosome on human gene.
  my $hum_chr = '';
  my @hum_genes = grep { $_->taxon_id == 9606 } $tree->leaves;
  if ( scalar(@hum_genes) > 0 ) {
    my $hum_gene = $hum_genes[0];
    my $gene     = $hum_gene->get_Gene;
    $hum_chr = $gene->slice->seq_region_name;
  }

  my $tags = {
    taxon_name => $taxon_name,
    taxon_id   => $taxon_id,
    taxon_mya  => $taxon_mya,
    bl_total   => $total_bl,
    bl_max     => $max_bl,
    human_chr  => $hum_chr
  };
  my $prefix = 'dupldiv';
  my $mapped_tags;
  map { $mapped_tags->{ $prefix . '_' . $_ } = $tags->{$_} } keys %{$tags};
  return $mapped_tags;
}

sub avg_distance {
  my $tree = shift;
  my $dist = 0;
  map { $dist += dist_to_root($_); } $tree->leaves;

  #map {$dist += dist_to_root($_);print dist_to_root($_)."\n";} $tree->leaves;
  $dist = $dist / scalar( $tree->leaves );

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
  map { $dist += $_->distance_to_parent } $tree->nodes;
  return sprintf "%.4f", $dist;
}

sub max_distance {
  my $tree = shift;
  my $max  = $tree->max_distance;
  return sprintf "%.4f", $max;
}

sub subtree_avg {
  my $tree   = shift;
  my $params = shift;

  $params->{'node_id'} = $tree->node_id;
  my $new_tree =
    Bio::EnsEMBL::Compara::ComparaUtils->get_tree_for_comparative_analysis( $tree->adaptor->db,
    $params );
  $new_tree->print_tree;
  return avg_distance($new_tree);
}

sub subtree_total {
  my $tree   = shift;
  my $params = shift;

  $params->{'node_id'} = $tree->node_id;
  my $new_tree =
    Bio::EnsEMBL::Compara::ComparaUtils->get_tree_for_comparative_analysis( $tree->adaptor->db,
    $params );
  return total_distance($new_tree);
}

sub subtree_leaves {
  my $tree   = shift;
  my $params = shift;

  $params->{'node_id'} = $tree->node_id;
  my $new_tree =
    Bio::EnsEMBL::Compara::ComparaUtils->get_tree_for_comparative_analysis( $tree->adaptor->db,
    $params );
  return scalar( $new_tree->leaves );
}

sub subtree_max {
  my $tree   = shift;
  my $params = shift;

  $params->{'node_id'} = $tree->node_id;
  my $new_tree =
    Bio::EnsEMBL::Compara::ComparaUtils->get_tree_for_comparative_analysis( $tree->adaptor->db,
    $params );
  return max_distance($new_tree);
}

sub run_r {
  my $class  = shift;
  my $rcmd   = shift;
  my $params = shift;

  my $temp_dir = "/tmp/eslr_rcmd";
  use File::Path;
  mkpath($temp_dir);

  my $temp_in = $temp_dir . "/temp_in.txt";

  #my $temp_out = $temp_dir . "/temp_out.txt";

  open( OUT, ">$temp_in" );
  print OUT $rcmd . "\n";
  close(OUT);

  my $vanilla = "--vanilla";
  $vanilla = "--slave" if ( $params->{'silent'} );

  my $r_cmd = $class->get_r_command($params);
  print "Using executable: '$r_cmd'\n";

  my $rc = system("$r_cmd $vanilla < $temp_in");
  die "R returned an error!" if ($rc);

  unlink($temp_in);
  unlink($temp_dir);
}

sub get_r_command {
  my $class = shift;
  my $params = shift || {};

  my $r_cmd = "R-2.10.0";
  if ( $params->{'farm'} ) {
    $r_cmd = "/software/bin/R-2.10.0";
  } elsif ( $params->{'bigmem'} ) {
    $r_cmd = "bsub -Is -R'select[mem>10000] rusage[mem=10000]' -M10000000 /software/R-2.9.0/bin/R ";
  } elsif ( $ENV{'USER'} =~ /gj1/ ) {
    $r_cmd = "/software/R-2.10.0/bin/R";
  } else {

  }
  return $r_cmd;
}

sub get_r_values {
  my $class    = shift;
  my $rcmd     = shift;
  my $temp_dir = shift;

  if ( !$temp_dir ) {
    $temp_dir = "/tmp/eslr_rcmd";
    use File::Path;
    mkpath($temp_dir);
  }

  use Digest::MD5 qw(md5_hex);
  my $digest = md5_hex( $rcmd . rand() );
  $digest = substr( $digest, 0, 10 );

  my $temp_in  = $temp_dir . "/temp_in_" . $digest . ".txt";
  my $temp_out = $temp_dir . "/temp_out" . $digest . ".txt";

  #  print "TEMP IN: $temp_in\n";
  open( OUT, ">$temp_in" );
  print OUT $rcmd . "\n";
  close(OUT);

  # cmd to run: /software/R-2.7.1/bin/R CMD BATCH $filename
  my $r_cmd = $class->get_r_command;
  my $rc    = system("$r_cmd --slave --vanilla < $temp_in > $temp_out");

  #my $rc = system("/software/R-2.7.1/bin/R --vanilla --slave < $temp_in ");

  open( IN, "$temp_out" );
  my @lines = <IN>;
  map {
    chomp $_;
    $_ =~ s/\s*\[\d+\]\s*//g
  } @lines;
  close(IN);

  if ($rc) {
    print join( "\n", @lines );
    die "R returned an error!";
  }

  unlink($temp_in);
  unlink($temp_out);

  #unlink($temp_dir);

  # Take each line and make sure it's split up by all whitespace.
  my @values;
  foreach my $line (@lines) {
    my @tokens = split /\s+/, $line;
    push @values, @tokens;
  }

  return @values;
}

## Please see file perltidy.ERR
sub mysql_array {
  my $class = shift;
  my $dbc   = shift;
  my $cmd   = shift;

  my $sth = $dbc->prepare($cmd);
  $sth->execute();

  my $array_ref = $sth->fetchall_arrayref( [0] );
  $sth->finish;
  my @vals = @{$array_ref};
  @vals = map { @{$_}[0] } @vals;    # Some weird mappings to unpack the numbers from the arrayrefs.
  return @vals;
}

sub get_short_name_for_parameter_set {
  my $class            = shift;
  my $parameter_set_id = shift;

  my $id  = $parameter_set_id;
  my $str = "";
  $str = "vertebrates"    if ( $id == 1 );
  $str = "no-2x"          if ( $id == 2 );
  $str = "only-2x"        if ( $id == 3 );
  $str = "no-seq-f"       if ( $id == 4 );
  $str = "no-aln-f"       if ( $id == 5 );
  $str = "no-f"           if ( $id == 6 );
  $str = "primates"       if ( $id == 7 );
  $str = "laurasiatheria" if ( $id == 9 );
  $str = "glires"         if ( $id == 8 );
  return ucfirst($str);
}

1;
