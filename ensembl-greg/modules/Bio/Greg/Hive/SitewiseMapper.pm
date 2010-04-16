package Bio::Greg::Hive::SitewiseMapper;

use strict;

use Cwd;
use Time::HiRes qw(sleep);

use Bio::EnsEMBL::Compara::Member;
use Bio::EnsEMBL::Compara::ComparaUtils;

use Bio::Greg::EslrUtils;

use base ('Bio::Greg::Hive::Process');

sub fetch_input {
  my $self = shift;

  ### DEFAULT PARAMETERS ###
  my $params = {};
  $params->{'do_mapping'}      = 1;
  $params->{'do_filter'}      = 1;
  $params->{'collect_pfam'}    = 1;
  $params->{'collect_uniprot'} = 1;
  $params->{'collect_exons'}   = 1;
  $params->{'collect_go'}      = 1;
  $params->{'pfam_taxon_ids'}  = '9606,10090';
  $params->{'go_taxon_ids'}    = '9606,10090';
  $params->{'genome_taxon_id'} = 9606;
  #########################

  # Fetch parameters from the two possible locations. Input_id takes precedence!

  $self->load_all_params($params);

}

sub run {
  my $self = shift;
  
  $self->{'start_time'} = time() * 1000;
  $self->compara_dba->dbc->disconnect_when_inactive(1);

  $self->param('parameter_set_id',0);

  if ($self->param('do_filter')) {
    print "Doing site-wise filter calculations...\n";
    $self->do_filter();
    print "  -> Finished doing site-wise filtering!\n";
  }
  
  if ( $self->param('do_mapping') ) {
    print "Mapping sitewise to genome...\n";
    $self->do_mapping();
    print "  -> Finished mapping values!\n";
  }

  if ( $self->param('collect_pfam') ) {
    print "Collecting Pfam annotations...\n";
    $self->collect_pfam();
    print "  -> Finished collecting Pfam!\n";
  }

  if ( $self->param('collect_uniprot') ) {
    print "Collecting UniProt annotations...\n";
    $self->collect_uniprot();
    print "  -> Finished collecting UniProt!\n";
  }

  if ( $self->param('collect_exons') ) {
    print "Collecting Exon annotations...\n";
    $self->collect_exons();
    print "  -> Finished collecting Exons!\n";
  }

  if ( $self->param('collect_go') ) {
    print "Collecting GO annotations...\n";
    $self->collect_go();
    print "  -> Finished collecting GO terms!\n";
  }

  $self->compara_dba->dbc->disconnect_when_inactive(0);
}

sub get_id {
  my $self = shift;

  return $self->param('node_id');
}

sub do_mapping {
  my $self = shift;

  my $node_id = $self->get_id;
  print "Mapping sitewise $node_id to genome...\n";

  my $mapping_taxon = $self->param('genome_taxon_id');
  my $mapped_sites = Bio::Greg::EslrUtils->mapSitewiseToGenome( $self->get_tree, $mapping_taxon );

  my @array_of_strings;
  foreach my $map ( @{$mapped_sites} ) {
    my @values =
      ( $node_id, 
	$self->param('parameter_set_id'),
	$map->{aln_position},
	$map->{member_id},
	'"'.$map->{chr}.'"',
	$map->{start},
	$map->{end});
    my $string = '(' . join( ',', @values ) . ')';
    #print $string."\n";
    push @array_of_strings, $string;
  }

  if ( scalar(@array_of_strings) > 0 ) {
    my $values_string = join( ",", @array_of_strings );
    my $cmd =
      "REPLACE INTO sitewise_genome (node_id,parameter_set_id,aln_position,member_id,chr_name,chr_start,chr_end) VALUES $values_string";
    my $sth = $self->compara_dba->dbc->prepare($cmd);
    $sth->execute;
    $sth->finish;
  }
}

sub collect_go {
  my $self = shift;

  my $tree = $self->get_tree;
  my @leaves = @{ $tree->get_all_leaves };

  my @taxon_ids;
  if ( defined $self->param('go_taxon_ids') ) {
    @taxon_ids = split( ",", $self->param('go_taxon_ids') );
  } else {
    @taxon_ids = (9606);
  }

  my %go_terms;
  foreach my $taxon_id (@taxon_ids) {
    my @taxon_leaves = grep { $_->taxon_id == $taxon_id } @leaves;

    foreach my $leaf (@taxon_leaves) {
      my $ts = $leaf->transcript;
      next if ( !defined $ts );

      my $db_entries = $ts->get_all_DBLinks;

      # Grep out all the GO xref entries.
      my @keepers = grep { $_->dbname eq "GO" } @{$db_entries};
      foreach my $db_e (@keepers) {
        my $iea = 0;
        foreach ( @{ $db_e->get_all_linkage_info } ) {
          $iea = 1 if ( $_->[0] eq 'IEA' );
        }
        if ( $self->param('go_ignore_iea') && $iea ) {

          # Don't insert!
        } else {
          $self->insert_go_term( $self->get_id, $leaf, $db_e );
          sleep(0.05);
        }
      }
    }
  }

  my @gos = keys(%go_terms);
  return \@gos;
}

sub insert_go_term {
  my $self = shift;
  my ( $node_id, $leaf, $db_e ) = @_;

  my $evidence = '';
  foreach ( @{ $db_e->get_all_linkage_info } ) {
    $evidence = $_->[0];
  }

  my $cmd =
    "REPLACE INTO go_terms (node_id,member_id,source_taxon,stable_id,go_term,evidence_code) values (?,?,?,?,?,?);";
  print "  "
    . join( " ", $leaf->dbID, $leaf->taxon_id, $leaf->stable_id, $db_e->display_id, $evidence )
    . "\n";
  my $sth = $self->compara_dba->dbc->prepare($cmd);
  $sth->execute( $node_id, $leaf->dbID, $leaf->taxon_id, $leaf->stable_id, $db_e->display_id,
    $evidence );
}

sub collect_uniprot {
  my $self = shift;

  my $tree = $self->get_tree;
  my $sa = $tree->get_SimpleAlign;
  my $pos_id_hash;

  my $orig_cwd = cwd();
  if ( $ENV{'USER'} =~ /gj1/ ) {
    chdir( $ENV{HOME} . "/src/greg-ensembl/projects/eslr/uniprot" );
  } else {
    chdir( $ENV{HOME} . "/lib/greg-ensembl/projects/eslr/uniprot" );
  }

  my $url = Bio::Greg::EslrUtils->urlFromConnection( $tree->adaptor->dbc );
  print "$url\n";
  print 'CWD:' . cwd() . "\n";

  my $cmd = qq^
    REPLACE INTO sitewise_tag (node_id,parameter_set_id,aln_position,tag,value,source) values(?,?,?,?,?,?);
  ^;
  my $sth = $self->compara_dba->dbc->prepare($cmd);

  foreach my $leaf ( $tree->leaves ) {
    next unless ( $leaf->taxon_id == 9606 || $leaf->taxon_id == 10090 );

    my $stable_id = $leaf->stable_id;

    my @output;
    my $rc;
    if ( $ENV{'USER'} =~ /gj1/ ) {
      my $proxy = "-Dhttp.proxyHost=wwwcache.sanger.ac.uk  -Dhttp.proxyPort=3128";
      open( JAVA,
        "java -Xmx512m $proxy  -cp uniprotjapi.properties -jar uniProtExtraction.jar $stable_id $url |"
      ) or $self->throw("Cannot run UniProt Collector!");
      @output = <JAVA>;
      $rc     = close(JAVA);
    } else {
      open( JAVA, "java -cp uniprotjapi.properties -jar uniProtExtraction.jar $stable_id $url |" )
        or $self->throw("Cannot run UniProt Collector!");
      @output = <JAVA>;
      $rc     = close(JAVA);
    }

    foreach my $line (@output) {
      chomp $line;
      my @tokens = split( "\t", $line );
      print $line. "\n";

      my $acc     = $tokens[0];
      my $seq_pos = $tokens[1];
      my $source  = $tokens[2];
      my $tag     = $tokens[3];
      my $value   = $tokens[4];
      my $residue = $tokens[5];

      next if ( !$seq_pos );

      my $node_id          = $self->get_id;
      my $parameter_set_id = $self->param('parameter_set_id');
      my $aln_position     = $sa->column_from_residue_number( $stable_id, $seq_pos );

      my ($seq) = $sa->each_seq_with_id($stable_id);
      my $actual_residue = substr( $seq->seq, $aln_position - 1, 1 );

      if ( $actual_residue ne $residue ) {
        print "Residues don't match: "
          . sprintf( "%s %s ensembl:%s pdb:%s", $node_id, $aln_position, $actual_residue,
          $residue );
        next;
      }

      $sth->execute( $node_id, $parameter_set_id, $aln_position, $tag, $value, $source );
    }
  }

  $sth->finish;

  chdir $orig_cwd;

}

sub collect_exons {
  my $self = shift;

  my $tree = $self->get_tree;
  my $sa = $tree->get_SimpleAlign;

  print $sa->length . "\n";

  use Time::HiRes qw(sleep);
  my $cmd =
    "REPLACE INTO sitewise_tag (node_id,aln_position,parameter_set_id,tag,value,source) VALUES (?,?,?,?,?,?)";
  my $sth = $tree->adaptor->prepare($cmd);

  foreach my $leaf ( $tree->leaves ) {
    next unless ( $leaf->taxon_id == 9606 );
    print $leaf->seq_length . "\n";
    print $leaf->stable_id . "\n";
    my $name  = $leaf->stable_id;
    my $tx    = $leaf->get_Transcript;
    my @exons = @{ $tx->get_all_Exons };
    my $i     = 0;
    foreach my $exon (@exons) {
      $i++;
      my $sequence = $leaf->sequence;
      my $exon_pep = '';
      eval { $exon_pep = $exon->peptide($tx)->seq; };
      next if ( $exon_pep eq '' );
      if ( $sequence =~ /($exon_pep)/g ) {
        my $pep_end   = pos($sequence);
        my $pep_start = $pep_end - length($exon_pep) + 1;
        my $len       = $pep_end - $pep_start;

        print $sequence. "\n";
        print ' ' x ( $pep_start - 1 ) . $exon_pep . "\n";

        my $exon_position = 'MIDDLE';
        $exon_position = 'FIRST' if ( $i == 1 );
        $exon_position = 'LAST'  if ( $i == scalar(@exons) );

        my $aln_start = $sa->column_from_residue_number( $leaf->stable_id, $pep_start );
        my $aln_end   = $sa->column_from_residue_number( $leaf->stable_id, $pep_end );

        foreach my $pos ( $aln_start .. $aln_end ) {
	  my $node_id = $self->get_id;
	  my $parameter_set_id = $self->param('parameter_set_id');
          # Store whether we're in a first, middle, or last exon.
          $sth->execute( $node_id, $pos, $parameter_set_id,
            'EXON', $exon_position, "EnsEMBL" );

          # Store the smallest distance to an exon junction.
          my $dist_left  = $aln_start - $pos;
          my $dist_right = $aln_end - $pos;

          my $dist = $dist_left;
          if ( $dist_right + $dist_left <= 0 ) {
            $dist = $dist_right;
          }

          $sth->execute( $node_id, $pos, $parameter_set_id,
            'SPLICE_DISTANCE', $dist, "EnsEMBL" );
        }

      } else {
        die("Exon $exon_pep did not match protein $sequence!\n");
      }
    }
  }
}

sub collect_pfam {
  my $self = shift;

  my @taxon_ids;
  my $pfam_taxon_ids = $self->param('pfam_taxon_ids');
  if ( defined $pfam_taxon_ids ) {
    @taxon_ids = split( ",", $pfam_taxon_ids );
  } else {
    @taxon_ids = (9606);
  }

  my $tree = $self->get_tree;
  my $sa = $tree->get_SimpleAlign;
  my $pos_id_hash;
  foreach my $leaf ( $tree->leaves ) {
    next unless ( grep { $leaf->taxon_id == $_ } @taxon_ids );

    print $leaf->stable_id . "\n";
    my $name = $leaf->stable_id;
    my $tx   = $leaf->get_Translation;
    next if ( !defined $tx );
    my @features = @{ $tx->get_all_ProteinFeatures('PFam') };
    foreach my $f (@features) {
      my $pf_id = $f->display_id;
      my $pf_lo = $f->hstart;       # Start and end in the hit (Pfam domain) coordinates.
      my $pf_hi = $f->hend;
      my $lo    = $f->start;
      my $hi    = $f->end;
      $hi = $tx->length if ( $hi > $tx->length );
      print "$pf_id  lo:$lo hi:$hi hstart:$pf_lo hend:$pf_hi\n";
      foreach my $i ( 0 .. ( $hi - $lo ) ) {
        my $pos     = $lo + $i;
        my $aln_col = $sa->column_from_residue_number( $name, $pos );
        my $obj     = {
          pf_pos => $pf_lo + $i,
          score  => $f->score
        };
        $pos_id_hash->{ $pf_id . "_" . $aln_col } = $obj;
      }
    }
  }

  my @array_of_strings = ();

  use Time::HiRes qw(sleep);
  foreach my $key ( sort keys %{$pos_id_hash} ) {
    my $id;
    my $pos;
    ( $id, $pos ) = split( "_", $key );
    my $obj    = $pos_id_hash->{$key};
    my $pf_pos = $obj->{'pf_pos'};
    my $score  = $obj->{'score'};
    printf( "%s %s %s %s\n", $id, $pos, $pf_pos, $score );
    my @values =
      ( $self->get_id, $pos, $self->param('parameter_set_id'), '"DOMAIN"', '"' . $id . '"', $score );
    my $string = '(' . join( ',', @values ) . ')';
    push @array_of_strings, $string;
  }

  if ( scalar(@array_of_strings) > 0 ) {
    my $values_string = join( ",", @array_of_strings );
    my $cmd =
      "REPLACE INTO sitewise_tag (node_id,aln_position,parameter_set_id,tag,value,source) VALUES $values_string";
    my $sth = $tree->adaptor->prepare($cmd);
    $sth->execute;
    $sth->finish;
  }
}


sub do_filter {
  my $self = shift;

  # Create a local params object to temporarly remove the subtree settings.
  my $params = $self->params;
  $params->{parameter_set_id} = 0;
  $params->{keep_species} = '';
  $params->{remove_species} = '';

  my $tree = $self->get_tree($params);
  my $sa = $tree->get_SimpleAlign;

  Bio::EnsEMBL::Compara::AlignUtils->pretty_print($sa,{length=>400});

  my $dba = $self->compara_dba;

  # Subroutines to return a list of taxon IDs with specific features.
  sub clade_taxon_ids {
    my $clade = shift || 1;
    my @genomes = Bio::EnsEMBL::Compara::ComparaUtils->get_genomes_within_clade($dba,$clade);
    my @taxon_ids = map {$_->taxon_id} @genomes;
    return @taxon_ids;
  }
  sub subtract {
    my $list_a = shift;
    my @remove_us = @_;
    my $hash;
    map {$hash->{$_}=1} @$list_a;
    foreach my $list_b (@remove_us) {
      map {delete $hash->{$_}} @$list_b;
    }
    return keys %$hash;
  }

  my @mammals_arr = clade_taxon_ids("Eutheria");
  my @primates_arr = clade_taxon_ids("Primates");
  my @glires_arr = clade_taxon_ids("Glires");
  my @laur_arr = clade_taxon_ids("Laurasiatheria");

  my @outgroup_arr = subtract(\@mammals_arr,\@primates_arr,\@glires_arr,\@laur_arr);

  print "primates: @primates_arr\n";

  my @array_of_strings;
  foreach my $aln_position (1 .. $sa->length) {
    my $primate_count = 0;
    my $glires_count = 0;
    my $laur_count = 0;
    my $outgroup_count = 0;
    my $total_count = 0;

    my $taxon_id_hash = {};
    
    my $slice = $sa->slice($aln_position,$aln_position);
    
    foreach my $leaf ($tree->leaves) {
      my $seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id($slice,$leaf->stable_id);
      if (!defined $seq || $seq->seq =~ m/[-Xx]/) {
	#print "-";
	next;
      }
      #next if ($seq->seq =~ m/[-Xx]/);
      #print $seq->seq;
      $taxon_id_hash->{$leaf->taxon_id} = 1;
    }

    $primate_count = scalar(grep {$taxon_id_hash->{$_}==1} @primates_arr);
    $glires_count = scalar(grep {$taxon_id_hash->{$_}==1} @glires_arr);
    $laur_count = scalar(grep {$taxon_id_hash->{$_}==1} @laur_arr);
    $outgroup_count = scalar(grep {$taxon_id_hash->{$_}==1} @outgroup_arr);
    $total_count = scalar(keys %$taxon_id_hash);

    my $filter_val = 0;
    $filter_val = 1 if ($total_count >= 3);
    $filter_val = 2 if ($total_count >= 6);
    # This corresponds to the filtering criteria used by Pollard et al. Gen Res 2010:
    $filter_val = 3 if ($primate_count >= 3 && $glires_count >= 2 && $outgroup_count >= 1);

    print "$aln_position $filter_val\n";

    my @values =
      ( $self->get_id, $aln_position, 0, '"FILTER"', $filter_val, 'NULL' );
    my $string = '(' . join( ',', @values ) . ')';
    push @array_of_strings, $string;
  }

  if ( scalar(@array_of_strings) > 0 ) {
    my $values_string = join( ",", @array_of_strings );
    my $cmd =
      "REPLACE INTO sitewise_tag (node_id,aln_position,parameter_set_id,tag,value,source) VALUES $values_string";
    #print $cmd."\n";
    my $sth = $tree->adaptor->prepare($cmd);
    $sth->execute;
    $sth->finish;
  }

}

sub create_plot {
  my $self = shift;

  my $base_dir = "/lustre/scratch103/ensembl/gj1/2xmammals_plots";
  my $node_id  = $self->get_id;

  use Digest::MD5 qw(md5_hex);
  my $digest       = md5_hex($node_id);
  my $small_digest = substr( $digest, 0, 10 );
  my $subdir       = substr( $digest, 0, 1 );
  my $sub_subdir   = substr( $digest, 1, 1 );

  print "$digest\n";
  print "DIRS: $subdir $sub_subdir\n";
  my $output_dir = $base_dir . "/" . $subdir . "/" . $sub_subdir;
  use File::Path qw(mkpath);
  mkpath( $output_dir, { mode => 0777 } );

  my $tree = $self->get_tree;
  my @human_peps = grep { $_->taxon_id == 9606 } $tree->leaves;
  foreach my $hum_pep (@human_peps) {
    my $stable_id   = $hum_pep->stable_id;
    my $output_file = $output_dir . "/${stable_id}.pdf";

    print "PLOT OUTPUT: $output_file\n";

    my $params = {
      parameter_set_id     => 0,
      remove_blank_columns => 1,
      mask_outside_subtree => 1,
      remove_subtree       => 0,
      mask_gblocks         => 0,
      sitewise_table       => 'sitewise_aln'
    };
    use Bio::Greg::EslrPlots;
    my $fresh_tree = $tree->copy;
    Bio::Greg::EslrPlots->plotTreeWithOmegas( $output_file, $self->params, $fresh_tree );
    $fresh_tree->release_tree;
  }

}

sub write_output {
  my $self = shift;

}

1;
