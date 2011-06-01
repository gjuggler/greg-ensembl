package Bio::Greg::Hive::SitewiseMapper;

use strict;

use Cwd;
use Time::HiRes qw(sleep);

use Bio::EnsEMBL::Compara::Member;
use Bio::EnsEMBL::Compara::ComparaUtils;

use Bio::Greg::EslrUtils;

use base ('Bio::Greg::Hive::Process');

sub param_defaults {
  return {
    
  };
}

sub run {
  my $self = shift;

  warn("This module is generally used by other runnables...");
}

sub get_id {
  my $self = shift;

  return $self->node_id;
}

sub human_or_mouse_member {
  my $self = shift;
  my $tree = shift;

  # Only use human or mouse genes.
  my ($leaf) = grep {$_->taxon_id == 9606} $tree->leaves;
  ($leaf) = grep {$_->taxon_id == 10090} $tree->leaves if (!defined $leaf);

  warn("No human or mouse gene!") unless (defined $leaf);

  return $leaf;
}

sub collect_uniprot {
  my $self = shift;

  my $tree = $self->get_tree;
  my $sa   = $tree->get_SimpleAlign;
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

      my $node_id          = $self->node_id;
      my $parameter_set_id = $self->parameter_set_id;
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
  my $tree = shift;
  my $pep_aln = shift;

  my $leaf = $self->human_or_mouse_member($tree);
  if (!defined $leaf) {
    warn("Failed to collect exon data -- no human or mouse genes!");
    return {};
  }

  my $site_objects = {};

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
      
#      print $sequence. "\n";
#      print ' ' x ( $pep_start - 1 ) . $exon_pep . "\n";
      
      my $exon_position = 'MIDDLE';
      $exon_position = 'FIRST' if ( $i == 1 );
      $exon_position = 'LAST'  if ( $i == scalar(@exons) );
      
      my $aln_start = $pep_aln->column_from_residue_number( $leaf->stable_id, $pep_start );
      my $aln_end   = $pep_aln->column_from_residue_number( $leaf->stable_id, $pep_end );
      
      foreach my $pos ( $aln_start .. $aln_end ) {
        # Store the smallest distance to an exon junction.
        my $dist_left  = $aln_start - $pos;
        my $dist_right = $aln_end - $pos;
        
        my $dist = $dist_left;
        if ( $dist_right + $dist_left <= 0 ) {
          $dist = $dist_right;
        }
        
        $site_objects->{$pos} = {
          exon_position => $exon_position,
          splice_distance => $dist
        };
      }
    } else {
      die("Exon $exon_pep did not match protein $sequence!\n");
    }
  }

  return $site_objects;
}

sub collect_pfam {
  my $self = shift;
  my $tree = shift;
  my $pep_aln = shift;

  print ("  collecting pfam data...\n") if ($self->debug);

  # Only take human or mouse protein for pfam reference.
  my $leaf = $self->human_or_mouse_member($tree);
  if (!defined $leaf) {
    warn("Failed to collect pfam data -- no human or mouse genes!");
    return {};
  }

  my $site_objects = {};

#  print $leaf->stable_id . "\n";
  my $name = $leaf->stable_id;
  my $tx   = $leaf->get_Transcript->translation;
  if ( !defined $tx ) {
    warn("Transcript has no translation!");
    return {};
  }

  my @features = @{ $tx->get_all_ProteinFeatures('PFam') };
  foreach my $f (@features) {
    my $pf_id = $f->display_id;
    my $pf_lo = $f->hstart;       # Start and end in the hit (Pfam domain) coordinates.
    my $pf_hi = $f->hend;
    my $lo    = $f->start;
    my $hi    = $f->end;
    $hi = $tx->length if ( $hi > $tx->length );
    print "  $pf_id  lo:$lo hi:$hi hstart:$pf_lo hend:$pf_hi\n" if ($self->debug);
    foreach my $i ( 0 .. ( $hi - $lo ) ) {
      my $pos     = $lo + $i;
      my $aln_col = $pep_aln->column_from_residue_number( $name, $pos );

      $site_objects->{$aln_col} = {
        pfam_domain => $pf_id,
        pfam_position => $pf_lo + $i,
        pfam_score  => $f->score
      };
    }
  }
  return $site_objects;
}

sub do_filter {
  my $self = shift;
  my $tree = shift;
  my $pep_aln = shift;

  my $dba = $self->compara_dba;

  # Subroutines to return a list of taxon IDs with specific features.
  sub clade_taxon_ids {
    my $clade = shift || 1;
    my @genomes = Bio::EnsEMBL::Compara::ComparaUtils->get_genomes_within_clade( $dba, $clade );
    my @taxon_ids = map { $_->taxon_id } @genomes;
    return @taxon_ids;
  }

  sub subtract {
    my $list_a    = shift;
    my @remove_us = @_;
    my $hash;
    map { $hash->{$_} = 1 } @$list_a;
    foreach my $list_b (@remove_us) {
      map { delete $hash->{$_} } @$list_b;
    }
    return keys %$hash;
  }

  my @mammals_arr  = clade_taxon_ids("Eutheria");
  my @primates_arr = clade_taxon_ids("Primates");
  my @glires_arr   = clade_taxon_ids("Glires");
  my @laur_arr     = clade_taxon_ids("Laurasiatheria");

  my @outgroup_arr = subtract( \@mammals_arr, \@primates_arr, \@glires_arr, \@laur_arr );

#  print "primates: @primates_arr\n";

  my $site_objects;

  foreach my $aln_position ( 1 .. $pep_aln->length ) {
    my $primate_count  = 0;
    my $glires_count   = 0;
    my $laur_count     = 0;
    my $outgroup_count = 0;
    my $total_count    = 0;

    my $taxon_id_hash = {};
    my $slice;
    eval { $slice = $pep_aln->slice( $aln_position, $aln_position, 1 ); };
    next if ($@);

    foreach my $leaf ( $tree->leaves ) {
      my $seq = Bio::EnsEMBL::Compara::AlignUtils->get_seq_with_id( $slice, $leaf->stable_id );
      if ( !defined $seq || $seq->seq =~ m/[-Xx]/ ) {
        #print "It's crap!";
        next;
      }

      #next if ($seq->seq =~ m/[-Xx]/);
      #print $seq->seq;
      $taxon_id_hash->{ $leaf->taxon_id } = 1;
    }

    $primate_count  = scalar( grep { $taxon_id_hash->{$_} == 1 } @primates_arr );
    $glires_count   = scalar( grep { $taxon_id_hash->{$_} == 1 } @glires_arr );
    $laur_count     = scalar( grep { $taxon_id_hash->{$_} == 1 } @laur_arr );
    $outgroup_count = scalar( grep { $taxon_id_hash->{$_} == 1 } @outgroup_arr );
    $total_count    = scalar( keys %$taxon_id_hash );

    my $filter_val = 0;
    $filter_val = 1 if ( $total_count >= 3 );
    $filter_val = 2 if ( $total_count >= 6 );

    # This corresponds to the filtering criteria used by Pollard et al. Gen Res 2010:
    $filter_val = 3 if ( $primate_count >= 3 && $glires_count >= 2 && $outgroup_count >= 1 );

    #print "$aln_position $filter_val\n";

    $site_objects->{$aln_position} = {
      filter_value => $filter_val
    };
  }

  return $site_objects;
}

sub write_output {
  my $self = shift;

}

1;
