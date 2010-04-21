package Bio::Greg::Hive::AlignmentScores;

use strict;
use Cwd;
use POSIX qw(ceil floor);
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::NestedSet;
use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::Process');

sub fetch_input {
  my ($self) = @_;

  my $params = {
    alignment_table       => 'protein_tree_member',
    alignment_score_table => 'protein_tree_member_score',
    alignment_scores_action => 'prank'    # Options: 'gblocks', 'prank', 'trimal', 'indelign'
  };

  $self->load_all_params($params);
  
  $self->compara_dba->dbc->disconnect_when_inactive(1);

  my $no_filter_param = $self->replace_params( $self->params, { alignment_score_filtering => 0 } );
  my ( $tree, $aln ) =
    Bio::EnsEMBL::Compara::ComparaUtils->get_tree_and_alignment( $self->compara_dba, $no_filter_param );
  $self->param('tree',$tree);
  $self->param('aln',$aln);
}

sub run {
  my $self = shift;

  my $tree = $self->param('tree');
  my $aln = $self->param('aln');

  my $params = $self->params;

  my $input_table = $params->{'alignment_table'};
  my $action      = $params->{'alignment_scores_action'};
  my $table       = $params->{'alignment_score_table'};

  if ( $action =~ m/gblocks/i ) {
    print " -> RUN GBLOCKS\n";
    my $score_hash = $self->run_gblocks( $tree, $aln, $params );
    $self->store_scores( $tree, $score_hash, $table );
  }

  if ( $action =~ 'prank' ) {
    print " -> RUN PRANK\n";
    $params->{prank_filtering_scheme} = $action;
    my $score_hash = $self->run_prank( $tree, $aln, $params );
    $self->store_scores( $tree, $score_hash, $table );
  }

  if ( $action =~ m/trimal/i ) {
    print " -> RUN TRIMAL\n";
    my $score_hash = $self->run_trimal( $tree, $aln, $params );
    $self->store_scores( $tree, $score_hash, $table );
  }

  if ( $action =~ m/indelign/i ) {
    print " -> RUN INDELIGN\n";
    eval {
      my $score_hash = $self->run_indelign( $tree, $aln, $params );
      $self->store_scores( $tree, $score_hash, $table );
    };
    if ($@) {
      print "Indelign error: $@\n";
    }
  }

  if ( $action =~ m/(coffee|score)/i ) {
    print " -> RUN TCOFFEE\n";
    my $score_hash = $self->run_tcoffee( $tree, $aln, $params );
    $self->store_scores( $tree, $score_hash, $table );
  }

}

sub store_scores {
  my $self         = shift;
  my $tree         = shift;
  my $score_hash   = shift;
  my $output_table = shift;

  $self->compara_dba->dbc->do("CREATE TABLE IF NOT EXISTS $output_table LIKE protein_tree_member_score");
  my $sth = $tree->adaptor->prepare(
    "REPLACE INTO $output_table (node_id,member_id,cigar_line) VALUES (?,?,?)");
  foreach my $leaf ( $tree->leaves ) {
    my $score_string = $score_hash->{ $leaf->stable_id };
    throw("No score string found when saving!") unless ( defined $score_string );
    $sth->execute( $leaf->node_id, $leaf->member_id, $score_string );
    printf "%20s %10s %s\n", $leaf->stable_id, $leaf->member_id, $score_string;
  }
  $sth->finish;
}

sub run_indelign {
  my $self   = shift;
  my $tree   = shift;
  my $aln    = shift;
  my $params = shift;

  my ( $i_obj, $d_obj, $ins_rate, $del_rate ) =
    Bio::EnsEMBL::Compara::AlignUtils->indelign( $aln, $tree, $params,
    $self->worker_temp_directory );
  my @ins = @$i_obj;
  my @del = @$d_obj;

  my $num_leaves         = scalar( $tree->leaves );
  my $tree_length        = Bio::EnsEMBL::Compara::TreeUtils->total_distance($tree);
  my $rate_sum           = $ins_rate + $del_rate;
  my $max_allowed_indels = $rate_sum * $tree_length * 1;
  $max_allowed_indels = floor($max_allowed_indels);
  $max_allowed_indels = 1 if ( $max_allowed_indels < 1 );

  my $blocks_string = "1" x $aln->length;
  for ( my $i = 0 ; $i < $aln->length ; $i++ ) {
    my $indel_sum = $ins[$i] + $del[$i];

    if ( $indel_sum > $max_allowed_indels ) {
      my $indel_excess = ( $indel_sum / $max_allowed_indels ) * 2;
      my $indel_score  = 9 - $indel_excess;
      $indel_score = 0 if ( $indel_score < 0 );
      substr $blocks_string, $i, 1, floor($indel_score);
    } else {
      substr $blocks_string, $i, 1, '9';
    }
  }

  my @column_scores = split( "", $blocks_string );

  #  print "@column_scores\n";
  my %scores_hash;
  foreach my $leaf ( $tree->leaves ) {
    my $aln_string = _apply_columns_to_leaf( \@column_scores, $leaf );
    print "$aln_string\n";
    $scores_hash{ $leaf->stable_id } = $aln_string;
  }
  return \%scores_hash;
}

sub run_tcoffee {
  my $self   = shift;
  my $tree   = shift;
  my $aln    = shift;
  my $params = shift;

  Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $aln, { length => 200 } );

  my $node_id = $tree->node_id;
  my $tmpdir  = $self->worker_temp_directory . "${node_id}/";
  mkdir($tmpdir);
  my $filename = "$tmpdir" . "tcoffee_aln_${node_id}.fasta";
  my $tmpfile  = Bio::AlignIO->new(
    -file   => ">$filename",
    -format => 'fasta'
  );
  $tmpfile->write_aln($aln);
  $tmpfile->close;

  my $prefix = "export HOME_4_TCOFFEE=\"$tmpdir\";";
  $prefix .= "export DIR_4_TCOFFEE=\"$tmpdir\";";
  $prefix .= "export TMP_4_TCOFFEE=\"$tmpdir\";";
  $prefix .= "export CACHE_4_TCOFFEE=\"$tmpdir\";";
  $prefix .= "export NO_ERROR_REPORT_4_TCOFFEE=1;";
  $prefix .= "export MAFFT_BINARIES=/nfs/users/nfs_g/gj1/bin/mafft-bins/binaries;"
    ;    # GJ 2008-11-04. What a hack!

  my $outfile = $filename . ".score_ascii";

  my $cmd = qq^t_coffee -mode=evaluate -infile=$filename -outfile=$outfile -output=score_ascii -multi_core no^;
  print $cmd. "\n";
  system( $prefix. $cmd );

  my $scores_file = $outfile;
  my %score_hash;
  if ( -e $scores_file ) {
    my $FH = IO::File->new();
    $FH->open($scores_file) || throw("Could not open tcoffee scores file!");
    <$FH>;    #skip header
    my $i = 0;
    while (<$FH>) {
      $i++;

      #      next if ($i < 7); # skip first 7 lines.
      next if ( $_ =~ /^\s+/ );    #skip lines that start with space
      if ( $_ =~ /:/ ) {

        #my ($id,$overall_score) = split(/:/,$_);
        #$id =~ s/^\s+|\s+$//g;
        #$overall_score =~ s/^\s+|\s+$//g;
        #print "___".$id."___".$overall_score."___\n";
        next;
      }
      chomp;
      my ( $id, $align ) = split;
      $_ = $id;
      $id =~ s^/.*^^;
      $score_hash{$id} = '' if ( !exists $score_hash{$id} );
      $score_hash{$id} .= $align;

      #print $id." ". $align."\n";
    }
    $FH->close;
  }

  foreach my $leaf ( $tree->leaves ) {
    my $id     = $leaf->stable_id;
    my $string = $score_hash{$id};

    throw("No score string found!") unless ( defined $string );

    $string =~ s/[^\d-]/9/g
      ; # Convert non-digits and non-dashes into 9s. This is necessary because t_coffee leaves some leftover letters.
        #print $string."\n";
        #print $leaf->alignment_string."\n";
        #exit(0);
    $score_hash{$id} = $string;
  }

  return \%score_hash;
}

sub run_gblocks {
  my $self   = shift;
  my $tree   = shift;
  my $aln    = shift;
  my $params = shift;

  my $defaults = {
    t  => 'p',    # Type of sequence (p=protein,c=codon,d=dna)
    b3 => '8',    # Max # of contiguous nonconserved positions
    b4 => '5',    # Minimum length of a block
    b5 => 'a',    # Allow gap positions (n=none, h=with half,a=all)
  };
  my $use_params = $self->replace_params( $defaults, $params );

  printf("Sitewise_dNdS::run_gblocks\n") if ( $self->debug );

  my $aln_length = $aln->length;
  my $tmpdir     = $self->worker_temp_directory;
  my $filename   = "$tmpdir" . "gblocks_aln.fasta";
  my $tmpfile    = Bio::AlignIO->new(
    -file   => ">$filename",
    -format => 'fasta'
  );
  $tmpfile->write_aln($aln);
  $tmpfile->close;

  my @leaves             = $tree->leaves;
  my $num_leaves         = scalar(@leaves);
  my $min_leaves_gblocks = int( ( $num_leaves + 1 ) / 2 + 0.5 );

  #Example command: Gblocks 2138.fasta -t=p -b2=20 -b3=50 -b4=3 -b5=a -p=s
  my $cmd = sprintf(
    "Gblocks %s -t=%s -b1=%s -b2=%s -b3=%s -b4=%s -b5=%s -p=s\n",
    $filename,           $use_params->{'t'},  $min_leaves_gblocks, $min_leaves_gblocks,
    $use_params->{'b3'}, $use_params->{'b4'}, $use_params->{'b5'}
  );
  print "GBLOCKS: $cmd\n";
  my $ret = system("$cmd");
  open FLANKS, "$filename-gb.txts" or die "$!\n";
  my $segments_string;
  while (<FLANKS>) {
    chomp $_;
    next unless ( $_ =~ /Flanks:/ );
    $segments_string = $_;
    last;
  }
  close FLANKS;
  $segments_string =~ s/Flanks\: //g;
  $segments_string =~ s/\s+$//g;

  print "SEGS: " . $segments_string . "\n";
  $_ = $segments_string;
  my @bs = /\[(.+?)\]/g;

  my $blocks_string = "0" x $aln->length;
  foreach my $b (@bs) {
    my ( $start, $end ) = split( " ", $b );
    for ( my $i = $start ; $i <= $end ; $i++ ) {
      substr $blocks_string, $i - 1, 1, '9';
    }
  }

  my @column_scores = split( "", $blocks_string );
  my %scores_hash;

  foreach my $leaf ( $tree->leaves ) {
    my $aln_string = _apply_columns_to_leaf( \@column_scores, $leaf );
    $scores_hash{ $leaf->stable_id } = $aln_string;
  }

  return \%scores_hash;
}

sub _apply_columns_to_leaf {
  my $columns_ref = shift;
  my $leaf        = shift;

  my @column_scores = @{$columns_ref};

  my $aln_str = $leaf->alignment_string;
  my @aln_chars = split( "", $aln_str );

  for ( my $i = 0 ; $i < scalar(@aln_chars) ; $i++ ) {
    if ( $aln_chars[$i] ne '-' ) {
      $aln_chars[$i] = $column_scores[$i] || '0';
    }
  }

  return join( "", @aln_chars );
}

sub run_trimal {
  my $self   = shift;
  my $tree   = shift;
  my $aln    = shift;
  my $params = shift;

  # Write temporary alignment.
  my $dir   = $self->worker_temp_directory;
  my $aln_f = $dir . "/aln_" . $tree->node_id . ".fa";
  Bio::EnsEMBL::Compara::AlignUtils->to_file( $aln, $aln_f );

  # Build a command for TrimAl.
  my $trim_params = '';
  $trim_params = $params->{'trimal_filtering_params'}
    if ( $params->{'trimal_filtering_params'} );    # Custom parameters if desired.
  my $cmd = "trimal -gt 0.5 -cons 30 -in $aln_f -out $aln_f -colnumbering $trim_params";
  print "CMD: $cmd\n";
  my $output = `$cmd`;

  #print "OUTPUT:". $output."\n";

  chomp $output;
  my @cons_cols = split( /[,]/, $output );
  print join( ",", @cons_cols ) . "\n";

  my @column_scores = 0 x $aln->length;
  foreach my $column (@cons_cols) {
    $column_scores[$column] = 9;
  }

  my %scores_hash;
  foreach my $leaf ( $tree->leaves ) {
    my $aln_string = _apply_columns_to_leaf( \@column_scores, $leaf );
    $scores_hash{ $leaf->stable_id } = $aln_string;
  }

  return \%scores_hash;
}

sub run_prank {
  my $self   = shift;
  my $tree   = shift;
  my $aln    = shift;
  my $params = shift;

  my $threshold = $params->{'prank_filtering_threshold'} || 7;
  my ( $scores, $blocks ) =
    Bio::EnsEMBL::Compara::AlignUtils->get_prank_filter_matrices( $aln, $tree, $params );

  return $scores;
}

1;
