package Bio::Greg::AlignmentScores;

use strict;
use Cwd;
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::NestedSet;
use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Hive::Process;

use Bio::Greg::ProcessUtils;

our @ISA = qw(Bio::EnsEMBL::Hive::Process Bio::Greg::ProcessUtils);

#
# Some global-ish variables.
#
my $dba;
my $pta;

# INPUT FILES / OBJECTS.
my $tree;
my $params;

# OUTPUT FILES / OBJECTS / STATES.
my %score_hash;

sub debug {1;}

sub fetch_input {
  my( $self) = @_;
  
  # Load up the Compara DBAdaptor.
  $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-DBCONN=>$self->db->dbc);
  $pta = $dba->get_ProteinTreeAdaptor;
  
  ### DEFAULT PARAMETERS ###
  $params = {
    alignment_table            => 'protein_tree_member',
    alignment_score_table => 'protein_tree_member_score',
    alignment_scores_action    => 'prank'             # Options: 'gblocks', 'prank', 'trimal', 'indelign'
    };
  
  # For aminof, codonf, and freqtype, see the SLR readme.txt for more info.

  #########################

  # Fetch parameters from the two possible locations. Input_id takes precedence!
  my $p_params = $self->get_params($self->parameters);
  my $i_params = $self->get_params($self->input_id);
  my $node_id = $i_params->{'protein_tree_id'};
  $node_id = $i_params->{'node_id'} if (!defined $node_id);
  my $t_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_tree_tags($dba,$node_id);
  
  $params = $self->replace_params($params,$p_params,$i_params,$t_params);

  $dba->dbc->disconnect_when_inactive(1);
}

sub run {
  my $self = shift;
  
  my $node_id = $params->{'node_id'};

  $tree = Bio::EnsEMBL::Compara::ComparaUtils->get_tree_for_comparative_analysis($dba,$params);
  $tree = $tree->minimize_tree;
  
  my $input_table = $params->{'alignment_table'};
  my $action = $params->{'alignment_scores_action'};

  if ($action =~ m/gblocks/i) {
    print " -> RUN GBLOCKS\n";
    my $score_hash = $self->run_gblocks($tree,$params);
    $self->store_scores($tree,$score_hash,$input_table."_".'gblocks');
  }

  if ($action =~ m/prank/i) {
    print " -> RUN PRANK\n";
    my $score_hash = $self->run_prank($tree,$params);
    $self->store_scores($tree,$score_hash,$input_table."_".'prank');
  }

  if ($action =~ m/trimal/i) {
    print " -> RUN TRIMAL\n";
    my $score_hash = $self->run_trimal($tree,$params);
    $self->store_scores($tree,$score_hash,$input_table."_".'trimal');
  }    

  if ($action =~ m/indelign/i) {
    print " -> RUN INDELIGN\n";
    eval {
      my $score_hash = $self->run_indelign($tree,$params);
      $self->store_scores($tree,$score_hash,$input_table."_".'indelign');
    };
    if ($@) {
      print "Indelible error: $@\n";
    }
  }

}

sub store_scores {
  my $self = shift;
  my $tree = shift;
  my $score_hash = shift;
  my $output_table = shift;

  $dba->dbc->do("CREATE TABLE IF NOT EXISTS $output_table LIKE protein_tree_member_score");
  my $sth = $tree->adaptor->prepare("INSERT INTO $output_table (node_id,member_id,cigar_line) VALUES (?,?,?) ON DUPLICATE KEY UPDATE cigar_line=?");
  foreach my $leaf ($tree->leaves) {
    my $score_string = $score_hash->{$leaf->stable_id};
    $sth->execute($leaf->node_id,
                  $leaf->member_id,
                  $score_string,
                  $score_string);
    print $score_string."\n";
  }
  $sth->finish;
}

sub run_indelign {
  my $self = shift;
  my $tree = shift;
  my $params = shift;

  my $aln = $tree->get_SimpleAlign();

  my ($i_obj,$d_obj,$ins_rate,$del_rate) = Bio::EnsEMBL::Compara::AlignUtils->indelign($aln,$tree,$params,$self->worker_temp_directory);
  my @ins = @$i_obj;
  my @del = @$d_obj;

  my $num_leaves = scalar ($tree->leaves);
  my $tree_length = Bio::EnsEMBL::Compara::TreeUtils->total_distance($tree);
  my $rate_sum = $ins_rate+$del_rate;

  my $blocks_string = "0" x $aln->length;  
  for (my $i=0; $i < $aln->length; $i++) {
    my $indel_sum = $ins[$i] + $del[$i];
    my $max_allowed_indels = $rate_sum*$tree_length*2;
    if ($indel_sum >= $max_allowed_indels) {
      substr $blocks_string,$i,1,'0';
    } else {
      substr $blocks_string,$i,1,'9';
    }
  }

  my @column_scores = split("",$blocks_string);  
  my %scores_hash;
  foreach my $leaf ($tree->leaves) {
    my $aln_string = _apply_columns_to_leaf(\@column_scores,$leaf);
    $scores_hash{$leaf->stable_id} = $aln_string;
  }
  return \%scores_hash;
}

sub run_gblocks {
  my $self = shift;
  my $tree = shift;
  my $params = shift;

  my $aln = $tree->get_SimpleAlign();

  my $defaults = {
    t    => 'p',        # Type of sequence (p=protein,c=codon,d=dna)
    b3   => '8',       # Max # of contiguous nonconserved positions
    b4   => '5',        # Minimum length of a block
    b5   => 'a',        # Allow gap positions (n=none, h=with half,a=all)
  };
  my $use_params = $self->replace_params($defaults,$params);

  printf("Sitewise_dNdS::run_gblocks\n") if($self->debug);

  my $aln_length = $aln->length;
  my $tmpdir = $self->worker_temp_directory;
  my $filename = "$tmpdir". "gblocks_aln.fasta";
  my $tmpfile = Bio::AlignIO->new
    (-file => ">$filename",
     -format => 'fasta');
  $tmpfile->write_aln($aln);
  $tmpfile->close;

  my @leaves = $tree->leaves;
  my $num_leaves = scalar(@leaves);
  my $min_leaves_gblocks = int(($num_leaves+1)/2 + 0.5);
  #Example command: Gblocks 2138.fasta -t=p -b2=20 -b3=50 -b4=3 -b5=a -p=s
  my $cmd = sprintf("Gblocks %s -t=%s -b1=%s -b2=%s -b3=%s -b4=%s -b5=%s -p=s\n",
		    $filename,
		    $use_params->{'t'},
		    $min_leaves_gblocks,
		    $min_leaves_gblocks,
		    $use_params->{'b3'},
		    $use_params->{'b4'},
		    $use_params->{'b5'});
  print "GBLOCKS: $cmd\n";
  my $ret = system("$cmd");
  open FLANKS, "$filename-gb.txts" or die "$!\n";
  my $segments_string;
  while (<FLANKS>) {
    chomp $_;
    next unless ($_ =~ /Flanks:/);
    $segments_string = $_;
    last;
  }
  close FLANKS;
  $segments_string =~ s/Flanks\: //g;
  $segments_string =~ s/\s+$//g;

  print "SEGS: ". $segments_string."\n";
  $_ = $segments_string;
  my @bs = /\[(.+?)\]/g;

  my $blocks_string = "0" x $aln->length;
  foreach my $b (@bs) {
    my ($start,$end) = split(" ",$b);
    for (my $i=$start; $i <= $end; $i++) {
      substr $blocks_string,$i-1,1,'9';
    }
  }

  my @column_scores = split("",$blocks_string);  
  my %scores_hash;

  foreach my $leaf ($tree->leaves) {
    my $aln_string = _apply_columns_to_leaf(\@column_scores,$leaf);
    $scores_hash{$leaf->stable_id} = $aln_string;
  }

  return \%scores_hash;
}

sub _apply_columns_to_leaf {
  my $columns_ref = shift;
  my $leaf = shift;

  my @column_scores = @{$columns_ref};

  my $aln_str = $leaf->alignment_string;
  my @aln_chars = split("",$aln_str);
  
  for (my $i=0; $i < scalar(@aln_chars); $i++) {
    if ($column_scores[$i] < 9 && $aln_chars[$i] ne '-') {
      $aln_chars[$i] = 0;
    } elsif ($aln_chars[$i] ne '-') {
      $aln_chars[$i] = 9;
    }
  }
  return join("",@aln_chars);
}

sub run_trimal {
  my $self = shift;
  my $tree = shift;
  my $params = shift;

  my $aln = $tree->get_SimpleAlign();

  # Write temporary alignment.
  my $dir = $self->worker_temp_directory;
  my $aln_f = $dir."/aln_".$tree->node_id.".fa";
  Bio::EnsEMBL::Compara::AlignUtils->to_file($aln,$aln_f);
  
  # Build a command for TrimAl.
  my $trim_params = '-cons 30 -gt 0.5 -w 2';
  $trim_params = $params->{'trimal_filtering_params'} if ($params->{'trimal_filtering_params'}); # Custom parameters if desired.
  my $cmd = "trimal -in $aln_f -out $aln_f -colnumbering $trim_params";
  my $output = `$cmd`;
  #print "OUTPUT:". $output."\n";
  
  chomp $output;
  my @cons_cols = split(/[,]/,$output);
  print join(",",@cons_cols)."\n";

  my @column_scores = 0 x $aln->length;
  foreach my $column (@cons_cols) {
    $column_scores[$column] = 9;
  }

  my %scores_hash;
  foreach my $leaf ($tree->leaves) {
    my $aln_string = _apply_columns_to_leaf(\@column_scores,$leaf);
    $scores_hash{$leaf->stable_id} = $aln_string;
  }

  return \%scores_hash;
}

sub run_prank {
  my $self = shift;
  my $tree = shift;
  my $params = shift;

  my $aln = $tree->get_SimpleAlign();
  my $threshold = $params->{'prank_filtering_threshold'} || 5;
  my ($scores,$blocks) = Bio::EnsEMBL::Compara::AlignUtils->get_prank_filter_matrices($aln,$tree,$params);

  return $scores;
}

1;
