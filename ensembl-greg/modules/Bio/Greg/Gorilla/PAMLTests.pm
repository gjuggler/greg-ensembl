package Bio::Greg::Gorilla::PAMLTests;

use strict;
use File::Path;
use Bio::Greg::Codeml;
use FreezeThaw qw(freeze thaw cmpStr safeFreeze cmpStrHard);
use Data::Types qw(:all);

use base (
  'Bio::Greg::StatsCollectionUtils',
  'Bio::Greg::Hive::Process',
);

sub param_defaults {
  return {
  };
}

sub fetch_input {
  my ($self) = @_;

  # Fetch parameters from all possible locations.
  $self->load_all_params();

  # Create tables if necessary.
  $self->create_table_from_params( $self->compara_dba, 'results',
    $self->results_table );
}

sub run {
  my $self = shift;

  $self->param('data_id', $self->param('lrt_lo').' to '.$self->param('lrt_hi'));

  my $aln = $self->_collect_aln;
}

sub _collect_aln {
  my $self = shift;

  my $f = $self->_save_file('codon_aln', 'fasta');
  print $f->{full_file}."\n";
  if (!-e $f->{full_file} || $self->param('force_recalc')) {
    my $lo = $self->param('lrt_lo');
    my $hi = $self->param('lrt_hi');

    print "  collecting aln from $lo to $hi...\n";

    my $query = qq^
SELECT %s AS string FROM 
  gj1_gorilla.seqs seqs JOIN gj1_gorilla.sites sites
  ON (seqs.data_id=sites.data_id AND seqs.aln_position=sites.aln_position)
  WHERE sites.lrt_stat*sign(sites.omega-.9999) > ? AND
        sites.lrt_stat*sign(sites.omega-.9999) < ?
  ^;
    
    my $aln = new Bio::SimpleAlign;
    foreach my $species ('h', 'c', 'g', 'o', 'm', 'r') {
      my $str = '';

      my $cur_query = sprintf($query, $species);
      my $sth = $self->dbc->prepare($cur_query);
      $sth->execute($lo, $hi);
      while ( my $obj = $sth->fetchrow_hashref ) {
        $str .= $obj->{string};
      }
      $sth->finish;

      my $bioseq = Bio::LocatableSeq->new( 
        -seq => $str,
        id => $species
        );
      $aln->add_seq($bioseq);
    }
    Bio::EnsEMBL::Compara::AlignUtils->to_file($aln, $f->{full_file});
    print "  done!\n";
  }
  print "  loading aln...\n";
  my $aln = Bio::EnsEMBL::Compara::AlignUtils->from_file($f->{full_file});
  print "  done!\n";

  $self->pretty_print($aln);
  $self->pretty_print($self->_tx_aln($aln));

  # Do any necessary filtering of the alignment - remove codons with gaps or Ns,
  # for example...
  $aln = Bio::EnsEMBL::Compara::AlignUtils->remove_regex_columns_in_threes($aln, "[-N]");

  $self->pretty_print($aln);
  $self->pretty_print($self->_tx_aln($aln));

  die "Alignment not a multiple of 3!" unless ($aln->length % 3 == 0);

  my $tree = Bio::EnsEMBL::Compara::TreeUtils->from_newick("(((((h,c),g),o),m),r);");
  print $tree->ascii."\n";

  $self->_run_paml($tree, $aln);

  return $aln;
}

sub _run_paml {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);

  my $res_f = $self->_save_file('paml_results', 'perlobj');
  if (!-e $res_f->{full_file} || $self->param('force_recalc')) {

    my $m0_params = {
      model => 0,
      fix_blength => 0,
      cleandata => 0,
      getSE => 1,
      Small_Diff => 1e-7
    };
    print "  running PAML m0...\n";
    my $m0 = Bio::Greg::Codeml->branch_model_likelihood( $treeI, $aln, $self->worker_temp_directory, $m0_params );

    my $m0_tree = Bio::Greg::Codeml->parse_codeml_results($m0->{lines});
    Bio::EnsEMBL::Compara::TreeUtils->transfer_branchlengths($m0_tree, $treeI);
    $treeI->root->branch_length(0.01); # Set small b.l. on the root node.

    print $treeI->newick."\n";

    my $branches_params = {
      model => 1,
      fix_blength => 0,
      method => 1,
      cleandata => 0,
      getSE => 1,
      Small_Diff => 1e-7
    };
    print "  running PAML branches...\n";
    my $branches = Bio::Greg::Codeml->branch_model_likelihood( $treeI, $aln, $self->worker_temp_directory, $branches_params);

    my $res_hash = {
      m0 => $m0->{lines},
      branches => $branches->{lines}
    };
    $self->frz($res_f->{full_file}, $res_hash);
  }
  print"  loading PAML results...\n";
  my $res_hash = $self->thw($res_f->{full_file});

  my $branches = $res_hash->{branches};
  my $m0 = $res_hash->{m0};

  my $branch_f = $self->_save_file('paml_branch', 'txt');
  open(OUT, ">".$branch_f->{full_file});
  print OUT join("\n", @{$branches}) . "\n";
  close(OUT);

  my $m0_f = $self->_save_file('paml_m0', 'txt');
  open(OUT, ">".$m0_f->{full_file});
  print OUT join("\n", @{$m0}) . "\n";
  close(OUT);

  my @omegas = Bio::Greg::Codeml->extract_omegas($m0);
  print "@omegas\n";
  my $branches_tree = Bio::Greg::Codeml->parse_codeml_results($branches);
  
  foreach my $node ($branches_tree->nodes) {
    print $node->enclosed_leaves_string."\n";
    $self->hash_print($node->get_tagvalue_hash);
  }

}

sub _tx_aln {
  my $self = shift;
  my $aln = shift;
  return Bio::EnsEMBL::Compara::AlignUtils->translate($aln);
}

sub results_table {
  return {
    data_id => 'int',
    unique_keys => 'data_id',
  };
}

sub _save_file {
  my $self = shift;
  my $filename_base = shift;
  my $ext = shift;

  my $id = $self->param('data_id');
  $id =~ s/[-]//g;
  $id =~ s/[\. ]/_/g;

  my $filename = "${id}_${filename_base}";

  my $subfolder = 'data';
  if ($self->param('aln_use_type') ne 'genomic_all') {
    $subfolder = 'data/' . $self->param('aln_use_type');
  }

  my $file_params = {
    id => $id,
    filename => $filename,
    extension => $ext,
    subfolder => $subfolder,
  };

  my $file_obj = $self->save_file($file_params);
  if (!defined $self->param('data_prefix')) {
    $self->store_param('data_prefix', $file_obj->{hash_folder});
  }

  return $file_obj;
}

1;
