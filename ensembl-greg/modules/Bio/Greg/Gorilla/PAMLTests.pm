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

  $self->{_table_structure} = $self->results_table;
}


sub write_output {
  my $self = shift;

  $self->create_table_from_params( $self->dbc, 'results', $self->{_table_structure});
  $self->store_params_in_table($self->dbc, 'results', $self->params);

}

sub run {
  my $self = shift;

  $self->param('data_id', $self->param('lrt_lo').' to '.$self->param('lrt_hi'));

  my $tree = Bio::EnsEMBL::Compara::TreeUtils->from_newick("(((((h,c),g),o),m),r);");
  print $tree->ascii."\n";

  my $filtered_aln = $self->_collect_aln('seqs');
  my $unfiltered_aln = $self->_collect_aln('seqs_nofilters');

  $self->pretty_print($self->_tx_aln($filtered_aln));
  $self->pretty_print($self->_tx_aln($unfiltered_aln));

  $self->store_param('aln_length', $filtered_aln->length);
  $self->store_param('aln_length_nofilt', $unfiltered_aln->length);

  $self->_run_paml($tree, $filtered_aln, '');
  $self->_run_paml($tree, $unfiltered_aln, 'nofilt');
}

sub _collect_aln {
  my $self = shift;
  my $table = shift;

  my $f = $self->_save_file("${table}_aln", 'fasta');
  print $f->{full_file}."\n";

  if (!-e $f->{full_file} || $self->param('force_recalc')) {
    my $lo = $self->param('lrt_lo');
    my $hi = $self->param('lrt_hi');

    print "  collecting aln from ${table} [$lo to $hi]...\n";

    my $query = qq^
SELECT %s AS string FROM 
  gj1_gorilla.${table} seqs JOIN gj1_gorilla.sites sites
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
  print "  loading aln from file...\n";
  my $aln = Bio::EnsEMBL::Compara::AlignUtils->from_file($f->{full_file});

  # Do any necessary filtering of the alignment - remove codons with gaps or Ns,
  # for example...
  $aln = Bio::EnsEMBL::Compara::AlignUtils->remove_regex_columns_in_threes($aln, "[-N]");

  die "Alignment not a multiple of 3!" unless ($aln->length % 3 == 0);

  return $aln;
}

sub _run_paml {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $prefix = shift;

  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);

  my $branch_f = $self->_save_file("paml_branch_${prefix}", 'txt');
  my $m0_f = $self->_save_file("paml_m0_${prefix}", 'txt');

  $prefix .= '_' if ($prefix ne '');

  if (!-e $m0_f->{full_file} || $self->param('force_recalc')) {
    my $m0_params = {
      model => 0,
      fix_blength => 0,
      cleandata => 0,
      getSE => 1,
      Small_Diff => 1e-7
    };
    print "  running PAML m0...\n";
    my $m0 = Bio::Greg::Codeml->branch_model_likelihood( $treeI, $aln, $self->worker_temp_directory, $m0_params );
    
    open(OUT, ">".$m0_f->{full_file});
    print OUT join("\n", @{$m0->{lines}}) . "\n";
    close(OUT);
  }

  print "  loading m0 results...\n";
  open(IN, $m0_f->{full_file});
  my @m0_lines = <IN>;
  close(IN);

  my ($dnds, $dnds_se) = Bio::Greg::Codeml->parse_m0_dnds(\@m0_lines);
  $self->store_param($prefix.'m0_dnds', $dnds);
  $self->store_param($prefix.'m0_dnds_se', $dnds_se);
  my $m0_tree = Bio::Greg::Codeml->parse_codeml_results(\@m0_lines);
  Bio::EnsEMBL::Compara::TreeUtils->transfer_branchlengths($m0_tree, $treeI);
  $treeI->root->branch_length(0.01); # Set small b.l. on the root node.
  print $treeI->newick."\n";
  
  if (!-e $branch_f->{full_file} || $self->param('force_recalc')) {
    my $branch_params = {
      model => 1,
      fix_blength => 1,
      method => 1,
      cleandata => 0,
      getSE => 1,
      Small_Diff => 1e-6
    };
    print "  running PAML branches...\n";
    my $branch = Bio::Greg::Codeml->branch_model_likelihood( $treeI, $aln, $self->worker_temp_directory, $branch_params );
    
    open(OUT, ">".$branch_f->{full_file});
    print OUT join("\n", @{$branch->{lines}}) . "\n";
    close(OUT);
  }

  print "  loading branch results...\n";
  open(IN, $branch_f->{full_file});
  my @branch_lines = <IN>;
  close(IN);
  my $branch_tree = Bio::Greg::Codeml->parse_codeml_results(\@branch_lines);
  
  foreach my $node ($branch_tree->nodes) {
    $self->hash_print($node->get_tagvalue_hash);
    my $enclosed = $node->enclosed_leaves_string('');
    $self->store_param($prefix.$enclosed.'_dnds', 0.0+$node->get_tag_value('dN/dS'));
    $self->store_param($prefix.$enclosed.'_dnds_se', 0.0+$node->get_tag_value('dN/dS_se'));
    $self->store_param($prefix.$enclosed.'_ds', 0.0+$node->get_tag_value('dS'));
  }

}

sub _tx_aln {
  my $self = shift;
  my $aln = shift;
  return Bio::EnsEMBL::Compara::AlignUtils->translate($aln);
}

sub results_table {
  my $p = {
    data_id => 'int',
    unique_keys => 'data_id'
  };
  return $p;
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
