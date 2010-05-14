package Bio::Greg::Hive::DumpAlignment;

use strict;
use Getopt::Long;
use IO::File;
use File::Basename;
use File::Path;

use Digest::MD5 qw(md5_hex);

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;

use base ('Bio::Greg::Hive::Process');

#
# Some global-ish variables.
#

sub fetch_input {
  my($self) = @_;
  
  ### DEFAULT PARAMETERS ###
  my $params = {
    output_folder => '',
    output_aln_aa => 1,
    output_aln_cdna => 1,
    output_tree => 1,
    hash_subfolders => 1,
  };
  
  #########################
  
  $self->load_all_params($params);
  
}

sub run {
  my $self = shift;
  
  my $tree = $self->get_tree();
  my $aln = $self->get_aln();
  ($aln) = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns($aln);
  my $cdna_aln = $self->get_cdna_aln();
  ($cdna_aln) = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns_in_threes($cdna_aln);
  
  my $output_folder = $self->param('output_folder') . '/' . $self->get_subfolder;
  mkpath([$output_folder]);
  die("No output folder set!") unless ($output_folder ne '');
  
  my $stable_id = $self->get_stable_id($tree);
  my $output_tree = "${output_folder}/${stable_id}.nh"; 
  my $output_aln_aa = "${output_folder}/${stable_id}_aa.fasta"; 
  my $output_aln_cdna = "${output_folder}/${stable_id}_cdna.fasta"; 

  # Output alignment.
  print "Outputting files ...\n";

  if ($self->param('output_aln_aa')) {
    print " -> $output_aln_aa\n";
    Bio::EnsEMBL::Compara::AlignUtils->pretty_print($aln, { length => 200 } );
    Bio::EnsEMBL::Compara::AlignUtils->to_file($aln,$output_aln_aa); # Write the alignment out to file.
  }

  if ($self->param('output_aln_cdna')) {
    print " -> $output_aln_cdna\n";
    Bio::EnsEMBL::Compara::AlignUtils->pretty_print($cdna_aln, { length => 200 } );
    Bio::EnsEMBL::Compara::AlignUtils->to_file($cdna_aln,$output_aln_cdna); # Write the alignment out to file.
  }
  
  if ($self->param('output_tree')) {
    print " -> $output_tree\n";    
    my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
    Bio::EnsEMBL::Compara::TreeUtils->to_file($treeI,$output_tree);
  }
}

sub get_stable_id {
  my $self = shift;
  my $tree = shift;

  my $stable_id;

  my @human_proteins = grep { $_->taxon_id == 9606 } $tree->leaves;
  if ( scalar @human_proteins > 0 ) {
    my $member = $human_proteins[0];
    $stable_id = $member->stable_id;
  } else {
    my @leaves = $tree->leaves;
    my $member = $leaves[0];
    $stable_id = $member->stable_id;
  }
  return $stable_id;
}

sub get_subfolder {
  my $self = shift;

  my $id = $self->data_id;
  
  if ($self->param('hash_subfolders')) {
    my $n = $self->param('hash_subfolders');
    my $lo = 0;
    my $hi = $n-1;

    my (@md5) = md5_hex($id) =~ /\G(..)/g;
    return join('/',@md5[$lo .. $hi]);
  } else {
    return '';
  }
}

1;
