package Bio::Greg::Slrsim::DumpSlrsimAlignment;

use strict;
use Getopt::Long;
use IO::File;
use File::Basename;
use File::Path;

use Digest::MD5 qw(md5_hex);

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;

use base ('Bio::Greg::Slrsim::SlrsimProcess');

sub fetch_input {
  my($self) = @_;
    
  $self->load_all_params;
}

sub run {
  my $self = shift;

  my $obj = 
    Bio::EnsEMBL::Compara::ComparaUtils->get_all_as_obj( $self->compara_dba, $self->params );
  
  my $aln = $obj->{aln};
  my $scores = $obj->{score_aln};
  my $tree = $obj->{tree};

  my $folder = $self->get_output_folder . '/alns';

  my $label = $self->param('slrsim_label');
  $label =~ s/[\(\)\/\\]//gi;
  my $exp = $self->param('experiment_name');
  my $rep = $self->param('slrsim_rep');
  my $filter = $self->param('alignment_score_threshold');
  
  print "REP: $rep\n";

  return unless ($rep == 1);

  print "SAVING!\n";

  throw("Folder doesn't exist!") unless (-d $folder);

  my $output_id = "${rep}_${exp}_${label}";
  if ($self->param('aln_output_id')) {
    $output_id = $self->param('aln_output_id');
  }

  my $base = "${folder}/${output_id}";

  print "$base\n";

  my $tree_f = "${base}.nh";
  Bio::EnsEMBL::Compara::TreeUtils->to_file($tree,$tree_f);
  my $aln_f = "${base}.fasta";
  Bio::EnsEMBL::Compara::AlignUtils->to_file($aln,$aln_f);

  # Output the params.
  my $params_f = "${base}.txt";
  my $out;
  open($out,">${params_f}");
  $self->hash_print($self->params,$out);
  close($out);
  
  # Output a scores string only for the most-filtered seqs.
  if ($filter == 9) {
    my $score_tree_f = "${base}_scores.nh";
    Bio::EnsEMBL::Compara::TreeUtils->to_file($tree,$score_tree_f);
    my $scores_f = "${base}_scores.fasta";
    Bio::EnsEMBL::Compara::AlignUtils->to_file($scores,$scores_f);
  }

}

1;
