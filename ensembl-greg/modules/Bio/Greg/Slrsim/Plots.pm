package Bio::Greg::Slrsim::Plots;

use strict;
use Cwd;
use POSIX qw(strftime mktime);
use File::Path;
use Bio::Greg::EslrUtils;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::Process');

sub fetch_input {
  my $self = shift;

  my $defaults = {};

  $self->load_all_params($defaults);
  $self->get_output_folder();
}

sub run {
  my $self = shift;

  my $experiment_name = $self->param('experiment_name');

  if ($experiment_name =~ m/filter_sweeps/i) {
    $self->output_params_file;
    eval {
      $self->filter_sweep_roc;
    };
  }

}

sub get_output_folder {
  my $self = shift;

  my $date_string = strftime("%Y-%m-%d",localtime);
  my $i=0;
  
  my $filename;
  do {
    $i++;
    $filename = "output/".$date_string."_".sprintf("%.2d",$i);
  } while (-e $filename);

  print "Output folder: $filename\n";
  $self->param('output_folder',$filename);
  mkpath([$filename]);
}

sub filter_sweep_roc {
  my $self = shift;

  my $filename = $self->param('output_folder') . '/filter_sweep_roc.png';

  my $rcmd = qq^
source("collect_slrsim.R")
data = get.all.data()

comb.roc = data.frame()
for (aln in sort(unique(data['alignment_name']))) {
  for (filt in sort(unique(data['filtering_name']))) {
    for (thresh in sort(unique(data['alignment_score_threshold']))) {
        sub = subset(data, alignment_name==aln & filtering_name==filt & alignment_score_threshold==thresh)
        print(paste(aln,filt,thresh,nrow(sub),sep="/"))
        if (nrow(sub)==0) {next}
        sub.roc = slr.roc(sub)
        sub.roc['filter'] = as.factor(filt)
        sub.roc['aln'] = as.factor(aln)
        sub.roc['thresh'] = as.factor(thresh)
        sub.roc['label'] = paste(aln,filt,thresh,sep="/")
        comb.roc = rbind(sub.roc,comb.roc)      
    }
  }
}

png(file="${filename}",width=600,height=600)
plot.roc(comb.roc)
dev.off()
^;
  print "$rcmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r($rcmd,$params);
}

sub output_params_file {
  my $self = shift;

  my $filename = $self->param('output_folder') . '/params.txt';
  print "$filename\n";
  my $out;
  open($out,">$filename");

  $self->hash_print($self->params,$out);
  
  close($out);
}

1;
