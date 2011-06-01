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

  $self->load_all_params();

  #$self->hash_print($self->params);
  $self->_output_folder();
}

sub run {
  my $self = shift;

  $self->output_params_file;
  $self->dump_sql;
  $self->dump_data;
#  $self->dump_alns;

  # Call the method corresponding to the current experiment.
  my $experiment_name = $self->param('experiment_name');
  $self->$experiment_name();
}

sub dump_sql {
  my $self     = shift;
  my $filename = $self->param('output_folder') . '/slrsim.sqldata';
  my $gzip     = $self->param('output_folder') . '/slrsim.sqldata.gz';

  my $dbc  = $self->compara_dba->dbc;
  my $u    = $dbc->username;
  my $p    = $dbc->password;
  my $h    = $dbc->host;
  my $port = $dbc->port;
  my $db   = $dbc->dbname;

  if (!-e $filename && !-e $gzip) {
    print "  mysqldumping...\n";
    my $cmd = qq^mysqldump -P$port -h$h -u$u -p$p $db > $filename;^;
    system($cmd);
  } else {
    print "  mysqldump already exists\n";
  }

  if (!-e $gzip) {
    print "  gzipping...\n";
    my $cmd = qq^gzip $filename;^;
    system($cmd);
    unlink($filename);
  } else {
    print "  gzip already exists\n";
  }
}

sub dump_alns {
  my $self = shift;

  my $folder = $self->get_output_folder;

  my $alns_gzip = "$folder/alns.tar.gz";
  if (!-e $alns_gzip) {
    my $cwd = cwd();
    chdir($folder);
    system('tar -cvvf alns.tar data; gzip alns.tar');
    chdir($cwd);
  }
}

sub dump_data {
  my $self = shift;

  my $sites_f = $self->param('output_folder') . '/sites.Rdata';
  my $genes_f = $self->param('output_folder') . '/genes.Rdata';
  my $collect_script   = $self->collect_script;
  my $functions_script   = $self->functions_script;

  my $dbname = $self->dbc->dbname;

  if ( !-e $sites_f || !-e $genes_f) {
    my $rcmd = qq^
dbname="$dbname"
source("${collect_script}")
source("${functions_script}")
sites <- get.all.data()
genes <- get.all.genes()
save(sites, genes, file="${sites_f}");
save(genes, file="${genes_f}");
^;
    #print "$rcmd\n";
    my $params = {};
    Bio::Greg::EslrUtils->run_r( $rcmd, $params );
  }

}

# Slrsim table: generic summary table output.
sub slrsim_table {
  my $self = shift;

  my $folder = $self->get_output_folder;
  my $file = "${folder}/sites.Rdata";
  my $functions = $self->base . "/projects/slrsim/slrsim.functions.R";
  my $table_file = "${folder}/table.csv";

my $rcmd = qq^
### Plots::slrsim_table

source("${functions}")
load("${file}")

sites[, 'slrsim_label'] <- paste( sites[, 'slrsim_label'], sites[, 'analysis'], sep='_')

paper.df <- ddply(sites, .(slrsim_label), paper.table)
paper.df <- format.numeric.df(paper.df, digits=3)
write.csv(paper.df, file="${table_file}", row.names=F, quote=F)
^;
  Bio::Greg::EslrUtils->run_r( $rcmd );

}

sub fig_zero {
  my $self = shift;

  my $folder = $self->get_output_folder;

  # Call slrsim_table to dump the table.
  my $table_file = "${folder}/table.csv";
  if (!-e $table_file) {
    $self->slrsim_table;
  }  

}


sub fig_one_a {
  my $self = shift;

  $self->fig_one;
}

sub fig_one_b {
  my $self = shift;

  $self->fig_one;
}

sub fig_one_c  {
  my $self = shift;

  $self->fig_one;
}

sub fig_one {
  my $self = shift;
  my $folder = $self->get_output_folder;

  # Call slrsim_table to dump the table.
  my $table_file = "${folder}/table.csv";
  if (!-e $table_file) {
    $self->slrsim_table;
  }

}

sub fig_two_a {
  my $self = shift;
  $self->fig_two;
}

sub fig_two_b {
  my $self = shift;
  $self->fig_two;
}

sub fig_two_c {
  my $self = shift;
  $self->fig_two;
}

sub fig_two_d {
  my $self = shift;
  $self->fig_two;
}

sub fig_two_e {
  my $self = shift;
  $self->fig_two;
}

sub fig_two_f {
  my $self = shift;
  $self->fig_two;
}

sub fig_two_h {
  my $self = shift;
  $self->fig_two;
}

sub fig_two_clustalw {
  my $self = shift;
  $self->fig_two;
}

sub fig_two_fsa {
  my $self = shift;
  $self->fig_two;
}

sub fig_three_a {
  my $self = shift;
  $self->fig_two;
}

sub fig_two {
  my $self = shift;

  my $folder = $self->get_output_folder;

  # Call slrsim_all to dump the data.
  my $sites_file = "${folder}/sites.Rdata";
  if (!-e $sites_file) {
    $self->slrsim_all;
  }

  # Call slrsim_table to dump the table.
  my $table_file = "${folder}/table.csv";
  if (!-e $table_file) {
    $self->slrsim_table;
  }
}

sub ari_indels {
  my $self = shift;

  my $folder = $self->get_output_folder;

  # Call slrsim_all to dump the data.
  my $sites_file = "${folder}/sites.Rdata";
  if (!-e $sites_file) {
    $self->slrsim_all;
  }

  # Call slrsim_table to dump the table.
  my $table_file = "${folder}/table.csv";
  if (!-e $table_file) {
    $self->slrsim_table;
  }
  
}

sub mammals_sim {
  my $self = shift;

  my $folder = $self->get_output_folder;

  # Call slrsim_all to dump the data.
  my $sites_file = "${folder}/sites.Rdata";
  if (!-e $sites_file) {
    $self->slrsim_all;
  }

  # Call slrsim_table to dump the table.
  my $table_file = "${folder}/table.csv";
  if (!-e $table_file) {
    $self->slrsim_table;
  }
  
}


# Filter table: return a table of results with multiple filtering results on each line.
sub filter_table {
  my $self = shift;

  my $folder = $self->get_output_folder;
  my $file = "${folder}/sites.Rdata";
  my $functions = $self->base . "/projects/slrsim/slrsim.functions.R";
  my $fdr_tpr = "${folder}/fdr_tpr_collated.csv";
  my $fpr_tpr = "${folder}/tpr_tpr_collated.csv";
  my $cor = "${folder}/cor_collated.csv";
  my $auc = "${folder}/auc_collated.csv";

my $rcmd = qq^
source("${functions}")

### Plots::filter_table

load("${file}")

x <- sites
x[, 'slrsim_label'] <- paste(x[, 'slrsim_tree_file'], x[,'alignment_name'], x[,'filtering_name'])
table.df <- ddply(x, .(slrsim_label, phylosim_insertrate, tree_mean_path), paper.table)
paper.df <- table.df

for (col in c('auc', 'fpr_tpr', 'fdr_tpr', 'cor')) {
  filter.df <- ddply(paper.df,
    .(phylosim_insertrate,tree_mean_path,tree,aligner),
    function(x) {
      filter.list <- dlply(x,.(filter),function(y) {y[1,col]})
      return(data.frame(x[1,],filter.list))
    })
  write.csv(filter.df,file=paste("${folder}/",col,"_collated.csv",sep=""),row.names=F,quote=F)
}

^;
  Bio::Greg::EslrUtils->run_r( $rcmd );  
}

sub output_params_file {
  my $self = shift;

  my $filename = $self->param('output_folder') . '/params.txt';

  if ( !-e $filename ) {
    print "$filename\n";
    my $out;
    open( $out, ">$filename" );
    $self->hash_print( $self->params, $out );
    close($out);
  }
}

sub _output_folder {
  my $self = shift;
  return $self->get_output_folder;
}

sub collect_script {
  my $self = shift;
  return $self->base . "/projects/slrsim/collect_slrsim.R";
}

sub functions_script {
  my $self = shift;
  return $self->base . "/projects/slrsim/slrsim.functions.R";
}

1;
