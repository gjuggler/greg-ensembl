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

  my $genes_f = $self->param('output_folder') . '/genes.Rdata';
  my $results_tbl = $self->param('output_folder') . '/table.csv'; 
  my $results_tbl_adj = $self->param('output_folder') . '/table_adj.csv'; 
  my $sites_f = $self->param('output_folder') . '/sites.Rdata';
  my $merged_f = $self->param('output_folder') . '/merged.Rdata'; 

  my $mysql_script = $self->base . "/scripts/mysql_functions.R";
  my $slrsim_script = $self->base . "/projects/slrsim/slrsim.functions.R";
  my $dbname = $self->dbc->dbname;
  
my $rcmd = qq^
source("${mysql_script}")
source("${slrsim_script}")
dbname="${dbname}"
con <- connect(dbname)

if (!file.exists("${genes_f}")) {
  genes <- dbGetQuery(con, "select * from genes;")
  save(genes, file="${genes_f}")
  rm(genes)
}

if (!file.exists("${results_tbl}")) {
  results <- dbGetQuery(con, "select * from slrsim_results;")
  write.csv(results, file="${results_tbl}", row.names=F)
  rm(results)
}

if (!file.exists("${results_tbl_adj}")) {
  results <- dbGetQuery(con, "select * from slrsim_results_adj;")
  write.csv(results, file="${results_tbl_adj}", row.names=F)
  rm(results)
}

if (!file.exists("${sites_f}")) {
  sites <- dbGetQuery(con, "select * from sites;")
  save(sites, file="${sites_f}")
  rm(sites)
}


if (!file.exists("${merged_f}")) {
  merged <- dbGetQuery(con, 
    "select label, tree, analysis, tree_length, ins_rate, aligner, filter, slrsim_rep, seq_position, lrt_stat, aln_dnds, true_dnds, true_type from merged;"
  )
  save(merged, file="${merged_f}")
  rm(merged)
}
^;
  Bio::Greg::EslrUtils->run_r( $rcmd);

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
