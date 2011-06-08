package Bio::Greg::Slrsim::CalculateResults;

use strict;
use Cwd;
use POSIX qw(strftime mktime);
use File::Path;
use Bio::Greg::EslrUtils;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::Process');

sub param_defaults {
  return {
    slrsim_label => 'NO LABEL'
  };
}

sub fetch_input {
  my $self = shift;

  $self->load_all_params();
}

sub run {
  my $self = shift;

  my $mysql_script = $self->base . "/scripts/mysql_functions.R";
  my $slrsim_script = $self->base . "/projects/slrsim/slrsim.functions.R";
  my $dbname = $self->dbc->dbname;
  
  my $results_table = "slrsim_results";
  my $merged_table = "merged";
  my $slrsim_label = $self->param('slrsim_label');

my $rcmd = qq^
source("${mysql_script}")
source("${slrsim_script}")
dbname="${dbname}"
con <- connect(dbname)

new.col.names <- list(
  c('slrsim_label', 'label'),
  c('slrsim_tree_file', 'tree'),
  c('slrsim_analysis_name', 'analysis'),
  c('slrsim_tree_mean_path', 'tree_length'),
  c('phylosim_insertrate', 'ins_rate'),
  c('alignment_name', 'aligner'),
  c('filtering_name', 'filter')
)
rename.cols <- function(data, old.new.list) {
  all.names <- names(data)
  for (old.new in old.new.list) {
    old <- old.new[1]
    new <- old.new[2]
    all.names[all.names == old] <- new
  }
  names(data) <- all.names
  return(data)
}

### Get genes.
sql <- sprintf("select * from genes where slrsim_label='%s'", "${slrsim_label}")
genes <- dbGetQuery(con, sql)
genes <- rename.cols(genes, new.col.names)
genes[, 'label'] <- as.factor(genes[, 'label'])
print(head(genes))

### Get sites.
nodes.s <- paste(unique(genes[, 'node_id']), collapse=',')
reps.s <- paste(unique(genes[, 'slrsim_rep']), collapse=',')
sql <- sprintf("select * FROM sites WHERE node_id IN (%s) AND slrsim_rep IN (%s)", nodes.s, reps.s)
sites <- dbGetQuery(con, sql)
if (max(sites[, 'lrt_stat'], na.rm=T) > 1.1) {
  sites[, 'lrt_stat'] <- sites[, 'lrt_stat'] * sign(sites[, 'aln_dnds'] - .9999)
}
print(head(sites))
print(paste(genes[1, 'label'], nrow(sites)))

merged <- merge(genes, sites, by=c('node_id', 'slrsim_rep'))
res.df <- paper.table(merged)

if(dbExistsTable(con, "${results_table}")) {
  dbWriteTable(con, "${results_table}", res.df, append = T)
} else {
  dbWriteTable(con, "${results_table}", res.df)
}

if(dbExistsTable(con, "${merged_table}")) {
  dbWriteTable(con, "${merged_table}", merged, append = T)
} else {
  dbWriteTable(con, "${merged_table}", merged)
}
^;
  Bio::Greg::EslrUtils->run_r( $rcmd );

}

1;
