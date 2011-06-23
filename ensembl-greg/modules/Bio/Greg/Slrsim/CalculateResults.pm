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
  my $adj_results_table = "slrsim_results_adj";
  my $merged_table = "merged";

  my $slrsim_label = $self->param('slrsim_label');
  my $slrsim_tree_mean_path = $self->param('slrsim_tree_mean_path');
  my $phylosim_insertrate = $self->param('phylosim_insertrate');

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

do.meta <- ("${slrsim_label}" == 'NO LABEL')

if (!do.meta) {
  sql <- sprintf("select * from genes where slrsim_label='%s'", "${slrsim_label}")
} else {
  print("  loading meta...")
  sql <- sprintf("select * from genes where slrsim_tree_mean_path between %.3f and %.3f AND phylosim_insertrate between %.3f and %.3f",
    ${slrsim_tree_mean_path}-0.01, ${slrsim_tree_mean_path}+0.01,
    ${phylosim_insertrate}-0.01, ${phylosim_insertrate}+0.01);
  print(sql)
}

genes <- dbGetQuery(con, sql)
genes <- rename.cols(genes, new.col.names)
genes[, 'label'] <- as.factor(genes[, 'label'])
print(head(genes))

### Get sites.
nodes.s <- paste(unique(genes[, 'node_id']), collapse=',')
reps.s <- paste(unique(genes[, 'slrsim_rep']), collapse=',')
sql <- sprintf("select * FROM sites WHERE node_id IN (%s) AND slrsim_rep IN (%s)", nodes.s, reps.s)
sites <- dbGetQuery(con, sql)
print(head(sites))
print(paste(genes[1, 'label'], nrow(sites)))

merged <- merge(genes, sites, by=c('node_id', 'slrsim_rep'))

if (do.meta) {
  print("  Collecting meta-results...")
  merged <- meta.sub(merged)
  merged[, 'label'] <- paste('meta', "${slrsim_tree_mean_path}", "${phylosim_insertrate}", sep='_')
  merged[, 'filter'] <- 'meta'
}

res.df <- paper.table(merged)
if(dbExistsTable(con, "${results_table}")) {
  dbUpdateVars(con, "${results_table}", res.df, 'label')
} else {
  dbWriteTable(con, "${results_table}", res.df, row.names=F)
  dbSendQuery(con, "ALTER TABLE ${results_table} ADD UNIQUE (label(64));")
}

adj.res.df <- adj.paper.table(merged)
row.names(adj.res.df) <- adj.res.df[, 'label']
if(dbExistsTable(con, "${adj_results_table}")) {
  dbUpdateVars(con, "${adj_results_table}", adj.res.df, 'label')
} else {
  dbWriteTable(con, "${adj_results_table}", adj.res.df, row.names=F)
  dbSendQuery(con, "ALTER TABLE ${adj_results_table} ADD UNIQUE (label(64));")
}

#if (!do.meta) {
#  if(dbExistsTable(con, "${merged_table}")) {
#    dbWriteTable(con, "${merged_table}", merged, append = T, row.names=F)
#  } else {
#    dbWriteTable(con, "${merged_table}", merged, row.names=F)
#    dbSendQuery(con, "ALTER TABLE ${merged_table} ADD UNIQUE (label(64),slrsim_rep,aln_position);")
#  }
#}
^;
  Bio::Greg::EslrUtils->run_r( $rcmd );

}

1;
