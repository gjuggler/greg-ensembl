package Bio::Greg::PrimateHIV::OutputPrimateHIVData;

use strict;
use Time::HiRes qw(sleep);

use POSIX qw(strftime mktime);
use Cwd;
use File::Path;
use File::Spec;

use Bio::EnsEMBL::Hive::Process;

use Bio::Greg::EslrUtils;

use base ('Bio::Greg::Hive::Process');

sub run {
  my $self = shift;

  $self->export_windows;
  $self->plot_manual_genes;

  # Tar and zip up all the data.
  my $dbname = $self->dbc->dbname;
  my $folder = $self->get_output_folder;
  my @dirs = File::Spec->splitdir( $folder );
  my $last_dir = $dirs[scalar(@dirs)-1];
  my $out_file = sprintf("${folder}/%s_%s.tar",$dbname,$last_dir);
  print "$out_file\n";
  if (!-e $out_file) {
    chdir($folder);
    `tar -cf $out_file *`;
    `gzip $out_file`;
  }
}

sub plot_manual_genes {
  my $self = shift;

  my $output_folder = $self->get_output_folder;
  my $sitewise_script = Bio::Greg::EslrUtils->baseDirectory . "/scripts/collect_sitewise.R";
  my $dbname = $self->hive_dba->dbc->dbname;

  my $pdf_folder = "${output_folder}/pdfs";
  mkpath($pdf_folder);

  my $cmd = qq^
setwd("${output_folder}/data")

library(phylosim)
source("~/src/greg-ensembl/projects/phylosim/Plots.R")

genes <- read.table(file="~/src/greg-ensembl/projects/primate_hiv/50_manual_genes.txt", stringsAsFactors=F)
genes <- genes[, 'V1']

load("${output_folder}/windows_all.Rdata")
genes.rows <- subset(data, gene_name %in% genes)
genes.rows <- genes.rows[!duplicated(genes.rows[, 'gene_name']), ]

for (i in 1:nrow(genes.rows)) {
  row <- genes.rows[i,]
  print(row[, 'gene_name'])
  #print(row)
 
  f.name = paste("${pdf_folder}", '/', row[,'gene_name'], '.pdf', sep='')
  pdf(file=f.name)

  sim <- PhyloSim()
  readAlignment(sim, row[, 'aln_file'])
  obj <- plotAlignment(sim, aln.do.plot=F, aln.plot.chars=T, aln.char.text.size=2, axis.text.size=3)
  p <- obj[['grob']]
  p <- p + opts(title=row[,'gene_name'])
  print(p)

  sim <- PhyloSim()
  readAlignment(sim, row[, 'pep_aln_file'])
  obj <- plotAlignment(sim, aln.do.plot=F, aln.plot.chars=T, aln.char.text.size=2, axis.text.size=3)
  p <- obj[['grob']]
  p <- p + opts(title=row[,'gene_name'])
  print(p)

  dev.off()

}
print(nrow(genes.rows))

^;
  Bio::Greg::EslrUtils->run_r( $cmd);
}

sub export_windows {
  my $self = shift;

  my $output_folder = $self->get_output_folder;
  my $sitewise_script = Bio::Greg::EslrUtils->baseDirectory . "/scripts/collect_sitewise.R";
  my $dbname = $self->hive_dba->dbc->dbname;

  my $cmd = qq^
output_folder = "${output_folder}";
dbname = "$dbname"

source("$sitewise_script",echo=T);
all.data <- data.frame()
table.name <- 'stats_windows'

cols <- 'parameter_set_id,aln_type,gene_name,n_leaves,n_pos_sites,n_sites,parameter_set_shortname,peptide_window_start,peptide_window_end,peptide_window_width,pval,stable_id_transcript'
query <- sprintf("SELECT * FROM %s",table.name)
cur.df <- get.vector(con,query,columns='all')
all.data <- rbind(all.data,cur.df)

print(nrow(all.data))

psets <- unique(all.data[,'parameter_set_shortname'])
aln.types <- unique(all.data[,'aln_type'])
print(psets)
print(aln.types)

if (length(psets) > 1 || length(aln.types) > 1) {
  for (pset in psets) {
    print(pset)
    for (aln.type in aln.types) {
      print(aln.type)
      data <- subset(all.data,parameter_set_shortname==pset & aln_type==aln.type)
      if (nrow(data)==0) {next}
      file.name <- paste(output_folder,"/windows_",pset,"_",aln.type,sep="")
      print(file.name)
      save(data,file=paste(file.name,".Rdata",sep=""))
      write.csv(data,file=paste(file.name,".csv",sep=""),row.names=F)
    }
  }
}

data <- all.data
save(data,file=paste(output_folder,"/windows_all.Rdata",sep=""))
^;
  print "$cmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r( $cmd, $params );
}

1;
