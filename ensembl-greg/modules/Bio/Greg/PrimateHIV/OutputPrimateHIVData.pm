package Bio::Greg::PrimateHIV::OutputPrimateHIVData;

use strict;
use Time::HiRes qw(sleep);

use POSIX qw(strftime mktime);
use Cwd;
use File::Path;

use Bio::EnsEMBL::Hive::Process;

use Bio::Greg::EslrUtils;

use base ('Bio::Greg::Hive::Process');

sub run {
  my $self = shift;

  my $base = '/nfs/users/nfs_g/gj1/scratch/primate_hiv/output';
  $self->get_output_folder($base);

  $self->export_windows;
}

sub export_windows {
  my $self = shift;

  my $genes_file = $self->get_output_folder . "/gene_windows.Rdata";
  my $genes_csv  = $self->get_output_folder . "/gene_windows.csv";

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
