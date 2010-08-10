package Bio::Greg::Iso::OutputIsoData;

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

  $self->get_output_folder;

  $self->export_genes();
}

sub export_genes {
  my $self = shift;

  my $genes_file = $self->get_output_folder . "/gene.csv";
  my $folder     = $self->get_output_folder;

  my $cmd = qq^
source("../../scripts/collect_sitewise.R");
genes <- get.genes(db="gj1_iso_57")
write.csv(genes,file="${genes_file}",row.names=F)
^;
  print "$cmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r( $cmd, $params );
}

sub get_output_folder {
  my $self = shift;

  if ( defined $self->param('output_folder') ) {
    return $self->param('output_folder');
  }

  my $date_string = strftime( "%Y-%m-%d", localtime );
  my $i = 0;

  my $filename;
  do {
    $i++;
    $filename = sprintf( "/nfs/users/nfs_g/gj1/scratch/iso/output/%s/%s_%.2d",
      $date_string, $date_string, $i );
  } while ( -e $filename );

  print "Output folder: $filename\n";
  $self->param( 'output_folder', $filename );

# We'll store this output folder in the meta table, so it will be re-used if this module is run again w/ the same database.
  $self->store_meta( { output_folder => $filename } );
  mkpath( [$filename] );
}

1;
