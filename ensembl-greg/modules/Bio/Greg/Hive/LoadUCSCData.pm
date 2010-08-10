package Bio::Greg::Hive::LoadUCSCData;

use strict;
use Time::HiRes qw(sleep);

use Cwd;
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::Process');

sub fetch_input {
  my $self = shift;

  ### DEFAULT PARAMETERS ###
  my $defaults = { collect_recomb_rate => 1 };
  ##########################

  $self->load_all_params($defaults);
}

sub run {
  my $self = shift;

  if ( $self->param('collect_recomb_rate') ) {
    collect_recomb_rate();
  }

}

sub collect_recomb_rate {
  my $self = shift;

  # Recombination rate was calculated for hg18; first we lift over to hg19.

# Download the hg18 to hg19 liftover.
#> wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
#> gunzip hg18ToHg19.over.chain.gz
#> mysql --skip-column-names --user=genome --host=genome-mysql.cse.ucsc.edu hg18 -e "select chrom,chromStart,chromEnd,concat_ws('_',chrom,chromStart,chromEnd) AS name from recombRate;" > recombRate.txt
#> /software/ensembl/compara/bin/liftOver recombRate.txt hg18ToHg19.over.chain recombRate_hg19.txt unmapped.txt

  my $cwd       = Bio::Greg::EslrUtils->baseDirectory . '/projects/2xmammals/data/ucsc/';
  my $chainFile = $cwd . 'hg18ToHg19.over.chain';
  my $outFile   = $cwd . 'recombRate_hg19.txt';

  open( IN, "$outfile" );
  while (<IN>) {
    my ( $chr, $start, $end, $name ) = split("\t");

  }
  close(IN);

}

sub write_output {

}

1;
