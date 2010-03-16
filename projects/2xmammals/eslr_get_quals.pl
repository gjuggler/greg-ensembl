#!/usr/bin/env perl
use warnings;
use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;
use File::Path;
use File::Basename;
use Bio::Greg::EslrUtils;

Bio::EnsEMBL::Registry->no_version_check(1);

my $url = 'mysql://ensadmin:ensembl@ens-research/gj1_2x_57';
GetOptions('url=s' => \$url);

my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-url => $url);
my $dbc = $dba->dbc;
my $dbh = $dbc->db_handle;
my $pta = $dba->get_ProteinTreeAdaptor();
my $mba = $dba->get_MemberAdaptor();
my $gba = $dba->get_GenomeDBAdaptor();

my $qual_base = "/nfs/users/nfs_g/gj1/scratch/2x_quality/assemblies/";

my $broad = "";
my $ucsc = "";

my $seq_map = {
  "Myotis lucifugus"                    => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/microbat/myoLuc1/assembly.bases.gz"
  };

my $qual_map = {
  "Aedes aegypti"                       => "",
  "Anolis carolinensis"                 => "http://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Anolis_carolinensis/bigZips/anoCar1.qual.qv.gz",
  "Anopheles gambiae"                   => "",
  "Bos taurus"                          => "http://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Bos_taurus/bigZips/bosTau4.quals.fa.gz",
  "Caenorhabditis elegans"              => "",
  "Canis familiaris"                    => "http://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Canis_familiaris/bigZips/canFam2.quals.fa.gz",
  "Cavia porcellus"                     => "http://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Cavia_porcellus/bigZips/cavPor3.quals.fa.gz",
  "Choloepus hoffmanni"                 => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/sloth/ChoHof1.0/assembly.quals.gz",
  "Ciona intestinalis"                  => "",
  "Danio rerio"                         => "",
  "Ciona savignyi"                      => "",
  "Danio rerio"                         => "",
  "Dasypus novemcinctus"                => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/armadillo/dasNov2/assembly.quals.gz",
  "Dipodomys ordii"                     => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/kangarooRat/Dipord1.0/assembly.quals.gz",
  "Drosophila melanogaster"             => "",
  "Echinops telfairi"                   => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/tenrec/echTel1/assembly.quals.gz",
  "Equus caballus"                      => "http://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Equus_caballus/bigZips/equCab2.quals.fa.gz",
  "Erinaceus europaeus"                 => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/hedgehog/eriEur1/assembly.quals.gz",
  "Felis catus"                         => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/cat/felCat3/assembly.quals.gz",
  "Gallus gallus"                       => "http://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Gallus_gallus/bigZips/galGal3.quals.fa.gz",
  "Gasterosteus aculeatus"              => "",
  "Gorilla gorilla"                     => "",
  "Homo sapiens"                        => "",
  "Loxodonta africana"                  => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/elephant/loxAfr2/assembly.quals.gz",
  "Macaca mulatta"                      => "",
  "Microcebus murinus"                  => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/mouseLemur/MicMur1.0/assembly.quals.gz",
  "Monodelphis domestica"               => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/monodelphis/monDom5/Monodelphis5.0.agp.chromosome.qual.gz",
  "Mus musculus"                        => "",
  "Myotis lucifugus"                    => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/microbat/myoLuc1/assembly.quals.gz",
  "Ochotona princeps"                   => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/pika/OchPri2.0/assembly.quals.gz",
  "Ornithorhynchus anatinus"            => "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Ornithorhynchus_anatinus/bigZips/ornAna1.quals.fa.gz",
  "Oryctolagus cuniculus"               => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/rabbit/oryCun1/assembly.quals.gz",
  "Oryzias latipes"                     => "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Oryzias_latipes/bigZips/oryLat2.quals.fa.gz",
  "Otolemur garnettii"                  => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/bushbaby/otoGar1/assembly.quals.gz",
  "Pan troglodytes"                     => "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Pan_troglodytes/bigZips/panTro2.quals.fa.gz",
  "Pongo pygmaeus"                      => "",
  "Procavia capensis"                   => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/rockHyrax/Procap1.0/assembly.quals.gz",
  "Pteropus vampyrus"                   => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/megabat/Ptevap1.0/assembly.quals.gz",
  "Rattus norvegicus"                   => "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Rattus_norvegicus/bigZips/rn4.quals.fa.gz",
  "Sorex araneus"                       => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/commonShrew/sorAra1/assembly.quals.gz",
  "Spermophilus tridecemlineatus"       => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/squirrel/speTri1/assembly.quals.gz",
  "Taeniopygia guttata"                 => "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Taeniopygia_guttata/bigZips/taeGut1.quals.fa.gz",
  "Takifugu rubripes"                   => "",
  "Tarsius syrichta"                    => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/tarsier/Tarsyr1.0/assembly.quals.gz",
  "Tetraodon nigroviridis"              => "",
  "Tupaia belangeri"                    => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/treeShrew/tupBel1/assembly.quals.gz",
  "Tursiops truncatus"                  => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/bottlenosedDolphin/Turtru1.0/assembly.quals.gz",
  "Vicugna pacos"                       => "ftp://ftp.broad.mit.edu/pub/assemblies/mammals/alpaca/VicPac1.0/assembly.quals.gz",
  "Xenopus tropicalis"                  => "",
};

my @gdbs = @{$gba->fetch_all()};
@gdbs = sort {$a->name cmp $b->name} @gdbs;
foreach my $gdb (@gdbs) {
  print $gdb->name."  ".$gdb->assembly()."\n";
  my $assembly = $gdb->assembly;
  my $species_name = $gdb->name;
  my $no_spaces = $species_name;
  $no_spaces =~ s/ /_/g;

  my $qual_file = $qual_base . "$no_spaces.quals.fa.gz";
  my $qual_uncompressed = $qual_base . "$no_spaces.quals.fa";
  if (!-e $qual_file && !-e $qual_uncompressed) {
    my $url = $qual_map->{$species_name};
    if (defined $url && $url ne "") {
      my $rc = system("wget --no-clobber -O $qual_file $url");
      die if ($rc);
    }
  }

  my $seq_file = $qual_base . "$no_spaces.bases.fa.gz";
  my $seq_uncompressed = $qual_base . "$no_spaces.bases.fa";
  if (!-e $seq_file && !-e $seq_uncompressed) {
    my $url = $seq_map->{$species_name};
    if (defined $url && $url ne '') {
      my $rc = system("wget --no-clobber -O $seq_file $url");
      die if ($rc);
    }
  }

}
