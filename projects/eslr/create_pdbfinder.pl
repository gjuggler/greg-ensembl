#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::Greg::EslrUtils;
use File::Path;
use File::Basename;

# A quick script to download and create a pdbfinder sqlite table.
# wget ftp://ftp.cmbi.ru.nl/pub/molbio/data/pdbfinder2/PDBFIND2.TXT.gz; gunzip PDBFIND2.TXT.gz;
# gunzip PDBFIND2.TXT.gz
# Make sure the PDBFIND2.TXT is in the current directory.

my $url = 'mysql://ensadmin:ensembl@compara2:3306/gj1_57';
GetOptions('url=s' => \$url);
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-url => $url);
my $dbc = $dba->dbc;

create_mysql();
#create_sqlite();

$dbc->do("BEGIN WORK;");
my $insert = $dbc->prepare("REPLACE INTO pdbfinder (id) VALUES (?);");
my $update = $dbc->prepare("UPDATE pdbfinder set seq=? WHERE id=?;");
my $update2 = $dbc->prepare("UPDATE pdbfinder set access=? WHERE id=?;");
my $update3 = $dbc->prepare("UPDATE pdbfinder set dssp=? WHERE id=?;");

open(IN,"PDBFIND2.TXT");
my $i=0;
my $id;
while (<IN>) {
  $i++;
#  last if ($i > 500);
  chomp $_;

  if ($_ =~ 'ID\s*:\s*(\S+)') {
    $id = $1;
    print "$id\n";
    $insert->execute($id);
  }
  if ($_ =~ 'Sequence\s*:\s*(\S+)') {
    $update->execute($1,$id);
  }
  if ($_ =~ 'Access\s*:\s*(\S+)') { 
    $update2->execute($1,$id);
  }
  if ($_ =~ 'DSSP\s*:\s*(\S+)') {
    $update3->execute($1,$id);
  }
}
close(IN);

$insert->finish;
$update->finish;
$update2->finish;
$update3->finish;

$dbc->do("END WORK;");

sub create_mysql {
  $dbc->do(qq^
CREATE TABLE IF NOT EXISTS pdbfinder (
  id      CHAR(8),
  seq     MEDIUMTEXT,
  access  MEDIUMTEXT,
  dssp    MEDIUMTEXT, 
  PRIMARY KEY (id)
);
^);
}

sub create_sqlite {
my $cmd = qq^
CREATE TABLE IF NOT EXISTS pdbfinder (
  id      TEXT UNIQUE,
  seq     TEXT,
  access  TEXT,
  dssp    TEXT
);
^;
`sqlite3 pdbfinder.db '$cmd'`;
}
