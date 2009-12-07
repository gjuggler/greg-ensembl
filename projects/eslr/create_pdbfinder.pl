#!/usr/bin/env perl
use warnings;
use strict;

# A quick script to download and create a pdbfinder sqlite table.
# wget ftp://ftp.cmbi.ru.nl/pub/molbio/data/pdbfinder2/PDBFIND2.TXT.gz; gunzip PDBFIND2.TXT.gz;
# gunzip PDBFIND2.TXT.gz
# Make sure the PDBFIND2.TXT is in the current directory.

my $url = "mysql://ensadmin:ensembl@compara2:3306/gj1_57";
GetOptions('url=s' => \$url);
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-url => $url);
my $dbc = $dba->dbc;

create_mysql();
#create_sqlite();

open(IN,"PDBFIND2.TXT");
my $pid = open(WRITEME, "| mysql -u") or die ("Couldn't fork $!\n");
my $i=0;
my $id;
print WRITEME "BEGIN TRANSACTION;\n";
while (<IN>) {
  $i++;
#  last if ($i > 500);
  chomp $_;

  if ($_ =~ 'ID\s*:\s*(\S+)') {
    $id = $1;
    print "$id\n";
    print WRITEME "INSERT OR IGNORE INTO pdbfinder (id) VALUES ('$id');\n";
  }
  if ($_ =~ 'Sequence\s*:\s*(\S+)') {
    print WRITEME "UPDATE pdbfinder SET seq='$1' WHERE id='$id';\n";
  }
  if ($_ =~ 'Access\s*:\s*(\S+)') { 
    print WRITEME "UPDATE pdbfinder SET access='$1' WHERE id='$id';\n";
  }
  if ($_ =~ 'DSSP\s*:\s*(\S+)') {
    print WRITEME "UPDATE pdbfinder SET dssp='$1' WHERE id='$id';\n";
  }
}
close(IN);
print WRITEME "END TRANSACTION;\n";
close(WRITEME);

sub create_mysql {
  $dbc->do(qq^
CREATE TABLE IF NOT EXISTS pdbfinder (
  id      CHAR(8),
  seq     MEDIUMTEXT,
  access  MEDIUMTEXT,
  dssp    MEDIUMTEXT, 
  PRIMARY KEY (id),
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
