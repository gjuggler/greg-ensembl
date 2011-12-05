#!/usr/bin/env perl

use warnings;
use strict;

use Bio::Greg::EslrUtils;
use File::Copy::Recursive qw(dircopy pathempty);

my $version = $ARGV[0];
die("Need to give an Ensembl API version number to use!") unless (defined $version);

### Step 1: make sure the ensembl and compara directories for the desired version exist.
my $base_folder = Bio::Greg::EslrUtils->baseDirectory;
my $target_ens_folder = "${base_folder}/ensembl-${version}/";
die("Target Ensembl folder ${target_ens_folder} doesn't exist!") unless (-d $target_ens_folder);
my $target_compara_folder = "${base_folder}/ensembl-compara-${version}/";
die("Target Compara folder ${target_compara_folder} doesn't exist!") unless (-d $target_compara_folder);
my $target_var_folder = "${base_folder}/ensembl-variation-${version}/";
die("Target Variation folder ${target_var_folder} doesn't exist!") unless (-d $target_compara_folder);

### Step 2: remove the existing Ensembl and Compara directories.
my $existing_ens_folder = "${base_folder}/ensembl/";
my $existing_compara_folder = "${base_folder}/ensembl-compara/";
my $existing_var_folder = "${base_folder}/ensembl-variation/";
print "  emptying existing directories...\n";
pathempty($existing_ens_folder) or die $!;
pathempty($existing_compara_folder) or die $!;
pathempty($existing_var_folder) or die $!;

print "  copying ensembl...\n";
dircopy($target_ens_folder, $existing_ens_folder) or die $!;
print "  copying ensembl-compara...\n";
dircopy($target_compara_folder, $existing_compara_folder) or die $!;
print "  copying ensembl-variation...\n";
dircopy($target_var_folder, $existing_var_folder) or die $!;

print "Now set up to use Ensembl v${version}!\n";
