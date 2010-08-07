#!/usr/bin/env perl
use warnings;
use strict;
use DBI;
use Bio::Greg::Hive::ComparaHiveLoaderUtils;

my $url = 'mysql://ensadmin:ensembl@ens-research:3306/gj1_subs';
my $clean = 1;

my $h = new Bio::Greg::Hive::ComparaHiveLoaderUtils;
$h->init($url);

# Clean up our mess.
if ($clean) {
  $h->clean_hive_tables;
}

# Define parameters (species sets, filtering options, etc).
collect_substitution_runs();
add_all_members();

sub add_all_members {
  Bio::EnsEMBL::Registry->load_registry_from_multiple_dbs({
                                                           -host => 'ens-livemirror',
                                                           -user => 'ensro'
                                                          });
  my $compara = Bio::EnsEMBL::Registry->get_DBAdaptor('multi','compara');
  my $mba = $compara->get_MemberAdaptor;
  my @members = @{$mba->fetch_all_by_source_taxon('ENSEMBLGENE',9606)};

  foreach my $member (@members) {
    print "$member\n";
    my $input_id = {
                    member_id => $member->member_id
                   };
    $h->add_job_to_analysis("CollectSubstitutionRuns",$input_id);
  }

}

sub collect_substitution_runs {
  my $logic_name = "CollectSubstitutionRuns";
  my $module = "Bio::Greg::Hive::CollectSubstitutionRuns";
  my $params = {
                ref_species => '9606',
                other_species => '9598,9593,9600,9544'
  };
  my $id = $h->create_analysis($logic_name,$module,$params,50,1);
}
