#!/usr/bin/env perl
print "hey!\n";
use warnings;
use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::Greg::Hive::ComparaHiveLoaderUtils;

#my ($url) = undef;
#GetOptions('url=s' => \$url);

my $url = 'mysql://ensadmin:ensembl@ens-research:3306/gj1_dnds_58';

my $clean = 1;

my $h = new Bio::Greg::Hive::ComparaHiveLoaderUtils;
$h->init($url);

# Clean up our mess.
if ($clean) {
  $h->clean_compara_analysis_tables;
  $h->clean_hive_tables;
}

# Create analyses.
node_sets();
tree_stats();
split_subtrees();
paml_omegas();
slr_omegas();
collect_stats();

# Connect the dots.
$h->connect_analysis("NodeSets","TreeStats",1);

$h->connect_analysis("SplitBySubtrees","PamlOmegas",1);
$h->connect_analysis("PamlOmegas","SlrOmegas",1);

$h->wait_for("PamlOmegas",["SplitBySubtrees","NodeSets","TreeStats"]);
#$h->wait_for("CollectStats",["TreeStats","NodeSets","PamlOmegas","SlrOmegas"]);

add_all_nodes();

sub add_all_nodes {
  my $cmd = "SELECT node_id FROM protein_tree_node WHERE parent_id=1;";
  my @nodes = _select_node_ids($cmd);
  my $logic_name = "NodeSets";
  _add_nodes_to_analysis($logic_name,\@nodes);

  my $logic_name = "SplitBySubtrees";
  $h->add_inputs_to_analysis($logic_name, [{}]);
}

sub node_sets {
  my $logic_name = "NodeSets";
  my $module = "Bio::Greg::Hive::NodeSets";
  my $params = {
    flow_node_set => 'Primates'
  };
  my $id = $h->create_analysis($logic_name,$module,$params,150,1);
}

sub tree_stats {
  my $logic_name = "TreeStats";
  my $module = "Bio::Greg::Hive::CollectTreeStats";
  my $params = {
  };
  my $id = $h->create_analysis($logic_name,$module,$params,150,1);
}

sub split_subtrees {
  my $logic_name = "SplitSubtrees";
  my $module = "Bio::Greg::Hive::SplitBySubtrees";
  my $params = {
    seed_species => '9606',
    out_to => 'Mammalia'
  };
  $h->create_analysis($logic_name,$module,$params,100,1);
}

sub paml_omegas {
  my $logic_name = "PamlOmegas";
  my $module = "Bio::Greg::Hive::GenomewideOmegas";
  my $params = {
    method => 'paml',
    foreground_species => '9606,9593'
  };
  $h->create_analysis($logic_name,$module,$params,100,1);
}

sub slr_omegas {
  my $logic_name = "SlrOmegas";
  my $module = "Bio::Greg::Hive::GenomewideOmegas";
  my $params = {
    method => 'slr'
  };
  $h->create_analysis($logic_name,$module,$params,100,1);
}

sub collect_stats {
  my $logic_name = "CollectStats";
  my $module = "Bio::Greg::Gorilla::CollectOmegaStats";
  my $params = {
  };
  $h->create_analysis($logic_name,$module,$params,50,1);
}

########*********########
#-------~~~~~~~~~-------#
########*********########

# Subroutines to return a list of taxon IDs with specific features.
sub clade_taxon_ids {
  my $clade = shift || 1;

  my $dba = $h->dba;
  my @genomes = Bio::EnsEMBL::Compara::ComparaUtils->get_genomes_within_clade($dba,$clade);
  my @taxon_ids = map {$_->taxon_id} @genomes;
  return @taxon_ids;
}
sub coverage_taxon_ids {
  my $coverage = shift;
  
  my $dba = $h->dba;
  my @output;
  my @all_gdb = Bio::EnsEMBL::Compara::ComparaUtils->get_genomes_within_clade($dba,1);
  foreach my $gdb (@all_gdb) {
    # This is finicky: we need to call the "db_adaptor" method to get the Bio::EnsEMBL::DBSQL::DBAdaptor object, and then the meta container.
    my $meta = $gdb->db_adaptor->get_MetaContainer;
    my $str = @{$meta->list_value_by_key('assembly.coverage_depth')}[0];
    push @output, $gdb->taxon_id if ($str eq $coverage);
  }  return @output;
}
sub subtract {
  my $list_a = shift;
  my @remove_us = @_;
  my $hash;
  map {$hash->{$_}=1} @$list_a;
  foreach my $list_b (@remove_us) {
    map {delete $hash->{$_}} @$list_b;
}
  return keys %$hash;
}
 
sub _select_node_ids {
  my $cmd = shift;
  if ( !defined $cmd ) {
    $cmd = "SELECT node_id FROM protein_tree_node WHERE parent_id=1 AND root_id=1";
  }
  my $dbc = $h->dbc;
  my $sth = $dbc->prepare($cmd);
  $sth->execute();

  my $array_ref = $sth->fetchall_arrayref( [0] );
  my @node_ids = @{$array_ref};
  @node_ids =
    map { @{$_}[0] } @node_ids;    # Some weird mappings to unpack the numbers from the arrayrefs.
  $sth->finish;
  return @node_ids;
}

sub _add_nodes_to_analysis {
  my $logic_name = shift;
  my $node_id_arrayref = shift;

  my @node_ids = @{$node_id_arrayref};
  foreach my $node_id (@node_ids) {
    my $input_id = {
      node_id => $node_id
    };
    $h->add_job_to_analysis($logic_name,$input_id);
  }
}
