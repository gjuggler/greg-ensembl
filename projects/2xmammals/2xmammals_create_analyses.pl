#!/usr/bin/env perl

use warnings;
use strict;
use DBI;
use Getopt::Long;
use File::Path;
use File::Basename;
use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::Greg::Hive::ComparaHiveLoaderUtils;

my ($url) = 'mysql://ensadmin:ensembl@ens-research/gj1_2x';
GetOptions('url=s' => \$url);

my $h = new Bio::Greg::Hive::ComparaHiveLoaderUtils;
$h->init($url);

my $clean = 1;

# Clean up our mess.
if ($clean) {
  $h->clean_compara_analysis_tables;
  $h->clean_hive_tables;
}

# Define parameters (species sets, filtering options, etc).
parameter_sets();

# Create analyses.
node_sets();
align();
sequence_quality();
split_by_parameter_set();
gene_omegas();
sitewise_omegas();
mapping();
collect_stats();
output_data();

# Connect the dots.
$h->connect_analysis("NodeSets","Align");
$h->connect_analysis("Align","SequenceQuality");
$h->connect_analysis("SequenceQuality","SplitByParameterSet");
$h->connect_analysis("SplitByParameterSet","GeneOmegas");
$h->connect_analysis("GeneOmegas","SitewiseOmegas");
$h->connect_analysis("SitewiseOmegas","CollectStats");

$h->connect_analysis("Align","Mapping");
$h->wait_for("CollectStats",["Mapping"]);
$h->wait_for("OutputTabularData",["CollectStats","SitewiseOmegas"]);

sub node_sets {
  my $logic_name = "NodeSets";
  my $module = "Bio::Greg::Hive::NodeSets";
  my $params = {
    flow_node_set => 'MammaleerPlusOutgroup'
  };
  my $analysis_id = $h->create_analysis($logic_name,$module,$params,50,1);

  # Add all root nodes to this analysis.
  $params = {};
  my $cmd = "SELECT node_id FROM protein_tree_node WHERE parent_id=1;";
  my @nodes = _select_node_ids($cmd);
  _add_nodes_to_analysis($analysis_id,$params,\@nodes);
}

sub parameter_sets {
  my $params;
  my $base_params={};

  # 2xmammals off-limits species.
  our @off = (
	     9593, # Gorilla
	     9600, # Orang
	     9483  # Marmoset
	     );
  
  # Subroutines to return a list of taxon IDs with specific features.
  sub clade_taxon_ids {
    my $clade = shift || 1;
    my @genomes = Bio::EnsEMBL::Compara::ComparaUtils->get_genomes_within_clade($h->dba,$clade);
    my @taxon_ids = map {$_->taxon_id} @genomes;
    
    return subtract(\@taxon_ids,\@off);
  }
  sub coverage_taxon_ids {
    my $coverage = shift;

    my @output;
    my @all_gdb = Bio::EnsEMBL::Compara::ComparaUtils->get_genomes_within_clade($h->dba,1);
    foreach my $gdb (@all_gdb) {
      # This is finicky: we need to call the "db_adaptor" method to get the Bio::EnsEMBL::DBSQL::DBAdaptor object, and then the meta container.
      my $meta = $gdb->db_adaptor->get_MetaContainer;
      my $str = @{$meta->list_value_by_key('assembly.coverage_depth')}[0];
      print "Coverage: $str\n";
      push @output, $gdb->taxon_id if ($str eq $coverage);
    }
    return subtract(\@output,\@off);
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

  my @all = clade_taxon_ids();
  my @mamms = clade_taxon_ids("Eutheria");
  my @primates_arr = clade_taxon_ids("Primates");
  my @glires_arr = clade_taxon_ids("Glires");

  my $hmrd = join(",",(9606, 10090, 10116, 9615));

  my $not_mammals = join(",",subtract(\@all,\@mamms));
  my $not_primates = join(",",subtract(\@mamms,\@primates_arr));
  my $not_glires = join(",",subtract(\@mamms,\@glires_arr));

  my $everything = join(",",clade_taxon_ids());
  my $mammals = join(",",clade_taxon_ids("Eutheria"));
  my $primates = join(",",clade_taxon_ids("Primates"));
  my $glires = join(",",clade_taxon_ids("Glires"));
  my $laurasiatheria = join(",",clade_taxon_ids("Laurasiatheria"));
  my $afrotheria = join(",",clade_taxon_ids("Afrotheria"));

  my $sauria = join(",",clade_taxon_ids("Sauria"));
  my $fishes = join(",",clade_taxon_ids("Clupeocephala"));
  
  # Get only hi-coverage genomes.
  my @hi_coverage_arr = coverage_taxon_ids("high");
  my @lo_coverage_arr = coverage_taxon_ids("low");

  my $hi_coverage = join(",",subtract(\@mamms,\@hi_coverage_arr));
  my $lo_coverage = join(",",subtract(\@mamms,\@lo_coverage_arr));

  $params = {
    parameter_set_name => "Mammals",
    parameter_set_shortname => 'm',
    keep_species => $mammals,
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  return;

  $params = {
    parameter_set_name => "Primates",
    parameter_set_shortname => 'p',
    keep_species => $primates
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "Glires",
    parameter_set_shortname => 'g',
    keep_species => $glires
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "Laurasiatheria",
    parameter_set_shortname => 'l',
    keep_species => $laurasiatheria
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "Sauria",
    parameter_set_shortname => 's',
    keep_species => $sauria,
    remove_species => $glires
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "Fishes",
    parameter_set_shortname => 'f',
    keep_species => $fishes
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "High Coverage Only",
    parameter_set_shortname => 'hi',
    keep_species => $hi_coverage
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "Low Coverage Only",
    parameter_set_shortname => 'lo',
    keep_species => $lo_coverage
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "No Primates",
    parameter_set_shortname => 'np',
    keep_species => $not_primates,
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "No Glires",
    parameter_set_shortname => 'ng',
    keep_species => $not_glires
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "Human/Mouse/Rat/Dog",
    parameter_set_shortname => 'hmrd',
    keep_species => $hmrd
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

}

sub align {
  my $logic_name = "Align";
  my $module = "Bio::Greg::Hive::Align";
  my $params = {
    parameter_set_id => 1,
    alignment_method => 'prank'
  };

  $h->create_analysis($logic_name,$module,$params,400,1);
}

sub split_by_parameter_set {
  my $logic_name = "SplitByParameterSet";
  my $module = "Bio::Greg::Hive::SplitByParameterSet";
  my $params = {
    flow_parameter_sets => '1'
  };
  $h->create_analysis($logic_name,$module,$params,100,1);
}

sub sequence_quality {
  my $analysis_id=103;
  my $logic_name = "SequenceQuality";
  my $module = "Bio::Greg::Hive::SequenceQualityLoader";
  my $params = {};
  
  $h->create_analysis($logic_name,$module,$params,50,1);
}

sub gene_omegas {
  my $logic_name = "GeneOmegas";
  my $module = "Bio::Greg::Hive::PhyloAnalysis";
  my $base_params = {
    sequence_quality_filtering => 1,
    analysis_action => 'hyphy_dnds'
    };
  $h->create_analysis($logic_name,$module,$base_params,500,1);
}

sub sitewise_omegas {
  my $logic_name = "SitewiseOmegas";
  my $module = "Bio::Greg::Hive::PhyloAnalysis";
  my $base_params = {
    sequence_quality_filtering => 1,
    analysis_action => 'slr'
    };
  $h->create_analysis($logic_name,$module,$base_params,500,1);
}

sub mapping {
  my $logic_name = "Mapping";
  my $module = "Bio::Greg::Hive::SitewiseMapper";
  my $params = {
  };
  $h->create_analysis($logic_name,$module,$params,50,1);
}


sub collect_stats {
  my $logic_name = "CollectStats";
  my $module = "Bio::Greg::Mammals::CollectMammalsStats";
  my $params = {
    sequence_quality_filtering => 1,
  };
  $h->create_analysis($logic_name,$module,$params,50,1);
}

sub output_data {
  my $logic_name = "OutputTabularData";
  my $module = "Bio::Greg::Mammals::OutputMammalsData";
  my $params = {
    sequence_quality_filtering => 1,
  };
  my $analysis_id = $h->create_analysis($logic_name,$module,$params,50,1);
  _add_nodes_to_analysis($analysis_id,{},[0]);  
}


########*********########
#-------~~~~~~~~~-------#
########*********########

sub _combine_hashes {
  my @hashes = @_;

  my $new_hash = {};
  foreach my $hash (@hashes) {
    foreach my $key (keys %{$hash}) {
      $new_hash->{$key} = $hash->{$key};
    }
  }
  return $new_hash;
}

our $param_set_counter;
sub _add_parameter_set {
  my $params = shift;

  my $dbc = $h->dba->dbc;

  $param_set_counter = 1 if (!$param_set_counter);
  my $parameter_set_id = $params->{'parameter_set_id'} || $param_set_counter++;
  $params->{'parameter_set_id'} = $parameter_set_id;

  if (exists $params->{'parameter_set_name'} ) {
    my $parameter_set_name = $params->{'parameter_set_name'};
    my $name_cmd = "REPLACE INTO parameter_set VALUES ('$parameter_set_id','name',\"$parameter_set_name\");";
    $dbc->do($name_cmd);
  }

  if (exists $params->{'parameter_set_shortname'} ) {
    my $parameter_set_shortname = $params->{'parameter_set_shortname'} || '';
    my $shortname_cmd = "REPLACE INTO parameter_set VALUES ('$parameter_set_id','parameter_set_shortname',\"$parameter_set_shortname\");";
    $dbc->do($shortname_cmd);
  }
  
  my $param_string = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
  my $cmd = "REPLACE INTO parameter_set VALUES ('$parameter_set_id','params',\"$param_string\");";
  $dbc->do($cmd);
}

sub _add_nodes_to_analysis {
  my $analysis_id = shift;
  my $params = shift || {};
  my $node_arrayref = shift;

  my $dbc = $h->dbc;

  my @node_ids = @{$node_arrayref};  
  my $sth = $dbc->prepare( qq{REPLACE INTO analysis_job SET
			       analysis_id=?,
			       input_id=?;}
			   );
  foreach my $node_id (@node_ids) {
    $params->{'node_id'} = $node_id;
    my $input_id = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
    $sth->execute($analysis_id,$input_id);
    my $analysis_job_id = $sth->{'mysql_insertid'};
    print "New AnalysisJob: $analysis_id  $analysis_job_id  $input_id\n";
  }
  $sth->finish;
}

sub _select_node_ids {
  my $cmd = shift;
  if (!defined $cmd) {
    $cmd = "SELECT node_id FROM protein_tree_node WHERE parent_id=1 AND root_id=1";
  }
  my $dbc = $h->dbc;
  my $sth = $dbc->prepare($cmd);
  $sth->execute();

  my $array_ref = $sth->fetchall_arrayref([0]);
  my @node_ids = @{$array_ref};
  @node_ids = map {@{$_}[0]} @node_ids;  # Some weird mappings to unpack the numbers from the arrayrefs.
  $sth->finish;
  return @node_ids;
}
