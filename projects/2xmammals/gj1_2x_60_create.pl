#!/usr/bin/env perl

use warnings;
use strict;
use Bio::Greg::Hive::ComparaHiveLoaderUtils;

my $url = 'mysql://ensadmin:ensembl@ens-research/gj1_2x_60';
my $h = new Bio::Greg::Hive::ComparaHiveLoaderUtils;
$h->init($url);

my $clean = 1;

# Clean up our mess.
if ($clean) {
  $h->clean_compara_analysis_tables;
  $h->clean_hive_tables;
  $h->truncate_tables([qw(meta genes_g genes_c sites_g sites_c)]);
}

# Define parameters (species sets, filtering options, etc).
parameter_sets();

# Create analyses.
node_sets();
split_by_parameter_set();
sitewise_mammals();
webpages();
output_data();

# Connect the dots.
$h->connect_analysis("NodeSets","SplitByParameterSet");
$h->connect_analysis("SplitByParameterSet","SitewiseMammals");
$h->connect_analysis("SitewiseMammals","Webpages");
$h->wait_for("OutputMammalsData",["NodeSets","SplitByParameterSet","SitewiseMammals"]);

# Add some trees.
add_all_nodes();

sub add_all_nodes {
  my $cmd = "SELECT node_id FROM protein_tree_node WHERE parent_id=1;";
  my @nodes = $h->select_node_ids($cmd);

  # TEMP: Only add the first N.
  #@nodes = @nodes[60..70];

  $h->add_nodes_to_analysis('NodeSets',\@nodes);
}

### Process definitions.
###

sub node_sets {
  my $logic_name = "NodeSets";
  my $module = "Bio::Greg::Hive::NodeSets";
  my $params = {
    flow_node_set => 'Mammals'
  };
  my $analysis_id = $h->create_analysis($logic_name,$module,$params,50,1);
}


sub split_by_parameter_set {
  my $logic_name = "SplitByParameterSet";
  my $module = "Bio::Greg::Hive::SplitByParameterSet";
  my $params = {
    flow_parameter_sets => 'all'
  };
  $h->create_analysis($logic_name,$module,$params,150,1);
}

sub sitewise_mammals {
  my $logic_name = "SitewiseMammals";
  my $module = "Bio::Greg::Mammals::SitewiseMammals";
  my $base_params = {
  };
  $h->create_analysis($logic_name,$module,$base_params,500,1);
}

sub Webpages {
  my $logic_name = "Webpages";
  my $module = "Bio::Greg::Mammals::Webpages";
  my $base_params = {
  };
  $h->create_analysis($logic_name,$module,$base_params,80,1);
}

sub output_data {
  my $logic_name = "OutputMammalsData";
  my $module = "Bio::Greg::Mammals::OutputMammalsData";
  my $params = {
  };
  $h->create_analysis($logic_name,$module,$params,1,1);
  $h->add_job_to_analysis($logic_name,{});
}

### Parameter sets.
###

sub parameter_sets {
  my $params;

  # 2xmammals off-limits species.
  our @off = (
	     9593, # Gorilla
	     9600, # Orang
	     9483  # Marmoset
	     );
  
  # Subroutines to return a list of taxon IDs with specific features.
  sub clade_taxon_ids {
    my $clade = shift || 1;
    my @genomes = Bio::EnsEMBL::Compara::ComparaUtils->get_genomes_within_clade($h->compara_dba,$clade);
    my @taxon_ids = map {$_->taxon_id} @genomes;
    
    return subtract(\@taxon_ids,\@off);
  }
  sub coverage_taxon_ids {
    my $coverage = shift;

    my @output;
    my @all_gdb = Bio::EnsEMBL::Compara::ComparaUtils->get_genomes_within_clade($h->compara_dba,1);
    foreach my $gdb (@all_gdb) {
      # This is finicky: we need to call the "db_adaptor" method to get the Bio::EnsEMBL::DBSQL::DBAdaptor object, and then the meta container.
      my $meta = $gdb->db_adaptor->get_MetaContainer;
      my $str = @{$meta->list_value_by_key('assembly.coverage_depth')}[0];
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

  foreach my $aln_type ('compara','genomic_all') {
    my $aln_short = 'm';
    $aln_short = 'c' if ($aln_type eq 'compara');
    $aln_short = 'g' if ($aln_type eq 'genomic_all');

    $params = {
      parameter_set_name => "Mammals".' '.$aln_type,
      parameter_set_shortname => 'm_'.$aln_short,
      quality_threshold => 30,
      keep_species => $mammals,
      genes_table => 'genes_'.$aln_short,
      sites_table => 'sites_'.$aln_short,
      aln_type => $aln_type,
      rename_sequences_by_taxon => 0
    };
    $params->{realign_with_prank} = 1 if ($aln_type eq 'compara');
    $h->add_parameter_set($params);
  }
}


### Misc. functions.
###

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
