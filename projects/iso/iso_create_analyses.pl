#!/usr/bin/env perl

use warnings;
use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::Greg::EslrUtils;
use File::Path;
use File::Basename;

my ($url) = 'mysql://ensadmin:ensembl@ens-research/gj1_iso_57';
GetOptions('url=s' => \$url);
my $clean = 1;
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-url => $url);
my $hive_dba = Bio::EnsEMBL::Hive::DBSQL::DBAdaptor->new(-url => $url);
my $dbc = $dba->dbc;
my $analysis_hash; # Hash to store mapping between analysis names and ID numbers.

# Clean up our mess.
clean_tables();

# Define parameters (species sets, filtering options, etc).
parameter_sets();

# Create analyses.
node_sets();
split_by_parameter_set();
dump_alignments();
collect_stats();
output_data();

# Connect the dots.
connect_analysis("NodeSets","SplitByParameterSet");
connect_analysis("SplitByParameterSet","DumpAlignments");
connect_analysis("DumpAlignments","CollectStats");

wait_for("OutputTabularData",["CollectStats"]);

sub node_sets {
  my $logic_name = "NodeSets";
  my $module = "Bio::Greg::Hive::NodeSets";
  my $params = {
    flow_node_set => 'MammalPlusOutgroup'
  };
  my $analysis_id = _create_analysis($logic_name,$module,$params,100,1);

  # Add all root nodes to this analysis.
  $params = {};
 
  # Add a subset for testing.
  my $lo = 0;
  foreach my $hi (0.3, 0.5, 1, 2, 4, 8, 16) {
    my $cmd = "SELECT human_protein FROM gj1_gor_57.stats_genes WHERE human_gene IS NOT NULL AND parameter_set_name='Mammals' AND tree_length between $lo and $hi limit 5;";
    my @nodes = _select_node_ids($cmd);
    print "$lo $hi: @nodes\n";
    _add_nodes_to_analysis($analysis_id,\@nodes,$params);
    $lo = $hi;
  }

}

sub split_by_parameter_set {
  my $logic_name = "SplitByParameterSet";
  my $module = "Bio::Greg::Hive::SplitByParameterSet";
  my $params = {
    flow_parameter_sets => 'all'
  };
  my $analysis_id = _create_analysis($logic_name,$module,$params,100,1);
}

sub dump_alignments {
  my $logic_name = "DumpAlignments";
  my $module = "Bio::Greg::Hive::DumpAlignment";
  my $params = {
    output_folder => '/nfs/users/nfs_g/gj1/scratch/iso/alignments',
    hash_subfolders => 1
  };
  _create_analysis($logic_name,$module,$params,20,1);
}

sub collect_stats {
  my $logic_name = "CollectStats";
  my $module = "Bio::Greg::Iso::CollectIsoStats";
  my $params = {
  };
  _create_analysis($logic_name,$module,$params,50,1);
}

sub output_data {
  my $logic_name = "OutputTabularData";
  my $module = "Bio::Greg::Iso::OutputIsoData";
  my $params = {
  };
  my $analysis_id = _create_analysis($logic_name,$module,$params,1,1);
  _add_nodes_to_analysis($analysis_id,{},[0]);  
}

sub parameter_sets {
  my $params;
  my $base_params={};

  # Subroutines to return a list of taxon IDs with specific features.
  sub clade_taxon_ids {
    my $clade = shift || 1;
    my @genomes = Bio::EnsEMBL::Compara::ComparaUtils->get_genomes_within_clade($dba,$clade);
    my @taxon_ids = map {$_->taxon_id} @genomes;
    return @taxon_ids;
  }
  sub coverage_taxon_ids {
    my $coverage = shift;

    my @output;
    my @all_gdb = Bio::EnsEMBL::Compara::ComparaUtils->get_genomes_within_clade($dba,1);
    foreach my $gdb (@all_gdb) {
      # This is finicky: we need to call the "db_adaptor" method to get the Bio::EnsEMBL::DBSQL::DBAdaptor object, and then the meta container.
      my $meta = $gdb->db_adaptor->get_MetaContainer;
      my $str = @{$meta->list_value_by_key('assembly.coverage_depth')}[0];
      push @output, $gdb->taxon_id if ($str eq $coverage);
    }
    return @output;
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

  my $everything = join(",",clade_taxon_ids());

  my $eutheria = join(",",clade_taxon_ids("Eutheria"));
  my $mammalia = join(",",clade_taxon_ids("Mammalia"));
  my $amniota = join(",",clade_taxon_ids("Amniota"));
  my $chordata = join(",",clade_taxon_ids("Chordata"));

  my $primates = join(",",clade_taxon_ids("Primates"));
  my $glires = join(",",clade_taxon_ids("Glires"));
  my $laurasiatheria = join(",",clade_taxon_ids("Laurasiatheria"));
  my $afrotheria = join(",",clade_taxon_ids("Afrotheria"));
  my $sauria = join(",",clade_taxon_ids("Sauria"));
  my $fishes = join(",",clade_taxon_ids("Clupeocephala"));

 
  $params = {
    parameter_set_name => "Eutheria",
    parameter_set_shortname => 'e',
    keep_species => $eutheria
  };
  _add_parameter_set($params);

  return;

  $params = {
    parameter_set_name => "Mammalia",
    parameter_set_shortname => 'm',
    keep_species => $mammalia
  };
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "Amniota",
    parameter_set_shortname => 'a',
    keep_species => $amniota
  };
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "Chordata",
    parameter_set_shortname => 'c',
    keep_species => $chordata
  };
  _add_parameter_set($params);

  
  
  $params = {
    parameter_set_name => "Primates",
    parameter_set_shortname => 'p',
    keep_species => $primates,
  };
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "Glires",
    parameter_set_shortname => 'g',
    keep_species => $glires
  };
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "Laurasiatheria",
    parameter_set_shortname => 'l',
    keep_species => $laurasiatheria
  };
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "Afrotheria",
    parameter_set_shortname => 'afr',
    keep_species => $afrotheria
  };
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "Sauria",
    parameter_set_shortname => 's',
    keep_species => $sauria
  };
  _add_parameter_set($params);

  $params = {
    parameter_set_name => "Fishes",
    parameter_set_shortname => 'f',
    keep_species => $fishes
  };
  _add_parameter_set($params);


}

our $param_set_counter;
sub _add_parameter_set {
  my $params = shift;

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


sub clean_tables {
  if ($clean) {    
    my @truncate_tables = qw^
      analysis analysis_job analysis_stats dataflow_rule hive
      parameter_set
      node_set_member node_set
      sitewise_omega sitewise_tag sitewise_genome
      go_terms
      counts_sites counts_genes
      stats_sites stats_genes
      ^;
    map {
      print "$_\n";
      eval {$dba->dbc->do("truncate table $_");}} @truncate_tables;
  }
}


########*********########
#-------~~~~~~~~~-------#
########*********########

our $analysis_counter = 0;
sub _create_analysis {
  my $logic_name = shift;
  my $module = shift;
  my $params = shift;
  my $hive_capacity = shift || 500;
  my $batch_size = shift || 1;

  my $analysis_id = ++$analysis_counter;
  
  $analysis_hash->{$logic_name} = $analysis_id;
  my $param_string = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
  my $cmd = qq{REPLACE INTO analysis SET
		 created=now(),
		 analysis_id=$analysis_id,
		 logic_name="$logic_name",
		 module="$module",
		 parameters="$param_string"
		 ;};
  $dbc->do($cmd);
  $cmd = qq{REPLACE INTO analysis_stats SET
	      analysis_id=$analysis_id,
	      hive_capacity=$hive_capacity,
	      batch_size=$batch_size
	      ;};
  $dbc->do($cmd);
  return $analysis_id;
}

sub connect_analysis {
  my $from_name = shift;
  my $to_name = shift;
  my $branch_code = shift;
  $branch_code = 1 unless (defined $branch_code);

  my $dataflow_rule_adaptor = $hive_dba->get_DataflowRuleAdaptor;
  my $analysis_adaptor = $hive_dba->get_AnalysisAdaptor;

  my $from_analysis = $analysis_adaptor->fetch_by_logic_name($from_name);
  my $to_analysis = $analysis_adaptor->fetch_by_logic_name($to_name);
  
  if($from_analysis and $to_analysis) {
    $dataflow_rule_adaptor->create_rule( $from_analysis, $to_analysis, $branch_code);
    warn "Created DataFlow rule: [$branch_code] $from_name -> $to_name\n";
  } else {
    die "Could not fetch analyses $from_analysis -> $to_analysis to create a dataflow rule";
  }
}

sub wait_for {
  my $waiting_name = shift;
  my $wait_for_list = shift;
  
  my $ctrl_rule_adaptor = $hive_dba->get_AnalysisCtrlRuleAdaptor;
  my $analysis_adaptor = $hive_dba->get_AnalysisAdaptor;

  my $waiting_analysis = $analysis_adaptor->fetch_by_logic_name($waiting_name);

  foreach my $wait_for_name (@$wait_for_list) {
    my $wait_for_analysis = $analysis_adaptor->fetch_by_logic_name($wait_for_name);

    if($waiting_analysis and $wait_for_analysis) {
      $ctrl_rule_adaptor->create_rule( $wait_for_analysis, $waiting_analysis);
      warn "Created Control rule: $waiting_name will wait for $wait_for_name\n";
    } else {
      die "Could not fetch $waiting_name -> $wait_for_name to create a control rule";
    }
  }
}

sub _add_genes_to_analysis {
  my $analysis_id = shift;
  my $gene_id_arrayref = shift;
  my $params = shift || {};

  my @node_ids = ();
  foreach my $gene_id (@$gene_id_arrayref) {
    my $member;

    my $ext_member = Bio::Greg::EslrUtils->find_member_by_external_id($dba,$gene_id);
    $member = $ext_member if (defined $ext_member);

    my $member = $dba->get_MemberAdaptor->fetch_by_source_stable_id(undef,$gene_id);
    my $gene_member = $member->gene_member;
    $member = $gene_member if (defined $gene_member);

    next if (!defined $member);

    my $tree = _root_node_for_stable_id($member->stable_id);

    push @node_ids, $tree->node_id;
  }
  
  _add_nodes_to_analysis($analysis_id,\@node_ids,$params);
}

sub _add_nodes_to_analysis {
  my $analysis_id = shift;
  my $node_arrayref = shift;
  my $params = shift || {};

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

sub _root_node_for_stable_id {
  my $stable_id = shift;
  my $member = $dba->get_MemberAdaptor->fetch_by_source_stable_id(undef,$node_id);
  my $gene_member = $member->gene_member;
  my $tree = _root_node_for_member($gene_member);
  return $root_node;
}

sub _root_node_for_member {
  my $member = shift;
  my $pta = $dba->get_ProteinTreeAdaptor;
  my $tree = $pta->fetch_by_gene_Member_root_id($member);
  return $tree;
}

sub _select_node_ids {
  my $cmd = shift;
  if (!defined $cmd) {
    $cmd = "SELECT node_id FROM protein_tree_node WHERE parent_id=1 AND root_id=1";
  }
  my $sth = $dbc->prepare($cmd);
  $sth->execute();

  my $array_ref = $sth->fetchall_arrayref([0]);
  my @node_ids = @{$array_ref};
  @node_ids = map {@{$_}[0]} @node_ids;  # Some weird mappings to unpack the numbers from the arrayrefs.
  $sth->finish;
  return @node_ids;
}
