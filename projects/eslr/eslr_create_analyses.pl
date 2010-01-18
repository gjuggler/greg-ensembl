#!/usr/bin/env perl

use warnings;
use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::Greg::EslrUtils;
use File::Path;
use File::Basename;

my $url = 'mysql://ensadmin:ensembl@compara2:3306/gj1_eslr';
my $clean = 1;
my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-url => $url);
my $dbc = $dba->dbc;
my $analysis_hash; # Hash to store mapping between analysis names and ID numbers.


###
# Most of the pipeline is defined with the following function calls.
###

parameter_sets();
node_sets();
omegas();
mapping();

connect_analysis("NodeSets","Omegas",1);
connect_analysis("Omegas","Mapping",1);

#print_database_tree();


sub clean_tables {
  if ($clean) {    
    my @truncate_tables = qw^
      analysis_job dataflow_rule hive
      sitewise_omega
      parameter_set
      node_set_member node_set
      ^;
    map {print "$_\n";eval {$dba->dbc->do("truncate table $_");}} @truncate_tables;
  }
}

sub node_sets {
  my $analysis_id=101;
  my $logic_name = "NodeSets";
  my $module = "Bio::Greg::NodeSets";
  my $params = {
    flow_node_set => 8
  };
  _create_analysis($analysis_id,$logic_name,$module,$params,30,1);

  $params = {};
  my $cmd = "SELECT node_id FROM protein_tree_node WHERE parent_id=1 AND root_id=1;";
  my @nodes = _select_node_ids($cmd);
  _add_nodes_to_analysis($analysis_id,$params,\@nodes);  
}

sub parameter_sets {
  my $params;
  my $base_params={};
  
  $params = {
    parameter_set_id => 1,
    parameter_set_name => "Mammals"
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_id => 2,
    parameter_set_name => "Primates",
    keep_species => "9606,9598,9544,9478,30611,30608"
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_id => 3,
    parameter_set_name => "Glires",
    keep_species => "10090,10116,43179,10020,10141,9986,9978"
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);


  $params = {
    parameter_set_id => 4,
    parameter_set_name => "Laurasiatheria",
    keep_species => "9365,42254,9796,59463,132908,30538,9739,9913,9615,9685"
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_id => 5,
    parameter_set_name => "Afrotheria",
    keep_species => "9813,9785,9731"
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_id => 6,
    parameter_set_name => "No 2x",
    keep_species => join(",",(9796, # horse
			      9913, # cow
			      9615, # dog
			      9606, # human
			      9600, # orang
			      9593, # gorilla
			      9598, # chimp
			      9544, # macaque
			      10090,# mouse
			      10116,# rat
			      10141,# guinea pig
			      # Chicken and tetraodon are hi-Q but are being ignored in all sets.
 			      ))
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_id => 7,
    parameter_set_name => "2x Only",
    keep_species => join(",",(9365,   # hedgehog
			      42254,  # shrew
			      59463,  # microbat
			      132908, # megabat
			      30538,  # alpaca
			      9739,   # dolphin
			      9685,   # cat
			      9478,   # tarsier
			      30611,  # bushbaby
			      30608,  # mouse lemur
			      43179,  # squirrel
			      10020,  # kangaroo rat
			      9986,   # rabbit
			      9978,   # pika
			      37347,  # tree shrew
			      9361,   # armadillo
			      9358,   # sloth
			      9813,   # hyrax
			      9785,   # elephant
			      9371    # tenrec
			      ))
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_id => 8,
    parameter_set_name => "No Primates",
    remove_species => join(",",(9606,   # human
				9598,   # chimp
				9544,   # macaque
				9478,   # tarsier
				30611,  # bushbaby
				30608,  # mouse lemur
				9593,   # gorilla
				9600,   # orang
				9483    # marmoset
				))
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

}

sub align {
  my $analysis_id = 102;
  my $logic_name = "Align";
  my $module = "Bio::EnsEMBL::Compara::RunnableDB::MCoffee";
  my $params = {
    method => 'cmcoffee'
  };

  _create_analysis($analysis_id,$logic_name,$module,$params,400,1);
}

sub sequence_quality {
  my $analysis_id=103;
  my $logic_name = "SequenceQuality";
  my $module = "Bio::Greg::SequenceQualityLoader";
  my $params = {};
  
  _create_analysis($analysis_id,$logic_name,$module,$params,100,1);
}

sub omegas {
  my $analysis_id = 105;
  my $logic_name = "Omegas";
  my $module = "Bio::EnsEMBL::Compara::RunnableDB::Sitewise_dNdS";
  my $base_params = {
    parameter_sets => "all",
    sequence_quality_filtering => 1,
    alignment_quality_filtering => 1,
    remove_species => join(",",(
				# Outside of Eutheria:
				9258, # platypus
  			        13616,# opossum
				9031, # chicken
                                9103, # turkey!!
				59729,# zebra finch
				28377,# anole lizard
				8364, # xenopus
				69293,# stickleback
				8090, # medaka
				99883,# tetraodon
				31033,# fugu
				7955, # zebrafish
				7719, # c. intestinalis
				51511,# c. savignyi
				7165, # anopheles
				7159, # aedes
				7227, # fruitfly
				6239, # c. elegans
				4932  # s. cerevisiae
				))
    };
  _create_analysis($analysis_id,$logic_name,$module,$base_params,400,1);

}

sub mapping {
  my $analysis_id=106;
  my $logic_name = "Mapping";
  my $module = "Bio::Greg::SitewiseMapper";
  my $params = {
    sitewise_table => 'sitewise_omega'
  };
  _create_analysis($analysis_id,$logic_name,$module,$params,20,1);
}

sub print_database_tree {
  my $nhx = Bio::EnsEMBL::Compara::ComparaUtils->get_genome_tree_nhx($dba,{labels => 'mnemonics',
									   images => 1});

  print $nhx."\n";
}

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


sub _add_parameter_set {
  my $params = shift;
  my $parameter_set_id = $params->{'parameter_set_id'};
  
  if (exists $params->{'parameter_set_name'} ) {
    my $parameter_set_name = $params->{'parameter_set_name'};
    delete $params->{'parameter_set_name'};
    my $name_cmd = "REPLACE INTO parameter_set VALUES ('$parameter_set_id','name',\"$parameter_set_name\");";
    $dbc->do($name_cmd);
  }
  
  my $param_string = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
  my $cmd = "REPLACE INTO parameter_set VALUES ('$parameter_set_id','params',\"$param_string\");";
  $dbc->do($cmd);
}


########*********########
#-------~~~~~~~~~-------#
########*********########

sub _create_analysis {
  my $analysis_id = shift;
  my $logic_name = shift;
  my $module = shift;
  my $params = shift;
  my $hive_capacity = shift || 500;
  my $batch_size = shift || 1;
  
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
  $cmd = qq{UPDATE analysis_stats SET
	      hive_capacity=$hive_capacity,
	      batch_size=$batch_size
	      WHERE analysis_id=$analysis_id
	      ;};
  $dbc->do($cmd);
}

sub connect_analysis {
  my $from_name = shift;
  my $to_name = shift;
  my $branch_code = shift;

  my $from_id = $analysis_hash->{$from_name};
  $branch_code = 1 unless (defined $branch_code);
  my $cmd = qq{REPLACE INTO dataflow_rule SET
		 from_analysis_id=$from_id,
		 to_analysis_url="$to_name",
		 branch_code=$branch_code
		 ;};
  $dbc->do($cmd);
}

sub _add_nodes_to_analysis {
  my $analysis_id = shift;
  my $params = shift || {};
  my $node_arrayref = shift;

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
  my $sth = $dbc->prepare($cmd);
  $sth->execute();

  my $array_ref = $sth->fetchall_arrayref([0]);
  my @node_ids = @{$array_ref};
  @node_ids = map {@{$_}[0]} @node_ids;  # Some weird mappings to unpack the numbers from the arrayrefs.
  $sth->finish;
  return @node_ids;
}
