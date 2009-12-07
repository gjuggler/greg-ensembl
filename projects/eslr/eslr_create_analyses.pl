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

Bio::EnsEMBL::Registry->no_version_check(1);

#my $url = Bio::Greg::EslrUtils->defaultMysqlURL();
my $url = "mysql://ensadmin:ensembl@compara2:3306/gj1_eslr";
GetOptions('url=s' => \$url);

my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-url => $url);
my $dbc = $dba->dbc;
my $dbh = $dbc->db_handle;
my $pta = $dba->get_ProteinTreeAdaptor();
my $mba = $dba->get_MemberAdaptor();

my $module = "Bio::EnsEMBL::Compara::RunnableDB::Sitewise_dNdS";

my $analysis_hash;

my $SKIP_ADDING_JOBS = 0;
my $LIMIT = "";
#my $LIMIT = "limit 10";
one_time_use();

node_sets();
parameter_sets();
align();
sequence_quality();
omegas();
mapping();
protein_features();

#output_trees();

#connect_analysis("NodeSets","SequenceQuality",1);
#connect_analysis("Align","SequenceQuality",1);
#connect_analysis("SequenceQuality","Omegas",1);
#connect_analysis("Omegas","Mapping",1);

my $ns = 1;

sub one_time_use {
  #my $member = $mba->fetch_by_source_stable_id(undef,'ENSMLUP00000013762');
  #print $member->stable_id."\n";
  #Bio::EnsEMBL::Compara::ComparaUtils->get_quality_string_for_member($member);
  #exit(0);

  my $reg = "Bio::EnsEMBL::Registry";
  $reg->load_registry_from_url('mysql://ensro@ens-livemirror/54:3306',0); 
  my $genome_db_adaptor = $dba->get_GenomeDBAdaptor();
  foreach my $gdb (@{$genome_db_adaptor->fetch_all}) {
    #print $gdb->name."\n";
    my $core_adaptor = $reg->get_DBAdaptor($gdb->name,"Core");
    next unless ($core_adaptor);
    my $locator = $core_adaptor->locator;
    $locator .= ";species=".$gdb->name;
    $locator .= ";disconnect_when_inactive=1";
    print $locator."\n";
    $gdb->locator($locator);
    $genome_db_adaptor->store($gdb);
    #print $core_adaptor->dbc->locator."\n";
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
  my $cmd = "SELECT node_id FROM protein_tree_node WHERE parent_id=1 AND root_id=1 $LIMIT;";
  my @nodes = _select_node_ids($cmd);
  $dbc->do("DELETE from analysis_job WHERE analysis_id=$analysis_id;");
  _add_nodes_to_analysis($analysis_id,$params,\@nodes);  
}

sub parameter_sets {
  my $params;
  
  $params = {
    parameter_set_id => 1,
    parameter_set_name => "Mammals"
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_id => 7,
    parameter_set_name => "Primates",
    keep_species => "9606,9598,9544,9478,30611,30608"
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_id => 8,
    parameter_set_name => "Glires",
    keep_species => "10090,10116,43179,10020,10141,9986,9978"
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);


  $params = {
    parameter_set_id => 9,
    parameter_set_name => "Laurasiatheria",
    keep_species => "9365,42254,9796,59463,132908,30538,9739,9913,9615,9685"
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_id => 2,
    parameter_set_name => "No 2x",
    keep_species => join(",",(9796, # horse
			      9913, # cow
			      9615, # dog
			      9606, # human
			      9598, # chimp
			      9544, # macaque
			      10090,# mouse
			      10116,# rat
			      10141,# guinea pig
			      # Orang, chicken and tetraodon are hi-Q but are being ignored in all sets.
			      ))
  };
  $params = _combine_hashes($base_params,$params);
  _add_parameter_set($params);

  $params = {
    parameter_set_id => 3,
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
    parameter_set_id => 6,
    parameter_set_name => "No Seq/Aln Filtering",
    sequence_quality_filtering => 0,
    alignment_quality_filtering => 0
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
  my $base_params = {
    parameter_sets => "1,2,3,4,5,6,7,8,9",
    sequence_quality_filtering => 1,
    alignment_quality_filtering => 1,
    remove_species => join(",",(9600, # orang not allowed.
				9593, # gorilla not allowed.
				9258, # platypus -- not allowed, or just ignored?
				# Outside of Eutheria:
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
    sitewise_table => 'sitewise_aln'
  };
  _create_analysis($analysis_id,$logic_name,$module,$params,20,1);
}

sub protein_features {
  my $analysis_id=107;
  my $logic_name = "ProteinFeatures";
  my $module = "Bio::EnsEMBL::Compara::RunnableDB::ProteinFeatureCollection";
  my $params = {
    action => 'go'
    };
  _create_analysis($analysis_id,$logic_name,$module,$params,50);
}


sub clean {
  #my $analysis_id = shift;
  #$dbc->do("DELETE FROM analysis_job WHERE analysis_id=$analysis_id;");
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

#  $cmd = qq{REPLACE INTO dataflow_rule SET
#	      from_analysis_id=$analysis_id,
#	      to_analysis_url="$logic_name",
#	      branch_code=2
#	      ;
#	  };
#  $dbc->do($cmd);
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
  return if ($SKIP_ADDING_JOBS);

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
