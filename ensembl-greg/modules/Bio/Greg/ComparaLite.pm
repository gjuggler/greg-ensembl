package Bio::Greg::ComparaLite;

use File::Basename;
use Bio::EnsEMBL::Compara::NestedSet;
use Bio::EnsEMBL::Compara::AlignUtils;
use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::TreeUtils;

# Define some constants.
my $DBA = "Bio::EnsEMBL::Compara::DBSQL::DBAdaptor";


sub load_trees_for_analysis {
  my $class = shift;
  my $dba = shift;
  my $tree_dir = shift;
  
  my $nsa = $dba->get_NestedSetAdaptor();
  my $mba = $dba->get_MemberAdaptor();
  my $pta = $dba->get_ProteinTreeAdaptor();

  my $root = Bio::EnsEMBL::Compara::NestedSet->new();
  $root->name("Trees loaded from $tree_dir");
  $pta->store($root);

  my @trees = <$tree_dir/*.{nh,nhx,txt}>;
  my @seqs = <$tree_dir/*.{fasta,fa}>;

  foreach my $tree_file (sort @trees) {
    my ($name,$path,$suffix) = fileparse($tree_file,(".nh",".nhx",".txt"));
    my @aln_files = grep {$_ =~ /$name/i} @seqs;
    my $aln_file = $aln_files[0];    
    $class->load_tree_for_analysis($dba,$root,$tree_file,$aln_file);
  }
}

sub load_tree_for_analysis {
  my $class = shift;
  my $dba = shift;
  my $parent_node = shift;
  my $tree_file = shift;
  my $aln_file = shift;

  my $nsa = $dba->get_NestedSetAdaptor();
  my $mba = $dba->get_MemberAdaptor();
  my $pta = $dba->get_ProteinTreeAdaptor();

  open(IN,"$tree_file");
  my $newick_str = join("",<IN>);
  # Clean up newlines.
  $newick_str =~ s/\n//g;
  close(IN);

  my ($name,$path,$suffix) = fileparse($tree_file,(".nh",".nhx",".txt"));
  printf "Loading $tree_file: %.50s\n", $newick_str;

  my $node = Bio::EnsEMBL::Compara::TreeUtils->from_newick($newick_str);

  foreach my $leaf ($node->leaves) {
    bless $leaf, "Bio::EnsEMBL::Compara::LocalMember";
    $leaf->stable_id($leaf->name);
    $leaf->source_name("$name");
    $mba->store($leaf);
  }

  # Store the node.
  $parent_node->add_child($node);
  #$node->parent($parent_node);
  $pta->store($node);
  $node->store_tag("input_file",$tree_file);

  if (defined $aln_file) {
    # Store the alignment.
    print "Loading $aln_file...\n";
    my $sa = Bio::EnsEMBL::Compara::AlignUtils->from_file($aln_file);
    
    Bio::EnsEMBL::Compara::ComparaUtils->store_SimpleAlign_into_table("protein_tree_member",$node,$sa);
  }

  return $node;
}


# TODO: Fix this guy up! (2009-05-26)
sub load_trees_for_simulation {
  my $class = shift;
  my $dba = shift;
  my $tree_dir = shift;
  my $params = shift;
  
  my $default_params = {
    tree_dir               => "trees",
    simulation_params      => "{'function'=>'constant','a'=>0,'length'=>'150','lambda'=>'{0,0}'}",
    replicates  => 1,
  };

  $params = {}; #Bio::EnsEMBL::Compara::ComparaUtils->

  # Load from a configuration file if available.
  use Bio::Greg::EslrUtils;
  $params = Bio::Greg::EslrUtils->loadConfiguration($params,$config) if (defined $config);
  my $base_sim_params = _string_to_hash($params->{'simulation_params'});

  my $nsa = $dba->get_NestedSetAdaptor();
  my $mba = $dba->get_MemberAdaptor();
  my $pta = $dba->get_ProteinTreeAdaptor();

  # The list of simulation replicates to use. These will be stored as tags in the XYZ_tag table,
  # and accessible using the $node->get_tagvalue("tag_key") method.
  my @simulation_sets = grep {$_ =~ /simset_/} keys %{$params};

  my $root = Bio::EnsEMBL::Compara::NestedSet->new();
  $root->name("Simulation root");
  $pta->store($root);
  $root->store_tag("simulation_params",$params->{'simulation_params'});
  
  my $tree_dir = $params->{'tree_dir'};
  my @files = <$tree_dir/*.{nh,nhx,txt}>;
  foreach my $file (sort @files) {   
    open(IN,"$file");
    my $newick_str = join("",<IN>);
    close(IN);

    foreach my $sim_set (@simulation_sets) {
      my $temp_sim_params = _copy_hash($base_sim_params);
      my $sim_param_str = $params->{$sim_set};
      my $sim_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_string($temp_sim_params,$sim_param_str);
      my $replicates = $sim_params->{'replicates'};

      foreach my $sim_rep (1 .. $replicates) {
	print "$file $sim_set $sim_rep\n";	
	my $node = Bio::EnsEMBL::Compara::TreeUtils->from_newick($newick_str);
	
	# Go through each leaf and store the member objects.
	foreach my $leaf ($node->leaves) {
	  $leaf->stable_id($leaf->name);
	  #$leaf->source_name("ENSEMBLGENE");
	  $leaf->source_name($sim_set."-".$sim_rep);
	  $mba->store($leaf);
	}

	# Store the node.
	$root->add_child($node);
	$pta->store($node);
	
	# Store the tags.
	
	$node->store_tag("input_file",$file);
	$node->store_tag("sim_set",$sim_set);
	$node->store_tag("sim_rep",$sim_rep);
	$node->store_tag("simulation_params",$sim_param_str);
      }
    }
  }
}

sub create_analysis {
  my $class = shift;
  my $dba = shift;
  my $analysis_id = shift;
  my $logic_name = shift;
  my $module = shift;
  my $params = shift;
  my $hive_capacity = shift || 500;
  my $batch_size = shift || 1;
  
  my $dbc = $dba->dbc;

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
  $cmd = qq{REPLACE INTO dataflow_rule SET
	      from_analysis_id=$analysis_id,
	      to_analysis_url="$logic_name",
	      branch_code=2
	      ;
	  };
  $dbc->do($cmd);
}



sub url_to_mysql_args {
  my $class = shift;
  my $url = shift;

  $url =~ m!mysql://(.*?):(.*?)@(.*?)/(.*)!;
  
  # Something like mysql -uensadmin -pensembl -hensdb-2-12 -P5106 -A gj1_compara_53
  # $1 = user, $2 = pw, $3 = server+port, $4 = db
  my ($u,$p,$h,$P,$db) = undef;
  $u = $1;
  $p = $2;
  $h = $3;
  $P = 3306;
  $db = $4;
  if ($h =~ m/(.*?):(.*)/) {
    $h = $1;
    $P = $2;
  }

  my $mysqlargs = sprintf("-u%s -p%s -h%s -P%s -A %s",
			  $u,
			  $p,
			  $h,
			  $P,
			  $db);
  return $mysqlargs;
}

sub url_to_hashref {
  my $class = shift;
  my $url = shift;

  my $args = $class->hive_url_to_mysql_args($url);

  $args =~ m/-u(.*) -p(.*) -h(.*) -P(.*) -A (.*)/;
  my $data = {user => $1,
	      password => $2,
	      host => $3,
	      port => $4,
	      database => $5};
  return $data;
}

# Creates a database $database in the given mysql server at $url.
# Returns the new full URL for accessing the database.
sub create_database {
  my $class = shift;
  my $url = shift;
  my $format_database = shift;

  $format_database = 0 unless (defined $format_database);

  my $obj = $class->hive_url_to_hashref($url);
  $url = sprintf("mysql://%s:%s@%s/",
		 $obj->{'user'},
		 $obj->{'password'},
		 $obj->{'host'}.":".$obj->{'port'}
		 );
  my $database = $obj->{'database'};
  #print "BASE URL: $url\n";
  my $temp_url = $url . "mysql";

  my $dba = $DBA->new(-url => $temp_url);
  eval {
    $dba->dbc->do("drop database if exists $database;") if ($format_database);
  };
  eval {
    $dba->dbc->do("create database if not exists $database;");
  };

}

# Creates the SQL tables necessary for running a Hive.
sub create_hive_tables {
  my $class = shift;
  my $url = shift; 		# URL of the desired MySQL location INCLUDING the database name.
  my $format_database = shift;	# Set to 1 if we should wipe the database clean.
  my $hive_sql_file = shift;

  open(IN,"$hive_sql_file");
  my $hive_sql_string = join("",<IN>);
  close(IN);

  $format_database = 1 unless (defined $format_database);

  my $dba = $DBA->new(-url => $url);
  # Create the tables.
  $class->run_sql_commands($dba,$hive_sql_string);
  # Format if necessary.
  if ($format_database) {
    print "Creating hive database: $url \n";
    my @hive_tables = qw(hive analysis analysis_stats analysis_stats_monitor analysis_ctrl_rule dataflow_rule
			 analysis_job analysis_job_file analysis_data analysis_description monitor
			 );
    map {eval{ $dba->dbc->do("truncate $_;")}} @hive_tables; # Truncate each table in turn.
  }
}

# Creates the SQL tables necessary for housing a ComparaLite database.
sub create_compara_tables {
  my $class = shift;
  my $url = shift;
  my $format_database = shift;
  my $compara_sql_file = shift;

  open(IN,"$compara_sql_file");
  my $compara_sql_string = join("",<IN>);
  close(IN);

  $format_database = 1 unless (defined $format_database);

  my $dba = $DBA->new(-url => $url);
  # Create the tables.
  $class->run_sql_commands($dba,$compara_sql_string);
  # Format if necessary.
  if ($format_database) {
    print "Creating Compara database: $url \n";
    my @hive_tables = qw(
			 member sequence protein_tree_member protein_tree_node sitewise_aln
			 );
    map {eval{ $dba->dbc->do("truncate $_;")}} @hive_tables; # Truncate each table in turn.
  }  
}

sub run_sql_commands {
  my $class = shift;
  my $dba = shift;
  my $string = shift;

  # Remove comment lines.
  $string =~ s/^#//g;

  # Split the string into separate commands.
  my @sql_array = split(";",$string);
  # Remove whitespace lines.
  @sql_array = grep {$_ !~ m/^\s+$/} @sql_array;
  # Use the map function to run each command in turn.
  map {eval {$dba->dbc->do($_.";")};} @sql_array; # Run the commands in turn.
}

# Create a string from a hashref.
sub hash_to_string {
  my $class = shift;
  my $hashref = shift;

  my %hash = %{$hashref};

  my $str = "{";
  my @keyvals = ();
  foreach my $key (sort keys %hash) {
    push(@keyvals, "'$key'" . "=>" . "'" . $hash{$key} . "'");
  }

  $str .= join(",",@keyvals);
  $str .= "}";
  return $str;
}

1;
