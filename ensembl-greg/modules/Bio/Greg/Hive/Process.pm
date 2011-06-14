package Bio::Greg::Hive::Process;

use strict;
use Bio::EnsEMBL::Utils::Argument;
use Bio::EnsEMBL::Utils::Exception;

use Data::UUID;
use Data::Types qw(:all);

use Digest::MD5 qw(md5_hex);

use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Hive::Utils 'stringify';    # import 'stringify()'

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

use File::Path qw(make_path);
use File::Basename;
use POSIX qw(strftime mktime);

use FreezeThaw qw(freeze thaw cmpStr safeFreeze cmpStrHard);

use Bio::Greg::EslrUtils;

use base ('Bio::EnsEMBL::Hive::Process');

sub get_tree {
  my $self   = shift;
  my $params = shift;

  $params = $self->params unless ( defined $params );

  my $tree =
    Bio::EnsEMBL::Compara::ComparaUtils->get_tree_for_comparative_analysis( $self->compara_dba,
    $params );

  if ($tree == -1) {
    $self->fail_and_die("Failed to get tree (maybe too few nodes)");
  }

  return $tree;
}

sub get_aln {
  my $self   = shift;
  my $params = shift;

  return $self->_get_aln( $params, 0 );
}

sub get_cdna_aln {
  my $self   = shift;
  my $params = shift;

  return $self->_get_aln( $params, 1 );
}

sub get_orig_aln {
  my $self = shift;

  $self->get_aln if ( !defined $self->param('full_aln_aa') );
  return $self->param('full_aln_aa');
}

sub get_orig_aln_position {
  my $self         = shift;
  my $aln          = shift;
  my $aln_position = shift;

  my $orig_aln = $self->get_orig_aln;
  return Bio::EnsEMBL::Compara::AlignUtils->map_alignment_position( $aln, $aln_position,
    $orig_aln );
}

sub _get_aln {
  my $self   = shift;
  my $params = shift;
  my $cdna   = shift;

  $params = $self->params unless ( defined $params );

  $params->{remove_blank_columns} = 1
    unless ( defined $params->{remove_blank_columns} );

  my $tree;
  if ( defined $params->{tree} ) {
    $tree = $params->{tree};
  } else {
    $tree = $self->get_tree($params);
  }

  my $aa_aln = $tree->get_SimpleAlign();
  my $cdna_aln = $tree->get_SimpleAlign( -cdna => 1, -hide_positions => 1 );
  my $aln =
    Bio::EnsEMBL::Compara::ComparaUtils->fetch_masked_alignment( $aa_aln, $cdna_aln, $tree, $params,
    $cdna );

  # Remove blank columns here. Don't forget to store the 'full' alignment as a parameter.
  if ( $params->{remove_blank_columns} ) {
    if ($cdna) {
      $self->param( 'full_aln_cdna', $aln );
      ($aln) = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns_in_threes($aln);
    } else {
      $self->param( 'full_aln_aa', $aln );
      ($aln) = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns($aln);
    }
  }

  # Here's where we'll split up alignments by slice.
  if ( $params->{alignment_slice} ) {
    my $slice_string = $params->{alignment_slice};
    my ( $lo, $hi ) = split( "-", $slice_string );

    #$lo = 1;
    #$hi = 10;
    if ($cdna) {
      $lo = ($lo) * 3 - 2;
      $hi = ($hi) * 3;
    }
    print "$lo $hi\n";
    print "Length: " . $aln->length . "\n";
    $aln = $aln->slice( $lo, $hi, 1 );
    print "Length: " . $aln->length . "\n";
  }

  return $aln;
}

sub fail {
  my $self    = shift;
  my $message = shift;

  $self->input_job->update_status('FAILED');
  warn($message);
}

sub fail_and_die {
  my $self    = shift;
  my $type = shift;
  my $message = shift;

  $self->create_table_from_params($self->hive_dba, 'failed_jobs', $self->failed_job_table_structure);

  # Create an entry in the failed job table structure.
  $self->param('analysis_id', $self->input_job->analysis_id);
  my $analysis = $self->hive_dba->get_AnalysisAdaptor->fetch_by_dbID($self->input_job->analysis_id);
  $self->param('logic_name', $analysis->logic_name);
  $self->param('job_id', $self->input_job->dbID);

  $self->param('fail_type', $type);
  $self->param('fail_msg', $message);
  $self->store_params_in_table($self->hive_dba, 'failed_jobs', $self->params);

  $self->input_job->update_status('FAILED');
  die($message);
}

sub failed_job_table_structure {
  return {
    job_id => 'int',
    analysis_id => 'int',
    logic_name => 'string',
    
    node_id => 'int',
    gene_name => 'string',
    gene_id => 'string',

    fail_type => 'string',
    fail_msg => 'string',

    unique_keys => 'job_id,analysis_id,node_id'
  };
}

sub check_tree_aln {
  my $self = shift;
  my $tree = shift;
  my $aln  = shift;

  if ( !defined $tree ) {
    $tree = $self->get_tree;
  }

  # Look at the tree size.
  if ( scalar $tree->leaves <= 1 ) {
    return -1;
  }

  if ( !defined $aln ) {
    $aln = $self->get_aln;
  }

  # Look at the alignment size.
  if ( scalar $aln->each_seq <= 1 || $aln->length == -1 ) {
    return -1;
  }

  # More complicated stuff -- look for alignments where only one sequence isn't entirely X'ed out.
  my $good_seq_count = 0;
  foreach my $seq ( $aln->each_seq ) {
    my $seq_str = $seq->seq;
    $good_seq_count++ if ( $seq_str =~ m/[^NX-]/i );
  }
  print "Good seq count: $good_seq_count\n";
  if ( $good_seq_count < 2 ) {
    return -1;
  }
  return 1;
}

sub pretty_print {
  my $self   = shift;
  my $aln    = shift;
  my $params = shift || {};

  if (ref($aln) =~ m/align/i) {
    Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $aln, $params );
  } elsif (ref($aln) =~ m/proteintree/i) {
    $aln->print_tree;
  }
}

sub load_registry {
  my $self = shift;

  Bio::EnsEMBL::Compara::ComparaUtils->load_registry();
}

sub compara_dba {
  my $self = shift;

  if ( !defined $self->{_compara_dba} ) {
    print " >> Getting new DBA!!!!\n" if ($self->debug);
    my $compara_dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new( -DBCONN => $self->db->dbc );
    eval {
      $compara_dba->dbc->do("select * from member limit 1;");
      1;
      }
      or do {
      print "ERROR: $@\n";
      print " >> No compara in hive DB -- falling back to ens-livemirror!!\n";
      $self->load_registry();
      $compara_dba = Bio::EnsEMBL::Registry->get_DBAdaptor( 'multi', 'compara' );
      
      printf " >> Using Compara DB at [%s/%s]\n",$compara_dba->dbc->host,$compara_dba->dbc->dbname;
      #      $compara_dba->disconnect_when_inactive(1);
      };
    $self->{_compara_dba} = $compara_dba;
  }
  return $self->{_compara_dba};
}

sub hive_dba {
  my $self = shift;

  if ( !defined $self->{_hive_dba} ) {
    #print " >>>> Getting new Hive DBA!!!!\n";
    my $hive_dba = Bio::EnsEMBL::Hive::DBSQL::DBAdaptor->new( -DBCONN => $self->db->dbc );
    $self->{_hive_dba} = $hive_dba;
  }
  return $self->{_hive_dba};
}

sub pta {
  my $self = shift;

  #  return $self->compara_dba->get_ProteinTreeAdaptor;

  if ( !defined $self->param('_pta') ) {
    #print " >>>> Getting new PTA!!!!\n";
    my $compara_dba = $self->compara_dba;
    my $pta         = $compara_dba->get_ProteinTreeAdaptor;
    $self->param( '_pta', $pta );
  }
  return $self->param('_pta');
}

sub mba {
  my $self = shift;

  #return $self->compara_dba->get_MemberAdaptor;

  if ( !defined $self->param('_mba') ) {
    #print " >>>> Getting new MBA!!!!\n";
    my $compara_dba = $self->compara_dba;
    my $pta         = $compara_dba->get_MemberAdaptor;
    $self->param( '_mba', $pta );
  }
  return $self->param('_mba');
}

sub db_handle {
  my $self  = shift;
  my $force = shift;

  warn('Shouldn\'t be accessing db_handle directly -- use $dbc->do and $dbc->prepare');

  if ( !defined $self->param('_compara_dbh') || $force == 1 ) {
    #print " >>>> Getting new Compara DBH!!!!\n";
    my $dbc = $self->dbc;
    my $dbh = $dbc->db_handle;
    $self->param( '_compara_dbh', $dbh );
    print $dbh. "\n";
  }

  my $dbh = $self->param('_compara_dbh');

  my $ping = $dbh->ping;
  #print "PING: $ping\n";

  #  if (!$ping) {
  #    $dbh = $dbh->clone;
  #    $self->param('_compara_dbh',$dbh);
  #  }
  return $dbh;
}

sub disconnect_when_inactive {
  my $self   = shift;
  my $discon = shift;

  $self->dbc->disconnect_when_inactive($discon);
}

sub within_hive {
  my $self = shift;

  return $self->db != undef;
}

sub dbc {
  my $self = shift;

  return $self->db->dbc;

  if ( !defined $self->param('_compara_dbc') ) {
    #print " >>>> Getting new Compara DBC!!!!\n";
    my $compara_dba = $self->compara_dba;
    my $dbc         = $compara_dba->dbc;

    # It's apparently important to turn off the inactive disconnect here, since we'll
    # be sharing this DBC throughout the lifetime of this Process.
    # TODO: Think of how to handle the case when a Process wants to let a connection
    # run idle...
    $dbc->disconnect_when_inactive(1);
    $self->param( '_compara_dbc', $dbc );
  }
  return $self->param('_compara_dbc');
}

sub node_id {
  my $self = shift;

  my $node_id = 0;
  $node_id = $self->param('node_id') if ( defined $self->param('node_id') );

  return $node_id;
}

sub data_id {
  my $self = shift;

  my $data_id = 0;
  $data_id = $self->param('data_id') if ( defined $self->param('data_id') );
  return $data_id;
}

sub stable_id {
  my $self = shift;
  my $tree = shift;
  my $s_id = $self->_get_pep->stable_id;
  return $s_id || 'unknown_gene';
}

sub gene_name {
  my $self = shift;
  my $tree = shift;
  my $name = $self->_get_pep->display_label;
  return $name || 'unknown_gene';
}

sub _get_pep {
  my $self = shift;
  
  my $pta = $self->compara_dba->get_ProteinTreeAdaptor;  
  my $tree = $pta->fetch_node_by_node_id($self->node_id);

  # Try human first.
  my $pep;
  ($pep) = grep {$_->taxon_id == 9606} $tree->leaves;
  if (!defined $pep) {
    ($pep) = $tree->leaves;
  }
  return $pep;
}

sub parameter_set_id {
  my $self = shift;

  my $id = 0;
  $id = $self->param('parameter_set_id')
    if ( defined $self->param('parameter_set_id') );
  return $id;
}

sub worker_temp_directory {
  my $self = shift;

  my $wtd = $self->SUPER::worker_temp_directory;

  if (!$self->within_hive) {
    my $pid = $$;
    $wtd = "/tmp/process_tmp_${pid}/";
    make_path($wtd);
  }

  chmod 0777, $wtd;
  return $wtd;
}

sub store_meta {
  my $self   = shift;
  my $params = shift;

  foreach my $key ( keys %$params ) {
    my $value = $params->{$key};
    $self->dbc->do("DELETE from meta where meta_key='$key';");
    $self->dbc->do("REPLACE into meta (meta_key,meta_value) VALUES ('$key','$value');");
  }
}

sub get_meta {
  my $self = shift;
  my $key = shift;

  my $sth = $self->dbc->prepare("SELECT * from meta where meta_key='$key' limit 1;");
  $sth->execute;
  while ( my $obj = $sth->fetchrow_hashref ) {
    $sth->finish;
    return $obj->{meta_value};
  }
}

sub get_parameter_sets {
  my $self = shift;
  my $dbname = shift || $self->dbc->dbname;

  my $sth = $self->dbc->prepare("SELECT * from ${dbname}.parameter_set;");
  $sth->execute;

  my $hash;
  while ( my $obj = $sth->fetchrow_hashref ) {
    my $id = $obj->{parameter_set_id};
    $hash->{$id} = {} unless ( defined $hash->{$id} );
    my $param_name  = $obj->{parameter_name};
    my $param_value = $obj->{parameter_value};
    $hash->{$id}->{$param_name} = $param_value;
    $hash->{$id}->{parameter_set_id} = $id;
  }
  $sth->finish;

  my @sets;
  foreach my $pset ( 1 .. scalar( keys %$hash ) ) {
    push @sets, $hash->{$pset};
  }
  return @sets;

}

sub new_data_id {
  my $self   = shift;
  my $params = shift;

  my $dbc = $self->dbc;
  $dbc->do("LOCK TABLES meta WRITE");

  my $sth = $dbc->prepare(
    "SELECT meta_value from meta WHERE meta_key='data_id_counter';");
  $sth->execute;
  my @row = $sth->fetchrow_array;
  $sth->finish;
  my $data_id = -1;
  if (@row) {
    $data_id = $row[0];
  }

  $data_id++;
  $dbc->do(
    "DELETE from meta where meta_key='data_id_counter'");
  $dbc->do("REPLACE into meta (meta_key,meta_value) VALUES ('data_id_counter',$data_id);");

  $self->param( 'data_id', $data_id );
  $params->{data_id} = $data_id;

  $dbc->do("UNLOCK TABLES");

  return $data_id;
}

sub param {
  my $self  = shift;
  my $param = shift;

  if (!$self->within_hive) {
    if (@_) {
      my $value = shift @_;
      print "Set local param: [$param] = [$value]\n"; 
      $self->{'_'.$param} = $value;
    }
    return $self->{'_'.$param};
  }

  my $param_value;
  my $orig_value = $self->SUPER::param($param);
  if (@_) {
    $self->SUPER::param( $param, shift @_ );
    $param_value = $self->SUPER::param($param);
    if ($self->debug && $param_value ne $orig_value) {
      print "Set $param => [$param_value]\n";
    }
  } else {
    $param_value = $self->SUPER::param($param);
  }
  
  return $param_value;
}

sub params {
  my $self = shift;

  # Make a copy!
  if (!$self->within_hive) {
    my $param_hash;
    foreach my $key (keys %{$self}) {
      next unless ($key =~ m/^_/);
      #print $key."\n";
      my $fixed_key = substr($key,1);
      $param_hash->{$fixed_key} = $self->{$key};
    }
    return $param_hash;
  }
  my $param_hash = $self->input_job->{_param_hash};
  my $new_params = {};
  foreach my $key ( keys %$param_hash ) {
    $new_params->{$key} = $param_hash->{$key};
  }
  return $new_params;
}

sub load_all_params {
  my $self = shift;

  my $node_id          = $self->param('node_id');
  my $parameter_set_id = $self->param('parameter_set_id');

  my $tree_tag_params = {};
  if ( defined $node_id ) {
    $tree_tag_params = $self->get_params_from_tree_tags( $self->compara_dba, $node_id );
    $self->set_params($tree_tag_params) if ($tree_tag_params);
  }

  my $param_set_params = {};
  if ( defined $parameter_set_id ) {
    $param_set_params = $self->get_params_from_param_set($parameter_set_id);
    $self->set_params($param_set_params) if ($param_set_params);
  }

  if ( !defined $self->param('data_id') ) {
    $self->param( 'data_id', $node_id );
  }
  if ( !defined $self->param('parameter_set_id') ) {
    $self->param( 'parameter_set_id', 0 );
  }
  if ( !defined $self->param('job_id')) {
    $self->param('job_id',$self->job_id);
  }

  if ($self->debug) {
    print "Bio::Greg::Hive::Process.pm - load all params\n";
    $self->hash_print( $self->params );
  }
}

sub store_param {
  my $self = shift;
  my $param = shift;
  my $value = shift;
  my $flag = shift;

  # First, store in the local param hashref.
  $self->param($param, $value);

  # Now, make sure this param is stored in our tables structure.
  my $str = $self->{_table_structure};
  if ($flag) {
    $str->{$param} = $flag;
  } else {
    # Default to string value.
    $str->{$param} = 'string';
    if (is_int($value)) {
      $str->{$param} = 'int';
    } elsif (is_float($value)) {
      $str->{$param} = 'float';
    }
  }
}

sub job_id {
  my $self = shift;
  return $self->input_job->dbID;
}

sub get_params_from_tree_tags {
  my $self    = shift;
  my $dba     = shift;
  my $node_id = shift;

  my $pta  = $dba->get_ProteinTreeAdaptor;

  print "node: $node_id\n";
  my $tree = $pta->fetch_node_by_node_id($node_id);

  my $tags = $tree->get_tagvalue_hash;

  my @param_objs;
  my $simple_tags;
  foreach my $tag ( keys %$tags ) {
    if ( $tag =~ /params/i ) {
      push @param_objs, $self->load_params_from_tag( $tree, $tag );
    } else {
      $simple_tags->{$tag} = $tags->{$tag};
    }
  }
  push @param_objs, $simple_tags;

  return $self->replace_params(@param_objs);
}

sub load_params_from_tag {
  my $self = shift;
  my $tree = shift;
  my $tag  = shift;

  my $tag_string = $tree->get_tagvalue($tag);
  return $self->string_to_hash($tag_string);
}

sub store_message {
  my $self = shift;
  my $msg = shift;

  $self->db()->get_JobMessageAdaptor()->register_message($self->input_job->dbID, $msg, 0 );
}

# Eval a hashref from a string.
sub string_to_hash {
  my $self   = shift;
  my $string = shift;

  my $hash = eval $string;
  return $hash if ( defined $hash );

  return {};
}

sub hash_print {
  my $class   = shift;
  my $hashref = shift;
  my $HANDLE  = shift;

  select($HANDLE) if ( defined $HANDLE );

  print "  {\n";
  foreach my $key ( sort keys %{$hashref} ) {
    printf "    %-40.40s => %-40s\n", $key, $hashref->{$key};
  }
  print "  }\n";

  select(STDOUT);
}

sub clone {
  my $self   = shift;
  my $params = shift;

  return $self->replace( $params, {} );
}

sub replace {
  my $self = shift;
  return $self->replace_params(@_);
}

sub replace_params {
  my $self   = shift;
  my @params = @_;

  my $final_params = shift @params;
  while (@params) {
    my $p = shift @params;
    if ( !ref $p ) {
      $p = eval($p);
    }

    my $new_params = {};
    foreach my $key ( keys %{$final_params} ) {
      $new_params->{$key} = $final_params->{$key};
    }
    foreach my $key ( keys %{$p} ) {
      $new_params->{$key} = $p->{$key};
    }
    $final_params = $new_params;
  }

  return $final_params;
}

sub get_params {
  my $self = shift;

  return $self->params;
}

sub set_params {
  my $self = shift;
  my $params = shift;

  my $debug_state = $self->debug;
  $self->debug(0);
  print("Set params [".$params."]\n");

  foreach my $key (keys %$params) {
    $self->param($key,$params->{$key});
  }

  $self->debug($debug_state);
}

sub get_params_from_string {
  my $self         = shift;
  my $param_string = shift;

  my $param_object;
  my $tmp_params = eval($param_string);
  foreach my $key ( keys %$tmp_params ) {
    $param_object->{$key} = $tmp_params->{$key};
  }
  return $param_object;
}

sub get_params_from_param_set {
  my $self         = shift;
  my $param_set_id = shift;

  my $dba = $self->hive_dba;
  my $dbc = $dba->dbc;

  throw "Undefined parameter_set_id!" unless ( defined $param_set_id );

  my $params;

  my $cmd =
    qq^SELECT parameter_value FROM parameter_set WHERE parameter_set_id=$param_set_id AND parameter_name="params";  ^;
  my $sth = $dbc->prepare($cmd);
  $sth->execute();
  my @row;
  while ( @row = $sth->fetchrow_array ) {
    $params = eval( $row[0] );
  }
  $sth->finish;

  my $cmd =
    qq^SELECT parameter_value FROM parameter_set WHERE parameter_set_id=$param_set_id AND parameter_name="name";  ^;
  $sth = $dbc->prepare($cmd);
  $sth->execute();
  my @row;
  while ( @row = $sth->fetchrow_array ) {
    $params->{parameter_set_name} = $row[0];
  }
  $sth->finish;

  return $params;
}

sub create_table_from_params {
  my $self       = shift;
  my $dba        = shift;
  my $table_name = shift;
  my $params     = shift;

  $self->dbc->disconnect_when_inactive(0);

  my $dbc = $self->dbc;

  # First, create the table with data_id and parameter_set columns.
  my $sth;

  eval {
  $sth = $dbc->do("select count(*) from $table_name where 1=0");
  };
  if ($@) {
    # Create the table shell.
    $dbc->do(
      qq^
      CREATE TABLE $table_name (
	data_id INT(10) NOT NULL
      ) engine=InnoDB;
      ^
      );
  }

  # Collect all existing fields into a hashref.
  $sth = $dbc->prepare("show fields from $table_name;");
  $sth->execute;
  my $fields;
  while ( my $obj = $sth->fetchrow_hashref ) {
    $fields->{$obj->{Field}} = 1;
  }
  $sth->finish;
  
  my $unique_keys = delete $params->{'unique_keys'};
  my $extra_keys  = delete $params->{'extra_keys'};
  foreach my $key ( sort keys %$params ) {
    if ($fields->{$key} == 1) {
      #print "  key $key exists!\n";
      next;
    } else {
      print "  creating column $key\n";
      $fields->{$key} = 1;
    }

    my $type = $params->{$key};
    
    my $not_null = 0;
    if ( $type =~ m/not null/gi ) {
      $type =~ s/\s*not null\s*//gi;
      $not_null = 1;
    }
    
    my $type_map = {
      'uuid'      => 'BINARY(16)',
      'timestamp' => 'DATETIME',
      'tinyint'   => 'TINYINT',
      'smallint'  => 'SMALLINT',
      'int'       => 'INT',
      'string'    => 'TEXT',
      'char1'     => 'CHAR(1)',
      'char2'     => 'CHAR(2)',
      'char4'     => 'CHAR(4)',
      'char8'     => 'CHAR(8)',
      'char16'    => 'CHAR(16)',
      'char32'    => 'CHAR(32)',
      'char64' => 'CHAR(64)',
      'float'     => 'FLOAT'
    };
    $type = $type_map->{$type};
    
    $type = $type . 'NOT NULL ' if ( $not_null == 1 );
    
    my $create_cmd = qq^ALTER TABLE $table_name ADD COLUMN `$key` $type^;
    $dbc->do($create_cmd);
  }

  $sth = $dbc->prepare("show indexes from $table_name;");
  $sth->execute;
  my $indexes;
  while ( my $obj = $sth->fetchrow_hashref ) {
    my $key_name = $obj->{Key_name};
    my $i = $obj->{Seq_in_index};
    my $colname = $obj->{Column_name};
    $indexes->{$key_name} = {} if (!defined $indexes->{$key_name});
    $indexes->{$key_name}->{$i} = $colname;
  }
  $sth->finish;
  
  if ($unique_keys) {
    my @toks = split(',', $unique_keys);
    my $found_match = 0;
    foreach my $key (keys %$indexes) {
      my $unq = $indexes->{$key};
      my $cur_match = 1;
      
      for (my $i=0; $i < scalar(@toks); $i++) {
        my $colname = $toks[$i];
        if ($unq->{$i+1} ne $colname) {
          $cur_match = 0;
        }
      }
      if ($cur_match == 1) {
        $found_match = 1;
      }
    }    
    if (!$found_match) {
      print "Creating UNIQUE $unique_keys\n";
      my $unique_cmd = qq^ALTER TABLE $table_name ADD UNIQUE ($unique_keys)^;
      $dbc->do($unique_cmd);
    } else {
      #print "  NOT creating unique $unique_keys\n";
    }
  }
  
  if ($extra_keys) {
    my @toks = split( ',', $extra_keys );
    
    for (my $i=0; $i < scalar(@toks); $i++) {
      my $colname = $toks[$i];
      my $match = 0;
      foreach my $key (keys %$indexes) {
        my $unq = $indexes->{$key};
        if ($unq->{1} eq $colname) {
          $match = 1;
        }
      }
      if ($match == 0) {
        print "Creating KEY $colname\n";
        my $key_cmd = qq^ALTER TABLE $table_name ADD KEY (`$colname`)^;
        $dbc->do($key_cmd);
      }
    }  
  }
}

sub output_rows_to_file {
  my $self = shift;
  my $rows_arrayref = shift;
  my $file = shift;

  my @rows = @{$rows_arrayref};

  my $first = $rows[0];

  my @keys = sort keys %$first;

  my $header = join("\t",@keys)."\n";

  open(OUT,">$file");
  print OUT $header;
  foreach my $row (@rows) {
    my @values;
    foreach my $key (@keys) {
      push @values, $row->{$key};
    }
    my $line = join("\t",@values)."\n";
    print OUT $line;
    #print "OUTPUT: $line";
  }
  close(OUT);
}

sub store_params_in_table {
  my $self       = shift;
  my $dbc        = shift;
  my $table_name = shift;
  my $params     = shift;

  my @fields;
  my $cache_key = '_fields_arrayref_' . $table_name;
  if ( !defined $self->{$cache_key} ) {
    my $sth = $dbc->prepare("SHOW FIELDS FROM $table_name");
    $sth->execute;
    my $hashref = $sth->fetchall_hashref( ['Field'] );
    $sth->finish;
    @fields = keys %$hashref;

    $self->{$cache_key} = \@fields;
  } else {
    @fields = @{ $self->{$cache_key} };
  }

  my $fields_string = '(`' . join( '`,`', @fields ) . '`)';
  my $question_marks_string = '(' . join( ',', ('?') x @fields ) . ')';  # Don't ask. It just works.
  my $sth2 =
    $dbc->prepare(qq^REPLACE INTO $table_name $fields_string VALUES $question_marks_string^);

  # Prepare an array of values.
  my @values = map {
    if ( defined $params->{$_} && $params->{$_} ne '' ) {
      $params->{$_};
    } elsif ( $_ eq 'parameter_set_id' ) {
      0;
    } else {
      undef;
    }
  } @fields;
  $sth2->execute(@values);
  $sth2->finish;
}

sub store_tag {
  my $self  = shift;
  my $tag   = shift;
  my $value = shift;

  my $node_id = $self->param('node_id') || 1;
  my $tree = $self->pta->fetch_node_by_node_id($node_id);

  print "  -> Storing tag [$tag] => [$value]\n" if ($self->debug);
  $tree->store_tag( $tag, $value );
}

sub base {
  my $self = shift;

  return Bio::Greg::EslrUtils->baseDirectory;
}

sub get_output_folder {
  my $self        = shift;

  # Double-check the meta table for a stored value.
  my $db_value = $self->get_meta('output_folder');
  if (defined $db_value) {
    $self->param('output_folder',$db_value);
  }

  if ( defined $self->param('output_folder') ) {
    my $folder = $self->param('output_folder');

    #print "Output folder: $folder\n" if ($self->debug);
    if ( !-e $folder ) {
      print
        "Warning: Output folder was already stored in the database, but the folder did not exist.\n";
      print " -> Creating folder: $folder\n";
      make_path( $folder );
    }
    return $folder;
  }

  my $output_base = Bio::Greg::EslrUtils->scratchDirectory.'/'.$self->hive_dba->dbc->dbname;

  #print "output_base: $output_base\n";

  my $date_string = strftime( "%Y-%m-%d", localtime );
  my $i = 0;

  my $filename;
  do {
    $i++;
    $filename = sprintf( "%s/%s_%.2d", $output_base, $date_string, $i );
  } while ( -e $filename );

  $self->store_meta( { output_folder => $filename } );

  make_path( $filename );
  $self->param( 'output_folder', $filename );
  return $filename;
}

sub save_aln {
  my $self   = shift;
  my $aln    = shift;
  my $params = shift;

  $params->{extension} = 'fasta';

  my $file_obj = $self->save_file($params);
  
  my $full_file = $file_obj->{full_file};

  Bio::EnsEMBL::Compara::AlignUtils->to_file( $aln, $full_file );

  return $file_obj;
}

sub save_file {
  my $self = shift;
  my $params = shift;

  my $hash_subfolders = 1;

  my $filename = "[unnamed_file]";
  $filename = $params->{filename} if ( defined $params->{filename} );

  my $id       = $filename;
  $id       = $params->{id}       if ( defined $params->{id} );

  my $subfolder = '';
  $subfolder = $params->{subfolder} if ( defined $params->{subfolder} );

  my $output_base = $self->get_output_folder;

  my $extension = $params->{extension} || 'txt';

  my $full_dir;
  my $full_file;
  my $rel_file;
  my $hash_folder = '';
  if ($hash_subfolders) {
    my $lo = 0;
    my $hi = 0;
    my (@md5) = md5_hex($id) =~ /\G(..)/g;
    my $hash_subfolder = join( '/', @md5[ $lo .. $hi ] );
    $hash_folder = $hash_subfolder;

    $full_dir = "$output_base/$subfolder/$hash_subfolder";
    $full_file = "$full_dir/$filename.$extension";
    $rel_file = "$hash_subfolder/$filename.$extension";
  } else {
    $full_dir = "$output_base/$subfolder";
    $full_file = "$full_dir/$filename.$extension";
    $rel_file = "$filename.$extension";
  }

  # If the file contains a directory, we need to make sure that's also created.
  my ($f,$d) = fileparse($full_file);
  make_path($d);

  return {
    full_file => $full_file,
    rel_file => $rel_file,
    full_dir => $full_dir,
    hash_folder => $hash_folder
  };
}

sub output_hash_tsv {
  my $self = shift;
  my $arry = shift;
  my $file = shift;

  my @rows = @{$arry};

  my $columns;
  map {
    foreach my $key (keys %$_) {
      $columns->{$key} = 1;
    }
  } @rows;

  open(OUT, ">$file");
  my @sorted_keys = sort( keys( %$columns ) );
  print OUT join("\t", @sorted_keys)."\n";
  #print join("\t", @sorted_keys)."\n";
  foreach my $row (@rows) {
    my @values = map { 
      my $v = $row->{$_};
      $v = 'NA' if (!defined $v);
      $v = 'NA' if ($v eq '');
      $v;
    } @sorted_keys;
    print OUT join("\t", @values)."\n";
    #print join("\t", @values)."\n";
  }
  close(OUT);
}

sub get_hashed_file {
  my $self      = shift;
  my $subfolder = shift;
  my $filename  = shift;

  my $output_base = $self->get_output_folder;

  my $lo = 0;
  my $hi = 1;
  my (@md5) = md5_hex($filename) =~ /\G(..)/g;
  my $hash_subfolder = join( '/', @md5[ $lo .. $hi ] );

  my $full_dir = "$output_base/$subfolder/$hash_subfolder";
  make_path( $full_dir );

  my $full_file = "$full_dir/$filename";
  print "$full_file\n" if ($self->debug);
  return $full_file;
}

sub get_r_dbc_string {
  my $self = shift;
  my $dbname = shift || $self->compara_dba->dbc->dbname;

  my $dbc  = $self->compara_dba->dbc;
  my $u    = $dbc->username;
  my $p    = $dbc->password;
  my $h    = $dbc->host;
  my $port = $dbc->port;

  return qq^
dbname='$dbname'
host = '$h'
port = $port
user = '$u'
password = '$p'
^;
}

sub cleanup_temp {
  my $self = shift;

  $self->worker->cleanup_worker_process_temp_directory;
  delete $self->worker->{'_tmp_dir'};
  $self->worker_temp_directory;
}

sub get_taxon_ids_from_keepers_list {
  my $self        = shift;
  my $list_string = shift;

  return Bio::EnsEMBL::Compara::ComparaUtils->get_taxon_ids_from_keepers_list( $self->compara_dba,
    $list_string );
}

sub DESTROY {
  my $self = shift;

  if (!$self->within_hive) {
    my $tmp = $self->worker_temp_directory;
    rmtree($tmp);
    #$self->cleanup_temp;
  }

  delete $self->{_param_hash};
  $self->SUPER::DESTROY if $self->can("SUPER::DESTROY");
}

sub thw {
  my $self = shift;
  my $file = shift;

  open(IN, $file);
  my @lines = <IN>;
  close(IN);
  my ($obj) = thaw(join('',@lines));
  return $obj;
}

sub frz {
  my $self = shift;
  my $file = shift;
  my $obj = shift;

  open(OUT,">$file");
  print OUT freeze($obj);
  close(OUT);
}

1;
