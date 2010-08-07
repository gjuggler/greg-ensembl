package Bio::Greg::Hive::Process;

use strict;
use Bio::EnsEMBL::Utils::Argument;
use Bio::EnsEMBL::Utils::Exception;

use Data::UUID;
use Digest::MD5 qw(md5_hex);

use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Hive::Utils 'stringify';  # import 'stringify()'

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

use File::Path;
use POSIX qw(strftime mktime);

use base ('Bio::EnsEMBL::Hive::ProcessWithParams');

sub get_tree {
  my $self = shift;
  my $params = shift;

  $params = $self->params unless (defined $params);
  
  my $tree = Bio::EnsEMBL::Compara::ComparaUtils->get_tree_for_comparative_analysis($self->compara_dba,$params);
  return $tree;
}

sub get_aln {
  my $self = shift;
  my $params = shift;

  return $self->_get_aln($params,0);
}

sub get_cdna_aln {
  my $self = shift;
  my $params = shift;

  return $self->_get_aln($params,1);
}

sub get_orig_aln {
  my $self = shift;

  $self->get_aln if (!defined $self->param('full_aln_aa'));
  return $self->param('full_aln_aa');
}

sub get_orig_aln_position {
  my $self = shift;
  my $aln = shift;
  my $aln_position = shift;

  my $orig_aln = $self->get_orig_aln;
  return Bio::EnsEMBL::Compara::AlignUtils->map_alignment_position($aln,$aln_position,$orig_aln);
}

sub _get_aln {
  my $self = shift;
  my $params = shift;
  my $cdna = shift;

  $params = $self->params unless (defined $params);

  $params->{remove_blank_columns} = 1 unless (defined $params->{remove_blank_columns});

  my $tree;
  if (defined $params->{tree}) {
    $tree = $params->{tree};
  } else {
    $tree = $self->get_tree($params);
  }

  my $aa_aln = $tree->get_SimpleAlign();
  my $cdna_aln = $tree->get_SimpleAlign(-cdna => 1, -hide_positions => 1);
  my $aln = Bio::EnsEMBL::Compara::ComparaUtils->fetch_masked_alignment($aa_aln,$cdna_aln,$tree,$params,$cdna);

  # Remove blank columns here. Don't forget to store the 'full' alignment as a parameter.
  if ($params->{remove_blank_columns}) {
    if ($cdna) {
      $self->param('full_aln_cdna',$aln);
      ($aln) = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns_in_threes($aln);
    } else {
      $self->param('full_aln_aa',$aln);
      ($aln) = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns($aln);
    }
  }

  # Here's where we'll split up alignments by slice.
  if ($params->{alignment_slice}) {
    my $slice_string = $params->{alignment_slice};
    my ($lo,$hi) = split("-",$slice_string);
    #$lo = 1;
    #$hi = 10;
    if ($cdna) {
      $lo = ($lo) * 3-2;
      $hi = ($hi) * 3;
    }
    print "$lo $hi\n";
    print "Length: ".$aln->length."\n";
    $aln = $aln->slice($lo,$hi,1);
    print "Length: ".$aln->length."\n";
  }
 
  return $aln;
}

sub fail {
  my $self = shift;
  my $message = shift;

  $self->input_job->update_status('FAILED');
  warn($message);
}

sub fail_and_die {
  my $self = shift;
  my $message = shift;

  $self->input_job->update_status('FAILED');
  throw($message);
}

sub check_tree_aln {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  if (!defined $tree) {
    $tree = $self->get_tree;
  }
  # Look at the tree size.
  if (scalar $tree->leaves <= 1) {
    return -1;
  }

  if (!defined $aln) {
    $aln = $self->get_aln;
  }
  # Look at the alignment size.
  if (scalar $aln->each_seq <= 1 || $aln->length == -1) {
    return -1;
  }

  # More complicated stuff -- look for alignments where only one sequence isn't entirely X'ed out.
  my $good_seq_count = 0;
  foreach my $seq ($aln->each_seq) {
    my $seq_str = $seq->seq;
    $good_seq_count++ if ($seq_str =~ m/[^NX-]/i);
  }
  print "Good seq count: $good_seq_count\n";
  if ($good_seq_count < 2) {
    return -1;
  }
  return 1;
}

sub pretty_print {
  my $self = shift;
  my $aln = shift;
  my $params = shift || {};

  Bio::EnsEMBL::Compara::AlignUtils->pretty_print($aln,$params);
}


sub compara_dba {
  my $self = shift;

  if (!defined $self->{_compara_dba}) {
    print " >> Getting new DBA!!!!\n";
    my $compara_dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new( -DBCONN => $self->db->dbc );
    eval {
      $compara_dba->do("select * from member limit 1;");
    };
    if ($@) {
      print " >> No compara in hive DB -- falling back to ens-livemirror!!\n";
      Bio::EnsEMBL::Registry->load_registry_from_multiple_dbs({
                                                               -host => 'ens-livemirror',
                                                               -user => 'ensro',
#                                                               -verbose => 1
                                                              });
      $compara_dba = Bio::EnsEMBL::Registry->get_DBAdaptor('multi','compara');
#      $compara_dba->disconnect_when_inactive(1);
    }
    $self->{_compara_dba} = $compara_dba;
  }
  return $self->{_compara_dba};
}

sub hive_dba {
  my $self = shift;

  if (!defined $self->{_hive_dba}) {
    print " >>>> Getting new Hive DBA!!!!\n";
    my $hive_dba = Bio::EnsEMBL::Hive::DBSQL::DBAdaptor->new( -DBCONN => $self->db->dbc );
    $self->{_hive_dba} = $hive_dba;
  }
  return $self->{_hive_dba};
}

sub pta {
  my $self = shift;

#  return $self->compara_dba->get_ProteinTreeAdaptor;

  if (!defined $self->param('_pta')) {
    print " >>>> Getting new PTA!!!!\n";
    my $compara_dba = $self->compara_dba;
    my $pta = $compara_dba->get_ProteinTreeAdaptor;
    $self->param('_pta',$pta);
  }
  return $self->param('_pta');
}

sub mba {
  my $self = shift;

  #return $self->compara_dba->get_MemberAdaptor;

  if (!defined $self->param('_mba')) {
    print " >>>> Getting new MBA!!!!\n";
    my $compara_dba = $self->compara_dba;
    my $pta = $compara_dba->get_MemberAdaptor;
    $self->param('_mba',$pta);
  }
  return $self->param('_mba');
}

sub db_handle {
  my $self = shift;
  my $force = shift;

  if (!defined $self->param('_compara_dbh') || $force == 1) {
    print " >>>> Getting new Compara DBH!!!!\n";
    my $dbc = $self->dbc;
    my $dbh = $dbc->db_handle;
    $self->param('_compara_dbh',$dbh);
    print $dbh."\n";
  }
  return $self->param('_compara_dbh');
}

sub dbc {
  my $self = shift;

  return $self->db->dbc;

  if (!defined $self->param('_compara_dbc')) {
    print " >>>> Getting new Compara DBC!!!!\n";
    my $compara_dba = $self->compara_dba;
    my $dbc = $compara_dba->dbc;
    print $dbc."\n";
    # It's apparently important to turn off the inactive disconnect here, since we'll
    # be sharing this DBC throughout the lifetime of this Process.
    # TODO: Think of how to handle the case when a Process wants to let a connection
    # run idle...
    $dbc->disconnect_when_inactive(1);
    $self->param('_compara_dbc',$dbc);
  }
  return $self->param('_compara_dbc');
}

sub node_id {
  my $self = shift;

  my $node_id = 0;
  $node_id = $self->param('node_id') if (defined $self->param('node_id'));

  return $node_id;
}

sub data_id {
  my $self = shift;

  my $data_id = 0;
  $data_id = $self->param('data_id') if (defined $self->param('data_id'));
  return $data_id;
}

sub parameter_set_id {
  my $self = shift;

  my $id = 0;
  $id = $self->param('parameter_set_id') if (defined $self->param('parameter_set_id'));
  return $id;
}

sub worker_temp_directory {
  my $self = shift;

  my $wtd = $self->SUPER::worker_temp_directory;

  chmod 0777, $wtd;

  #print "ANALYSIS ID:". $self->worker->analysis->dbID."\n";
  #$wtd = $wtd . $self->data_id."/".$self->worker->analysis->dbID."/";
  #mkpath([$wtd]);
  #print "TEMP: $wtd\n";
  return $wtd;
}


sub breadcrumb_param {
  my $self = shift;
  my $param = shift;

  my $breadcrumb = $self->param('breadcrumb') || '';
  #print "bc: $breadcrumb\n";

  if ($breadcrumb ne '') {
    my @crumbs = split("\\.",$breadcrumb);
    my $key_prefix = shift @crumbs;
    foreach my $crumb (@crumbs) {
      $key_prefix .= '.' . $crumb;
      
      #print "Searching for param $param with breadcrumb $key_prefix...\n";
      my $key = $key_prefix.".".$param;
      my $value = $self->{_param_hash}->{$key};
      if (defined $value) {
	#print "Got it! $value\n";
	return $value;
      }
    }
  }
  return undef;
}

sub add_breadcrumb {
  my $self = shift;
  my $params = shift;
  my $crumb = shift;

  my $breadcrumb = $params->{breadcrumb} || 'bcrmb';

  $breadcrumb .= '.' unless ($breadcrumb eq '');
  $breadcrumb .= $crumb;
  
  $params->{breadcrumb} = $breadcrumb;
}

sub store_meta {
  my $self = shift;
  my $params = shift;

  foreach my $key (keys %$params) {
    my $value = $params->{$key};
    $self->dbc->do("DELETE from meta where meta_key='$key';");
    $self->dbc->do("REPLACE into meta (meta_key,meta_value) VALUES ('$key','$value');");
  }
}

sub get_parameter_sets {
  my $self = shift;
  my $dbname = shift || $self->dbc->dbname;

  my $sth = $self->dbc->prepare("SELECT * from ${dbname}.parameter_set;");
  $sth->execute;
  
  my $hash;
  while (my $obj = $sth->fetchrow_hashref) {
    my $id = $obj->{parameter_set_id};
    $hash->{$id} = {} unless (defined $hash->{$id});
    my $param_name = $obj->{parameter_name};
    my $param_value = $obj->{parameter_value};
    $hash->{$id}->{$param_name} = $param_value;
    $hash->{$id}->{parameter_set_id} = $id;
  }

  my @sets;
  foreach my $pset (1 .. scalar(keys %$hash)) {
    push @sets,$hash->{$pset};
  }
  return @sets;
}

sub new_data_id {
  my $self = shift;
  my $params = shift;

  my $ug = new Data::UUID;
  my $uuid = $ug->create();
  

  my $dbc = $self->dbc;
  $dbc->do("LOCK TABLES protein_tree_tag WRITE");

  my $sth = $self->dbc->prepare("SELECT value from protein_tree_tag WHERE node_id=0 AND tag='data_id_counter';");
  $sth->execute;
  my @row = $sth->fetchrow_array;
  $sth->finish;
  my $data_id=1;
  if (@row) {
    $data_id = $row[0];
  }

  $data_id++;
  $self->dbc->do("REPLACE into protein_tree_tag (node_id,tag,value) VALUES (0,'data_id_counter',$data_id);");
  $self->add_breadcrumb($params,$data_id);

  $self->param('data_id',$data_id);
  $params->{data_id} = $data_id;

  $dbc->do("UNLOCK TABLES");

  return $data_id;
}

sub param {
  my $self = shift;
  my $param = shift;

  if (@_) {
    $self->SUPER::param($param,shift @_);
  }

  my $param_value = $self->SUPER::param($param);

  #print "$param value: $param_value\n" if ($param =~ m/slr/);
  if (!defined $param_value && $param ne 'breadcrumb') {
    $param_value = $self->breadcrumb_param($param);
  }
  return $param_value;
}

sub params {
  my $self = shift;  
  $self->param_init();

  # Make a copy!
  my $param_hash = $self->{_param_hash};
  my $new_params = {};
  foreach my $key (keys %$param_hash) {
    $new_params->{$key} = $param_hash->{$key};
  }
  return $new_params;
}

sub load_all_params {
  my $self = shift;
  my $default_params = shift || {};

  $self->param_init(%$default_params);

  my $node_id  = $self->param('node_id');
  my $parameter_set_id  = $self->param('parameter_set_id');

  my $tree_tag_params = {};
  if (defined $node_id) {
    $tree_tag_params = $self->get_params_from_tree_tags( $self->compara_dba, $node_id ) || {};
  }

  my $param_set_params = {};
  if (defined $parameter_set_id) {
    $param_set_params = $self->get_params_from_param_set($parameter_set_id) || {};
  }
  
  my $old_param_hash = $self->{'_param_hash'};
  $self->{'_param_hash'} = { %$old_param_hash, %$param_set_params, %$tree_tag_params };

  if (!defined $self->param('data_id')) {
    $self->param('data_id',$node_id);
  }
  if (!defined $self->param('parameter_set_id')) {
    $self->param('parameter_set_id',0);
  }

  print "Bio::Greg::Hive::Process.pm - load all params\n";
  $self->hash_print($self->params);

}

sub get_params_from_tree_tags {
  my $self    = shift;
  my $dba     = shift;
  my $node_id = shift;

  my $pta  = $dba->get_ProteinTreeAdaptor;
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
  my $HANDLE = shift;

  select($HANDLE) if (defined $HANDLE);

  print  "  {\n";
  foreach my $key ( sort keys %{$hashref} ) {
    printf "    %-40.40s => %-40s\n", $key, $hashref->{$key};
  }
  print "  }\n";

  select(STDOUT);
}

sub clone {
  my $self = shift;
  my $params = shift;

  return $self->replace($params,{});
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
  my $self = shift;
  my $param_set_id = shift;

  my $dba = $self->hive_dba;
  my $dbc = $self->dbc;

  throw "Undefined parameter_set_id!" unless (defined $param_set_id);

  my $params;

  my $cmd = qq^SELECT parameter_value FROM parameter_set WHERE parameter_set_id=$param_set_id AND parameter_name="params";  ^;
  my $sth = $dbc->prepare($cmd);
  $sth->execute();
  my @row;
  while (@row = $sth->fetchrow_array) {
    $params = eval($row[0]);
  }
  $sth->finish;

  my $cmd = qq^SELECT parameter_value FROM parameter_set WHERE parameter_set_id=$param_set_id AND parameter_name="name";  ^;
  $sth = $dbc->prepare($cmd);
  $sth->execute();
  my @row;
  while (@row = $sth->fetchrow_array) {
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

  my $dbc = $self->dbc;

  # First, create the table with data_id and parameter_set columns.
  eval {
    $dbc->do(
      qq^
	     CREATE TABLE $table_name (
				       data_id INT(10) NOT NULL,
				       parameter_set_id TINYINT(3) NOT NULL DEFAULT 0
				       ) engine=InnoDB;
	     ^
    );
  };
  if ($@) {
    print "TABLE $table_name EXISTS!!\n";
    return;
  }
  print "Creating new table $table_name ...\n";

  my $unique_keys = delete $params->{'unique_keys'};
  my $extra_keys   = delete $params->{'extra_keys'};
  eval {

    foreach my $key ( sort keys %$params ) {
      next if ($key eq 'data_id' || $key eq 'parameter_set_id');
      print "Creating column $key\n";
      my $type = $params->{$key};

      my $not_null = 0;
      if ($type =~ m/not null/gi) {
	$type =~ s/\s*not null\s*//gi;
	$not_null = 1;
      }

      my $type_map = {
	'uuid' => 'BINARY(16)',
	'timestamp' => 'DATETIME',
        'tinyint' => 'TINYINT',
        'smallint'    => 'SMALLINT',
        'int'    => 'INT',
        'string' => 'TEXT',
        'char8' => 'CHAR(8)',
        'char16' => 'CHAR(16)',
        'char32' => 'CHAR(32)',
        'float'  => 'FLOAT'
      };
      $type = $type_map->{$type};
      
      $type = $type .'NOT NULL ' if ($not_null == 1);

      my $create_cmd = qq^ALTER TABLE $table_name ADD COLUMN `$key` $type^;
      $dbc->do($create_cmd);
    }

    my $unique_cmd = ""; #qq^ALTER TABLE $table_name ADD UNIQUE (data_id)^;
    if ($unique_keys) {
      print "Creating UNIQUE $unique_keys\n";
      $unique_cmd = qq^ALTER TABLE $table_name ADD UNIQUE ($unique_keys)^;
      $dbc->do($unique_cmd);
    }

    if ($extra_keys) {
      my @keys = split(',',$extra_keys);
      foreach my $key (@keys) {
	print "Creating KEY $key\n";
	my $key_cmd = qq^ALTER TABLE $table_name ADD KEY (`$key`)^;
	$dbc->do($key_cmd);
      }
    }
  };

}

sub store_params_in_table {
  my $self       = shift;
  my $dbc        = shift;
  my $table_name = shift;
  my $params     = shift;

  my @fields;
  my $cache_key = '_fields_arrayref_'.$table_name;
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
    } elsif ($_ eq 'parameter_set_id') {
      0;
} else {
      undef;
    }
  } @fields;
  $sth2->execute(@values);
  $sth2->finish;
}


sub store_tag {
  my $self = shift;
  my $tag = shift;
  my $value = shift;

  my $node_id = $self->param('node_id') || 1;
  my $tree = $self->pta->fetch_node_by_node_id($node_id);

  if ($self->param('breadcrumb')) {
    $tag = $self->param('breadcrumb') . "." . $tag;
  }

  print "  -> Storing tag [$tag] => [$value]\n";
  $tree->store_tag($tag,$value);
}

sub base {
  my $self = shift;

  return Bio::Greg::EslrUtils->baseDirectory;
}

sub get_output_folder {
  my $self = shift;
  my $output_base = shift;

  if (defined $self->param('output_folder')) {
    my $folder = $self->param('output_folder');
    #print "Output folder: $folder\n";
    if (!-e $folder) {
      print "Warning: Output folder was already stored in the database, but the folder did not exist.\n";
      print " -> Creating folder: $folder\n";
      mkpath([$folder]);
    }
    return $folder;
  }

  die("No output base folder specified!") unless ($output_base);
    
  my $date_string = strftime("%Y-%m-%d",localtime);
  my $i=0;

  my $filename;
  do {
    $i++;
    $filename = sprintf("%s/%s/%s_%.2d",$output_base,$date_string,$date_string,$i);
  } while (-e $filename);
  
  print "Output folder: $filename\n";
  mkpath([$filename]);
  $self->param('output_folder',$filename);
  
  # We'll store this output folder in the meta table, so it will be re-used if this module is run again w/ the same database.
  $self->store_meta({output_folder => $filename});
  mkpath([$filename]);

  return $filename;
}

sub save_aln {
  my $self = shift;
  my $aln = shift;
  my $params = shift;
  
  my @seqs = $aln->each_seq;
  my $first = @seqs[0];
  my $filename = $first->id;
  my $id = $filename;
  $filename = $params->{filename} if (defined $params->{filename});
  $id = $params->{id} if (defined $params->{id});

  my $subfolder = '';
  $subfolder = $params->{subfolder} if (defined $params->{subfolder});
  
  my $output_base = $self->get_output_folder;

  my $lo = 0;
  my $hi = 0;
  my (@md5) = md5_hex($id) =~ /\G(..)/g;
  my $hash_subfolder = join('/',@md5[$lo .. $hi]);

  my $full_dir = "$output_base/$subfolder/$hash_subfolder";
  mkpath([$full_dir]);

  my $full_file = "$full_dir/$filename.fasta";
  print $full_file."\n";

  my $rel_file = "$hash_subfolder/$filename.fasta";

  # TODO: Output the aln.
  Bio::EnsEMBL::Compara::AlignUtils->to_file($aln,$full_file);

  return $rel_file;
}

sub get_hashed_file {
  my $self = shift;
  my $subfolder = shift;
  my $filename = shift;

  my $output_base = $self->get_output_folder;

  my $lo = 0;
  my $hi = 1;
  my (@md5) = md5_hex($filename) =~ /\G(..)/g;
  my $hash_subfolder = join('/',@md5[$lo .. $hi]);
  
  my $full_dir = "$output_base/$subfolder/$hash_subfolder";
  mkpath([$full_dir]);

  my $full_file = "$full_dir/$filename";
  print "$full_file\n";
  return $full_file;
}


sub get_r_dbc_string {
  my $self = shift;
  my $dbname = shift || $self->compara_dba->dbc->dbname;

  my $dbc = $self->compara_dba->dbc;
  my $u = $dbc->username;
  my $p = $dbc->password;
  my $h = $dbc->host;
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
  my $self = shift;
  my $list_string = shift;

  return Bio::EnsEMBL::Compara::ComparaUtils->get_taxon_ids_from_keepers_list($self->compara_dba,$list_string);
}

sub DESTROY {
    my $self = shift;

    $self->cleanup_temp;

    delete $self->{_param_hash};
    $self->SUPER::DESTROY if $self->can("SUPER::DESTROY");
}


1;
