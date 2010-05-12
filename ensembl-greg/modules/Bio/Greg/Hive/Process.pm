package Bio::Greg::Hive::Process;

use strict;
use Bio::EnsEMBL::Utils::Argument;
use Bio::EnsEMBL::Utils::Exception;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Hive::Utils 'stringify';  # import 'stringify()'

use File::Path;

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

sub _get_aln {
  my $self = shift;
  my $params = shift;
  my $cdna = shift;

  $params = $self->params unless (defined $params);
  my $tree = $self->get_tree;

  my $aa_aln = $tree->get_SimpleAlign();
  my $cdna_aln = $tree->get_SimpleAlign(-cdna => 1, -hide_positions => 1);
  my $aln = Bio::EnsEMBL::Compara::ComparaUtils->fetch_masked_alignment($aa_aln,$cdna_aln,$tree,$params,$cdna);

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
    $aln = $aln->slice($lo,$hi);
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
  if ($good_seq_count < 2) {
    return -1;
  }
  return 1;
}


sub compara_dba {
  my $self = shift;

  if (!defined $self->param('compara_dba')) {
    my $compara_dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new( -DBCONN => $self->db->dbc );
    $self->{_compara_dba} = $compara_dba;
  }
  return $self->{_compara_dba};
}

sub db_handle {
  my $self = shift;

  my $dba = $self->compara_dba();
  return $dba->dbc->db_handle();
}

sub dbc {
  my $self = shift;

  return $self->compara_dba->dbc;
}

sub data_id {
  my $self = shift;

  my $id = $self->param('node_id');
  $id = $self->param('data_id') if (defined $self->param('data_id'));
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

sub parameter_set_id {
  my $self = shift;

  my $id = 0;
  $id = $self->param('parameter_set_id') if (defined $self->param('parameter_set_id'));
  return $id;
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
    $self->db_handle->do("DELETE from meta where meta_key='$key';");
    $self->db_handle->do("REPLACE into meta (meta_key,meta_value) VALUES ('$key','$value');");
  }
}

sub new_data_id {
  my $self = shift;
  my $params = shift;

  my $sth = $self->db_handle->prepare("SELECT value from protein_tree_tag WHERE node_id=0 AND tag='data_id_counter';");
  $sth->execute;
  my @row = $sth->fetchrow_array;
  $sth->finish;
  my $data_id=1;
  if (@row) {
    $data_id = $row[0];
  }

  $data_id++;
  $self->db_handle->do("REPLACE into protein_tree_tag (node_id,tag,value) VALUES (0,'data_id_counter',$data_id);");
  $self->add_breadcrumb($params,$data_id);

  $self->param('data_id',$data_id);
  #$params->{data_id} = $data_id;

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
#  $self->hash_print($self->params);

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

  my $dba = $self->compara_dba;
  my $dbh = $self->db_handle;

  throw "Undefined parameter_set_id!" unless (defined $param_set_id);

  my $params;

  my $cmd = qq^SELECT parameter_value FROM parameter_set WHERE parameter_set_id=$param_set_id AND parameter_name="params";  ^;
  my $sth = $dbh->prepare($cmd);
  $sth->execute();
  my @row;
  while (@row = $sth->fetchrow_array) {
    $params = eval($row[0]);
  }
  $sth->finish;

  my $cmd = qq^SELECT parameter_value FROM parameter_set WHERE parameter_set_id=$param_set_id AND parameter_name="name";  ^;
  $sth = $dbh->prepare($cmd);
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

  my $dbh = $self->db_handle;

  # First, create the table with node_id and parameter_set columns.
  eval {
    $dbh->do(
      qq^
	     CREATE TABLE $table_name (
				       data_id INT(10) NOT NULL,
				       parameter_set_id TINYINT(3) NOT NULL DEFAULT 0
				       )
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
      print "Creating key $key\n";
      my $type = $params->{$key};

      my $type_map = {
        'tinyint' => 'TINYINT(3)',
        'int'    => 'INT',
        'string' => 'TEXT',
        'float'  => 'FLOAT'
      };
      $type = $type_map->{$type};

      my $create_cmd = qq^ALTER TABLE $table_name ADD COLUMN `$key` $type^;
      $dbh->do($create_cmd);
    }

    my $unique_cmd = ""; #qq^ALTER TABLE $table_name ADD UNIQUE (data_id)^;
    if ($unique_keys) {
      $unique_cmd = qq^ALTER TABLE $table_name ADD UNIQUE ($unique_keys)^;
      $dbh->do($unique_cmd);
    }

    if ($extra_keys) {
      my @keys = split(',',$extra_keys);
      foreach my $key (@keys) {
	my $key_cmd = qq^ALTER TABLE $table_name ADD KEY (`$key`)^;
	$dbh->do($key_cmd);
      }
    }
  };

}

sub store_params_in_table {
  my $self       = shift;
  my $dba        = shift;
  my $table_name = shift;
  my $params     = shift;

  my $dbh = $self->db_handle;

  my @fields;
  my $cache_key = '_fields_arrayref_'.$table_name;
  if ( !defined $self->{$cache_key} ) {
    my $sth = $dbh->prepare("SHOW FIELDS FROM $table_name");
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
    $dbh->prepare(qq^REPLACE INTO $table_name $fields_string VALUES $question_marks_string^);

  # Prepare an array of values.
  my @values = map {
    if ( defined $params->{$_} && $params->{$_} ne '' ) {
      $params->{$_};
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

  my $node_id = $self->param('node_id');
  my $tree = $self->compara_dba->get_ProteinTreeAdaptor->fetch_node_by_node_id($node_id);

  if ($self->param('breadcrumb')) {
    $tag = $self->param('breadcrumb') . "." . $tag;
  }

  print "  -> Storing tag [$tag] => [$value]\n";
  $tree->store_tag($tag,$value);
}

sub get_output_folder {
  my $self = shift;
  my $output_base = shift;

  die("No output base folder specified!") unless ($output_base);
  
  if (defined $self->param('output_folder')) {
    return $self->param('output_folder');
  }
  
  my $date_string = strftime("%Y-%m-%d",localtime);
  my $i=0;

  my $filename;
  do {
    $i++;
    $filename = sprintf("%s/%s/%s_%.2d",$output_base,$date_string,$date_string,$i);
  } while (-e $filename);
  
  print "Output folder: $filename\n";
  $self->param('output_folder',$filename);
  
  # We'll store this output folder in the meta table, so it will be re-used if this module is run again w/ the same database.
  $self->store_meta({output_folder => $filename});
  mkpath([$filename]);

sub cleanup_temp {
  my $self = shift;

  $self->worker->cleanup_worker_process_temp_directory;
  delete $self->worker->{'_tmp_dir'};
  $self->worker_temp_directory;
}

sub DESTROY {
    my $self = shift;

    $self->cleanup_temp;

    delete $self->{_param_hash};
    $self->SUPER::DESTROY if $self->can("SUPER::DESTROY");
}


1;
