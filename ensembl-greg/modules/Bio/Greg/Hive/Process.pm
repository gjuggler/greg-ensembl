package Bio::Greg::Hive::Process;

use strict;
use Bio::EnsEMBL::Utils::Argument;
use Bio::EnsEMBL::Utils::Exception;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

use base ('Bio::EnsEMBL::Hive::ProcessWithParams');

sub get_tree {
  my $self = shift;
  my $params = shift;

  $params = $self->params unless (defined $params);
  
  my $tree = Bio::EnsEMBL::Compara::ComparaUtils->get_tree_for_comparative_analysis($self->compara_dba,$params);
  return $tree;
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
    $tree_tag_params = $self->get_params_from_tree_tags( $self->compara_dba, $node_id );
  }

  my $param_set_params = {};
  if (defined $parameter_set_id) {
    $param_set_params = $self->get_params_from_param_set($parameter_set_id);
  }
  
  my $old_param_hash = $self->{'_param_hash'};
  $self->{'_param_hash'} = { %$old_param_hash, %$param_set_params, %$tree_tag_params };

  print "Bio::Greg::Hive::Process.pm - load all params\n";
  $self->hash_print($self->{'_param_hash'});
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

  print "  {\n";
  foreach my $key ( sort keys %{$hashref} ) {
    printf( "    %-40.40s => %-40s\n", $key, $hashref->{$key} );
  }
  print "  }\n";
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
				       node_id INT(10) NOT NULL,
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
  my $extra_key   = delete $params->{'extra_key'};
  eval {

    foreach my $key ( sort keys %$params ) {
      my $type = $params->{$key};

      my $type_map = {
        'int'    => 'INT',
        'string' => 'TEXT',
        'float'  => 'FLOAT'
      };
      $type = $type_map->{$type};

      my $create_cmd = qq^ALTER TABLE $table_name ADD COLUMN `$key` $type^;
      $dbh->do($create_cmd);
    }

    my $unique_cmd = qq^ALTER TABLE $table_name ADD UNIQUE (node_id,parameter_set_id)^;
    if ($unique_keys) {
      $unique_cmd = qq^ALTER TABLE $table_name ADD UNIQUE ($unique_keys)^;
    }
    $dbh->do($unique_cmd);

    if ($extra_key) {
      my $key_cmd = qq^ALTER TABLE $table_name ADD KEY ($extra_key)^;
      $dbh->do($unique_cmd);
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

  print "  -> Storing tag $tag -> $value\n";
  $tree->store_tag($tag,$value);
}


1;
