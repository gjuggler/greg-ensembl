package Bio::Greg::ProcessUtils;

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;

use Bio::EnsEMBL::Hive;

our @ISA = qw();

sub load_params_from_tree_tags {
  my $self = shift;
  my $dba = shift;
  my $node_id = shift;

  my $pta = $dba->get_ProteinTreeAdaptor;
  my $tree = $pta->fetch_node_by_node_id($node_id);

  my $tags = $tree->get_tagvalue_hash;

  my @param_objs;
  my $simple_tags;
  foreach my $tag (keys %$tags) {
    if ($tag =~ /params/i) {
      push @param_objs,$self->load_params_from_tag($tree,$tag);
    } else {
      $simple_tags->{$tag} = $tags->{$tag};
    }
  }
  push @param_objs,$simple_tags;

  return $self->replace_params(@param_objs);
}

sub load_params_from_tag {
  my $self = shift;
  my $tree = shift;
  my $tag = shift;

  my $tag_string = $tree->get_tagvalue($tag);
  return $self->string_to_hash($tag_string);
}

# Eval a hashref from a string.
sub string_to_hash {
  my $self = shift;
  my $string = shift;

  my $hash = eval $string;
  return $hash if (defined $hash);

  return {};
}

sub hash_print {
  my $class = shift;
  my $hashref = shift;

  print "{\n";
  foreach my $key (sort keys %{$hashref}) {
    printf("    %-40.40s => %-40s\n",$key,$hashref->{$key});
  }
  print "}\n";
}

sub replace_params {
  my $self = shift;  
  my @params = @_;

  my $final_params = shift @params;
  while (@params) {
    my $p = shift @params;
    if (!ref $p) {
      $p = eval($p);
    }
    
    my $new_params = {};
    foreach my $key (keys %{$final_params}) {
      $new_params->{$key} = $final_params->{$key};
    }
    foreach my $key (keys %{$p}) {
      $new_params->{$key} = $p->{$key};
    }
    $final_params = $new_params;
  }

  return $final_params;
}

sub get_params {
  my $self = shift;
  my $param_string = shift;
  
  my $param_object;
  my $tmp_params = eval($param_string);  
  foreach my $key (keys %$tmp_params) {
    $param_object->{$key} = $tmp_params->{$key};
  }
  return $param_object;
}


sub create_table_from_params {
  my $self = shift;
  my $dba = shift;
  my $table_name = shift;
  my $params = shift;
  
  my $dbh = $dba->dbc->db_handle;

  # First, create the table with node_id and parameter_set columns.
  eval {
    $dbh->do(qq^
	     CREATE TABLE $table_name (
				       node_id INT(10),
				       parameter_set_id TINYINT(3)
				       )
	     ^);
  };
  if ($@) {
    print "TABLE $table_name EXISTS!!\n";
    return;
  }
  print "Creating new table $table_name ...\n";

  my $unique_keys = delete $params->{'unique_keys'};
  eval {
    
    foreach my $key (sort keys %$params) {
      my $type = $params->{$key};
      
      my $type_map = {
	'int' => 'INT',
	'string' => 'TEXT',
	'float' => 'FLOAT'
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
  };

}

sub store_params_in_table {
  my $self = shift;
  my $dba = shift;
  my $table_name = shift;
  my $params = shift;

  my $dbh = $dba->dbc->db_handle;

  my @fields;
  if (!defined $params->{_fields_arrayref}) {
    my $sth = $dbh->prepare("SHOW FIELDS FROM $table_name");
    $sth->execute;
    my $hashref = $sth->fetchall_hashref(['Field']);
    $sth->finish;
    @fields = keys %$hashref;
  
    $params->{_fields_arrayref} = \@fields;
  } else {
    @fields = @{$params->{_fields_arrayref}};
  }

  my $fields_string = '(`'.join('`,`',@fields).'`)';
  my $question_marks_string = '('. join(',',('?') x @fields) . ')';  # Don't ask. It just works.
  my $sth2 = $dbh->prepare(qq^REPLACE INTO $table_name $fields_string VALUES $question_marks_string^);
  
  # Prepare an array of values.
  my @values = map {
    if (defined $params->{$_} && $params->{$_} ne '') {
      $params->{$_};
    } else {
      undef;
    }
  } @fields;
  $sth2->execute(@values);
  $sth2->finish;
}

1;
