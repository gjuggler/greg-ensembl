package Bio::Greg::ProcessUtils;

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::ComparaUtils;

use Bio::EnsEMBL::Hive;

our @ISA = qw();

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

  # Check if the table already exists.
  my $sth = $dbh->table_info("%",'',$table_name);
  if ($sth->fetchrow_arrayref) {
    # table already exists!
    return;
  }

  # First, create the table with node_id and parameter_set columns.
  $dbh->do(qq^
	   CREATE TABLE $table_name (
	     node_id INT(10),
	     parameter_set_id TINYINT(3),
	     UNIQUE KEY (node_id,parameter_set_id)
	   )
	   ^);

  foreach my $key (keys %$params) {
    my $type = $params->{$key};

    my $type_map = {
      'int' => 'INT',
      'string' => 'TEXT',
      'float' => 'FLOAT'
    };
    $type = $type_map->{$type};

    my $create_cmd = qq^ALTER TABLE $table_name ADD COLUMN $key $type^;
    $dbh->do($create_cmd);
  }

}

sub store_params_in_table {
  my $self = shift;
  my $dba = shift;
  my $table_name = shift;
  my $params = shift;

  my $dbh = $dba->dbc->db_handle;

  my $sth = $dbh->prepare("SHOW FIELDS FROM $table_name");
  $sth->execute;
  my $hashref = $sth->fetchall_hashref(['Field']);
  my @fields = keys %$hashref;

  my $fields_string = '('.join(',',@fields).')';
  my $question_marks_string = '('. join(',',('?') x scalar(@fields)) . ')';  # Don't ask. It just works.
  
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
  
}

1;
