package Bio::Greg::Hive::ComparaHiveLoaderUtils;

use warnings;
use strict;

use base ('Bio::Greg::Hive::HiveLoaderUtils');

our $param_set_counter = 1;

sub new {
  my ( $class, @args ) = @_;
  my $self = bless {}, $class;

  Bio::EnsEMBL::Compara::ComparaUtils->load_registry();

  return $self;
}

sub add_inputs_to_analysis {
  my $self                 = shift;
  my $analysis_name        = shift;
  my $job_hashref_arrayref = shift;

  my @input_ids = @$job_hashref_arrayref;
  foreach my $input_id (@input_ids) {
    $self->add_job_to_analysis( $analysis_name, $input_id );
  }
}

sub add_genes_to_analysis {
  my $self             = shift;
  my $analysis_name    = shift;
  my $gene_id_arrayref = shift;

  my $compara_dba = $self->compara_dba;
  my $mba = $compara_dba->get_MemberAdaptor;
  my $pta = $compara_dba->get_ProteinTreeAdaptor;

  my @node_ids = ();
  foreach my $gene_id (@$gene_id_arrayref) {
    chomp $gene_id;
    my $member;
    my $gene_member;

    print "Trying to add $gene_id...\n";

    if ($gene_id !~ m/ensp/gi) {
      print "Finding gene by external ID...\n";
      my $ext_member = $self->_find_member_by_external_id($gene_id);
      $member = $ext_member if ( defined $ext_member );
      $gene_member = $member->gene_member if ( defined $member );
    }

    if ( !defined $member ) {
      print "Finding gene by stable id...\n";
      $member = $mba->fetch_by_source_stable_id( undef, $gene_id );
      $gene_member = $member->gene_member if ( defined $member );
    }

    if ( !defined $member ) {
      print STDERR " -> $gene_id not found!\n";
      next;
    }

    my $tree = $pta->fetch_by_Member_root_id($member,0);
    $tree = $pta->fetch_by_gene_Member_root_id($member,0) if (!defined $tree);
    $tree = $pta->fetch_by_Member_root_id($gene_member,0) if (!defined $tree && defined $gene_member);
    $tree = $pta->fetch_by_gene_Member_root_id($gene_member,0) if (!defined $tree && defined $gene_member);
    if (!defined $tree) {
      print STDERR " -> $gene_id tree not found!\n";
      next;
    }
    $self->add_job_to_analysis( $analysis_name, {
      node_id => $tree->node_id,
      gene_id => $gene_id,
      flow_node_set => 'one2one'
                                } );
    $tree->release_tree;
  }
}

sub add_nodes_to_analysis {
  my $self             = shift;
  my $analysis_name    = shift;
  my $node_id_arrayref = shift;

  foreach my $node_id (@$node_id_arrayref) {
    $self->add_job_to_analysis( $analysis_name, { node_id => $node_id } );
  }
}

sub compara_dba {
  my $self = shift;
  my $version = shift;

  if (!defined $self->{_compara_dba}) {
    warn("No Compara version specified -- loading from livemirror!") if (!defined $version);
    Bio::EnsEMBL::Compara::ComparaUtils->load_registry;
    my $compara_dba = Bio::EnsEMBL::Registry->get_DBAdaptor( 'multi', 'compara' );
    $self->{_compara_dba} = $compara_dba;
  }
  return $self->{_compara_dba};
}

sub select_node_ids {
  my $self = shift;
  my $cmd  = shift;

  my $dba = $self->compara_dba;
  if ( !defined $cmd ) {
    $cmd = "SELECT node_id FROM protein_tree_node WHERE parent_id=1 AND root_id=1";
  }
  my $sth = $dba->dbc->prepare($cmd);
  $sth->execute();

  my $array_ref = $sth->fetchall_arrayref( [0] );
  my @node_ids = @{$array_ref};
  @node_ids =
    map { @{$_}[0] } @node_ids;    # Some weird mappings to unpack the numbers from the arrayrefs.
  $sth->finish;
  return @node_ids;
}

sub _find_member_by_external_id {
  my $self = shift;
  my $id   = shift;

  my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "gene" );
  my $member_adaptor = $self->compara_dba->get_MemberAdaptor;

  my @genes = ();
  push @genes, @{ $gene_adaptor->fetch_all_by_external_name($id) };
  if ( scalar @genes == 0 ) {
    my $stable_id_gene = $gene_adaptor->fetch_by_stable_id($id);
    push @genes, $stable_id_gene if ( defined $stable_id_gene );
  }

  foreach my $gene (@genes) {
    print $gene->stable_id ."\n";
    my $stable_id = $gene->stable_id;

    my $member = $member_adaptor->fetch_by_source_stable_id( undef, $stable_id );

    #print $member."\n";
    return $member if ( defined $member );
  }
  return undef;
}

sub add_parameter_set {
  my $self   = shift;
  my $params = shift;

  my $dbc = $self->dbc;

  my $parameter_set_id = $param_set_counter++;
  $params->{'parameter_set_id'} = $parameter_set_id;

  if ( exists $params->{'parameter_set_name'} ) {
    my $parameter_set_name = $params->{'parameter_set_name'};
    my $name_cmd =
      "REPLACE INTO parameter_set VALUES ('$parameter_set_id','name',\"$parameter_set_name\");";
    $dbc->do($name_cmd);
  }

  if ( exists $params->{'parameter_set_shortname'} ) {
    my $parameter_set_shortname = $params->{'parameter_set_shortname'} || '';
    my $shortname_cmd =
      "REPLACE INTO parameter_set VALUES ('$parameter_set_id','parameter_set_shortname',\"$parameter_set_shortname\");";
    $dbc->do($shortname_cmd);
  }

  my $param_string = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
  my $cmd = "REPLACE INTO parameter_set VALUES ('$parameter_set_id','params',\"$param_string\");";
  $dbc->do($cmd);
}

sub clean_compara_analysis_tables {
  my $self            = shift;
  my $dba             = $self->dba;
  my @truncate_tables = qw^
    sitewise_omega sitewise_tag sitewise_genome
    parameter_set
    node_set_member node_set
    go_terms
    stats_sites stats_genes
    sites genes
    ^;
  map {
    print "$_\n";
    eval { $dba->dbc->do("truncate table $_") };
  } @truncate_tables;
  print "Dropping sites...\n";
  eval { $dba->dbc->do("drop table stats_sites") };
  print "Dropping genes...\n";
  eval { $dba->dbc->do("drop table stats_genes") };
  print "Dropping brcmb tags...\n";
  eval { $dba->dbc->do("delete from protein_tree_tag where tag like 'bcrmb%';") };

  eval { $dba->dbc->do("delete from meta where meta_key='output_folder';") };

}

sub clean_compara_tree_tables {
  my $self            = shift;
  my $dba             = $self->dba;
  my @truncate_tables = qw^
    protein_tree_node protein_tree_tag protein_tree_member
    member
    sequence sequence_cds sequence_quality
    ^;
  map {
    print "$_\n";
    eval { $dba->dbc->do("truncate table $_"); }
  } @truncate_tables;

}

1;
