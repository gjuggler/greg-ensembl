package Bio::Greg::Hive::ComparaHiveLoaderUtils;

use warnings;
use strict;

use base ('Bio::Greg::Hive::HiveLoaderUtils');

our $param_set_counter = 1;

sub add_genes_to_analysis {
  my $self = shift;
  my $analysis_name = shift;
  my $gene_id_arrayref = shift;

  my @node_ids = ();
  foreach my $gene_id (@$gene_id_arrayref) {
    my $member;

    my $ext_member = $self->_find_member_by_external_id($gene_id);
    $member = $ext_member if (defined $ext_member);

    if (!defined $member) {
      $member = $self->dba->get_MemberAdaptor->fetch_by_source_stable_id(undef,$gene_id);
      my $gene_member = $member->gene_member if (defined $member);
      $member = $gene_member if (defined $gene_member);
    }

    if (!defined $member) {
      print STDERR " -> $gene_id not found!\n";
      next;
    }

    my $pta = $self->dba->get_ProteinTreeAdaptor;
    my $tree = $pta->fetch_by_gene_Member_root_id($member);
    $self->add_job_to_analysis($analysis_name,{node_id => $tree->node_id});
  }
}

sub _find_member_by_external_id {
  my $self = shift;
  my $id = shift;

  Bio::EnsEMBL::Registry->load_registry_from_multiple_dbs(
    {
      -host => 'ensembldb.ensembl.org',
      -user => 'anonymous'
    });
  my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor("human","core","gene");
  my $member_adaptor = $self->dba->get_MemberAdaptor;
  
  my @genes = @{ $gene_adaptor->fetch_all_by_external_name($id)};
  push @genes, $gene_adaptor->fetch_by_stable_id($id) if (scalar @genes == 0);
  
  foreach my $gene (@genes) {
    my $stable_id = $gene->stable_id;
    #print $stable_id."\n";
    my $member = $member_adaptor->fetch_by_source_stable_id(undef,$stable_id);
    #print $member."\n";
    return $member if (defined $member);
  }
  return undef;
}

sub add_parameter_set {
  my $self = shift;
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
  my $self = shift;
  my $dba = $self->dba;
  my @truncate_tables = qw^
      sitewise_omega sitewise_tag sitewise_genome
      parameter_set
      node_set_member node_set
      go_terms
      stats_sites stats_genes
      ^;
  map {
    print "$_\n";
    eval { $dba->dbc->do("truncate table $_")};
  } @truncate_tables;
  eval { $dba->dbc->do("drop table stats_sites")};
  eval { $dba->dbc->do("drop table stats_genes")};
  eval {$dba->dbc->do("delete from protein_tree_tag where tag like 'bcrmb%';")};
}

sub clean_compara_tree_tables {
 my $self = shift;
  my $dba = $self->dba;
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