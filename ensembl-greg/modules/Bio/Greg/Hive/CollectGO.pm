package Bio::Greg::Hive::CollectGO;

use strict;

use Cwd;
use Time::HiRes qw(sleep);

use Bio::EnsEMBL::Compara::Member;
use Bio::EnsEMBL::Compara::ComparaUtils;

use Bio::Greg::EslrUtils;

use base ('Bio::Greg::Hive::Process');

sub go_table_def {
  my $self = shift;

  return {
    data_id => 'int',
    member_id => 'int',
    taxon_id => 'smallint',
    protein_id => 'char32',
    gene_id => 'char32',
    go_term => 'char16',
    dbname => 'char16',
    evidence_code => 'char8',
    unique_keys => 'member_id,go_term,dbname'
  };
}

sub fetch_input {
  my $self = shift;

  ### DEFAULT PARAMETERS ###
  my $params = {};
  $params->{'go_taxon_ids'}    = '9606,10090,9593';
  $params->{'go_table'} = 'go_terms';
  #########################

  $self->load_all_params($params);
  $self->create_table_from_params($self->compara_dba, 'go_terms', $self->go_table_def);
}

sub run {
  my $self = shift;

  $self->collect_go;
}

sub collect_go {
  my $self = shift;

  my $tree = $self->get_tree;
  my @leaves = @{ $tree->get_all_leaves };

  my @taxon_ids;
  if ( defined $self->param('go_taxon_ids') ) {
    @taxon_ids = split( ",", $self->param('go_taxon_ids') );
  } else {
    @taxon_ids = (9606);
  }

  my %go_terms;
  foreach my $taxon_id (@taxon_ids) {
    my @taxon_leaves = grep { $_->taxon_id == $taxon_id } @leaves;

    foreach my $leaf (@taxon_leaves) {
      my $ts = $leaf->transcript;
      next if ( !defined $ts );

      print $leaf->stable_id."\n";
      my $db_entries = $ts->get_all_DBLinks;

      # Grep out all the GO xref entries.
      my @keepers = grep { $_->dbname =~ m/go/i || $_->dbname =~ m/goslim_goa/i } @{$db_entries};
      foreach my $db_e (@keepers) {
	$self->insert_go_term( $leaf, $db_e );
      }
    }
  }
}


sub insert_go_term {
  my $self = shift;
  my ( $leaf, $db_e ) = @_;

  my $evidence = '';
    if ($db_e->isa('Bio::EnsEMBL::GoXref')) {
      $evidence = join ', ', @{$db_e->get_all_linkage_types}; 
    }

  my $dbname = $db_e->dbname;

  my $values = $self->replace($self->get_params,{
    data_id => $self->data_id,
    member_id => $leaf->dbID,
    taxon_id => $leaf->taxon_id,
    gene_id => $leaf->gene_member->stable_id,
    protein_id => $leaf->stable_id,
    go_term => $db_e->display_id,
    dbname => $db_e->dbname,
    evidence_code => $evidence
				     });

  $self->store_params_in_table( $self->db_handle, $self->param('go_table'), $values);
}

1;
