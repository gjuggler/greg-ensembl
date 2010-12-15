package Bio::Greg::Hive::CollectGO;

use strict;

use Cwd;
use Time::HiRes qw(sleep);

use Bio::EnsEMBL::Compara::Member;
use Bio::EnsEMBL::Compara::ComparaUtils;

use Bio::Greg::EslrUtils;

use base ('Bio::Greg::Hive::Process', 'Bio::Greg::Gorilla::CollectDuplicationStats');

sub go_table_def {
  my $self = shift;

  return {
    data_id           => 'int',
    member_id         => 'int',
    taxon_id          => 'smallint',
    protein_id        => 'char32',
    gene_id           => 'char32',
    go_term           => 'char16',
    name              => 'string',
    ontology          => 'char8',
    namespace         => 'char32',
    subset            => 'char16',
    evidence_code     => 'char8',
    ancestral_mapping => 'smallint',
    extra_keys        => 'protein_id,taxon_id,go_term',
    unique_keys       => 'member_id,go_term,subset,ancestral_mapping'
  };
}

sub param_defaults {
  my $self = shift;
  ### DEFAULT PARAMETERS ###
  my $params = {    
    collect_duplication_species => '9606,9598,9593,9600',
    genes_table                 => 'stats_dups',
  };

  $params->{'go_taxon_ids'} = '9606,10090,9593';
  $params->{'go_table'}     = 'go_terms';
  $params->{'go_subsets'}   = 'goslim_goa,GO';
  #########################
  return $params;
}

sub fetch_input {
  my $self = shift;

  $self->load_all_params();
#  $self->create_table_from_params( $self->compara_dba, 'go_terms', $self->go_table_def );

  # Get a GO term adaptor and a gene adaptor (for human).
#  my $go_dba = Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Ontology', 'GOTerm' );
#  die if ( !defined $go_dba );
#  $self->param( 'go_dba', $go_dba );
}

sub run {
  my $self = shift;

#  $self->collect_go;

  # Remove me: temporary collection of gorilla gene stats...
  my $gene_stats = $self->get_gene_stats_def;
  $self->create_table_from_params( $self->compara_dba, 'stats_dups', $gene_stats );  
  $self->get_gene_data($self->param('node_id'),$self->param('parameter_set_id'));

}

sub collect_go {
  my $self = shift;

  my $tree   = $self->get_tree;
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

      print $leaf->stable_id . "\n";
      my $db_entries = $ts->get_all_DBLinks;

      # Grep out all the GO xref entries.
      my $allowed_subsets;
      map { $allowed_subsets->{$_} = 1 } split( ',', $self->param('go_subsets') );
      my @keepers = grep { $allowed_subsets->{ $_->dbname } == 1 } @{$db_entries};
      foreach my $db_e (@keepers) {

        #print $db_e->dbname."\n";
        $self->insert_go_term( $leaf, $db_e );
      }
    }
  }
}

sub insert_go_term {
  my $self = shift;
  my ( $leaf, $db_e ) = @_;

  my $go_dba = $self->param('go_dba');
  $self->throw("Error: GO DBA not defined") if ( !defined $go_dba );

  # Fetch the GO term objects from Ensembl's GO database.
  my $this_term = $go_dba->fetch_by_accession( $db_e->display_id );

  #print $this_term->name."\n";

  my $allowed_subsets;
  map { $allowed_subsets->{$_} = 1 } split( ',', $self->param('go_subsets') );
  my @subsets = @{ $this_term->subsets };
  @subsets = grep { $allowed_subsets->{$_} == 1 } @subsets;
  push @subsets, 'go' if ( $db_e->dbname eq 'GO' );

  foreach my $subset (@subsets) {

    # Collect all ancestral terms as well.
    my @all_terms;
    if ( $subset ne 'go' ) {
      @all_terms = @{ $go_dba->fetch_all_by_descendant_term( $this_term, $subset ) };
    } else {
      @all_terms = @{ $go_dba->fetch_all_by_descendant_term($this_term) };
    }

    # Remove any "copies" of the descendent term.
    @all_terms = grep { $_->accession ne $this_term->accession } @all_terms;

    # Add the original / descendent term to the list.
    push @all_terms, $this_term;

    foreach my $term (@all_terms) {
      my $is_ancestral = 1;
      $is_ancestral = 0 if ( $term == $this_term );

      my $evidence = '';
      print " " . $term->name . "\n";
      if ( $db_e->isa('Bio::EnsEMBL::GoXref') ) {
        $evidence = join ', ', @{ $db_e->get_all_linkage_types };
      } else {
      }

      my $values = {
        data_id    => $self->data_id,
        member_id  => $leaf->dbID,
        taxon_id   => $leaf->taxon_id,
        gene_id    => $leaf->gene_member->stable_id,
        protein_id => $leaf->stable_id,

        go_term   => $term->accession,
        namespace => $term->namespace,
        ontology  => $term->ontology,
        name      => $term->name,

        subset            => $subset,
        evidence_code     => $evidence,
        ancestral_mapping => $is_ancestral
      };

      $self->store_params_in_table( $self->dbc, $self->param('go_table'), $values );
    }
  }
}

1;
