package Bio::Greg::Hive::SplitBySubtrees;

use strict;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::Process');

sub fetch_input {
  my $self = shift;

  ### DEFAULT PARAMETERS ###
  my $defaults = {
    seed_species => 9606,
    out_to       => 'Mammalia'
  };
  ##########################

  $self->load_all_params($defaults);
}

sub run {
  my $self = shift;

  my $out_to = $self->param('out_to');
  my $tax_tree =
    Bio::EnsEMBL::Compara::ComparaUtils->get_genome_taxonomy_below_level( $self->compara_dba,
    $out_to );

  my $seed_leaf;
  foreach my $leaf ( $tax_tree->leaves ) {
    $seed_leaf = $leaf if ( $leaf->taxon_id == $self->param('seed_species') );
  }

  # "Zoom out" from humans by taking
  my $cur_node = $seed_leaf->parent;
  while ( defined $cur_node ) {
    $self->flow_taxon_subtree($cur_node);
    $cur_node = $cur_node->parent;
  }
}

sub flow_taxon_subtree {
  my $self = shift;
  my $node = shift;

  print "Flowing " . $node->newick_format . "\n";

  my @taxon_ids = map { $_->taxon_id } $node->leaves;

  my $keep_taxon_string = join( ',', @taxon_ids );

  my $added_params = { keep_species => $keep_taxon_string };

  my $input_params = $self->_parse_string('input_id');
  my $output_params = { %$input_params, %$added_params };

  $self->new_data_id($output_params);

  my ($job_id) = @{ $self->dataflow_output_id( $output_params, 1 ) };
  print "  -> Fanned job: $job_id \n";
}

1;
