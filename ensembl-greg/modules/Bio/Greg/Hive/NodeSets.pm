package Bio::Greg::NodeSetsB;

use strict;
use Time::HiRes qw(time gettimeofday tv_interval);
use Cwd;
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::NestedSet;
use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Hive::Process;
use Time::HiRes qw(sleep);
use Bio::EnsEMBL::Registry;
use Bio::Greg::Hive::Process;

use base ('Bio::Greg::Hive::Process');

#
# Some global-ish variables.
#
my $dba;
my $pta;

# INPUT FILES / OBJECTS.
my $tree;
my $params;

# OUTPUT FILES / OBJECTS / STATES.
my %node_set_hash;

sub fetch_input {
  my ($self) = @_;

  # Load up the Compara DBAdaptor.
  $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new( -DBCONN => $self->db->dbc );
  $pta = $dba->get_ProteinTreeAdaptor;

  ### DEFAULT PARAMETERS ###
  $params = {
    flow_node_set => 'Primates',
    flow_parent_and_children => 0,
    debug => 0
  };
  ##########################

  $self->load_all_params;

  $self->param('tree',$self->get_tree);
}

sub run {
  my $self = shift;

  my $tree = $self->param('tree');

  #print $tree->nhx_format('display_label') . "\n";

  $self->tag_root_nodes( $tree, "Primates" );
  $self->tag_root_nodes( $tree, "Glires" );
  $self->tag_root_nodes( $tree, "Laurasiatheria" );
  $self->tag_root_nodes( $tree, "Mammals" );
  $self->tag_root_nodes( $tree, "MammalPlusOutgroup" );
  $self->tag_root_nodes( $tree, "MammalPlusTwoOutgroups" );

  $self->tag_root_nodes( $tree, "Fish" );
  $self->tag_root_nodes( $tree, "Sauria" );

  $self->tag_nodes_with_clade_coverage( $tree, "Primates" );
  $self->tag_nodes_with_clade_coverage( $tree, "Glires" );
  $self->tag_nodes_with_clade_coverage( $tree, "Laurasiatheria" );
  $self->tag_nodes_with_clade_coverage( $tree, "Sauria" );
  $self->tag_nodes_with_clade_coverage( $tree, "Clupeocephala" );

}

sub write_output {
  my $self = shift;

  my $tree = $self->param('tree');
  
  if (defined $params->{flow_node_set}) {
    $self->autoflow_inputjob(0);
    my $flow_set = $params->{flow_node_set};
    
    foreach my $node ($tree->nodes) {
      next if ($node->is_leaf);
      
      my $id = $node->node_id;
      if ($node->has_tag("cc_root_".$flow_set)) {
	print " -> Flowing node $id\n";

        my $output_id = { node_id => $id };
        $self->dataflow_output_id( $output_id, 1 );
        if ( $params->{flow_parent_and_children} ) {
          my $i = 0;
          foreach my $child ( @{ $node->children } ) {
            my $output_id = { 
	      node_id => $child->node_id,
	      node_set_parent_id => $id, node_set_child_number => $i++ 
	      };
	    $self->dataflow_output_id( $output_id, 1 );
	    print "  --> Flowing child $output_id\n";
	  }
	}
      }
    }
  }
}

sub tag_root_nodes {
  my $self        = shift;
  my $tree        = shift;
  my $method_name = shift;

  my $base_p = {
    tree             => $tree,
    subtree_function => \&does_parent_have_clade_children,
    cc_min_size      => 4
  };

  my $params;

  if ( $method_name eq 'Primates' ) {
    $params = $self->replace_params( $base_p, { cc_Primates => 0.3, } );
  } elsif ( $method_name eq 'Glires' ) {
    $params = $self->replace_params( $base_p, { cc_Glires => 0.3, } );
  } elsif ( $method_name eq 'Laurasiatheria' ) {
    $params = $self->replace_params( $base_p, { cc_Laurasiatheria => 0.3, } );
  } elsif ( $method_name eq 'Mammals' ) {
    $params = $self->replace_params(
      $base_p, {
        cc_Primates       => 0.1,
        cc_Glires         => 0.1,
        cc_Laurasiatheria => 0.1,
      }
    );
  } elsif ( $method_name eq 'MammalPlusOutgroup' ) {
    $params = $self->replace_params(
      $base_p, {
        cc_Eutheria => 0.1,
        cc_any      => [ 'Sauria', 'Clupeocephala', 'Ciona', 'Marsupialia' ]
      }
    );
  } elsif ( $method_name eq 'MammalPlusTwoOutgroups' ) {
    $params = $self->replace_params(
      $base_p, {
        cc_Eutheria      => 0.1,
        cc_Sauria        => 0.01,
        cc_Clupeocephala => 0.01
      }
    );
  } elsif ( $method_name eq 'Fish' ) {
    $params = $self->replace_params( $base_p, { cc_Clupeocephala => 0.2 } );
  } elsif ( $method_name eq 'Sauria' ) {
    $params = $self->replace_params( $base_p, { cc_Sauria => 0.2 } );
  }

  my $tree = $params->{tree};

  #my $subtree_function = $params->{subtree_function};
  my @root_nodes = $self->get_smallest_subtrees_from_node( $tree, $params );

  foreach my $node (@root_nodes) {
    printf " -> %s %d %d\n", $method_name, scalar $node->leaves, $node->node_id;
    $node->store_tag( "cc_root_" . $method_name, 1 );
  }
}

sub tag_nodes_with_clade_coverage {
  my $self  = shift;
  my $tree  = shift;
  my $clade = shift;

  print "  -> Tagging clade coverage for clade: $clade...\n";

  foreach my $node ( $tree->nodes ) {
    next if ( $node->is_leaf );

    my $coverage_fraction = $self->clade_coverage_for_node( $node, $clade );

    if ( $coverage_fraction > 0 ) {

      #printf "%.20s %.3f\n",$node->newick_format,$coverage_fraction;
      $node->store_tag( "cc_$clade", sprintf( "%.3f", $coverage_fraction ) );
    }
  }
}

sub clade_coverage_for_node {
  my $self  = shift;
  my $node  = shift;
  my $clade = shift;

  my $key = "_taxon_ids_" . $clade;
  if ( !defined $self->{$key} ) {
    my @genomes = $self->get_genomes_within_clade( $dba, $clade );
    my @taxon_ids = map { $_->taxon_id } @genomes;
    $self->{$key} = \@taxon_ids;
  }
  my @taxon_ids = @{ $self->{$key} };

  my $taxon_hash;
  map { $taxon_hash->{ $_->taxon_id } = 1 } $node->leaves;

  my $total     = scalar @taxon_ids;
  my $clade_sum = 0;
  foreach my $genome_taxon (@taxon_ids) {
    $clade_sum++ if ( $taxon_hash->{$genome_taxon} == 1 );
  }

  my $coverage_fraction = $clade_sum / $total;
  return $coverage_fraction;
}

sub get_smallest_subtrees_from_node {
  my $self   = shift;
  my $node   = shift;
  my $params = shift;

  my $inclusion_function = $params->{subtree_function};

  my @children  = @{ $node->children };
  my @node_list = ();
  foreach my $child (@children) {
    push @node_list, $self->get_smallest_subtrees_from_node( $child, $params );
  }

  # If none of the children matched, then push ourselves onto the list if applicable.
  if ( scalar(@node_list) == 0 ) {
    if ( $inclusion_function->( $self, $node, $params ) ) {
      push @node_list, $node;
    }
  }

  return @node_list;
}

sub does_parent_have_clade_children {
  my $self   = shift;
  my $node   = shift;
  my $params = shift;

  return $self->generic_parent_has_good_children( $node, \&does_tree_have_clade_coverage, $params );
}

sub does_tree_have_clade_coverage {
  my $self   = shift;
  my $tree   = shift;
  my $params = shift;

  foreach my $key ( keys %$params ) {
    next unless ( $key =~ /cc_/ );
    my $value = $params->{$key};

    $key =~ s/cc_//;

    # the 'cc_any' tag signifies we'll take anything from the following clades.
    if ( $key eq 'any' ) {
      my $exists_any = 0;
      foreach my $clade ( @{$value} ) {
        my $coverage = $self->clade_coverage_for_node( $tree, $clade, $params );
        $exists_any = 1 if ( $coverage > 0 );
      }
      return 0 unless ($exists_any);
    } elsif ( $key eq 'min_size' ) {
      return 0 unless ( scalar( $tree->leaves ) >= $value );
    } else {
      my $coverage = $self->clade_coverage_for_node( $tree, $key, $params );
      return 0 unless ( $coverage >= $value );
    }
  }
  
  return 1;
}

sub generic_parent_has_good_children {
  my $self               = shift;
  my $node               = shift;
  my $inclusion_function = shift;
  my $params             = shift;

  # 1. check that this tree has superfamily coverage.
  # 2. check that our sister node has superfamily coverage.
  # 3. if (1) and (2) are met, return true.

  my $tree             = $node;
  my $parent           = $node->parent;
  my @parents_children = @{ $parent->children };
  if ( $parent->node_id == 1 ) {
    my $value = $inclusion_function->($self, $node, $params );
    #print "Root node. Values is $value\n";
    return $value;
  }

  my $sister;
  foreach my $ch (@parents_children) {
    $sister = $ch if ( $ch->node_id != $tree->node_id );
  }

  if ( $inclusion_function->( $self, $node, $params )
    && $inclusion_function->( $self, $sister, $params ) ) {
    return 1;
  }
  return 0;
}

sub DESTROY {
  my $self = shift;
  $tree->release_tree if ($tree);
  $tree = undef;
  $self->SUPER::DESTROY if $self->can("SUPER::DESTROY");
}

# Returns the NCBI taxnomy of Ensembl genomes below a given taxonomic clade.
sub get_genome_taxonomy_below_level {
  my $self          = shift;
  my $dba           = shift;
  my $root_taxon_id = shift || 'Fungi/Metazoa group';
  my $verbose       = shift || 0;

  my @gdbs = $self->get_all_genomes($dba);

  my @ncbi_ids = map { $_->taxon_id } @gdbs;

  my $taxon_a = $dba->get_NCBITaxonAdaptor;

  # Try first with a taxon label.
  my $root = $taxon_a->fetch_node_by_name($root_taxon_id);
  if ( !defined $root ) {
    $root = $taxon_a->fetch_node_by_taxon_id($root_taxon_id);
  }

  # Collect all genome_db leaves, plus their internal lineages, into an array.
  my %keepers;
  $keepers{ $root->node_id } = $root;
  foreach my $gdb (@gdbs) {
    my $tx    = $gdb->taxon;
    my $tx_id = $tx->taxon_id;

    if ( !$self->has_ancestor_node_id( $tx, $root ) ) {

      #print "not below tax level!!!\n";
      next;
    } else {
      print "Okay: " . $tx->name . "\n" if ($verbose);
    }

    my $node = $tx;
    while ( defined $node ) {
      last if ( !$self->has_ancestor_node_id( $node, $root ) );
      print $node->name . " " if ($verbose);
      $keepers{ $node->node_id } = $node;
      $node = $node->parent;
    }
    print "\n" if ($verbose);
  }

  my @nodes = values %keepers;
  print "Size: " . scalar(@nodes) . "\n" if ($verbose);

  #$taxon_a->clear_cache;
  my $new_tree = $taxon_a->_build_tree_from_nodes( \@nodes );

  # The call to "copy" seems to be important here...
  $new_tree = $new_tree->copy->minimize_tree;
  return $new_tree;
}

sub get_genomes_within_clade {
  my $self  = shift;
  my $dba   = shift;
  my $clade = shift || 1;

  my @gdbs = $self->get_all_genomes($dba);
  my $species_tree = $self->get_genome_taxonomy_below_level( $dba, $clade );

  my @genomes;
  foreach my $gdb (@gdbs) {
    my $leaf = $species_tree->find_node_by_node_id( $gdb->taxon->taxon_id );
    push @genomes, $gdb if ($leaf);
  }

  return @genomes;
}

sub get_all_genomes {
  my $self = shift;
  my $dba  = shift;

  my $gda     = $dba->get_GenomeDBAdaptor();
  my $all_dbs = $gda->fetch_all();
  my @all_genomes;
  foreach my $db (@$all_dbs) {
    push @all_genomes, $db if ( $db->taxon );
  }
  return @all_genomes;
}

# Returns whether $node has an ancestor with the same node_id as $ancestor.
sub has_ancestor_node_id {
  my $self     = shift;
  my $node     = shift;
  my $ancestor = shift;
  $node->throw("[$ancestor] must be a Bio::EnsEMBL::Compara::NestedSet object")
    unless ( $ancestor and $ancestor->isa("Bio::EnsEMBL::Compara::NestedSet") );
  $node = $node->parent;
  while ($node) {
    return 1 if ( $node->node_id == $ancestor->node_id );
    $node = $node->parent;
  }
  return 0;
}

1;
