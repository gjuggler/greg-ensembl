package Bio::Greg::NodeSets;

use strict;
use Time::HiRes qw(time gettimeofday tv_interval);
use Cwd;
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::NestedSet;
use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Hive::Process;
use Time::HiRes qw(sleep);
use Bio::EnsEMBL::Registry;
use Bio::Greg::ProcessUtils;

our @ISA = qw(Bio::EnsEMBL::Hive::Process Bio::Greg::ProcessUtils);

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

  $pta->local_mode(-1);

  ### DEFAULT PARAMETERS ###
  $params = {
    flow_node_set => 0,
    debug         => 0
  };
  ##########################

  # Fetch parameters from the two possible locations. Input_id takes precedence!
  my $p_params = $self->get_params( $self->parameters );
  my $i_params = $self->get_params( $self->input_id );
  my $node_id  = $i_params->{'protein_tree_id'} || $i_params->{'node_id'};
  my $t_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_tree_tags( $dba, $node_id );

  $params = $self->replace_params( $params, $p_params, $i_params, $t_params );
  Bio::EnsEMBL::Compara::ComparaUtils->hash_print($params);

  $tree = Bio::EnsEMBL::Compara::ComparaUtils->get_tree_for_comparative_analysis( $dba, $params );
  $tree = $tree->minimize_tree;
}

sub run {
  my $self = shift;

  print $tree->newick_format . "\n";

  my $base = {
    tree          => $tree,
    node_set_hash => \%node_set_hash,
    node_set_id   => 1
  };
  my $p;

  print "Adding human gene subtrees...\n";
  $p = $self->replace_params(
    $base, {
      node_set_id      => 1,
      subtree_function => \&isa_human_gene_subtree
    }
  );
  $self->add_all_subtrees_to_hash($p);

  print "Adding small dups...\n";
  $p = $self->replace_params( $base, { node_set_id => 2 } );
  $self->add_small_balanced_dups_to_hash($p);

  print "Adding large dups...\n";
  $p = $self->replace_params( $base, { node_set_id => 3 } );
  $self->add_large_balanced_dups_to_hash($p);

  print "Adding overlapping dups...\n";
  $p = $self->replace_params( $base, { node_set_id => 4 } );
  $self->add_overlapping_dups_to_hash($p);

  print "Adding small subtrees...\n";
  $p = $self->replace_params(
    $base, {
      node_set_id      => 5,
      subtree_function => \&is_tree_small
    }
  );
  $self->add_all_subtrees_to_hash($p);

  print "Adding med subtrees...\n";
  $p = $self->replace_params(
    $base, {
      node_set_id      => 6,
      subtree_function => \&is_tree_med
    }
  );
  $self->add_all_subtrees_to_hash($p);

  print "Adding large subtrees...\n";
  $p = $self->replace_params(
    $base, {
      node_set_id      => 7,
      subtree_function => \&is_tree_large
    }
  );
  $self->add_all_subtrees_to_hash($p);

  print "Adding superfamily subtrees...\n";
  $p = $self->replace_params(
    $base, {
      node_set_id      => 8,
      subtree_function => \&does_parent_have_superfamily_children
    }
  );
  $self->add_all_subtrees_to_hash($p);

  print "Adding Tim superfamily subtrees...\n";
  $p = $self->replace_params(
    $base, {
      node_set_id      => 9,
      subtree_function => \&does_parent_have_tim_children
    }
  );
  $self->add_all_subtrees_to_hash($p);

  print "Adding Primate subtrees...\n";
  $p = $self->replace_params(
    $base, {
      node_set_id      => 10,
      subtree_function => \&does_parent_have_primate_children,
      num_primates     => 3,
      num_glires       => 0
    }
  );
  $self->add_all_subtrees_to_hash($p);

  $p = $self->replace_params(
    $base, {
      node_set_id      => 11,
      subtree_function => \&does_parent_have_primate_children,
      num_primates     => 4,
      num_glires       => 0
    }
  );
  $self->add_all_subtrees_to_hash($p);

  $p = $self->replace_params(
    $base, {
      node_set_id      => 12,
      subtree_function => \&does_parent_have_primate_children,
      num_primates     => 5,
      num_glires       => 0
    }
  );
  $self->add_all_subtrees_to_hash($p);

  $p = $self->replace_params(
    $base, {
      node_set_id      => 13,
      subtree_function => \&does_parent_have_primate_children,
      num_primates     => 4,
      num_glires       => 1
    }
  );
  $self->add_all_subtrees_to_hash($p);

  my $sql = "INSERT IGNORE INTO node_set (node_set_id,name,allows_overlap) values (?,?,?)";
  my $sth = $dba->dbc->prepare($sql);
  $sth->execute( 1,  "Human Gene Subtrees",          0 );
  $sth->execute( 2,  "Small Balanced Duplications",  0 );
  $sth->execute( 3,  "Large Balanced Duplications",  0 );
  $sth->execute( 4,  "Overlapping Duplications",     1 );
  $sth->execute( 5,  "Small Subtrees",               0 );
  $sth->execute( 6,  "Med Subtrees",                 0 );
  $sth->execute( 7,  "Large Subtrees",               0 );
  $sth->execute( 8,  "Superfamily Subtrees",         0 );
  $sth->execute( 9,  "Tim Subtrees",                 0 );
  $sth->execute( 10, "Primates n=3",                 0 );
  $sth->execute( 11, "Primates n=4",                 0 );
  $sth->execute( 12, "Primates n=5",                 0 );
  $sth->execute( 13, "Primates n=4 plus Glires n=1", 0 );
  $sth->finish;
}

sub add_all_subtrees_to_hash {
  my $self   = shift;
  my $params = shift;

  my $tree             = $params->{tree};
  my $hash             = $params->{node_set_hash};
  my $node_set_id      = $params->{node_set_id};
  my $subtree_function = $params->{subtree_function};

  my @root_node_ids = $self->get_smallest_subtrees_from_node( $tree, $subtree_function, $params );

  $hash->{$node_set_id} = \@root_node_ids;
}

sub add_small_balanced_dups_to_hash {
  my $self   = shift;
  my $params = shift;

  my $tree        = $params->{tree};
  my $hash        = $params->{node_set_hash};
  my $node_set_id = $params->{node_set_id};

  my $subtree_function = \&isa_balanced_good_duplication;
  my @keeper_ids = $self->get_smallest_subtrees_from_node( $tree, $subtree_function, $params );

  $hash->{$node_set_id} = \@keeper_ids;
}

sub add_large_balanced_dups_to_hash {
  my $self   = shift;
  my $params = shift;

  my $tree        = $params->{tree};
  my $hash        = $params->{node_set_hash};
  my $node_set_id = $params->{node_set_id};

  my $subtree_function = \&isa_balanced_good_duplication;
  my @keeper_ids = $self->get_largest_subtrees_from_node( $tree, $subtree_function, $params );

  $hash->{$node_set_id} = \@keeper_ids;
}

sub add_overlapping_dups_to_hash {
  my $self   = shift;
  my $params = shift;

  my $tree        = $params->{tree};
  my $hash        = $params->{node_set_hash};
  my $node_set_id = $params->{node_set_id};

  my @keeper_ids = ();
  foreach my $node ( $tree->nodes ) {

    #    print "Added ".$node->node_id." to list!\n";
    push @keeper_ids, $node->node_id if ( is_balanced_dup_node($node) );
  }

  $hash->{$node_set_id} = \@keeper_ids;
}

# A generic method for gathering a set of small, non-overlapping subtrees.

sub get_smallest_subtrees_from_node {
  my $self               = shift;
  my $node               = shift;
  my $inclusion_function = shift;
  my $params             = shift;

  my @children  = @{ $node->children };
  my @node_list = ();
  foreach my $child (@children) {
    push @node_list, $self->get_smallest_subtrees_from_node( $child, $inclusion_function, $params );
  }

  # If none of the children matched, then push ourselves onto the list if applicable.
  if ( scalar(@node_list) == 0 ) {
    if ( $inclusion_function->( $node, $params ) ) {
      push @node_list, $node->node_id;
    }
  }

  return @node_list;
}

# A generic method for gathering a set of large, non-overlapping subtrees.
sub get_largest_subtrees_from_node {
  my $self               = shift;
  my $node               = shift;
  my $inclusion_function = shift;

  if ( $inclusion_function->($node) ) {

    # Return this node immediately if it satisfies the conditions.
    print "Added " . $node->node_id . " to list!\n";
    return ( $node->node_id );
  }

  # Recurse through child nodes to see if they satisfy the conditions.
  my @children  = @{ $node->children };
  my @node_list = ();
  foreach my $child (@children) {
    push @node_list, $self->get_largest_subtrees_from_node( $child, $inclusion_function );
  }
  return @node_list;
}

sub isa_human_gene_subtree {
  my $node = shift;

  return 1 if ( is_tree_big_enough($node) && is_contains_human_genes($node) );
  return 0;
}

sub isa_balanced_good_duplication {
  my $node = shift;

  return 1 if ( is_balanced_dup_node( $node, 6 ) );
  return 0;
}

sub is_duplication {
  my $node = shift;

  my $is_dup = 0;
  my $val    = $node->get_tagvalue("Duplication");
  $is_dup = 1 if ( defined($val) && $val ne '' && $val > 0 );

  my $dubious = 0;
  $dubious = 1 if ( $node->get_tagvalue("dubious_duplication") );

  return 1 if ( $is_dup && !$dubious );
  return 0;
}

sub is_contains_human_genes {
  my $node      = shift;
  my @leaves    = $node->leaves;
  my @hum_genes = grep { $_->taxon_id == 9606 } @leaves;
  return 1 if ( scalar(@hum_genes) > 0 );
  return 0;
}

sub is_tree_big_enough {
  my $node = shift;
  return is_tree_small($node);
}

sub is_tree_small {
  my $node = shift;
  return ( scalar( $node->leaves ) >= 6 );
}

sub is_tree_med {
  my $node = shift;
  return ( scalar( $node->leaves ) >= 10 );
}

sub is_tree_large {
  my $node = shift;
  return ( scalar( $node->leaves ) >= 20 );
}

sub does_parent_have_superfamily_children {
  my $node = shift;
  return generic_parent_has_good_children( $node, \&does_tree_have_superfamily_coverage );
}

sub does_parent_have_primate_children {
  my $node   = shift;
  my $params = shift;
  return generic_parent_has_good_children( $node, \&does_tree_have_primate_coverage, $params );
}

sub does_parent_have_primate_children_b {
  my $node = shift;
  return generic_parent_has_good_children( $node, \&does_tree_have_primate_coverage, 4 );
}

sub does_parent_have_tim_children {
  my $node = shift;
  return generic_parent_has_good_children( $node, \&does_tree_have_tim_family_coverage );
}

sub generic_parent_has_good_children {
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
    return $inclusion_function->( $node, $params );
  }

  my $sister;
  foreach my $ch (@parents_children) {
    $sister = $ch if ( $ch->node_id != $tree->node_id );
  }

  if ( $inclusion_function->( $node, $params ) && $inclusion_function->( $sister, $params ) ) {
    return 1;
  }
  return 0;
}

sub does_tree_have_primate_coverage {
  my $tree   = shift;
  my $params = shift;

  my $n        = $params->{num_primates};
  my $n_glires = $params->{num_glires};

  # 1. require coverage in N different primates.
  my $adaptor = $tree->adaptor;
  my $ncbi    = $adaptor->db->get_NCBITaxonAdaptor();

  if ( !exists $params->{_primate_taxonomy} ) {
    $params->{_primate_taxonomy} =
      Bio::EnsEMBL::Compara::ComparaUtils->get_genome_taxonomy_below_level( $adaptor->db, 9443 );
  }
  my $primates = $params->{_primate_taxonomy};
  my %primate_hash;
  map { $primate_hash{ $_->taxon_id } = 1 } $primates->leaves;

  if ( !exists $params->{_glires_taxonomy} ) {
    $params->{_glires_taxonomy} =
      Bio::EnsEMBL::Compara::ComparaUtils->get_genome_taxonomy_below_level( $adaptor->db, 314147 );
  }
  my $glires = $params->{_glires_taxonomy};
  my %glires_hash;
  map { $glires_hash{ $_->taxon_id } = 1 } $glires->leaves;

  my %primates_found;
  my @leaves = $tree->leaves;
  foreach my $leaf ( $tree->leaves ) {
    $primates_found{ $leaf->taxon_id } = 1 if ( exists $primate_hash{ $leaf->taxon_id } );
  }

  my %glires_found;
  my @leaves = $tree->leaves;
  foreach my $leaf ( $tree->leaves ) {
    $glires_found{ $leaf->taxon_id } = 1 if ( exists $glires_hash{ $leaf->taxon_id } );
  }

  my $num_primates = scalar keys %primates_found;
  my $num_glires   = scalar keys %glires_found;

  return 1 if ( $num_primates >= $n && $num_glires >= $n_glires );
  return 0;
}

sub does_tree_have_superfamily_coverage {
  my $tree = shift;

# 1. require 2 of 4 mammalian families to be represented (Laurasiatheria, Primates, Glires, Afrotheria)
# 2. require 1 of 3 outgroups to be represented (Opossum, Chicken, Tetraodon)

  my %mamm_hash = (
    '314145' => 'Laurasiatheria',
    '9443'   => 'Primates',
    '314147' => 'Glires',
    '311790' => 'Afrotheria'
  );
  my %outgroup_hash = (
    '9031'  => 'Chicken',
    '13616' => 'Opossum',
    '99883' => 'Tetraodon'
  );

  my %fams_found;
  my %outgroups_found;

  my @nodes = $tree->nodes;
  foreach my $node (@nodes) {
    my $tax_id = $node->taxon_id;
    $fams_found{$tax_id}      = 1 if ( exists $mamm_hash{$tax_id} );
    $outgroups_found{$tax_id} = 1 if ( exists $outgroup_hash{$tax_id} );
  }

  my $num_mamm_fams = scalar( keys %fams_found );
  my $num_outgroups = scalar( keys %outgroups_found );

#  printf "Node ID: %d  Mammal coverage: %d/%d  Outgroups: %d/%d\n",$tree->node_id,$num_mamm_fams,4,$num_outgroups,3;
  return 1 if ( $num_mamm_fams >= 2 && $num_outgroups >= 1 );
  return 0;
}

sub does_tree_have_tim_family_coverage {
  my $tree = shift;

# 1. require 2 of 4 mammalian families to be represented (Laurasiatheria, Primates, Glires, Afrotheria)

  my %mamm_hash = (
    '314145' => 'Laurasiatheria',
    '9443'   => 'Primates',
    '314147' => 'Glires',
    '311790' => 'Afrotheria'
  );
  my %fams_found;

  my @nodes = $tree->nodes;
  foreach my $node (@nodes) {
    my $tax_id = $node->taxon_id;
    $fams_found{$tax_id} = 1 if ( exists $mamm_hash{$tax_id} );
  }

  my $num_mamm_fams = scalar( keys %fams_found );

  #printf "Node ID: %d  Mammal coverage: %d/%d\n",$tree->node_id,$num_mamm_fams,4;
  return 1 if ( $num_mamm_fams >= 2 );
  return 0;
}

sub is_tree_tim_worthy {
  my $tree    = shift;
  my $node_id = shift;

  my $subtree = $tree->find_node_by_node_id($node_id);

  my @nodes  = $subtree->nodes;
  my @leaves = $subtree->leaves;

  return 0 if ( scalar(@leaves) < 8 );

  my $dup_count = 0;
  foreach my $node (@nodes) {
    $dup_count++ if ( is_duplication($node) );
    return 0 if ( is_duplication($node) && scalar( $node->leaves ) > 2 );
  }

  #  $subtree->print_tree;
  return 1 if ( $dup_count <= 4 );
  return 0;
}

sub is_balanced_dup_node {
  my $node         = shift;
  my $minimum_size = shift;

  my $dup     = $node->get_tagvalue("Duplication");
  my $dubious = $node->get_tagvalue("dubious_duplication");
  if ( $dup && !$dubious ) {

    # Ensure that the subtrees are balanced.
    my @children = @{ $node->children };
    if ( scalar(@children) == 2 ) {

      #print "two children!\n";
      my $a = $children[0];
      my $b = $children[1];

      my @a_leaves = $a->leaves;
      my @b_leaves = $b->leaves;

      my $ratio = scalar(@a_leaves) / scalar(@b_leaves);
      $ratio = 1 / $ratio if ( $ratio > 1 );

      #print "Left-right ratio: $ratio\n";

      my $n_a = scalar @a_leaves;
      my $n_b = scalar @b_leaves;

      if ( $ratio > 0.5
        && $n_a >= $minimum_size
        && $n_b > $minimum_size ) {

        # Only now can we add the dup to the list.
        return 1;
      }
    }
  }
  return 0;
}

sub write_output {
  my $self = shift;
  $self->autoflow_inputjob(0);

  my @value_string_array;

  foreach my $key ( sort { $a <=> $b } keys %node_set_hash ) {
    print "NODE SET ID: $key\n";
    my @ids = @{ $node_set_hash{$key} };
    foreach my $id (@ids) {
      print "  -> Inserting node set member $key -> $id\n";
      push @value_string_array, "($key,$id)";

      my $tree = $pta->fetch_node_by_node_id($id);

      #      print " ->". $tree->newick_format."\n";
      if ( $params->{'debug'} ) {
        foreach my $leaf ( $tree->leaves ) {
          if ( $leaf->taxon_id == 9606 ) {
            print "  " . $leaf->stable_id . "\n";
          }
        }
      }

      if ( $key == $params->{'flow_node_set'} ) {
        my $output_id = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string( { node_id => $id } );
        $self->dataflow_output_id( $output_id, 1 );
        print "  -> Flowing node $output_id\n";
        if ( $params->{'flow_parent_and_children'} ) {
          foreach my $child ( @{ $tree->children } ) {
            my $output_id =
              Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string( { node_id => $child->node_id } );
            $self->dataflow_output_id( $output_id, 1 );
            print "  -> Flowing child $output_id\n";
          }
        }
      }
    }
  }

  my $value_string = join( ",", @value_string_array );

  if ( scalar(@value_string_array) > 0 ) {
    my $sql = "REPLACE INTO node_set_member (node_set_id,node_id) values $value_string";
    $dba->dbc->do($sql);
  }
}

sub DESTROY {
  my $self = shift;
  $tree->release_tree if ($tree);
  $tree = undef;
  $self->SUPER::DESTROY if $self->can("SUPER::DESTROY");
}

1;
