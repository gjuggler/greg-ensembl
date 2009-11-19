package Bio::EnsEMBL::Compara::TreeUtils;

use Bio::TreeIO;
use Bio::EnsEMBL::Compara::LocalMember;

#
# A grab bag of useful methods for tree manipulations.
#

my $TREE = "Bio::EnsEMBL::Compara::TreeUtils";
my $ALN = "Bio::EnsEMBL::Compara::AlignUtils";
my $COMPARA = "Bio::EnsEMBL::Compara::ComparaUtils";

my $TREEI = "Bio::Tree::TreeI";
my $NSET = "Bio::EnsEMBL::Compara::NestedSet";
my $PTREE = "Bio::EnsEMBL::Compara::ProteinTree";

# Returns whether $node has an ancestor with the same node_id as $ancestor.
sub has_ancestor_node_id {
  my $class = shift;
  my $node = shift;
  my $ancestor = shift;
  $node->throw("[$ancestor] must be a Bio::EnsEMBL::Compara::NestedSet object")
    unless ($ancestor and $ancestor->isa("Bio::EnsEMBL::Compara::NestedSet"));
  $node = $node->parent;
  while($node) {
    return 1 if ($node->node_id == $ancestor->node_id);
    $node = $node->parent;
  }
  return 0;
}

# Returns a Bio::Tree::TreeI object given a Bio::EnsEMBL::Compara::NestedSet object.
sub to_treeI {
  my $class = shift;
  my $tree = shift;

  my $newick = "";
  if (!ref $tree) {
    if (-e $tree) {
      my $treeI = $class->treeI_from_file($tree);
      return $treeI;
    } else {
      $newick = $tree;
    }
  } elsif ($tree->isa($TREEI)) {
    return $tree;
  } elsif ($tree->isa($NSET)) {
    $newick = $tree->newick_format();
  }
  
  open(my $fake_fh, "+<", \$newick);
  my $treein = new Bio::TreeIO
    (-fh => $fake_fh,
     #-verbose => 1,
     -format => 'newick');
  my $treeI = $treein->next_tree;
  $treein->close;
  return $treeI;
}

sub treeI_from_newick {
  my $class = shift;
  my $newick = shift;

  $newick .= ";" if ($newick !~ /;/);
  $newick =~ s/\s//g;

  # Parse out and store the NHX annotations... because bioperl-live sucks.
  my $annotation;
  while ($newick =~ m/[,\(\)]+(.*?)(:.*?)?\[&&nhx:(.*?)\]/gi) {
    my $hash;
    my @keyvals = split(/:/,$3);
    foreach my $keyval (@keyvals) {
      my ($key,$val) = split(/=/,$keyval);
      $hash->{$key} = $val;
    }
    $annotation->{$1} = $hash if (length $1 > 1);
  }
  
  # Now remove the NHX annotations from the string.
  $newick =~ s/\[&&nhx.*?\]//ig;
  
  # Load the tree using TreeIO.
  open(my $fake_fh, "+<", \$newick);
  my $treein = new Bio::TreeIO(-fh => $fake_fh, -format => 'newick');
  my $treeI = $treein->next_tree;
  $treein->close();

  # Put annotations back into nodes.
  foreach my $node ($treeI->get_nodes) {
    my $hash = $annotation->{$node->id};
    if (defined $hash) {
      foreach my $key (keys %{$hash}) {
	$node->add_tag_value($key,$hash->{$key});
      }
    }
  }
  
  # Test outputting the tree to a string.
  #my $new_newick;
  #open(my $fake_fh,">",\$new_newick);
  #my $treeout = new Bio::TreeIO(-fh => $fake_fh, -format => 'newick');
  #$treeout->write_tree($treeI);
  #print "new_new:\n" . $new_newick."\n";

  return $treeI;
}

sub treeI_from_file {
  my $class = shift;
  my $file = shift;
  
  my $treein = new Bio::TreeIO(-file=>$file, -format => 'newick');
  my $treeI = $treein->next_tree;
  $treein->close();
  return $treeI;
}

# Creates a Bio::EnsEMBL::Compara::ProteinTree object from a Newick string.
sub from_newick {
  my $class = shift;
  my $newick = shift;

  return $class->from_treeI($class->treeI_from_newick($newick));
}

# Creates a Bio::EnsEMBL::Compara::NestedSet from a file.
sub from_file {
  my $class = shift;
  my $file = shift;
  
  open(IN,$file);
  my @lines = <IN>;
  close(IN);
  my $newick = join("",@lines);
  
  use Bio::EnsEMBL::Compara::Graph::NewickParser;
  return Bio::EnsEMBL::Compara::Graph::NewickParser::parse_newick_into_tree($newick);  
}

# Creates a ProteinTree from a Bio::Tree::TreeI object.
sub from_treeI {
  my $class = shift;
  my $treeI = shift;
  
  my $rootI = $treeI->get_root_node;
  my $node = new Bio::EnsEMBL::Compara::NestedSet;
  $node->init;
  
  # Kick off the recursive, parallel node adding.
  $class->add_nodeI_to_node($node,$rootI);
  
  return $node;
}

# Recursive helper for new_from_TreeI.
sub add_nodeI_to_node {
  my $class = shift;
  my $node = shift; # Our node object (Compara)
  my $nodeI = shift; # Our nodeI object (BioPerl)
  
  foreach my $c ($nodeI->each_Descendent) {
    my $child = ref($node)->new;
    
    my $name = $c->id || "";
    $name =~ s/^\s+//;
    $name =~ s/\s+$//;
    
    if ($c->is_Leaf) {
      $child = Bio::EnsEMBL::Compara::LocalMember->new();
      $child->init;
      $child->stable_id($name);
      $child->source_name("");
    }
    
    # Set name.
    $child->name($name) if (length $name > 1);
    
    # Set branch length.
    $node->add_child($child,$c->branch_length);

    #print "Adding node $name ".$c->branch_length."\n";
    
    # Add the tags.
    my $tags = $c->{'_tags'};
    foreach my $tag (keys %{$tags}) {
      $tags->{$tag} = $tags->{$tag}[0]
    }
    $child->{'_tags'} = $tags;

    # Recurse.
    $class->add_nodeI_to_node($child,$c);
  }
  $node->{'_children_loaded'} = 1;
}


# Scales all the branches in a tree by a given factor.
# NOTE: Returns a COPY of the tree!
# @created GJ 2008-12-12
sub scale {
  my $class = shift;
  my $tree = shift;
  my $scale = shift;

  if ($tree->isa("Bio::Tree::TreeI")) {
    # Create a copy of the tree.
    my $newick = $class->to_newick($tree);
    $tree = $class->to_treeI($newick);

    for my $node ($tree->get_nodes) {
      my $bl = $node->branch_length;
      $bl = 0 unless (defined $bl);
      $node->branch_length($bl*$scale);
    }
  } elsif ($tree->isa("Bio::EnsEMBL::Compara::NestedSet")) {
    foreach my $node (@{$tree->get_all_nodes}) {
      my $bl = $node->distance_to_parent;
      $bl = 0 unless (defined $bl);
      $node->distance_to_parent($bl*$scale);
    }
  }
  return $tree;
}

sub scale_max_to {
  my $class = shift;
  my $tree = shift;
  my $new_max = shift;

  my $max_dist = $tree->max_distance;
  my $scale_factor = $new_max / $max_dist;
  return $class->scale($tree,$scale_factor);
}

# Scales a tree to a certain total branch length.
# NOTE: Tree is modified in-place!
# @created GJ 2009-01-09
sub scale_to {
  my $class = shift;
  my $tree = shift;
  my $new_total = shift;

  my $scale_factor;
  if ($tree->isa("Bio::Tree::TreeI")) {
    my $total_dist = $class->max_distance($tree);
    $scale_factor = $new_total / $total_dist;
  } elsif ($tree->isa("Bio::EnsEMBL::Compara::NestedSet")) {
#    my $total_dist = $tree->max_distance;
#    my $total_dist = $class->max_distance_unrooted($tree);
    my $total_dist = $class->total_distance($tree);
    $scale_factor = $new_total / $total_dist;
  }
  return $class->scale($tree,$scale_factor);
}

sub total_distance {
  my $class = shift;
  my $tree = shift;

  my $sum = 0;
  foreach my $n ($tree->nodes) {
    $sum += $n->distance_to_parent;
  }
  return $sum;
}

sub max_distance_unrooted {
  my $class = shift;
  my $tree = shift;

  my $max_dist = 0;
  foreach my $a ($tree->leaves) {
    foreach my $b( $tree->leaves) {
      my $sum = $a->distance_to_node($b);
      $max_dist = $sum if ($sum > $max_dist);
    }
  }
  return $max_dist;
}

# Returns the maximum root-to-tip distance for a Bio::Tree:TreeI object.
# @created GJ 2009-01-09
sub max_distance {
  my $class = shift;
  my $tree = shift;

  die ("$tree not a $TREEI!") unless ($tree->isa($TREEI));
  
  my $max_dist = 0;
  foreach my $node ($tree->get_leaf_nodes()) {
    my $dist = $class->distance_to_root($tree,$node);
    $max_dist = $dist if ($dist > $max_dist);
  }
  return $max_dist;
}

# GJ 2009-01-09 : Distance from a given node to the root.
# Returns the distance from the given node to the root of the given TreeI.
sub distance_to_root {
  my $class = shift;
  my $tree = shift;
  my $node = shift;

  die ("$tree not a $TREEI!") unless ($tree->isa($TREEI));

  my $branch_length_sum = 0;
  while (defined $node) {
    $branch_length_sum += $node->branch_length;
    $node = $node->ancestor();
  }
  return $branch_length_sum;
}

# Convert a TreeI or ProteinTree object ot a Newick string.
# @created GJ 2009-01-09
sub to_newick {
  my $class = shift;
  my $tree = shift;

  my $ref = ref $tree;
  if ($ref =~ /Bio::Tree/i) {
    use Bio::TreeIO;
    my $string = "";
    open(my $fake_fh, "+>", \$string);
    my $out = Bio::TreeIO->new(-fh => $fake_fh,
			       -format => "newick");
    $out->write_tree($tree);
    $out->close();
    return $string;
  } elsif ($ref =~ /Bio::EnsEMBL/i) {
    return $tree->newick_format();
  } elsif (-e $tree) {
    my $treeI = $class->treeI_from_file($tree);
    return $class->to_newick($treeI);
  } else {
    return $tree;
  }
}

sub to_file {
  my $class = shift;
  my $tree = shift;
  my $out_file = shift;
  my $params = shift;

  my $newick = $class->to_newick($tree);
  if ($params->{'node_ids'}) {
    $newick = $tree->newick_format('int_node_id');
  }

  open(OUT,">".$out_file);
  print OUT $newick;
  print $newick."\n";
  close(OUT);
  return $out_file;
}

# Deletes the lineage leading to $del_me.
sub delete_lineage {
  my $class = shift;
  my $tree = shift;
  my $del_me = shift;

  die ("$tree not a $NSET!") unless ($tree->isa($NSET));    

  my $parent = $del_me->parent;
  while ($parent) {
    my $num_children = scalar @{$parent->children};
    if ($num_children > 1) {
      $tree->remove_nodes([$del_me]);
      return $tree;
    } elsif ($num_children == 1) {
      $tree->remove_nodes([$del_me]);
      $del_me = $parent;
      $parent = $del_me->parent;
    }
  }
  return $tree;
}

sub remove_members_by_member_id {
  my $class = shift;
  my $tree = shift;
  my $node_list = shift;

  my @node_ids = split(",",$node_list);
  my $ids_hash;
  map {$ids_hash{$_}=1} @node_ids;

  print "Removing specified nodes from tree: @node_ids \n";
  print "  Before:" . scalar($tree->leaves) . " leaves\n";

  foreach my $node ($tree->nodes) {
    if (exists $ids_hash{$node->node_id}) {
      $TREE->delete_lineage($tree,$node);
      $tree->minimize_tree();
    }
  }
  print "  After:" . scalar($tree->leaves) . " leaves\n";
  return $tree;
}

sub remove_members_by_taxon_id {
  my $class = shift;
  my $tree = shift;
  my $species_arrayref = shift;

  my @tax_ids = @{$species_arrayref};
  # Turn the arrayref into a hash.
  my %tax_hash;
  map {$tax_hash{$_}=1} @tax_ids;

  print "Pruning leaves with taxon_ids: ". join(",",sort(@tax_ids))."\n";
  print "  Before: " . scalar($tree->leaves) . "\n";

  foreach my $leaf ($tree->leaves) {
    if (exists $tax_hash{$leaf->taxon_id}) {
      $class->delete_lineage($tree,$leaf);
      $tree->minimize_tree();
    }
  }

  print "  After: " . scalar($tree->leaves) . "\n";
  return $tree;
}

# Returns a list of all taxon IDs within a given tree. 
# @created GJ 2009-07-29
sub get_species_in_tree {
  my $class = shift;
  my $tree = shift;
  my @leaves = $tree->leaves;
  my %taxon_hash;
  map {$taxon_hash{$_->taxon_id}=1} @leaves;
  return keys %taxon_hash;
}

# Deletes all leaves that don't lie within the specified taxonomic clade.
# Useful for pruning subtrees to be within a certain taxonomic range.
# NOTE: this function does NOT copy the tree before pruning (i.e., make a copy BEFORE calling!)
# @created GJ 2009-01-12
sub prune_leaves_outside_taxon {
  my $class = shift;
  my $tree = shift;
  my $taxon_id = shift;

  print "Pruning leaves outside taxon $taxon_id.\n";
  print "  Before: " . scalar(@{$tree->get_all_leaves}) . " leaves.\n";

  my $db = $tree->adaptor->db;
  my $ncbi_a = $db->get_NCBITaxonAdaptor;

  sub is_leaf_within_taxon {
    my $taxon_id = shift;
    my $leaf = shift;
    my $db = shift;

    my $leaf_taxon = $leaf->taxon_id;
    my $sth = $db->dbc->prepare(qq^
				SELECT * from ncbi_taxa_node p, ncbi_taxa_node c
				WHERE p.taxon_id=$taxon_id AND c.taxon_id=$leaf_taxon AND
				c.left_index BETWEEN p.left_index AND p.right_index;
				^);
    $sth->execute();
    if (defined $sth->fetchrow_arrayref) {
      $sth->finish;
      return 1;
    } else {
      $sth->finish;
      return 0;
    }
  }

  my @remove_me;
  foreach my $leaf ($tree->leaves) {
    next if (is_leaf_within_taxon($node_id,$leaf,$db));
    
    print "  -> Deleting leaf: " . $leaf->stable_id . "  " . $leaf->taxon->binomial . "\n";
    $class->delete_lineage($tree,$leaf);
  }

  print "  After: " . scalar(@{$tree->get_all_leaves}) . " leaves.\n";

  return $tree;
}

# Extracts a subtree given a ProteinTree and an arrayref of node_ids.
# @updated GJ 2009-01-14 : Smarter version, uses the inner method NestedSetAdaptor->_build_tree_from_nodes.
# @created GJ 2009-01-12
sub extract_subtree_from_leaves {
  my $class = shift;
  my $tree = shift;
  my $node_ids = shift;	# Array ref of node_ids.

  die("Object not a NestedSet!") unless ($tree->isa("Bio::EnsEMBL::Compara::NestedSet"));

  my $copy = $tree->copy;
  my @keepers = @{$node_ids};
  my @all = @{$copy->get_all_nodes};

  # Add all ancestors of kept nodes to the keep list.
  my @all_keepers = ();
  foreach my $keeper (@keepers) {
    my $node = $copy->find_node_by_node_id($keeper);
    if (!defined $node) {
      #$node = $copy->find_node_by_name($keeper);
      if (!defined $node) {
	print "Node $keeper not found in tree $tree\n";
	next;
      }
    }

    push @all_keepers, $keeper;

    my $parent = $node->parent;
    while (defined $parent) {
      push @all_keepers, $parent->node_id;
      $parent = $parent->parent;
    }
  }

  my @remove_me = ();
  foreach my $node (@all) {
    push @remove_me, $node unless (grep {$node->node_id == $_} @all_keepers);
  }
  $copy->remove_nodes(\@remove_me);
  return $copy;
}

sub get_minimum_ancestor_from_leaves {
  my $self = shift;
  my $leaf_list = shift;

  my @leaves = @{$leaf_list};

  my $ancestor = shift @leaves;
  while (scalar @leaves > 0) {
    my $node = shift @leaves;
    $ancestor = $ancestor->find_first_shared_ancestor($node);
  }
  return $ancestor;
}

sub build_tree_from_nodes {
  my $self = shift;
  my $node_list = shift;

  #first hash all the nodes by id for fast access
  my %node_hash;
  foreach my $node (@{$node_list}) {
    $node->no_autoload_children;
    $node_hash{$node->node_id} = $node;
  }
  
  #next add children to their parents
  my $root = undef;
  foreach my $node (@{$node_list}) {
    my $parent = $node_hash{$node->_parent_id};
    if($parent) { $parent->add_child($node); } 
    else { $root = $node; }
  }
  return $root;
}



1; # Keep perl happy.
