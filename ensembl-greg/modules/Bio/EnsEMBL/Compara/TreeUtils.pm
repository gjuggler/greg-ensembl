package Bio::EnsEMBL::Compara::TreeUtils;

use strict;
use warnings;

use Bio::TreeIO;
use Bio::EnsEMBL::Compara::LocalMember;
use File::Path;

#
# A grab bag of useful methods for tree manipulations.
#

my $TREE = "Bio::EnsEMBL::Compara::TreeUtils";
my $ALN = "Bio::EnsEMBL::Compara::AlignUtils";
my $COMPARA = "Bio::EnsEMBL::Compara::ComparaUtils";

my $TREEI = "Bio::Tree::TreeI";
my $NSET = "Bio::EnsEMBL::Compara::NestedSet";
my $PTREE = "Bio::EnsEMBL::Compara::ProteinTree";

sub species_count {
  my $class = shift;
  my $tree = shift;

  my $species_hash;
  map {$species_hash->{$_->taxon_id} = 1} $tree->leaves;

  return scalar(keys %$species_hash);
}

sub is_polytomy {
  my $class = shift;
  my $tree = shift;
  
  return 1 if ($tree->get_child_count > 2);
  return 0;
}

sub copy_tree {
  my $self = shift;
  my $tree = shift;

  my $node_hash;
  my $root;

  if ($tree->isa("Bio::EnsEMBL::Compara::NestedSet")) {
    foreach my $node ($tree->nodes) {
      my $cls = ref($node);
      $node_hash->{$node} = new $cls;
    }
    
    foreach my $node ($tree->nodes) {
      my $new_n = $node_hash->{$node};
      
      # Need to copy each tag value over.
      $node->_load_tags if ($node->can('_load_tags'));
      $new_n->{_tags} = {};
      my $tags = $node->{_tags};
      foreach my $key (keys %$tags) {
        $new_n->{$key} = $tags->{$key};
      }

      # Add node ID, name, and distance.
      $new_n->node_id($node->node_id);
      $new_n->name($node->name);
      $new_n->distance_to_parent($node->distance_to_parent);
      if ($node->isa("Bio::EnsEMBL::Compara::Member")) {
        $new_n->source_name($node->source_name);
        $new_n->stable_id($node->stable_id);
        $new_n->dbID($node->dbID);
        $new_n->sequence_id($node->sequence_id);
        if ($node->adaptor && $node->genome_db_id) {
          $new_n->genome_db($node->genome_db);
          $new_n->genome_db_id($node->genome_db_id);
        }
        $new_n->taxon_id($node->taxon_id) if ($node->can('taxon_id'));
        $new_n->{_tags}->{taxon_id} = $node->taxon_id if ($node->can('taxon_id'));
        #$new_n->{core_transcript} = $node->{core_transcript};
        $new_n->{core_gene} = $node->{core_gene};
        $new_n->gene_member($node->gene_member);
        $new_n->gene_member_id($node->gene_member_id);
        $new_n->{_adaptor} = $node->{_adaptor};
      } elsif ($node->isa("Bio::EnsEMBL::Compara::ProteinTree")) {
        # You'd think we should just use $node->taxon_id to set the new taxon ID here,
        # but Compara keeps us on our toes by instead requiring us to manually set a tagvalue.
        # W.T.F.
        #$new_n->taxon_id($node->taxon_id) if ($node->can('taxon_id') && $new_n->can('taxon_id'));
        #$new_n->store_tag('taxon_id',$node->taxon_id) if ($node->can('taxon_id')); 
      }
      
      if (defined $node->parent) {
        my $parent = $node_hash->{$node->parent};
        if (defined $parent) {
          $parent->add_child($new_n);
        } else {
          $root = $new_n;
        }
      } else {
        $root = $new_n;
      }
      
      # This is kind of important -- if we don't call no_autoload_children,
      # then the NCBITaxonomy tree will try to load subspecies into the NestedSet
      # object upon copying. Not nice!!!
      $new_n->no_autoload_children;
      
      $new_n->adaptor($node->adaptor);
    }
  } else {
    die("copy_tree not implemented for non- NestedSets!");
  }

  return $root;
}

sub robinson_foulds_dist {
  # Setup steps:
  # 1) Download http://hashrf.googlecode.com/files/hashrf-6.0.1.tgz . Extract to directory.
  # 2) >./configure
  # 3) >make
  # 4) >cp hashrf ~/bin

  my $class = shift;
  my $tree_a = shift;
  my $tree_b = shift;
  my $params = shift;

  my $temp_dir = $params->{temp_dir};
  if (!defined $temp_dir) {
    $temp_dir = '/tmp/rfdist/';
    rmtree([$temp_dir]);
    mkpath([$temp_dir]);
  }

  my $string = "";
  $string .= $TREE->to_newick($tree_a);
  $string .= "\n";
  $string .= $TREE->to_newick($tree_b);
  $string .= "\n";

  my $file_in = $temp_dir."trees.tre";
  open(OUT,">$file_in");
  print OUT $string;
  close(OUT);
  
  my $file_out = $temp_dir."result.rf";
  my $cmd = "hashrf $file_in 0 -o $file_out ";
  if ($params->{weighted_dist}) {
    $cmd .= " -w";
  }

  print $cmd."\n";
  system($cmd);

  open(IN,"$file_out");
  while (<IN>) {
    print $_;
  }
  close(IN);
  
}

sub k_tree_dist {
  # 1) Get http://molevol.cmima.csic.es/castresana/Ktreedist/Ktreedist_v1.tar.gz
  # 2) > cp Ktreedist.pl ~/bin

  my $class = shift;
  my $tree_a = shift;
  my $tree_b = shift;
  my $params = shift;

  my $temp_dir = $params->{temp_dir};
  if (!defined $temp_dir) {
    $temp_dir = '/tmp/ktreedist/';
    rmtree([$temp_dir]);
    mkpath([$temp_dir]);
  }

  my $file_a = $temp_dir."tree_a.nh";
  my $file_b = $temp_dir."tree_b.nh";

  $class->to_file($tree_a,$file_a);
  $class->to_file($tree_b,$file_b);

  my $cmd = "Ktreedist.pl -rt $file_a -ct $file_b -r";
  system($cmd);
  
}

# Root a tree at its midpoint.
sub midpoint_root {
  
  # TODO!!
}

sub unroot {
  my $class = shift;
  my $tree = shift;

  my $n_at_root = scalar @{$tree->children};
  return $tree if ($n_at_root == 3); # No need to unroot.
  return $tree if (scalar $tree->leaves == 2);

  if ($n_at_root == 2) {
    my $new_root_node;
    my $moving_node;
    my ($child_a,$child_b) = @{$tree->children};

    if (scalar @{$child_a->children} == 2) {
      $new_root_node = $child_a;
      $moving_node = $child_b;
    } elsif (scalar @{$child_b->children} == 2) {
      $new_root_node = $child_b;
      $moving_node = $child_a;
    } else {
      $tree->throw("Error unrooting tree!");
    }
    
    my $new_dist = $moving_node->distance_to_parent + $new_root_node->distance_to_parent;
    $moving_node->disavow_parent;
    $new_root_node->add_child($moving_node);
    $moving_node->distance_to_parent($new_dist);
    return $new_root_node;
  } else {
    $tree->throw("Error unrooting tree!");
  }
}

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
  my $format = shift;

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
    $newick = $tree->newick_format($format);
    $newick .= ';' unless ($newick =~ m/;/);
    #warn("Temporary Newick for NSET to TREEI conversion: [$newick]\n");
  }

  my $treein = new Bio::TreeIO
    (-string => $newick,
     -format => 'newick');
  my $treeI = $treein->next_tree;
  $treein->close;
  return $treeI;
}

sub treeI_from_newick {
  my $class = shift;
  my $newick = shift;

  $newick .= ";" if ($newick !~ /;/);
#  $newick =~ s/\s//g;

  # Parse out and store the NHX annotations... because bioperl-live sucks.
#  my $annotation;
#  while ($newick =~ m/[,\(\)]+(.*?)(:.*?)?\[&&nhx:(.*?)\]/gi) {
#    my $hash;
#    my @keyvals = split(/:/,$3);
#    foreach my $keyval (@keyvals) {
#      my ($key,$val) = split(/=/,$keyval);
#      $hash->{$key} = $val;
#    }
#    $annotation->{$1} = $hash if (length $1 > 1);
#  }
  
  # Now remove the NHX annotations from the string.
#  $newick =~ s/\[&&nhx.*?\]//ig;
  
  # Load the tree using TreeIO.
  open(my $fake_fh, "+<", \$newick);
  my $treein = new Bio::TreeIO(-fh => $fake_fh, -format => 'newick');
  my $treeI = $treein->next_tree;
  $treein->close();

  # Put annotations back into nodes.
#  foreach my $node ($treeI->get_nodes) {
#    my $hash = $annotation->{$node->id};
#    if (defined $hash) {
#      foreach my $key (keys %{$hash}) {
#	$node->add_tag_value($key,$hash->{$key});
#      }
#    }
#  }
  
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

  my $tree = $class->from_treeI($class->treeI_from_newick($newick));
  return $tree;
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

  # Kick off the recursive, parallel node adding.
  _add_nodeI_to_node($node,$rootI);
  
  return bless($node,'Bio::EnsEMBL::Compara::NestedSet');
}

# Recursive helper for new_from_TreeI.
sub _add_nodeI_to_node {
  my $node = shift; # Our node object (Compara)
  my $nodeI = shift; # Our nodeI object (BioPerl)
  my $node_id_counter = shift; # Node ID counter

  foreach my $c ($nodeI->each_Descendent) {
    my $child = ref($node)->new;
    
    my $name = $c->id || "";
    $name =~ s/^\s+//;
    $name =~ s/\s+$//;
    
    if ($c->is_Leaf) {
      $child = Bio::EnsEMBL::Compara::LocalMember->new();
      $child->stable_id($name);
      $child->source_name("");
    }

    $child->node_id($node_id_counter++);  # unless ($c->is_Leaf);
    
    # Set name.
    $child->name($name);
    $child->store_tag("name",$name);
    
    # Set branch length.
    $node->add_child($child,$c->branch_length);

    # Add the tags.
    my $tags = $c->{'_tags'};
    foreach my $tag (keys %{$tags}) {
      $tags->{$tag} = $tags->{$tag}[0]
    }
    $child->{'_tags'} = $tags;

    # Recurse.
    _add_nodeI_to_node($child, $c, $node_id_counter);
  }
}

# Scales all the branches in a tree by a given factor.
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

sub mean_path {
  my $class = shift;
  my $tree  = shift;

  my $dist  = 0;
  map { $dist += $class->dist_to_root($_); } $tree->leaves;
  $dist = $dist / scalar( $tree->leaves );
  return sprintf "%.3f", $dist;
}

sub dist_to_root {
  my $class = shift;
  my $leaf  = shift;

  my $d = $leaf->distance_to_parent;
  my $p = $leaf->parent;
  while ($p) {
    $d += $p->distance_to_parent;
    $p = $p->parent;
  }
  return $d;
}

sub scale_mean_to {
  my $class = shift;
  my $tree = shift;
  my $new_mean = shift;

  my $mean_dist = $class->mean_path($tree);
  my $scale_factor = $new_mean / $mean_dist;
  return $class->scale($tree,$scale_factor);
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

  my $treeI = $tree;
  if (!$tree->isa($TREEI)) {
    $treeI = $class->to_treeI($tree);
  }
  
  my $max_dist = 0;
  foreach my $node ($treeI->get_leaf_nodes()) {
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

  my $treeI = $tree;
  if (!$tree->isa($TREEI)) {
    $treeI = $class->to_treeI($tree);
  }

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
  my $params = shift || {};

  my $ref = ref $tree;
  if ($ref =~ /Bio::Tree/i) {
    use Bio::TreeIO;
    my $string = "";
    open(my $fake_fh, "+>", \$string);
    my $out = Bio::TreeIO->new(-fh => $fake_fh,
			       -format => "newick");
    $out->write_tree($tree);
    $out->close();
    $string =~ s/\n//g; # Remove newlines.
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

  die("Need to provide a tree and file!") unless ($tree && $out_file);

  my $newick = $class->to_newick($tree);
  if ($params->{'node_ids'}) {
    $newick = $tree->newick_format('int_node_id');
  } elsif (defined $params->{nhx_format} && $params->{nhx_format} == 1) {
    $newick = $tree->nhx_format;
  }

  open(OUT,">".$out_file);
  print OUT $newick;
  #print $newick."\n";
  close(OUT);
  return $out_file;
}

sub remove_elbows {
  my $class = shift;
  my $tree = shift;

  my @del_nodes = ();
  foreach my $node ($tree->nodes) {
    my @children = @{$node->children};
    if (scalar @children == 1 && defined $node->parent) {
      print "Deleting: ".$node->name."\n";
      push @del_nodes, $node;
    }

    if (defined $node->parent && !defined $node->parent->parent) {
      my @p_children = @{$node->parent->children};
      print scalar @p_children."\n";
      if (scalar @p_children == 1) {
	print "DELETING PARENT\n";
#	push @del_nodes, $node;
      }
    }
  }

  foreach my $del_me (@del_nodes) {
    $tree->remove_nodes([$del_me]);
  }

  return $tree;
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

sub remove_members_by_node_id {
  my $class = shift;
  my $tree = shift;
  my $node_list = shift;

  my @node_ids = @{$node_list};
  my %ids_hash;
  map {$ids_hash{$_}=1} @node_ids;

  #print "Removing specified nodes from tree: @node_ids \n";
  #print "  Before:" . scalar($tree->leaves) . " leaves\n";

  foreach my $node ($tree->nodes) {
    if (exists $ids_hash{$node->node_id}) {
      $TREE->delete_lineage($tree,$node);
      $tree->minimize_tree();
    }
  }
    #print "  After:" . scalar($tree->leaves) . " leaves\n";
  return $tree;
}

sub keep_members_by_method_call {
  my $class = shift;
  my $tree = shift;
  my $value_arrayref = shift;
  my $method = shift || 'taxon_id';

  my %value_hash;
  map {$value_hash{$_}=1} @{$value_arrayref};

  # Take the complement of the "delete me" set and extract a subtree.
  my @keep_me;
  foreach my $leaf ($tree->leaves) {
#    print $leaf->$method()."\n";
    if (exists $value_hash{$leaf->$method()}) {
      push @keep_me, $leaf;
    } else {
    }
  }
  if (scalar @keep_me > 0) {
    $tree = $class->extract_subtree_from_leaf_objects($tree,\@keep_me);
  } else {
    $tree = new $tree;
  }
  return $tree;

}



sub remove_members_by_method_call {
  my $class = shift;
  my $tree = shift;
  my $species_arrayref = shift;
  my $method = shift || 'taxon_id';

  my @tax_ids = @{$species_arrayref};
  # Turn the arrayref into a hash.
  my %tax_hash;
  map {$tax_hash{$_}=1} @tax_ids;

#  printf "Pruning leaves by $method: %-40.40s ...\n", join(",",sort(@tax_ids));
#  print " > Before: " . scalar($tree->leaves) . "\n";

  # Take the complement of the "delete me" set and extract a subtree.
  my @keep_me;
  foreach my $leaf ($tree->leaves) {
    if (!exists $tax_hash{$leaf->$method()}) {
      #print ">". $leaf->name." ".$leaf->node_id."\n";
      push @keep_me, $leaf->node_id;
    } else {
      #print "Exists: ".$leaf->$method()."\n";
    }
  }
  if (scalar @keep_me > 0) {
    $tree = $class->extract_subtree_from_leaves($tree,\@keep_me);
#    print " > After: " . scalar($tree->leaves) . "\n";
  }
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

#  print "Pruning leaves with taxon_ids: ". join(",",sort(@tax_ids))."\n";
#  print "  Before: " . scalar($tree->leaves) . "\n";

  # Take the complement of the "delete me" set and extract a subtree.
  my @keep_me;
  foreach my $leaf ($tree->leaves) {
    if (!exists $tax_hash{$leaf->taxon_id}) {
      push @keep_me, $leaf;
      #$class->delete_lineage($tree,$leaf);
      #$tree->minimize_tree();
    }
  }

  $tree = $class->extract_subtree_from_leaf_objects($tree,\@keep_me);
  
#  print "  After: " . scalar($tree->leaves) . "\n";
  return $tree;
}

sub get_leaves_for_species {
  my $class = shift;
  my $tree = shift;
  my $taxon_ids_obj = shift;

  my @taxon_ids = @{$taxon_ids_obj};

  my @leaves = ();
  foreach my $leaf ($tree->leaves) {
    if (grep {$leaf->taxon_id == $_} @taxon_ids) {
#      printf "%s %s\n",$leaf->taxon_id,$leaf->stable_id;
      push @leaves, $leaf;
    }
  }
  return @leaves;
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
    my $tx_id = shift;
    my $leaf = shift;
    my $db = shift;

    my $leaf_taxon = $leaf->taxon_id;
    my $sth = $db->dbc->prepare(qq^
				SELECT * from ncbi_taxa_node p, ncbi_taxa_node c
				WHERE p.taxon_id=$tx_id AND c.taxon_id=$leaf_taxon AND
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
    next if (is_leaf_within_taxon($taxon_id,$leaf,$db));
    
    print "  -> Deleting leaf: " . $leaf->stable_id . "  " . $leaf->taxon->binomial . "\n";
    $class->delete_lineage($tree,$leaf);
  }

  print "  After: " . scalar(@{$tree->get_all_leaves}) . " leaves.\n";

  return $tree;
}

# Returns the subtree of the input tree defined by sequnces which are present in the alignment.
sub extract_subtree_from_aln {
  my $class = shift;
  my $tree = shift;
  my $aln = shift;

  my $keep_name_hash;
  map {$keep_name_hash->{$_->id} = 1} $aln->each_seq;

  my @keep_node_ids;
  foreach my $leaf ($tree->leaves) {
    if ($keep_name_hash->{$leaf->name}) {
      push @keep_node_ids, $leaf->node_id;
    }
  }

  return $class->extract_subtree_from_leaves($tree, \@keep_node_ids);
}

sub extract_subtree_from_names {
  my $class = shift;
  my $tree = shift;
  my $name_arrayref = shift;
  my $return_copy = shift;

  foreach my $node ($tree->nodes) {
    if (!defined $node->node_id || $node->node_id eq '') {
      print "No node_id for ".$node->name."\n";
    } else {
#      print "  ".$node->node_id."\n";
    }
  }

  # Just go through and get the node_ids corresponding to the name input.
  my @names = @$name_arrayref;
  my @node_ids;
  my @nodes = $tree->nodes;
  foreach my $name (@names) {
    my ($match) = grep {$_->name eq $name} @nodes;
    if ($match) {
      push @node_ids, $match->node_id;
    } else {
      print $tree->newick_format."\n";
      die("Can't find node with name $name\n");
    }
  }
  #print "keeping @node_ids\n";
  return $class->extract_subtree_from_leaves($tree, \@node_ids, $return_copy);
}

# Extracts a subtree given a ProteinTree and an arrayref of node_ids.
# @updated GJ 2010-03-27 : Overhaul for Gorilla project, lots of stuff was broken here.
# @updated GJ 2009-01-14 : Smarter version, uses the inner method NestedSetAdaptor->_build_tree_from_nodes.
# @created GJ 2009-01-12
sub extract_subtree_from_leaf_objects {
  my $class = shift;
  my $tree = shift;
  my $leaf_objs = shift;

  die("Object not a NestedSet!") unless ($tree->isa("Bio::EnsEMBL::Compara::NestedSet"));

  my @keepers = @{$leaf_objs};

  # Add all ancestors of kept nodes to the keep list.
  my %keepers_hash;
  foreach my $keeper (@keepers) {
    $keepers_hash{$keeper} = 1;

    my $parent = $keeper->parent;
    while (defined $parent) {
      $keepers_hash{$parent} = 1;
      $parent = $parent->parent;
    }
  }

  # Remove all nodes NOT in the keepers hash.
  my @remove_me = ();
  foreach my $node ($tree->nodes) {
    if ($keepers_hash{$node}) {
      #print $node->node_id."\n";
    } else {
      push @remove_me, $node;
    }
  }
  
  $tree = $tree->remove_nodes(\@remove_me);
  #$tree = $tree->minimize_tree;
  return $tree;
}

sub extract_subtree_from_node_ids {
  my $class = shift;
  my $tree = shift;
  my $node_id_arrayref = shift;
  
  return $class->extract_subtree_from_leaves($tree, $node_id_arrayref);
}



# Extracts a subtree given a ProteinTree and an arrayref of node_ids.
# @updated GJ 2010-03-27 : Overhaul for Gorilla project, lots of stuff was broken here.
# @updated GJ 2009-01-14 : Smarter version, uses the inner method NestedSetAdaptor->_build_tree_from_nodes.
# @created GJ 2009-01-12
sub extract_subtree_from_leaves {
  my $class = shift;
  my $tree = shift;
  my $node_ids = shift;	# Array ref of node_ids.
  my $return_copy = shift;

  $return_copy = 1 unless (defined $return_copy);

  die("Object not a NestedSet!") unless ($tree->isa("Bio::EnsEMBL::Compara::NestedSet"));

  my @keepers = @{$node_ids};

  if (scalar @keepers == 1) {
    # Special case: a single keeper node.
    my $keep_id = $keepers[0];
    my $node = $tree->find_node_by_node_id($keep_id);
    $node = $tree->find_leaf_by_name($keep_id) unless (defined $node);
    return $node;
  }

  # Add all ancestors of kept nodes to the keep list.
  my %keepers_node_id_hash;
  foreach my $keep_id (@keepers) {
    my $node = $tree->find_node_by_node_id($keep_id);
    $node = $tree->find_leaf_by_name($keep_id) unless (defined $node);
    if (!defined $node) {
      printf "Node [%s] not found in tree [%s]!\n", $keep_id, $tree->newick_format;
      return undef;
    }

    $keepers_node_id_hash{$node->node_id} = 1;

    my $parent = $node->parent;
    while (defined $parent) {
      $keepers_node_id_hash{$parent->node_id} = 1;
      $parent = $parent->parent;
    }
  }

  # Remove all nodes NOT in the keepers hash.
  my @remove_me = ();
  foreach my $node ($tree->nodes) {
    if ($keepers_node_id_hash{$node->node_id}) {
      #print $node->node_id."\n";
    } else {
      push @remove_me, $node;
    }
  }

  if ($return_copy) {
    my $copy = $class->copy_tree($tree);
    
    my @nodes = $copy->nodes; # Store a cache of all nodes to use for grepping in the sub.
    
    # Get the list of all nodes to remove by finding them in the copied tree. We use the node_id
    # as the connection between the original and copied tree.
    my @copy_remove_me = map {
      my $target;
      my $node_id = $_->node_id;
      ($target) = grep {$_->node_id == $node_id} @nodes;
      die("Node ".$_->node_id." ".$_->name." not found in copied tree!") unless (defined $target);
      $target;
    } @remove_me;
    
    $copy = $copy->remove_nodes(\@copy_remove_me);
    $copy = $copy->minimize_tree if (defined $copy);
    return $copy;
  } else {
    $tree = $tree->remove_nodes(\@remove_me);
    $tree = $tree->minimize_tree if (defined $tree);
    return $tree;
  }
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

sub has_one_to_one_orthology {
  my $class = shift;
  my $tree = shift;
  my $taxon_id_arrayref = shift;

  my @species_list = @$taxon_id_arrayref;

  my @keeper_leaves = $class->get_leaves_for_species($tree,\@species_list);
  my @keeper_ids = map {$_->node_id} @keeper_leaves;
  my $pruned_tree = $class->extract_subtree_from_leaves($tree,\@keeper_ids);

#  print $pruned_tree->newick_format."\n";
#  print map {$_->stable_id." "} @keeper_leaves;
#  print "\n";

  # Test whether it's a good 1-1-1 orthology.
  my $is_good_tree = 1;
  my %keeper_hash;
  map {$keeper_hash{$_->taxon_id}=1} @keeper_leaves;
  
  # Not good if we're missing any species.
  map {$is_good_tree = 0 if (!defined $keeper_hash{$_})} @species_list;

  # Not good if the species and leaf counts don't match -- this means we have a duplication somewhere.
  $is_good_tree = 0 if ($#keeper_leaves != $#species_list);  

  if (!$is_good_tree) {
    return 0;
  } else {
    return 1;
  }
}

sub subtree {
  my $class = shift;
  my $tree = shift;
  my $stable_id_arrayref = shift;

  $tree = $tree->copy;
  my $subtree = $class->keep_members_by_method_call($tree,$stable_id_arrayref,'stable_id');
  return $subtree;
}

sub translate_ids {
  my $class = shift;
  my $tree = shift;
  my $map = shift;
  my $params = shift;

  $tree = $class->copy_tree($tree); # Important!!
  
  my $ensure_unique = 1;
  $ensure_unique = $params->{ensure_unique} if (defined $params->{ensure_unique});

  my $used_ids;
  foreach my $node ($tree->nodes) {
    my $name = $node->name;
    my $node_id = $node->node_id;

    my $new_name;
    $new_name = $map->{$node_id} if (defined $node_id && !defined $new_name); # Defined by node ID.
    $new_name = $map->{$name} if (defined $name && !defined $new_name); # Defined by node name.
    
    if (defined $new_name) {
      if ($ensure_unique) {
        while ($used_ids->{$new_name}) {
          $new_name =~ m/_(\d+)$/;
          my $num = $1;
          #print "ID in use: [$new_id]\n";
          $new_name =~ s/_\d+$//;
          my $new_int = $num + 1;
          $new_name .= "_$new_int";
          #print "going to use [$new_id]\n";
        }
      }

      $node->name($new_name);
      $used_ids->{$new_name} = 1;
    }
  }

  return $tree;
}

sub name_to_stable_id {
  my $class = shift;
  my $tree = shift;

  my $copy = $class->copy_tree($tree);

  map {$_->stable_id($_->name)} $copy->leaves;
  return $copy;
}

sub pretty_print {
  my $class = shift;
  my $tree = shift;
  my $params = shift;

  my $treeI = $class->to_treeI($tree);
  if ($treeI->can('ascii')) {
    print $treeI->ascii(0,0,0);
  }
}

sub transfer_branchlengths {
  my $class = shift;
  my $source = shift;
  my $target = shift;

  my @s_nodes = $source->nodes;
  my @t_nodes = $target->nodes;
  foreach my $t_node (@t_nodes) {
    # Find the equivalent source node.
    my ($s_node) = grep {$_->enclosed_leaves_string eq $t_node->enclosed_leaves_string} @s_nodes;
    die("No source node!") unless (defined $s_node);
    
    $t_node->branch_length($s_node->branch_length);
  }
}

1; # Keep perl happy.
