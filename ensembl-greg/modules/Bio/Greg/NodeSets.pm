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

our @ISA = qw(Bio::EnsEMBL::Hive::Process);

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
  my( $self) = @_;
  

  # Load up the Compara DBAdaptor.
  if ($self->dba) {
    $dba = $self->dba;
  } else {
    $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-DBCONN=>$self->db->dbc);
  }
  $pta = $dba->get_ProteinTreeAdaptor;
  
  ### DEFAULT PARAMETERS ###
  $params = {
    flow_node_set => 1,
    debug => 0
  };
  
  # For aminof, codonf, and freqtype, see the SLR readme.txt for more info.

  #########################

  print "TEMP: ".$self->worker_temp_directory."\n";
  print "PARAMS: ".$self->parameters."\n";
  print "INPUT_ID: ".$self->input_id."\n";

  # Fetch parameters from the two possible locations. Input_id takes precedence!

  $params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_string($params,$self->parameters);
  $params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_string($params,$self->input_id);

  #########################

  $self->check_if_exit_cleanly;

  # Deal with alternate table inputs.
  $pta->protein_tree_member($params->{'input_table'}) if (defined $params->{'input_table'});

  #
  # Load the tree. This logic has been delegated to ComparaUtils.
  #
  $tree = $pta->fetch_node_by_node_id($params->{'node_id'});
  if (!defined $tree) {
    $tree = $pta->fetch_node_by_node_id($params->{'protein_tree_id'});
  }
}

sub run {
  my $self = shift;

  #print $tree->newick_format."\n";
  $tree->print_tree;

  print "Adding human gene subtrees...\n";
  $self->add_human_gene_subtrees_to_hash($tree,\%node_set_hash,1);

  print "Adding small dups...\n";
  $self->add_small_balanced_dups_to_hash($tree,\%node_set_hash,2);

  print "Adding large dups...\n";
  $self->add_large_balanced_dups_to_hash($tree,\%node_set_hash,3);

  print "Adding overlapping dups...\n";
  $self->add_overlapping_dups_to_hash($tree,\%node_set_hash,4);

  print "Adding small subtrees...\n";
  $self->add_all_subtrees_to_hash($tree,\%node_set_hash,5,"small");

  print "Adding med subtrees...\n";
  $self->add_all_subtrees_to_hash($tree,\%node_set_hash,6,"med");

  print "Adding large subtrees...\n";
  $self->add_all_subtrees_to_hash($tree,\%node_set_hash,7,"large");

  print "Adding superfamily subtrees...\n";
  $self->add_all_subtrees_to_hash($tree,\%node_set_hash,8,"superfamily");

  print "Adding Tim superfamily subtrees...\n";
  $self->add_all_subtrees_to_hash($tree,\%node_set_hash,9,"superfamily_tim",\&is_tree_tim_worthy);

  my $sql = "INSERT IGNORE INTO node_set (node_set_id,name,allows_overlap) values (?,?,?)";
  my $sth = $dba->dbc->prepare($sql);
  $sth->execute(1,"Human Gene Subtrees",0);
  $sth->execute(2,"Small Balanced Duplications",0);
  $sth->execute(3,"Large Balanced Duplications",0);
  $sth->execute(4,"Overlapping Duplications",1);
  $sth->execute(5,"Small Subtrees",0);
  $sth->execute(6,"Med Subtrees",0);
  $sth->execute(7,"Large Subtrees",0);
  $sth->execute(8,"Superfamily Subtrees",0);
  $sth->execute(9,"Tim Subtrees",0);
  $sth->finish;
}

sub add_all_subtrees_to_hash {
  my $self = shift;
  my $tree = shift;
  my $hash = shift;
  my $node_set_id = shift;
  my $size = shift;
  
  my $ref = \&is_tree_small;
  $ref = \&is_tree_med if ($size eq "med");
  $ref = \&is_tree_large if ($size eq "large");
  $ref = \&does_parent_have_superfamily_children if ($size eq "superfamily");
  $ref = \&does_parent_have_tim_children if ($size eq 'superfamily_tim');
  my @root_node_ids = $self->get_smallest_subtrees_from_node($tree,$ref);
  
  $hash->{$node_set_id} = \@root_node_ids;
}

sub add_human_gene_subtrees_to_hash {
  my $self = shift;
  my $tree = shift;
  my $hash = shift;
  my $node_set_id = shift;
  
  my $ref = \&isa_human_gene_subtree;
  my @root_node_ids = $self->get_smallest_subtrees_from_node($tree,$ref);
  
  $hash->{$node_set_id} = \@root_node_ids;
}

sub add_small_balanced_dups_to_hash {
  my $self = shift;
  my $tree = shift;
  my $hash = shift;
  my $node_set_id = shift;
  
  my $ref = \&isa_balanced_good_duplication;
  my @keeper_ids = $self->get_smallest_subtrees_from_node($tree,$ref);

  $hash->{$node_set_id} = \@keeper_ids;
}

sub add_large_balanced_dups_to_hash {
  my $self = shift;
  my $tree = shift;
  my $hash = shift;
  my $node_set_id = shift;
  
  my $ref = \&isa_balanced_good_duplication;
  my @keeper_ids = $self->get_largest_subtrees_from_node($tree,$ref);

  $hash->{$node_set_id} = \@keeper_ids;
}

sub add_overlapping_dups_to_hash {
  my $self = shift;
  my $tree = shift;
  my $hash = shift;
  my $node_set_id = shift;
  
  my @keeper_ids = ();
  foreach my $node ($tree->nodes) {
    print "Added ".$node->node_id." to list!\n";
    push @keeper_ids,$node->node_id if (is_balanced_dup_node($node));
  }

  $hash->{$node_set_id} = \@keeper_ids;
}


# A generic method for gathering a set of small, non-overlapping subtrees.
sub get_smallest_subtrees_from_node {
  my $self = shift;
  my $node = shift;
  my $inclusion_function = shift;

  my @children = @{$node->children};
  my @node_list = ();
  foreach my $child (@children) {
    push @node_list, $self->get_smallest_subtrees_from_node($child,$inclusion_function);
  }

  # If none of the children matched, then push ourselves onto the list if applicable.
  if (scalar(@node_list) == 0) {
    if ($inclusion_function->($node)) {
      push @node_list, $node->node_id;
      print "Added ".$node->node_id." to list!\n";
    }
  }

  return @node_list;
}

# A generic method for gathering a set of large, non-overlapping subtrees.
sub get_largest_subtrees_from_node {
  my $self = shift;
  my $node = shift;
  my $inclusion_function = shift;

  if ($inclusion_function->($node)) {
    # Return this node immediately if it satisfies the conditions.
    print "Added ".$node->node_id." to list!\n";
    return ($node->node_id);
  }

  # Recurse through child nodes to see if they satisfy the conditions.
  my @children = @{$node->children};
  my @node_list = ();
  foreach my $child (@children) {
    push @node_list, $self->get_largest_subtrees_from_node($child,$inclusion_function);
  }
  return @node_list;
}

sub isa_human_gene_subtree {
  my $node = shift;

  return 1 if (is_tree_big_enough($node) && is_contains_human_genes($node));
  return 0;
}

sub isa_balanced_good_duplication {
  my $node = shift;

  return 1 if (is_balanced_dup_node($node,6));
  return 0;
}

sub is_duplication {
  my $node = shift;

  my $is_dup = 0;
  my $val = $node->get_tagvalue("Duplication");
  $is_dup = 1 if (defined($val) && $val ne '' && $val > 0);
  
  my $dubious = 0;
  $dubious = 1 if ($node->get_tagvalue("dubious_duplication"));

  return 1 if ($is_dup && !$dubious);
  return 0;
};

sub is_contains_human_genes {
  my $node = shift;
  my @leaves = $node->leaves;
  my @hum_genes = grep {$_->taxon_id == 9606} @leaves;
  return 1 if (scalar(@hum_genes) > 0);
  return 0;
}

sub is_tree_big_enough {
  my $node = shift;
  return is_tree_small($node);
}

sub is_tree_small {
  my $node = shift;
  return (scalar($node->leaves) >= 6);
}

sub is_tree_med {
  my $node = shift;
  return (scalar($node->leaves) >= 10);
}

sub is_tree_large {
  my $node = shift;
  return (scalar($node->leaves) >= 20);
}

sub does_parent_have_superfamily_children {
  my $node = shift;
  return generic_parent_has_good_children($node,\&does_tree_have_superfamily_coverage);
}

sub does_parent_have_tim_children {
  my $node = shift;
  return generic_parent_has_good_children($node,\&does_tree_have_tim_family_coverage);
}

sub generic_parent_has_good_children {
  my $node = shift;
  my $inclusion_function = shift;

  # 1. check that this tree has superfamily coverage.
  # 2. check that our sister node has superfamily coverage.
  # 3. if (1) and (2) are met, return true.

  my $tree = $node;
  my $parent = $node->parent;
  my @parents_children = @{$parent->children};
  #print $parent->node_id."\n";
  #print length(@parents_children)."\n";
  if ($parent->node_id == 1) {
    return $inclusion_function->($node);
  }

  my $sister;
  foreach my $ch (@parents_children) {
    $sister = $ch if ($ch->node_id != $tree->node_id);
  }
  
  if ($inclusion_function->($node) && $inclusion_function->($sister)) {
    return 1;
  }
  return 0;
}

sub does_tree_have_superfamily_coverage {
  my $tree = shift;

  # 1. require 2 of 4 mammalian families to be represented (Laurasiatheria, Primates, Glires, Afrotheria)
  # 2. require 1 of 3 outgroups to be represented (Opossum, Chicken, Tetraodon)

  my %mamm_hash = ('314145' => 'Laurasiatheria',
		   '9443' => 'Primates',
		   '314147' => 'Glires',
		   '311790' => 'Afrotheria');
  my %outgroup_hash = ('9031' => 'Chicken',
		       '13616' => 'Opossum',
		       '99883' => 'Tetraodon');

  my %fams_found;
  my %outgroups_found;

  my @nodes = $tree->nodes;
  foreach my $node (@nodes) {
    my $tax_id = $node->taxon_id;
    $fams_found{$tax_id}=1 if (exists $mamm_hash{$tax_id});
    $outgroups_found{$tax_id}=1 if (exists $outgroup_hash{$tax_id});
  }

  my $num_mamm_fams = scalar(keys %fams_found);
  my $num_outgroups = scalar(keys %outgroups_found);

  printf "Node ID: %d  Mammal coverage: %d/%d  Outgroups: %d/%d\n",$tree->node_id,$num_mamm_fams,4,$num_outgroups,3;
  return 1 if ($num_mamm_fams >= 2 && $num_outgroups >= 1);
  return 0;
}

sub does_tree_have_tim_family_coverage {
  my $tree = shift;

  # 1. require 2 of 4 mammalian families to be represented (Laurasiatheria, Primates, Glires, Afrotheria)

  my %mamm_hash = ('314145' => 'Laurasiatheria',
		   '9443' => 'Primates',
		   '314147' => 'Glires',
		   '311790' => 'Afrotheria');
  my %fams_found;

  my @nodes = $tree->nodes;
  foreach my $node (@nodes) {
    my $tax_id = $node->taxon_id;
    $fams_found{$tax_id}=1 if (exists $mamm_hash{$tax_id});
  }

  my $num_mamm_fams = scalar(keys %fams_found);

  #printf "Node ID: %d  Mammal coverage: %d/%d\n",$tree->node_id,$num_mamm_fams,4;
  return 1 if ($num_mamm_fams >= 2);
  return 0;
}

sub is_tree_tim_worthy {
  my $tree = shift;
  my $node_id = shift;
  
  my $subtree = $tree->find_node_by_node_id($node_id);

  my @nodes = $subtree->nodes;
  my @leaves = $subtree->leaves;

  return 0 if (scalar(@leaves) < 8);
  
  my $dup_count = 0;
  foreach my $node (@nodes) {
    $dup_count++ if (is_duplication($node));
    return 0 if (is_duplication($node) && scalar($node->leaves) > 2);
  }
#  $subtree->print_tree;
  return 1 if ($dup_count <= 4);
  return 0;
}


sub is_balanced_dup_node {
  my $node = shift;
  my $minimum_size = shift;

  my $dup = $node->get_tagvalue("Duplication");
  my $dubious = $node->get_tagvalue("dubious_duplication");
  if ($dup && !$dubious) {
    # Ensure that the subtrees are balanced.
    my @children = @{$node->children};
    if (scalar(@children) == 2) {
      #print "two children!\n";
      my $a = $children[0];
      my $b = $children[1];
      
      my @a_leaves = $a->leaves;
      my @b_leaves = $b->leaves;
      
      my $ratio = scalar(@a_leaves) / scalar(@b_leaves);
      $ratio = 1 / $ratio if ($ratio > 1);
      #print "Left-right ratio: $ratio\n";
      
      my $n_a = scalar @a_leaves;
      my $n_b = scalar @b_leaves;
      
      if ($ratio > 0.5 &&
	  $n_a >= $minimum_size &&
	  $n_b > $minimum_size ) {
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

  foreach my $key (sort keys %node_set_hash) {
    my $sql = "INSERT IGNORE INTO node_set_member (node_set_id,node_id) values (?,?)";
    my $sth = $dba->dbc->prepare($sql);
    print "NODE SET ID: $key\n";
    my @ids = @{$node_set_hash{$key}};
    foreach my $id (@ids) {
      print "  -> Inserting node set member $key -> $id\n";
      $sth->execute($key,$id);
      sleep(0.1);
      my $tree = $pta->fetch_node_by_node_id($id);
      if ($params->{'debug'}) {
	foreach my $leaf ($tree->leaves) {
	  if ($leaf->taxon_id == 9606) {
	    print "  ".$leaf->stable_id."\n";
	    #print "  -> ".$leaf->gene->external_name."\n";
	  }
	}
      }

      if ($key == $params->{'flow_node_set'}) {
	my $output_id = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string({node_id => $id});
	$self->dataflow_output_id($output_id,1);
	print "  -> Flowing node $output_id\n";
	if ($params->{'flow_parent_and_children'}) {
	  foreach my $child (@{$tree->children}) {
	    my $output_id = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string({node_id => $child->node_id});
	    $self->dataflow_output_id($output_id,1);
	    print "  -> Flowing child $output_id\n";
	  }
	}
      }
    }
    $sth->finish;
  }
}

sub DESTROY {
    my $self = shift;
    $tree->release_tree if ($tree);
    $tree = undef;
    $self->SUPER::DESTROY if $self->can("SUPER::DESTROY");
}

1;
