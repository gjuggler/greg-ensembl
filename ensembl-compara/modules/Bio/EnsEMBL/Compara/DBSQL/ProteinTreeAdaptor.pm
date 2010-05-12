=head1 NAME

ProteinTreeAdaptor - DESCRIPTION of Object

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONTACT

  Contact Jessica Severin on implemetation/design detail: jessica@ebi.ac.uk
  Contact Ewan Birney on EnsEMBL in general: birney@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Compara::DBSQL::ProteinTreeAdaptor;

use strict;
use Bio::EnsEMBL::Compara::ProteinTree;
use Bio::EnsEMBL::Compara::AlignedMember;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Compara::DBSQL::NestedSetAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor;
our @ISA = qw(Bio::EnsEMBL::Compara::DBSQL::NestedSetAdaptor);

###########################
# FETCH methods
###########################


sub store_tree_in_treeset {
    my $self = shift;
    my $tree = shift;
    my $treeset_name = shift;

    my $treeset = $self->fetch_treeset_by_name($treeset_name);

    throw("Treeset $treeset_name not found!\n") unless (defined $treeset);

    $treeset->add_child($tree,0);
    $self->store($treeset);
}

sub delete_all_from_treeset {
    my $self = shift;
    my $treeset_name = shift;
    my $one = shift;
    my $two = shift;
    my $three = shift;

    unless ($one == 1 and $two == 2 and $three == 3) {
	throw("Did you really want to delete ALL trees from $treeset_name?\n");
    }
    my $treeset = $self->fetch_treeset_by_name($treeset_name);

    if ($treeset) {
	foreach my $node ($treeset->each_child) {
	    #print "Deleting node: ".$node->node_id."\n";
	    $treeset->delete_node_and_under($node);
	}
    }

}

=head2 fetch_by_Member_root_id

  Arg[1]     : Bio::EnsEMBL::Compara::Member
  Arg[2]     : [optional] int clusterset_id (def. 1)
  Example    : $protein_tree = $proteintree_adaptor->fetch_by_Member_root_id($member);

  Description: Fetches from the database the protein_tree that contains the
               member. If you give it a clusterset id of 0 this will cause
               the search span across all known clustersets.
  Returntype : Bio::EnsEMBL::Compara::ProteinTree
  Exceptions :
  Caller     :

=cut


sub fetch_by_Member_root_id {
  my ($self, $member, $clusterset_id) = @_;
  $clusterset_id = 1 if ! defined $clusterset_id;

  my $root_id = $self->gene_member_id_is_in_tree($member->member_id);
  return undef unless (defined $root_id);
  my $aligned_member = $self->fetch_AlignedMember_by_member_id_root_id
    (
     $member->get_canonical_peptide_Member->member_id,
     $clusterset_id);
  return undef unless (defined $aligned_member);
  my $node = $aligned_member->subroot;
  return undef unless (defined $node);
  my $protein_tree = $self->fetch_node_by_node_id($node->node_id);

  return $protein_tree;
}


sub fetch_by_gene_Member_root_id {
  my ($self, $member, $clusterset_id) = @_;
  $clusterset_id = 1 if ! defined $clusterset_id;

  my $root_id = $self->gene_member_id_is_in_tree($member->member_id);
  return undef unless (defined $root_id);
  my $protein_tree = $self->fetch_node_by_node_id($root_id);

  return $protein_tree;
}


=head2 fetch_AlignedMember_by_member_id_root_id

  Arg[1]     : int member_id of a peptide member (longest translation)
  Arg[2]     : [optional] int clusterset_id (def. 0)
  Example    :

      my $aligned_member = $proteintree_adaptor->
                            fetch_AlignedMember_by_member_id_root_id
                            (
                             $member->get_canonical_peptide_Member->member_id
                            );

  Description: Fetches from the database the protein_tree that contains the member_id
  Returntype : Bio::EnsEMBL::Compara::AlignedMember
  Exceptions :
  Caller     :

=cut


sub fetch_AlignedMember_by_member_id_root_id {
  my ($self, $member_id, $clusterset_id) = @_;

  my $constraint = "WHERE tm.member_id = $member_id and m.member_id = $member_id";
  $constraint .= " AND t.clusterset_id = $clusterset_id" if($clusterset_id and $clusterset_id>0);
  my $final_clause = "order by tm.node_id desc";
  $self->final_clause($final_clause);
  my ($node) = @{$self->_generic_fetch($constraint)};
  return $node;
}

# This is experimental -- use at your own risk
sub fetch_first_shared_ancestor_indexed {
  my $self = shift;
  my $node1 = shift;
  my $node2 = shift;

  my $root_id1 = $node1->_root_id; # This depends on the new root_id field in the schema
  my $root_id2 = $node2->_root_id;

  return undef unless ($root_id1 eq $root_id2);

  my $left_node_id1 = $node1->left_index;
  my $left_node_id2 = $node2->left_index;

  my $right_node_id1 = $node1->right_index;
  my $right_node_id2 = $node2->right_index;

  my $min_left;
  $min_left = $left_node_id1 if ($left_node_id1 < $left_node_id2);
  $min_left = $left_node_id2 if ($left_node_id2 < $left_node_id1);

  my $max_right;
  $max_right = $right_node_id1 if ($right_node_id1 > $right_node_id2);
  $max_right = $right_node_id2 if ($right_node_id2 > $right_node_id1);

  my $constraint = "WHERE t.root_id=$root_id1 AND left_index < $min_left";
  $constraint .= " AND right_index > $max_right";
  $constraint .= " ORDER BY (right_index-left_index) LIMIT 1";

  my $ancestor = $self->_generic_fetch($constraint)->[0];

  return $ancestor;
}

=head2 fetch_AlignedMember_by_member_id_mlssID

  Arg[1]     : int member_id of a peptide member (longest translation)
  Arg[2]     : [optional] int clusterset_id (def. 0)
  Example    :

      my $aligned_member = $proteintree_adaptor->
                            fetch_AlignedMember_by_member_id_mlssID
                            (
                             $member->get_longest_peptide_Member->member_id, $mlssID
                            );

  Description: Fetches from the database the protein_tree that contains the member_id
  Returntype : Bio::EnsEMBL::Compara::AlignedMember
  Exceptions :
  Caller     :

=cut


sub fetch_AlignedMember_by_member_id_mlssID {
  my ($self, $member_id, $mlss_id) = @_;
    
  my $constraint = "WHERE tm.member_id = $member_id and m.member_id = $member_id";
  $constraint .= " AND tm.method_link_species_set_id = $mlss_id" if($mlss_id and $mlss_id>0);
  my ($node) = @{$self->_generic_fetch($constraint)};
  return $node;
}

sub gene_member_id_is_in_tree {
  my ($self, $member_id) = @_;

  my $sth = $self->prepare("SELECT ptm1.root_id FROM member m1, protein_tree_member ptm1 WHERE ptm1.member_id=m1.member_id AND m1.gene_member_id=? LIMIT 1");
  $sth->execute($member_id);
  my($root_id) = $sth->fetchrow_array;

  if (defined($root_id)) {
    return $root_id;
  } else {
    return undef;
  }
}

sub fetch_all_AlignedMembers_by_root_id {
  my ($self, $root_id) = @_;

  my $constraint = "WHERE tm.root_id = $root_id";
  my $nodes = $self->_generic_fetch($constraint);
  return $nodes;

}

###########################
# STORE methods
###########################

sub store_all {
  my ($self,$node) = @_;

  my $mba = $self->db->get_MemberAdaptor();

  unless($node->isa('Bio::EnsEMBL::Compara::ProteinTree')) {
    throw("set arg must be a [Bio::EnsEMBL::Compara::ProteinTree] not a $node");
  }
  
  # Store each leaf if it's not stored yet.
  foreach my $leaf ($node->leaves) {
    $mba->store($leaf);
  }

  # Store the node structure.
  $self->store($node);

}

sub store {
  my ($self, $node) = @_;

  unless($node->isa('Bio::EnsEMBL::Compara::NestedSet')) {
    throw("set arg must be a [Bio::EnsEMBL::Compara::NestedSet] not a $node");
  }
  
  $self->store_node($node);
  
  # recursively do all the children
  foreach my $child_node (@{$node->sorted_children}) {  
    $self->store($child_node);
  }
  
  return $node->node_id;
}

sub store_node {
  my ($self, $node) = @_;

  unless($node->isa('Bio::EnsEMBL::Compara::NestedSet')) {
    throw("set arg must be a [Bio::EnsEMBL::Compara::NestedSet] not a $node");
  }
 
  if($node->adaptor and 
     $node->adaptor->isa('Bio::EnsEMBL::Compara::DBSQL::NestedSetAdaptor') and
     $node->adaptor eq $self) 
  {
      #already stored so just update
      #print "Node ".$node->node_id." already stored, just updating!\n";
      return $self->update_node($node);
  }

  my $parent_id = 0;
  my $root_id = 0;

  if(defined $node->parent) {
      $parent_id = $node->parent->node_id ;
      $root_id = $node->root->node_id;
  } else {
      #print "NO parent for ".$node->newick_format()."\n";
      $parent_id = 0;
      $root_id = 0;
  }
  #printf("INSERT NODE parent_id = %d, root_id = %d\n", $parent_id, $root_id);
  
  my $ptn = $self->protein_tree_node;
  my $sth = $self->prepare("INSERT INTO $ptn
                             (parent_id,
                              root_id,
                              left_index,
                              right_index,
                              distance_to_parent)  VALUES (?,?,?,?,?)");
  $sth->execute($parent_id, $root_id, $node->left_index, $node->right_index, $node->distance_to_parent);

  $node->node_id( $sth->{'mysql_insertid'} );
#  printf("  new node_id %d\n", $node->node_id);
  $node->adaptor($self);
  $sth->finish;

  if($node->isa('Bio::EnsEMBL::Compara::AlignedMember')) {
      my $ptm = $self->protein_tree_member;
      $sth = $self->prepare("INSERT INTO $ptm
                               (node_id,
                                member_id,
                                method_link_species_set_id,
                                cigar_line)  VALUES (?,?,?,?)");
      $sth->execute($node->node_id, $node->member_id, $node->method_link_species_set_id, $node->cigar_line);
      $sth->finish;

      #printf "INSERT protein_tree_member node_id = %d, member_id = %d\n", $node->node_id,$node->member_id;
  }
  
  return $node->node_id;
}


sub update_node {
  my ($self, $node) = @_;

  unless($node->isa('Bio::EnsEMBL::Compara::NestedSet')) {
    throw("set arg must be a [Bio::EnsEMBL::Compara::NestedSet] not a $node");
  }
  
  my $parent_id = 0;
  if($node->parent) {
    $parent_id = $node->parent->node_id ;
  }
  my $ptn = $self->protein_tree_node;
  my $sth = $self->prepare("UPDATE $ptn SET
                              parent_id=?,
                              left_index=?,
                              right_index=?,
                              distance_to_parent=? 
                            WHERE node_id=?");
  $sth->execute($parent_id, $node->left_index, $node->right_index, 
                $node->distance_to_parent, $node->node_id);

  $node->adaptor($self);
  $sth->finish;

  if($node->isa('Bio::EnsEMBL::Compara::AlignedMember')) {
    my $ptm = $self->protein_tree_member;
    my $sql = "UPDATE $ptm SET ". 
              "cigar_line='". $node->cigar_line . "'";
    $sql .= ", cigar_start=" . $node->cigar_start if($node->cigar_start);              
    $sql .= ", cigar_end=" . $node->cigar_end if($node->cigar_end);              
    $sql .= ", method_link_species_set_id=" . $node->method_link_species_set_id if($node->method_link_species_set_id);              
    $sql .= " WHERE node_id=". $node->node_id;
    $self->dbc->do($sql);
  }

}


sub merge_nodes {
  my ($self, $node1, $node2) = @_;

  unless($node1->isa('Bio::EnsEMBL::Compara::NestedSet')) {
    throw("set arg must be a [Bio::EnsEMBL::Compara::NestedSet] not a $node1");
  }
  
  # printf("MERGE children from parent %d => %d\n", $node2->node_id, $node1->node_id);
  
  my $sth = $self->prepare("UPDATE ".$self->protein_tree_node." SET
                              parent_id=?,
			                     WHERE parent_id=?");
  $sth->execute($node1->node_id, $node2->node_id);
  $sth->finish;
  
  $sth = $self->prepare("DELETE from ".$self->protein_tree_node." WHERE node_id=?");
  $sth->execute($node2->node_id);
  $sth->finish;
}


sub delete_node {
  my $self = shift;
  my $node = shift;
  
  my $node_id = $node->node_id;
  #print("delete node $node_id\n");
  $self->dbc->do("UPDATE ".$self->protein_tree_node." dn, ".$self->protein_tree_node." n SET ". 
            "n.parent_id = dn.parent_id WHERE n.parent_id=dn.node_id AND dn.node_id=$node_id");
  $self->dbc->do("DELETE from ".$self->protein_tree_node." WHERE node_id = $node_id");
  $self->dbc->do("DELETE from ".$self->protein_tree_tag." WHERE node_id = $node_id");
  $self->dbc->do("DELETE from ".$self->protein_tree_member." WHERE node_id = $node_id");
}


sub delete_nodes_not_in_tree
{
  my $self = shift;
  my $tree = shift;

  unless($tree->isa('Bio::EnsEMBL::Compara::NestedSet')) {
    throw("set arg must be a [Bio::EnsEMBL::Compara::NestedSet] not a $tree");
  }
  #print("delete_nodes_not_present under ", $tree->node_id, "\n");
  my $dbtree = $self->fetch_node_by_node_id($tree->node_id);
  my @all_db_nodes = $dbtree->get_all_subnodes;
  foreach my $dbnode (@all_db_nodes) {
    next if($tree->find_node_by_node_id($dbnode->node_id));
    $self->delete_node($dbnode);
  }
  $dbtree->release_tree;
}

sub delete_node_and_under {
  my $self = shift;
  my $node = shift;

  my @all_subnodes = $node->get_all_subnodes;
  foreach my $subnode (@all_subnodes) {
    my $subnode_id = $subnode->node_id;
    $self->dbc->do("DELETE from ".$self->protein_tree_node." WHERE node_id = $subnode_id");
    $self->dbc->do("DELETE from ".$self->protein_tree_tag." WHERE node_id = $subnode_id");
    $self->dbc->do("DELETE from ".$self->protein_tree_member." WHERE node_id = $subnode_id");
  }
  my $node_id = $node->node_id;
  $self->dbc->do("DELETE from ".$self->protein_tree_node." WHERE node_id = $node_id");
  $self->dbc->do("DELETE from ".$self->protein_tree_tag." WHERE node_id = $node_id");
  $self->dbc->do("DELETE from ".$self->protein_tree_member." WHERE node_id = $node_id");
}

###################################
#
# tagging 
#
###################################

sub _load_tagvalues {
  my $self = shift;
  my $node = shift;
  
  unless($node->isa('Bio::EnsEMBL::Compara::NestedSet')) {
    throw("set arg must be a [Bio::EnsEMBL::Compara::NestedSet] not a $node");
  }

  my $sth = $self->prepare("SELECT tag,value from ".$self->protein_tree_tag." where node_id=?");
  $sth->execute($node->node_id);  
  while (my ($tag, $value) = $sth->fetchrow_array()) {
    $node->add_tag($tag,$value);
  }
  $sth->finish;
}


sub _store_tagvalue {
  my $self = shift;
  my $node_id = shift;
  my $tag = shift;
  my $value = shift;
  
  $value="" unless(defined($value));

  my $sql = "INSERT ignore into ".$self->protein_tree_tag." (node_id,tag) values ($node_id,\"$tag\")";
  #print("$sql\n");
  $self->dbc->do($sql);

  $sql = "UPDATE ".$self->protein_tree_tag." set value=\"$value\" where node_id=$node_id and tag=\"$tag\"";
  #print("$sql\n");
  $self->dbc->do($sql);
}


##################################
#
# subclass override methods
#
##################################

sub columns {
  my $self = shift;

  my @cols;
  @cols = @{Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor->columns};
  return ['t.node_id',
          't.parent_id',
          't.root_id',
          't.left_index',
          't.right_index',
          't.distance_to_parent',
          
          'tm.cigar_line',
          'tm.cigar_start',
          'tm.cigar_end',

	  @cols
          ];
}

sub tables {
  my $self = shift;
  return [[$self->protein_tree_node, 't']];
}

sub left_join_clause {
    my $self = shift;
  return "left join ".$self->protein_tree_member." tm on t.node_id = tm.node_id left join member m on tm.member_id = m.member_id";
}

sub default_where_clause {
  return "";
}


sub create_instance_from_rowhash {
  my $self = shift;
  my $rowhash = shift;
  
  my $node;  
  if($rowhash->{'member_id'}) {
    $node = new Bio::EnsEMBL::Compara::AlignedMember;    
  } else {
    $node = new Bio::EnsEMBL::Compara::ProteinTree;
  }
  $self->init_instance_from_rowhash($node, $rowhash);
  return $node;
}

sub _objs_from_sth {
  my ($self, $sth) = @_;
  my $node_list = [];
  
  while(my $rowhash = $sth->fetchrow_hashref) {
    my $node = $self->create_instance_from_rowhash($rowhash);        
    push @$node_list, $node;
  }
  
  return $node_list;
}


sub init_instance_from_rowhash {
  my $self = shift;
  my $node = shift;
  my $rowhash = shift;
  
  #SUPER is NestedSetAdaptor
  $self->SUPER::init_instance_from_rowhash($node, $rowhash);

  if ($rowhash->{'member_id'}) {
     Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor->init_instance_from_rowhash($node, $rowhash);
	
	$node->cigar_line($rowhash->{'cigar_line'});
	#cigar_start and cigar_end does not need to be set.
	#$node->cigar_start($rowhash->{'cigar_start'});
	#$node->cigar_end($rowhash->{'cigar_end'});
  }
  # print("  create node : ", $node, " : "); $node->print_node;
  
  $node->adaptor($self);
  
  return $node;
}

sub reset_tables {
    my $self = shift;

    # Create tables if they don't exist, to avoid errors / warnings.
    $self->create_tables;

    my @tables = ($self->protein_tree_node,$self->protein_tree_tag,
		  $self->protein_tree_member, $self->protein_tree_omegas,
		  $self->protein_tree_score);
    foreach my $table (@tables) {
	warn "ProteinTreeMember.pm resetting TABLE $table";
	$self->dbc->do("truncate $table;");
    }
}

sub create_tables {
    my $self = shift;

    $self->db->get_MemberAdaptor->create_tables;
    $self->db->get_SequenceAdaptor->create_tables;

    #print "Creating tree/aln/tag tables: ".$self->table_base."_XYZ  ...\n";

    my $protein_tree_node = $self->protein_tree_node;
    my $protein_tree_tag = $self->protein_tree_tag;
    my $protein_tree_member = $self->protein_tree_member;
    my $protein_tree_score = $self->protein_tree_score;
    my $omegas = $self->protein_tree_omegas;

    my $cmd;

#
# Create the member and sequence tables.
#
    $cmd = <<"STRING";
CREATE TABLE IF NOT EXISTS $protein_tree_node (
  node_id                         int(10) unsigned NOT NULL auto_increment, # unique internal id
  parent_id                       int(10) unsigned NOT NULL,
  root_id                         int(10) unsigned NOT NULL,
  left_index                      int(10) NOT NULL,
  right_index                     int(10) NOT NULL,
  distance_to_parent              double default 1.0 NOT NULL,

  PRIMARY KEY (node_id),
  KEY (parent_id),
  KEY (root_id),
  KEY (left_index),
  KEY (right_index)
    ) COLLATE=latin1_swedish_ci;
STRING
    
    $self->dbc->do($cmd);

    $cmd = <<"STRING";
CREATE TABLE IF NOT EXISTS $protein_tree_member (
  node_id                     int(10) unsigned NOT NULL,
  member_id                   int(10) unsigned NOT NULL,
  method_link_species_set_id  int(10) unsigned NOT NULL,
  cigar_line                  mediumtext,
  cigar_start                 int(10),
  cigar_end                   int(10),

  FOREIGN KEY (node_id) REFERENCES $protein_tree_node(node_id),

  UNIQUE (node_id),
  KEY (member_id)
) COLLATE=latin1_swedish_ci;
STRING

    $self->dbc->do($cmd);


    $cmd = <<"STRING";
CREATE TABLE IF NOT EXISTS $protein_tree_score (
  node_id                     int(10) unsigned NOT NULL,
  member_id                   int(10) unsigned NOT NULL,
  method_link_species_set_id  int(10) unsigned NOT NULL,
  cigar_line                  mediumtext,
  cigar_start                 int(10),
  cigar_end                   int(10),

  FOREIGN KEY (node_id) REFERENCES $protein_tree_node(node_id),

  UNIQUE (node_id),
  KEY (member_id)
) COLLATE=latin1_swedish_ci;
STRING

    $self->dbc->do($cmd);


    $cmd = <<"STRING";
CREATE TABLE IF NOT EXISTS $protein_tree_tag (
  node_id                int(10) unsigned NOT NULL,
  tag                    varchar(50),
  value                  mediumtext,

  FOREIGN KEY (node_id) REFERENCES $protein_tree_node(node_id),

  UNIQUE tag_node_id (node_id, tag),
  KEY (node_id),
  KEY (tag)
) COLLATE=latin1_swedish_ci;
STRING

    $self->dbc->do($cmd);

#
# Create the omega values table.
#
    $cmd = <<"STRING";
CREATE TABLE IF NOT EXISTS $omegas (
  sitewise_id                 int(10) unsigned NOT NULL auto_increment, # unique internal id
  aln_position                int(10) unsigned NOT NULL,
  node_id                     int(10) unsigned NOT NULL,
  tree_node_id                int(10) unsigned NOT NULL,
  ncod                        int(10) unsigned,
  omega                       float(10,5),
  omega_lower                 float(10,5),
  omega_upper                 float(10,5),
  optimal                     float(10,5),
  threshold_on_branch_ds      float(10,5),
  type                        ENUM('na','all_gaps','constant','default','negative1','negative2','negative3','negative4','positive1','positive2','positive3','positive4',
                              'synonymous','random','single_char') NOT NULL,

  FOREIGN KEY (node_id) REFERENCES $protein_tree_node(node_id),

  UNIQUE aln_position_node_id_ds (aln_position,node_id,threshold_on_branch_ds),
  PRIMARY KEY (sitewise_id),
  KEY (tree_node_id),
  KEY (node_id)
) COLLATE=latin1_swedish_ci;
STRING
    
    $self->dbc->do($cmd);
}


# Get / set the table base name.
sub table_base {
    my $self = shift;
    $self->{'_table_base'} = shift if(@_);

    $self->{'_table_base'} = "protein_tree" if (!defined $self->{'_table_base'});

    return $self->{'_table_base'};
}

sub protein_tree_member {
    my $self = shift;
    $self->{'_protein_tree_member'} = shift if(@_);

    # If it's been explicitly set, use the stored value.
    return $self->{'_protein_tree_member'} if (defined $self->{'_protein_tree_member'});

    # Otherwise, build it up from the table base.
    return $self->table_base . "_member";
}

sub protein_tree_score {
    my $self = shift;
    $self->{'_protein_tree_score'} = shift if(@_);

    # If it's been explicitly set, use the stored value.
    return $self->{'_protein_tree_score'} if (defined $self->{'_protein_tree_score'});

    # Otherwise, build it up from the protein_tree_member table name.
    return $self->protein_tree_member . "_score";
}

sub protein_tree_omegas {
    my $self = shift;
    $self->{'_protein_tree_omegas'} = shift if(@_);

    # If it's been explicitly set, use the stored value.
    return $self->{'_protein_tree_omegas'} if (defined $self->{'_protein_tree_omegas'});

    # Otherwise, build it up from the protein_tree_member table name.
    #return $self->table_base . "_omegas";
    return "sitewise_aln";
}

sub protein_tree_tag {
    my $self = shift;
    $self->{'_protein_tree_tag'} = shift if(@_);

    # If it's been explicitly set, use the stored value.
    return $self->{'_protein_tree_tag'} if (defined $self->{'_protein_tree_tag'});

    # Otherwise, build it up from the table base.
    return $self->table_base . "_tag";
}

sub protein_tree_node {
    my $self = shift;
    $self->{'_protein_tree_node'} = shift if(@_);

    # If it's been explicitly set, use the stored value.
    return $self->{'_protein_tree_node'} if (defined $self->{'_protein_tree_node'});

    # Otherwise, build it up from the table base.
    return $self->table_base . "_node";
}


##########################################################
#
# explicit method forwarding to MemberAdaptor
#
##########################################################

sub _fetch_sequence_by_id {
  my $self = shift;
  return $self->db->get_MemberAdaptor->_fetch_sequence_by_id(@_);
}

sub fetch_gene_for_peptide_member_id { 
  my $self = shift;
  return $self->db->get_MemberAdaptor->fetch_gene_for_peptide_member_id(@_);
}

sub fetch_peptides_for_gene_member_id {
  my $self = shift;
  return $self->db->get_MemberAdaptor->fetch_peptides_for_gene_member_id(@_);
}

sub fetch_longest_peptide_member_for_gene_member_id {
  my $self = shift;
  return $self->db->get_MemberAdaptor->fetch_longest_peptide_member_for_gene_member_id(@_);
}


1;
