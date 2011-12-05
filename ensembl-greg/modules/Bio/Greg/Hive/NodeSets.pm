package Bio::Greg::Hive::NodeSets;

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

sub param_defaults {
  my $self = shift;
  my $params = {
    flow_node_set            => 'Primates',
    flow_parent_and_children => 0,
    debug                    => 0,
    keep_at_least_root       => 0,
    merge_by_gene_names      => 0,
    create_orthologs_table => 1
  };
  return $params;
}

sub fetch_input {
  my ($self) = @_;
  $self->load_all_params();

  $self->create_table_from_params( $self->compara_dba, 'trees',
                                   $self->trees_table );

}

sub run {
  my $self = shift;

  $self->get_output_folder;

  my $tree = $self->get_tree;

  $self->param('orig_leaf_count', scalar($tree->leaves));
  $self->param('tree', $tree);

  return if (scalar($tree->leaves) < 2);
  
  my $genome_tree = Bio::EnsEMBL::Compara::ComparaUtils->get_genome_tree_subset($self->compara_dba,$self->params);

  my $human_gene_count = scalar(grep {$_->taxon_id==9606} $tree->leaves);
  $self->param('orig_human_gene_count',$human_gene_count);
  $self->param('genome_species_count',scalar($genome_tree->leaves));

  print "Tagging root nodes...\n";
  my @sets = $self->set_names;
  foreach my $set (@sets) {
    $self->tag_root_nodes($tree, $set);
  }

  $self->_mark_ortholog_trees($tree, 9606, 'Human Orthologs');
  $self->_mark_ortholog_trees($tree, 10090, 'Mouse Orthologs');
  $self->_mark_ortholog_trees($tree, 7955, 'Zebrafish Orthologs');
  $self->_mark_ortholog_trees($tree, 7227, 'Drosophila Orthologs');
}

sub set_names {
  my $self = shift;

  return (
    'Primates', 'Glires', 'Laurasiatheria',
    'Eukaryota', 'Euteleostomi', 'Amniota', 'Eutheria',
    'Clupeocephala', 'Sauria',
    'MammalSubgroups',
    'MammalSubgroupsPlusOutgroup',
    'Ensembl Roots'
    );
}

sub write_output {
  my $self = shift;

  if ($self->param('create_orthologs_table')) {
    $self->_store_in_table($self->param('tree'));
    return;
  }

  my $tree = $self->param('tree');

  my $flowed_human_gene_count = 0;

  if (defined $self->param('flow_node_set') && $self->param('flow_node_set') eq 'one2one') {
    # The node set method is one2one orthologs. Ensure that a gene ID is set, and return
    # to let the autoflow bring us to the next Process.
    # When one2one is set, the get_tree_for_comparative_analysis method in ComparaUtils
    # will restrict the tree to one-to-one orthologs of the member defined by gene_id.
    my $gene_id = $self->param('gene_id');
    my $node_id = $self->param('node_id');
    my $params = {
      flow_node_set => 'one2one',
      gene_id => $gene_id,
      node_id => $node_id
    };
    my ($output_job_id) = @{ $self->dataflow_output_id( $params, 1 ) };
    my $param_string = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
    print " -> Flowed node $node_id [$param_string] [job id: $output_job_id]\n";

    return;
  }

  if ( defined $self->param('flow_node_set') ) {
    $self->input_job->autoflow(0);
    my $flow_set = $self->param('flow_node_set');

    foreach my $node ( $tree->nodes ) {
      next if ( $node->is_leaf );
      my $params = {
        orig_node_id => $self->param('node_id'),
        node_id => $node->node_id,
        node_set_leaf_count => scalar($node->leaves)
      };
      my $id = $node->node_id;

      if ( $self->param('cc_root_'.$flow_set.'_'.$node->node_id) == 1) {

        # Don't flow this node if this job has a specific target gene ID.
        if (defined $self->param('gene_id')) {
          my $gene_id = $self->param('gene_id');
          print "Looking for gene $gene_id...\n";
          next unless ($self->tree_contains_a_gene_with_name($node,$gene_id));
          print "Found it!\n";
          $params->{gene_id} = $gene_id;
        }

        $self->new_data_id($params);

        my $param_string = Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params);
        
        my ($output_job_id) = @{ $self->dataflow_output_id( $params, 1 ) };
        print " -> Flowed node $id [$param_string] [job id: $output_job_id]\n";

        $flowed_human_gene_count += scalar(grep {$_->taxon_id==9606} $node->leaves);

        if ( $self->param('flow_parent_and_children') ) {
          my $i = 0;
          foreach my $child ( @{ $node->children } ) {
            my $output_id = $self->replace($params,{
              node_id               => $child->node_id,
              orig_node_id          => $self->param('node_id'),
              node_set_parent_id    => $id,
              node_set_child_number => $i++
            });
            my ($output_job_id) = @{ $self->dataflow_output_id( $output_id, 1 ) };
            print "  --> Flowed child $output_job_id\n";
          }
        }
      }
    }
  }

  $self->param('flowed_human_gene_count',$flowed_human_gene_count);
}

sub _mark_ortholog_trees {
  my $self = shift;
  my $tree = shift;
  my $ref_taxon_id = shift;
  my $set_name = shift;

  my @ref_members = grep {$_->taxon_id == $ref_taxon_id} $tree->leaves;

#  print $tree->ascii."\n";

  foreach my $ref (@ref_members) {
    print "Getting ortholog tree for ".$ref->stable_id."\n";
    my $ortholog_tree = Bio::EnsEMBL::Compara::ComparaUtils->get_one_to_one_ortholog_tree($self->compara_dba, $ref, 'ortholog.*');
    next if (!defined $ortholog_tree);

    #print $ortholog_tree->ascii."\n";
    print $ortholog_tree . "  ". $ortholog_tree->node_id." " . $ref->stable_id."\n";
    $self->param('cc_root_'.$set_name.'_'.$ortholog_tree->node_id, 1);
    $self->_store_tree_in_table($ortholog_tree, $set_name, $ref);
  }
}

sub _store_in_table {
  my $self = shift;
  my $tree = shift;

  my @set_names = $self->set_names;

  foreach my $node ( $tree->nodes ) {
    next if ( $node->is_leaf );
    foreach my $set (@set_names) {
      my $key = 'cc_root_'.$set.'_'.$node->node_id;
      if ($self->param('cc_root_'.$set.'_'.$node->node_id)) {
        print "$key\n";
        $self->_store_tree_in_table($node, $set, undef);
      }
    }
  }
  
}

sub species_taxid {
  my $self = shift;
  
  return {
    human => 9606,
    chimp => 9598,
    gorilla => 9593,
    orang => 9601,
    rhesus => 9544,
    gibbon => 61853,
    marmoset => 9483,
    tarsier => 9478,
    bushbaby => 30611,
    mouse_lemur => 30608,
    tree_shrew => 37347,

    mouse => 10090,
    rat => 10116,
    kangaroo_rat => 10020,
    squirrel => 43179,
    guinea_pig => 10141,
    pika => 9978,
    rabbit => 9986,
    
    hedgehog => 9365,
    shrew => 42254,
    microbat => 59463,
    megabat => 132908,
    dog => 9615,
    panda => 9646,
    cat => 9685,
    pig => 9823,
    alpaca => 30538,
    dolphin => 9739,
    cow => 9913,
    horse => 9796,

    tenrec => 9371,
    elephant => 9785,
    hyrax => 9813,

    sloth => 9358,
    armadillo => 9361,

    opossum => 13616,
    wallaby => 9315,

    platypus => 9258,

    chicken => 9031,
    turkey => 9103,
    zebrafinch => 59729,

    lizard => 28377,

    xenopus => 8364,

    stickleback => 69293,
    medaka => 8090,
    fugu => 31033,
    tetraodon => 99883,
    zebrafish => 7955,

    c_intestinalis => 7719,
    c_savignyi => 51511,

    drosophila => 7227,
    c_elegans => 6239,
    yeast => 4932
  };
}

sub _store_tree_in_table {
  my $self = shift;
  my $subtree = shift;
  my $set_name = shift;
  my $ref_member = shift;

  my $p = $self->params;

  $p->{method_id} = $set_name;

  if (defined $ref_member) {
    $p->{tree_id} = $ref_member->stable_id;
  } else {
    $p->{tree_id} = $subtree->node_id;
  }

  my $taxids;
  map {$taxids->{$_->taxon_id} = 1} $subtree->leaves;
  $p->{species_count} = scalar(keys %$taxids);
  $p->{leaf_count} = scalar($subtree->leaves);
  
  # Don't allow branches with length > 2.
  my $max_bl = 2;
  map {$_->branch_length($max_bl) if ($_->branch_length > $max_bl)} $subtree->nodes;
  my $treei = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($subtree);
  $p->{tree_length} = $treei->total_branch_length;
  $p->{tree_mpl} = $treei->root->mean_path_length;

  # Gene names
  my @good_genes = grep {$_->taxon_id == 9606 || $_->taxon_id == 10090 || $_->taxon_id == 7227} $subtree->leaves;
  my @gene_names = map {$_->get_Gene->external_name} @good_genes;
  my $gene_hash;
  map {$gene_hash->{lc($_)} = 1} @gene_names;
  delete $gene_hash->{''};
  my $name_str = join(', ', sort keys %$gene_hash);
  $p->{gene_names} = $name_str;

  my $species_taxid = $self->species_taxid;
  
  my @leaves = $subtree->leaves;
  foreach my $species (keys %$species_taxid) {
    my $taxid = $species_taxid->{$species};
    my $count = scalar(grep {$_->taxon_id == $taxid} @leaves);
    $p->{$species.'_count'} = $count;
  }

  $self->store_params_in_table($self->dbc, 'trees', $p);
}

sub tree_contains_a_gene_with_name {
  my $self = shift;
  my $tree = shift;
  my $name = shift;

  my @members = $tree->leaves;
  foreach my $member (@members) {
    # Note: the 'gene_name' can be either ENSG, ENST, ENSP, or gene-name like string.
    if ($member->stable_id eq $name) {
      return 1;
    }
    if ($member->gene_member->stable_id eq $name) {
      return 1;
    }
    my $gene = $member->get_Gene;
    if (!$gene) {
      next;
    }
    my $member_name = $gene->external_name;
    if ($member_name eq $name) {
      return 1;
    }
  }
  return 0;
}

sub tag_root_nodes {
  my $self        = shift;
  my $tree        = shift;
  my $method_name = shift;

  my $base_p = {
    tree             => $tree,
    subtree_function => \&does_parent_have_clade_children,
    min_size      => 2,
    max_size => ($self->param('genome_species_count') * 3)
  };

  my $params = $base_p;

  if ($method_name eq 'MammalSubgroups') {
    $params = $self->replace_params(
      $base_p, {
        cc_Laurasiatheria => 0.1,
        cc_Glires => 0.1,
        cc_Primates => 0.1
      }
    );
  } elsif ( $method_name eq 'MammalSubgroupsPlusOutgroup' ) {
    $params = $self->replace_params(
      $base_p, {
        cc_Primates       => 0.1,
        cc_Glires         => 0.1,
        cc_Laurasiatheria => 0.1,
        any      => [ 'Sauria', 'Clupeocephala', 'Ciona', 'Marsupialia' ]
      }
    );
  } elsif ( $method_name eq 'MammalPlusOutgroup' ) {
    $params = $self->replace_params(
      $base_p, {
        cc_Eutheria => 0.1,
        any      => [ 'Sauria', 'Clupeocephala', 'Ciona', 'Marsupialia' ]
      }
    );

  } else {
    $params = $self->replace_params( $base_p, {
      'cc_'.$method_name => 0.6,
      min_size => 2,
      max_size => 500
                                     });
  }

  #my $subtree_function = $params->{subtree_function};
  my @root_nodes;
  if ($method_name eq 'Ensembl Roots') {
    @root_nodes = ($tree);
    delete $params->{max_size};
  } else {
    @root_nodes = $self->get_smallest_subtrees_from_node( $tree, $params );
  }

  # Limit the max size of trees.
  if (defined $params->{max_size}) {
    @root_nodes = grep {scalar($_->leaves) <= $params->{max_size} } @root_nodes;
  }

  foreach my $node (@root_nodes) {
    print "$node\n";
    printf " -> %s %d %d\n", $method_name, scalar $node->leaves, $node->node_id;
    $self->param('cc_root_'.$method_name.'_'.$node->node_id,1);
    #$node->store_tag( "cc_root_" . $method_name, 1 );
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
      #$node->store_tag( "cc_$clade", sprintf( "%.3f", $coverage_fraction ) );
    }
  }
}

sub clade_coverage_for_node {
  my $self  = shift;
  my $node  = shift;
  my $clade = shift;

  my $key = "_taxon_ids_" . $clade;
  if ( !defined $self->{$key} ) {
    my @genomes = $self->get_genomes_within_clade( $self->compara_dba, $clade );
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
    my $value = $params->{$key};

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
    } elsif ($key =~ m/cc_/) {
      $key =~ s/cc_//;
      my $coverage = $self->clade_coverage_for_node( $tree, $key, $params );
      #print "$key $coverage $value\n" if ($coverage > 0);
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

  # 1. check that this tree satisfies inclusion function
  # 2. check that our sister node satisfies inclusion function.
  # 3. if (1) and (2) are met, return true.

  my $tree             = $node;
  my $parent           = $node->parent;
  my @parents_children = @{ $parent->children };
  if ( $parent->node_id == 1 ) {
    my $value = $inclusion_function->( $self, $node, $params );

    #print "Root node. Values is $value\n";
    if ( $self->param('keep_at_least_root') == 1 ) {
      return 1;
    } else {
      return $value;
    }
  }

  my $sister;
  foreach my $ch (@parents_children) {
    $sister = $ch if ( $ch->node_id != $tree->node_id );
  }

  if ( $inclusion_function->( $self, $node, $params ) == 1
    && $inclusion_function->( $self, $sister, $params ) == 1 ) {
    if ( $self->param('merge_by_gene_names') == 1 ) {
      return 1 if ( $self->no_similar_gene_names( $node, $sister ) == 1 );
    } else {
      return 1;
    }
  }
  return 0;
}

sub no_similar_gene_names {
  my $self   = shift;
  my $node_a = shift;
  my $node_b = shift;

  my $names_a;
  map {
    my $gene = $_->get_Gene;
    if ( defined $gene ) {
      $names_a->{ lc( $gene->external_name ) } = 1;
    }
  } $node_a->leaves;

  my $names_b;
  map {
    my $gene = $_->get_Gene;
    if ( defined $gene ) {
      $names_b->{ lc( $gene->external_name ) } = 1;
    }
  } $node_b->leaves;

  foreach my $key ( keys %$names_a ) {

    #print "$key\n";
    if ( $names_b->{$key} == 1 ) {
      return 0;
    }
  }
  return 1;
}

sub DESTROY {
  my $self = shift;

  my $tree = $self->get_tree;
  $tree->release_tree if ($tree);
  $tree = undef;
  $self->SUPER::DESTROY if $self->can("SUPER::DESTROY");
}

# Returns the NCBI taxnomy of Ensembl genomes below a given taxonomic clade.
sub get_genome_taxonomy_below_level {
  my $self          = shift;
  my $dba           = shift;
  my $root_taxon_id = shift || 33154;
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
  my $species_tree = Bio::EnsEMBL::Compara::ComparaUtils->get_genome_taxonomy_below_level( $dba, $clade );

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

sub trees_table {
  my $self = shift;

  my $p = {
    data_id => 'int',
    tree_id => 'char32',
    method_id => 'char32',

    tree_length => 'float',
    tree_mpl => 'float',
    
    species_count => 'int',
    leaf_count => 'int',

    gene_names => 'string',

    unique_keys => 'data_id,method_id,tree_id'
  };

  my $species_taxid = $self->species_taxid;
  foreach my $species (sort keys %$species_taxid) {
    $p->{$species.'_count'} = 'int';
  }

  return $p;
}

1;
