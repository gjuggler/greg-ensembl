### Runnable for PAML to detect positive selection in protein family alignments.
###

package Bio::EnsEMBL::Compara::RunnableDB::PAML;

use strict;
use Time::HiRes qw(time gettimeofday tv_interval);
use Cwd;
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::NestedSet;
use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Hive::Process;
use Bio::Tools::Run::Phylo::PAML::Codeml;

our @ISA = qw(Bio::EnsEMBL::Hive::Process);

#
# Some global-ish variables.
#
my $dba;
my $pta;

# INPUT FILES / OBJECTS.
my $tree;
my $input_cdna;
my $input_aa;
my $params;
my $tags;

# OUTPUT FILES / OBJECTS / STATES.
my $codeml; # The Bio::Tools::Run::... object.

sub debug {1;}

sub fetch_input {
  my( $self) = @_;
  
  # Load up the Compara DBAdaptor.
  $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-DBCONN=>$self->db->dbc);
  $pta = $dba->get_ProteinTreeAdaptor;
  
  ### DEFAULT PARAMETERS ###
  $params = {
    input_table            => 'protein_tree_member',
    output_table           => 'sitewise_aln',

    # Alignment Masking Parameters
    parameter_set_id       => 0,                    # This should be set to something sensible.
    mask_type              => 'auto',               # Options: 'none', 'threshold', 'auto'
    mask_param             => 5,
    remove_type            => 'auto',               # Options: 'none', 'branch length', 'auto'
    remove_param           => 0.8,

    # PAML Stuff.
    action                 => 'lrt',                # The action to perform with PAML.
                                                      # 'sitewise'   - Infer sitewise dN/dS ratios using the given model and BEB.
                                                      # 'reoptimize' - Reoptimize branch lengths using the given model.
                                                      # 'lrt'        - Perform a LRT using two models.
    model                  => 7,                    # The model to use.
    model_b                => 8,
                                                      # 0 - 8, the various M[0-8] models allowed in PAML.
                                                      # N.B.: 'model_b' is only relevant for the LRT action.

  };
  
  #########################
  
  # Fetch parameters from the two possible locations. Input_id takes precedence!
  $params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_string($params,$self->parameters);
  $params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_string($params,$self->input_id);

  #########################

#    $self->check_if_exit_cleanly;

  # Deal with alternate table inputs by changing the table used by our ProteinTreeAdaptor.
  $pta->protein_tree_member($params->{'input_table'});

  # Load the tree.
  # The node ID can come from either "protein_tree_id" or "node_id" in the parameters or input_id.
  if (defined $params->{'protein_tree_id'}) {
    $tree = $pta->fetch_node_by_node_id($params->{'protein_tree_id'});
  } elsif (defined $params->{'node_id'}) {
    $tree = $pta->fetch_node_by_node_id($params->{'node_id'});
  } else {
    throw("No protein tree input ID!\n");
  }

  foreach my $leaf ($tree->leaves) {
    #print $leaf->alignment_string."\n";
  }

  # Some last-minute adjustments based on retry counts or somesuch.

  # Think of reasons why we want to fail the job.

  ### EARLY EXIT IF THE TREE AINT RIGHT ###
}

sub run {
  my $self = shift;
  $self->{'start_time'} = time()*1000;

  # Get the cdna alignment and tree for input.
  my $cdna_aln = $tree->get_SimpleAlign(-cdna => 1);
  my $aa_aln = $tree->get_SimpleAlign();
  $input_cdna = Bio::EnsEMBL::Compara::ComparaUtils->fetch_masked_alignment($cdna_aln,$tree,$params,1);
  $input_aa = Bio::EnsEMBL::Compara::ComparaUtils->fetch_masked_alignment($aa_aln,$tree,$params,0);
  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);

  if ($params->{'action'} eq "reoptimize") {

    ### Reoptimize branch lengths.
    print $tree->newick_format."\n";
    my $new_treeI = Bio::Tools::Run::Phylo::PAML::Codeml->get_m0_tree($treeI,$input_cdna);
    my $new_pt = Bio::EnsEMBL::Compara::TreeUtils->from_treeI($new_treeI);
    print $new_pt->newick_format."\n";

    # Store new branch lengths back in the original protein_tree_node table.
    foreach my $leaf ($tree->leaves) {
      my $new_leaf = $new_pt->find_leaf_by_name($leaf->name);
      $leaf->distance_to_parent($new_leaf->distance_to_parent);
      $leaf->store;
      print "  -> Branch lengths updated!\n";
    }
    return;
  } elsif ($params->{'action'} eq "lrt") {

    ### Perform a likelihood ratio test between two models.
    my $model_a = $params->{'model'};
    my $model_b = $params->{'model_b'};
    
    my ($twice_lnL,$codeml_a,$codeml_b) = Bio::Tools::Run::Phylo::PAML::Codeml->NSsites_ratio_test($treeI,$cdna_aln,$model_a,$model_b);
    my $test_label = sprintf("PAML LRT M%s-M%s",$model_a,$model_b);
    $tags->{$test_label} = $twice_lnL;
    Bio::EnsEMBL::Compara::ComparaUtils->store_tags($tree,$tags);
    return;
  } else {

    ### Perform a BEB sitewise analysis of omegas.
    # Scale by a factor of 3 (PAML wants branch lengths in subst. / codon but Compara has them in subst. / site)
    #Bio::EnsEMBL::Compara::TreeUtils->scale($treeI,3);
    
    my $params = {
      NSsites => $params->{'model'}
    };
    
    $codeml = Bio::Tools::Run::Phylo::PAML::Codeml->new( -params => $params, 
							 -alignment => $input_cdna, 
							 -tree => $treeI);
    $codeml->run();
  }
  
  $self->{'end_time'} = time()*1000;
}


sub write_output {
  my $self = shift;
  
  if ($params->{'action'} ne "sitewise") {
    Bio::EnsEMBL::Compara::ComparaUtils->store_tags($tree,$tags);
      return;
    }
  
  my $table = $params->{'output_table'};
  
  my $aln_length = $input_aa->length;
  my $root_id = $tree->node_id;
  my $tree_node_id = 0;
  my $parameter_set_id = 0;
  
  $tree_node_id = $tree->subroot if ($tree->subroot);
  $parameter_set_id = $params->{'parameter_set_id'} if (defined $params->{'parameter_set_id'});
  
  my ($naive,$bayes,$se) = $codeml->extract_empirical_bayes();
  
  foreach my $site ( 1 .. $aln_length) {
    my $nongaps = Bio::EnsEMBL::Compara::AlignUtils->get_nongaps_at_column($input_aa,$site);
    
    my $omega = $naive->{$site};
    my $se = $se->{$site};
    $omega = $bayes->{$site} if ($bayes);
    
    next if (!defined $omega);

    my $lower = $omega - $se;
    my $upper = $omega + $se;

    my $sth = $tree->adaptor->prepare
      ("REPLACE INTO $table
                           (aln_position,
                            parameter_set_id,
                            node_id,
                            tree_node_id,
                            omega,
                            omega_lower,
                            omega_upper,
                            ncod
                            ) VALUES (?,?,?,?,?,?,?,?)");
    print "SITE: $site OMEGA: $omega\n";
    $sth->execute($site,
		  $parameter_set_id,
		  $root_id,
		  $tree_node_id,
		  $omega,
		  $lower,
		  $upper,
		  $nongaps
		  );
    $sth->finish();
  }

  # Store the log-likelihood as a protein tree tag.
  my $lnL = $codeml->extract_lnL();
  my $tag_name = "PAML lnL M" . $params->{'model'};
  $tags->{$tag_name} = $lnL;

  Bio::EnsEMBL::Compara::ComparaUtils->store_tags($tree,$tags);
}

1;
