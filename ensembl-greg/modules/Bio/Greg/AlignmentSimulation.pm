#
# You may distribute this module under the same terms as perl itself
#

package Bio::Greg::AlignmentSimulation;

use strict;
use Getopt::Long;
use IO::File;
use File::Basename;
use File::Path;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::BaseAlignFeature;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::PeptideAlignFeatureAdaptor;
use Bio::EnsEMBL::Compara::Member;
use Bio::EnsEMBL::Compara::ProteinTree;
use Bio::EnsEMBL::Compara::NestedSet;
use Bio::EnsEMBL::Compara::ComparaUtils;

use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Hive::Process;

use Bio::AlignIO;

our @ISA = qw(Bio::EnsEMBL::Hive::Process);

my $dba;
my $pta;

my $params;
my $tags;

my $tree;
my $final_cdna;
my @sitewise_omegas;

sub fetch_input {
  my $self = shift;

  # Load up the Compara DBAdaptor.
  $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-DBCONN=>$self->db->dbc);
  $pta = $dba->get_ProteinTreeAdaptor;

  ### DEFAULT PARAMETERS ###
  $params = {
    input_table_base             => 'protein_tree',
    output_table_base            => 'protein_tree',

    # Dawg Parameters.
    length                       => 300,
    lambda                       => "{0.02, 0.02}",           # {in,del}
    dawg_mult                    => .33,

    # Evolver Parameters.
    function                     => 'constant',   # Options: 'constant' (a) , 'uniform' (a,b), 'gamma' (a,b), 'beta' (p,q)
#    a                            => 0.2,
#    b                            => .5,
#    p                            => .5,
#    q                            => .5,
#    num_bins                     => 5,
  };
  
  #########################
  
  # Fetch parameters from the two possible locations. Input_id takes precedence!
  # (this utility method is from AlignmentProcess.pm)
  $params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_string($params,$self->parameters);
  $params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_string($params,$self->input_id);

  #########################

  #$self->check_if_exit_cleanly;

  # Load the tree.
  my $node_id;
  $node_id = $params->{'protein_tree_id'};
  $node_id = $params->{'node_id'} if (!defined $node_id);
  $pta->table_base($params->{'input_table_base'});
  $tree = $pta->fetch_node_by_node_id($node_id);

  print "Simulation tree: ".$tree->newick_format."\n";

  # Fetch sim params from (a) the simulation root and (b) the tree root.
  my $sim_root = $tree->parent;
  my $sp_str = $sim_root->get_tagvalue('simulation_params');
  $params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_string($params,$sp_str);
  $sp_str = $tree->get_tagvalue('simulation_params');
  $params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_string($params,$sp_str);

  print Bio::EnsEMBL::Compara::ComparaUtils->hash_to_string($params)."\n";

  # Some last-minute adjustments based on retry counts or somesuch.
  # Or think of reasons why we want to fail the job.
  throw("No protein tree!") unless (defined $tree);

}

sub run {
  my $self = shift;
  $self->check_if_exit_cleanly;
  $self->{'start_time'} = time()*1000;

  $final_cdna = $self->sim_simulate_alignment_from_tree($tree,$params);
  
  $self->{'end_time'} = time()*1000;
}


sub write_output {
  my $self = shift;
  my $final_aa = $final_cdna->translate();

  Bio::EnsEMBL::Compara::AlignUtils->pretty_print($final_aa);

  # Output the alignment and omegas.
  $pta->table_base($params->{'output_table_base'});
  my $out_table = $pta->protein_tree_member;
  Bio::EnsEMBL::Compara::ComparaUtils->store_SimpleAlign_into_table($out_table,$tree,$final_aa,$final_cdna);
  my $omega_out = $pta->protein_tree_omegas;
  $self->_store_sitewise_omegas($omega_out,\@sitewise_omegas,$final_aa);
}


###
###
###

sub sim_simulate_alignment_from_tree {
  my $self = shift;
  my $tree = shift;
  my $sim_params = shift;

  my $treeI;
  if (!ref $tree) {
    $treeI = Bio::EnsEMBL::Compara::TreeUtils->treeI_from_newick($tree);
  } elsif ($tree->isa("Bio::Tree::TreeI")) {
    $treeI = $tree;
  } elsif ($tree->isa("Bio::EnsEMBL::Compara::NestedSet")) {
    $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
  } else {
    print "asdf!\n";
  }

  # TODO (someday): simulate "chunks" with separate subst / indel parameters, to approximate
  # the domain / loop structures seen in real proteins. This could be an interesting project
  # in its own right -- maybe use an HMM to describe the domain structure of proteins??

  # GJ 2009-01-15 : TODO (soon): allow root sequences to be specified, so eventually we can start
  # "interjecting" simulated sequences into Ensembl gene trees.
  # This could be interesting in two ways:
  # 1) It may show us the maximum precision if we were to sequence all possible species
  # 2) We could perform a "omega bootstrap", to re-simulate evolution at a certain node and see
  #   how strong the evidence for the current omega value is...

  my $dawg_params = {
    'length' => 300,
    'lambda' => "{0.02, 0.02}",           # {in,del}
    temp_dir => "/tmp/dawg",
    dawg_mult => .33
    };
  my $omega_params = {
    num_bins         => 5,
    function         => 'gamma',
    a                => .5,
    b                => .5,
    p                => .5,
    q                => .5
    };
  my $evolver_params = {
    # Nothing yet!
  };

  foreach my $param (keys %{$omega_params}) {
    $omega_params->{$param} = $sim_params->{$param} if (defined $sim_params->{$param});
  }
  foreach my $param (keys %{$dawg_params}) {
    $dawg_params->{$param} = $sim_params->{$param} if (defined $sim_params->{$param});
  }
  
  # Scale the tree by a given length.
  if (defined $sim_params->{'tree_length'}) {
    $treeI = Bio::EnsEMBL::Compara::TreeUtils->scale_to($treeI,$sim_params->{'tree_length'});
  }
  
  # Simulate the alignment shape using Dawg.
  $dawg_params->{'Length'} = $dawg_params->{'length'};
  $dawg_params->{'Lambda'} = $dawg_params->{'lambda'};
  my $temp_dir = $dawg_params->{'temp_dir'}; # Extract the temp dir.
  my $dawg_mult = $dawg_params->{'dawg_mult'};
  delete $dawg_params->{'lambda'};
  delete $dawg_params->{'length'};
  delete $dawg_params->{'temp_dir'};
  delete $dawg_params->{'dawg_mult'};
  my $dawg_tree = Bio::EnsEMBL::Compara::TreeUtils->scale($treeI,$dawg_mult);
  my $dawg_cdna = $self->sim_simulate_dawg($dawg_tree,$temp_dir,$dawg_params);
  $dawg_cdna = $dawg_cdna->sort_by_tree($treeI);

  # Update the evolver simulation length.
  $evolver_params->{'nuclsites'} = $dawg_cdna->length/3;
  $omega_params->{'nuclsites'} = $dawg_cdna->length/3;
  
  # Sample omegas from the given distribution.
  my @omegas = @{$self->sim_sample_omegas($treeI,$omega_params)};
  # Run Evolver.
  $evolver_params->{'omegas'} = \@omegas;
  my $evolver_cdna = $self->sim_simulate_evolver($treeI,$evolver_params);
  $evolver_cdna = $evolver_cdna->sort_by_tree($treeI);
  
  # Stamp the dawg alignment onto the evolver codons.
  my $stamped_cdna = $self->sim_stamp_shape_onto_seqs($dawg_cdna,$evolver_cdna);
  
  # Fill up the sitewise_omegas with the omega values.
  for (my $i=0; $i < $stamped_cdna->length/3; $i++) {
    $sitewise_omegas[$i] = {aln_position => $i+1, omega_lower => $omegas[$i], omega_upper => $omegas[$i], omega => $omegas[$i], node_id => $tree->node_id};
  }
  
  return $stamped_cdna;
}


sub sim_simulate_dawg {
  my ($self,$tree,$temp_dir,$dawg_params) = @_;

  mkdir($temp_dir);
  
  my $ctrl_file = $temp_dir . "/dawg.ctrl";
  my $dawg_output = $temp_dir . "/dawg_output.fasta";

  my $newick = Bio::EnsEMBL::Compara::TreeUtils->to_newick($tree);

  # The indel params are taken from http://scit.us/projects/dawg/wiki/DawgFaq.
  my $params_to_use = {
    'Tree'   => $newick,
    'File'   => qq("$dawg_output"),
    'Format' => qq("Fasta"),
    'Width'  => 3,		              # Simulate codon blocks.
    'GapModel' => qq("PL"),
    'GapParams' => "{1.33, 100}",
    'TreeScale' => 1                  # I think we need to scale down to 1/3 to account for the x3 block-width. TODO: Check on this...
    };

  foreach my $key (keys %$dawg_params) {
    $params_to_use->{$key} = $dawg_params->{$key};
  }

  open(OUT,">$ctrl_file");
  foreach my $key (sort keys %{$params_to_use}) {
    print OUT $key . "=" . $params_to_use->{$key} . "\n";
  }
  close(OUT);

  my $rc = system("dawg $ctrl_file");
  if ($rc != 0) {
#    $self->DESTROY();
    die("Dawg simulation failed: $!\n");
  }

  # Read the output into a SimpleAlign object and return it.
  my $alignio = Bio::AlignIO->new(-file => $dawg_output,
				  -format => "fasta");
  my $aln = $alignio->next_aln();

  print "DAWG CDNA_length: " . $aln->length. "   AA_length: " . $aln->length/3 . "\n";
  return $aln;
}


sub sim_sample_omegas {
  my $self = shift;
  my $tree = shift;
  my $fp = shift; # Evolver Params.

  my $function = $fp->{'function'};
  my $num_sites = $fp->{'nuclsites'};

  # The array in which we'll store omega values.
  my @omegas = ();

  if ($function eq "constant") {
    @omegas = ($fp->{'a'}) x $num_sites;
  } else {
    @omegas = $self->sim_get_bins_from_R($function,$fp);
  }

  return \@omegas;  # return as arrayref.
}


# Calls R to extract the equiprobable bins for a given distribution.
sub sim_get_bins_from_R {
  my $self = shift;
  my $function = shift;
  my $params = shift;

  my $temp_in = "temp.txt";
  my $temp_out = "temp_out.txt";

  my $num_sites = $params->{'nuclsites'};
  my $k = $params->{'num_bins'};
  $k = 10 unless (defined $k);

  my $r_cmd = qq^
    dd = discretize_distribution(n=$k)
    means = dd\$means
    vals = sample(means,$num_sites,replace=TRUE)
    write(vals,file="$temp_out",ncolumns=1)
    ^;
  
  my $r_pre = "";
  if ($function eq "uniform") {
    my $a = $params->{'a'};
    my $b = $params->{'b'};
    
    $r_pre = qq^
      lo = $a;
    hi = $b;
    p1 = function(x) {
      return(punif(x,min=lo,max=hi));
    }
    d1 = function(x) {
      return(dunif(x,min=lo,max=hi));
    }
    q1 = function(x) {
      return(qunif(x,min=lo,max=hi));
    }
    ^;
  } elsif ($function eq "gamma") {
    my $a = $params->{'a'};
    my $b = $params->{'b'};

    $r_pre = qq^
      shape = $a;
      rate = $b;
      p1 = function(x) {
        return(pgamma(x,shape=shape,rate=rate));
      }
      d1 = function(x) {
        return(dgamma(x,shape=shape,rate=rate));
      }
      q1 = function(x) {
        return(qgamma(x,shape=shape,rate=rate));
      }
    ^;
  } elsif ($function eq "beta") {
    my $p = $params->{'p'};
    my $q = $params->{'q'};

    $r_pre = qq^
      shape1 = $p;
      shape2 = $q;
      p1 = function(x) {
        return(pbeta(x,shape1=shape1,shape2=shape2));
      }
      d1 = function(x) {
        return(dbeta(x,shape1=shape1,shape2=shape2));
      }
      q1 = function(x) {
        return(qbeta(x,shape1=shape1,shape2=shape2));
      }
    ^;
  }

  my $r_prefix = qq^
    $r_pre
    discretize_distribution = function(n=4,p=p1,d=d1,q=q1) {
      num_categories = n
	quantile_bounds = seq(from=0,to=1,length.out=num_categories+1)
	x_bounds = q(quantile_bounds)
	
	if(is.infinite(x_bounds[1])) {
	  x_bounds[1] = -100
	}
      if(is.infinite(x_bounds[num_categories+1])) {
	x_bounds[num_categories+1] = 100
      }
      
      mean_indices = c()
	for (i in 1:(length(x_bounds)-1)) {
	  lo_cum = p(x_bounds[i])
	  hi_cum = p(x_bounds[i+1])
	  
	  x = seq(from=x_bounds[i],to=x_bounds[i+1],by=0.001)
	  cumulative = p(x)
	  dens = d(x)
	  
	  mid_cum = (hi_cum + lo_cum) / 2
	  mid_index = sum(cumulative <= mid_cum)
	  mid_x = x[mid_index]
	  mean_indices[i] = mid_x
	}
      ret = list(bounds=x_bounds,means=mean_indices)
      return(ret)
    }
    test = function() {
      dd = discretize_distribution(n=5)
      mean_indices = dd\$means
      x_bounds = dd\$bounds
      
      x = seq(from=-10,to=10,by=0.01)
      plot(x,d(x),ylim=c(0,4),xlim=c(0,5),type='l')
      abline(v=x_bounds,col='black');
      abline(v=mean_indices,col=rgb(1,0,0,0.25))
    }
  ^;
  
  $r_cmd = $r_prefix . "\n" . $r_cmd;
  
  open(OUT,">$temp_in");
  print OUT $r_cmd."\n";
  close(OUT);
  
# cmd to run: /software/R-2.7.1/bin/R CMD BATCH $filename
  my $rc = system("/software/R-2.7.1/bin/R --vanilla --slave < $temp_in ");
  die "R returned an error!" if ($rc);
  
  open(IN,"$temp_out");
  my @lines = <IN>;
  map {chomp} @lines; 
  close(IN);
  
  unlink($temp_in);
  unlink($temp_out);
  
  return @lines;
}

sub sim_simulate_evolver {
  my ($self,$treeI,$evolver_params) = @_;
  use Bio::Tools::Run::Phylo::PAML::Evolver;
  
  my $num_codons = $evolver_params->{'nuclsites'};
  my @omegas = @{$evolver_params->{'omegas'}};
  delete $evolver_params->{'omegas'};
  
  my $full_aln = Bio::SimpleAlign->new();
  
  # Find each unique omega value.
  @omegas = map {sprintf("%.3f",$_)} @omegas;  # Round values to a certain precision.
  my %hash = map {$_ => 1} @omegas;
  my @keys = sort {$a<=>$b} keys(%hash);

  my $hash_of_hashrefs;
  foreach my $omega (@keys) {
    my $n=0;
    map {$n++ if ($_ == $omega) } @omegas;

    print "Omega: $omega  Codon count: $n \n";

    my $evolver = Bio::Tools::Run::Phylo::PAML::Evolver->new();
    foreach my $key (keys %$evolver_params) {
      $evolver->set_parameter($key,$evolver_params->{$key});
    }
    $evolver->set_parameter("omega",$omega);
    $evolver->set_parameter("nuclsites",$n);

    $evolver->tree($treeI);

    #$evolver->verbose(1);
    #$evolver->save_tempfiles(1);
    #print $evolver->tempdir();

    $evolver->prepare;
    $evolver->run;	
    my $aln = $evolver->alignment;

    # Go through the input omega list, and fill in the appropriate alignment columns with codons from this run.
    my $j=0;    # J holds the position in Evolver's output that we're currently grabbing from.
    for (my $i=0; $i < scalar(@omegas); $i++) { # go through all omega values.
      if ($omegas[$i] == $omega) {
	# Grab out the cdna slice for this position, and increment the aln placeholder $j.
	my $slice = $aln->slice(($j*3)+1,($j*3)+3);
	$j++;
	foreach my $seq ($slice->each_seq) {
	  # Put the current codon into the correct position of the codon array.
	  $hash_of_hashrefs->{$seq->id} = {} unless (defined $hash_of_hashrefs->{$seq->id});
	  my $asdf = $hash_of_hashrefs->{$seq->id};
	  $asdf->{$i} = $seq->seq;
	}
      }
    }
  }

  foreach my $key (keys %{$hash_of_hashrefs}) {
    my $hash_ref = $hash_of_hashrefs->{$key};
    
    my @xs = sort {$a <=> $b} keys(%{$hash_ref});
    my @aas = map {$hash_ref->{$_}} @xs;

    my $seq = join("",@aas);
    my $seqI = new Bio::LocatableSeq(-seq => $seq,
				     -id => $key
				     );
    $full_aln->add_seq($seqI);
  }
  return $full_aln;
}

sub _store_sitewise_omegas {
  my $self = shift;
  my $output_table = shift;
  my $sitewise_ref = shift;
  my $sa = shift;

  my @blocks = ();
  my $node_id=$tree->node_id;
  my $parameter_set_id = 0;
  $parameter_set_id = $params->{'parameter_set_id'} if (defined $params->{'parameter_set_id'});

  # Insert new omegas into the node.
  my $sth = $self->dbc->prepare("REPLACE INTO $output_table (aln_position,node_id,parameter_set_id,omega,omega_lower,omega_upper,ncod) values (?,?,?,?,?,?,?);");
  foreach my $hr (@{$sitewise_ref}) {
    my $ncod = Bio::EnsEMBL::Compara::AlignUtils->get_nongaps_at_column($sa,$hr->{aln_position});

    $sth->execute($hr->{aln_position},
		  $hr->{node_id},
		  $parameter_set_id,
		  $hr->{omega},
		  $hr->{omega_lower},
		  $hr->{omega_upper},
		  $ncod);
  }
  $sth->finish();
}


sub sim_stamp_shape_onto_seqs {
  my $self = shift;
  my $gap_aln = shift;
  my $seq_aln = shift;

  my $result_aln = new $seq_aln;

  my @gaps = $gap_aln->each_seq;
  my @seqs = $seq_aln->each_seq;
  # Looping through all sequences in both alignments...
  foreach my $gap (@gaps) {
    my $seq = shift @seqs;

    # The sequence order should be equal (after sorting by tree above).
    if ($gap->id ne $seq->id) {
      print "Stamper: gap aln and sequence aln don't match!  " . $gap->id."  ".$seq->id."\n";
    }
    my $gap_str = $gap->seq;
    my $seq_str = $seq->seq;

    #print "GAP> ".$gap_str."\n";
    #print "SEQ> ".$seq_str."\n";
    while ($gap_str =~ m/(-)+/g) {
      my $gap_end = pos($gap_str);
      my $gap_len = length($&);
      substr($seq_str,$gap_end-$gap_len,$gap_len) = $&;
    }
    #print "FIN> ".$seq_str."\n";

    # Create a new sequence and set the ID and seq attributes.
    my $new_seq = new $seq;
    $new_seq->id($seq->id);
    $new_seq->seq($seq_str);

    # Add the "stamped" sequence to the alignment.
    $result_aln->add_seq($new_seq);
  }
  return $result_aln;
}

1; # keep perl happy.
