package Bio::Greg::Hive::PhyloSim;

use strict;
use File::Basename;
use File::Path;
use Time::HiRes;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::Member;
use Bio::EnsEMBL::Compara::ProteinTree;
use Bio::EnsEMBL::Compara::NestedSet;
use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Hive::Process;
use Bio::Greg::EslrUtils;

use Bio::AlignIO;

use String::CRC32;

use base ('Bio::Greg::Hive::Process');

sub param_defaults {
  my $self = shift;
  return {
    ### Omega Distribution Parameters ###
    # For lognormal model.
    #phylosim_meanlog => -4.079,
    #phylosim_sdlog => 1.23,

    # for M3 models.
    #phylosim_p0 => 0.386,
    #phylosim_p1 => 0.535,
    #phylosim_p2 => 0.079,
    #phylosim_w0 => 0.018,
    #phylosim_w1 => 0.304,
    #phylosim_w2 => 1.691,
    # etc. etc. if you want to continue the M3 model.

    # for M8 models.
    #phylosim_p0 => 0.9432,
    #phylosim_p => 0.572,
    #phylosim_q => 2.172,
    #phylosim_w => 2.081,
    ### End Omega Distribution Parameters ###

    phylosim_num_bins           => 50,
    phylosim_omega_distribution => 'M3',

    # We can start to simulate based on domains using this type of format:
    #    phylosim_domains => [
    #      {
    #        length => 30,
    #        insertrate => 0.05,
    #        deleterate => 0.05
    #      },
    #      {
    #        length => 30,
    #        insertrate => 0.001,
    #        deleterate => 0.001,
    #      },
    #      {
    #        length => 30,
    #        insertrate => 0.1,
    #        deleterate => 0.1
    #      }
    #      ],

    phylosim_simulation_program => 'indelible',
    phylosim_seq_length         => 500,
    phylosim_kappa              => 4,
    phylosim_insertrate         => 0.05,
    phylosim_deleterate         => 0.05,
    phylosim_insertmodel        => 'NB 0.2 2',
    phylosim_deletemodel        => 'NB 0.2 2',
  };
}

sub fetch_input {
  my $self = shift;

  $self->load_all_params;
}

sub run {
  my $self = shift;

  $self->param( 'alignment_table', '' );
  print "Simulation tree: " . $self->get_tree->newick_format . "\n";

  my $aln = $self->simulate_alignment($self->get_tree);

  my $length;
  foreach my $seq ($aln->each_seq) {
    $length = $seq->length unless (defined $length);
    die("Phylosim aln not all the same length!") unless ($length == $seq->length);
  }

  $self->param( 'aln', $aln );
}

sub simulate_alignment {
  my $self = shift;
  my $tree = shift;

  my $params = Bio::Greg::Hive::PhyloSim::param_defaults;
  $params = $self->replace($params, $self->params);
  $self->set_params($params);

  my $obj;
  if ( $self->param('phylosim_simulation_program') eq 'indelible' ) {
    $obj = $self->simulate_alignment_indelible( $tree, $self->params );
  } elsif ( $self->param('phylosim_simulation_program') eq 'phylosim' ) {
    $obj = $self->simulate_alignment_phylosim( $tree, $self->params );
  }
  
  return $obj;
}

sub write_output {
  my $self = shift;

  my $final_cdna = $self->param('aln');

  Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $final_cdna, { length => 150, full => 1 } );
  my $final_aa = Bio::EnsEMBL::Compara::AlignUtils->translate($final_cdna);

  my $out_table = "protein_tree_member";
  #print "STORING ALIGNMENT\n";
  Bio::EnsEMBL::Compara::ComparaUtils->store_SimpleAlign_into_table( $out_table, $self->get_tree,
    $final_aa, $final_cdna );
  #print "STORING OMEGAS\n";
  $self->_store_sitewise_omegas( "sitewise_omega", $self->param('sitewise_omegas'),
    $final_aa, $self->params );

  #Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $final_aa, { length => 150, full => 1 } );

  my $new_tree = $self->get_tree();
  my $sa = $new_tree->get_SimpleAlign();
  if (!Bio::EnsEMBL::Compara::ComparaUtils->aligned_members_equal_length($new_tree)) {
      Bio::EnsEMBL::Compara::AlignUtils->pretty_print($new_tree->get_SimpleAlign(),{full=>1,width=>150});
      throw("PhyloSim: alignment strings not same length!");
  }  

  $new_tree->release_tree;
}

sub DESTROY {
  my $self = shift;
  
  $self->param('aln','');
  delete $self->input_job->{_param_hash};

  $self->SUPER::DESTROY if $self->can("SUPER::DESTROY");
}


sub get_indelible_multiplier {
  my $self = shift;
  my $tree = shift;
  my $params = shift;

  my @bins  = $self->get_equally_spaced_bins($params);
  my @probs = $self->get_equally_spaced_probs($params);

  my $mean = 0;
  for (my $i=0; $i < scalar(@bins); $i++) {
    my $cur_bin = $bins[$i];
    my $cur_prob = $probs[$i];

    $mean += $cur_bin * $cur_prob;
  }

  my $average_omega = $mean;
  my $proportion_synonymous = 0.3; # proportion syn. sites = 0.3 when kappa=4

  my $mult = 3*($average_omega*(1-$proportion_synonymous) + $proportion_synonymous);
  print "MEAN W: $mean  INDELIBLE TREE MULTIPLIER: $mult\n";
  return $mult;
}

sub simulate_alignment_indelible {
  my $self   = shift;
  my $tree   = shift;
  my $params = shift;

  my $mult = $self->get_indelible_multiplier($tree,$params);
  #print "Before scaling: [".$tree->newick_format."]\n";
  $tree = Bio::EnsEMBL::Compara::TreeUtils->scale($tree,$mult);
  #print "After scaling: [".$tree->newick_format."]\n";

  my $models_trees_partitions = '';

  my $domains_object = $params->{phylosim_domains};
  if ( defined $domains_object ) {

    # Use the domain-parameter list to create a partitioned simulation.
    $models_trees_partitions = $self->domains_to_indelible_setup( $domains_object, $tree, $params );
  } else {
    my $length    = $params->{phylosim_seq_length};
    my $ins_rate  = $params->{phylosim_ins_rate};
    my $del_rate  = $params->{phylosim_del_rate};
    my $ins_model = $params->{phylosim_insertmodel};
    my $del_model = $params->{phylosim_deletemodel};

    my $domain = [ {
        length      => $length,
        insertrate  => $ins_rate,
        deleterate  => $del_rate,
        insertmodel => $ins_model,
        deletemodel => $del_model
      }
    ];
    $models_trees_partitions = $self->domains_to_indelible_setup( $domain, $tree, $params );
  }

  my $tmp_dir = $self->worker_temp_directory;

  my $output_f = $tmp_dir . "sim.txt";
  my $ctrl_f   = $tmp_dir . "control.txt";

  my $randomseed_string = $self->param('random_seed');
  if (!defined $randomseed_string) {
    $randomseed_string = $models_trees_partitions;
    $randomseed_string .= $params->{slrsim_rep};
    $self->param('random_seed', $randomseed_string);
  }
  my $randomseed = crc32($randomseed_string);

  my $ctrl_str = qq^
[TYPE] CODON 1
[SETTINGS]
  [output] FASTA
  [printrates] TRUE
  [randomseed] $randomseed

$models_trees_partitions

[EVOLVE] partition 1 $output_f
  ^;
  print $ctrl_str. "\n";
  open( OUT, ">$ctrl_f" );
  print OUT $ctrl_str;
  close(OUT);

  use Cwd;
  my $cwd = getcwd;
  chdir($tmp_dir);
  my @output = `indelible $ctrl_f\r\n`;

  my $aln_f = $output_f . "_TRUE.fas";

  #  if (!-e $aln_f) {
  #    open(LOG,'LOG.txt');
  #    while (<LOG>) {
  #      chomp;
  #      print "$_\n";
  #    }
  #    close(LOG);
  #  }

  #print "output:\n<<<<@output>>>>\n";

  my $aln = Bio::EnsEMBL::Compara::AlignUtils->from_file($aln_f);

  my $orig_pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($aln);
  my ($pep_aln, $new_to_old, $old_to_new) = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns($orig_pep_aln);
  $aln = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns_in_threes($aln);

  my $sitewise_hash;
  my $class_to_omega = $self->param('class_to_omega');

  # Collect the omega values.
  my $rate_f = $output_f . "_RATES.txt";
  open( IN, "$rate_f" );
  while (<IN>) {
    if (/(\d+)\t(\d+)\t(\d+)/) {
      my $site = int($1);

      # Important -- we removed gaps from the resulting alignment, so we need to map the original
      # Indelible output position (which is what we just read from the file) to that column's new
      # position in the gap-removed alignment!!
      $site = $old_to_new->{$site};

      my $nongaps = Bio::EnsEMBL::Compara::AlignUtils->get_nongaps_at_column( $pep_aln, $site );
      my $omg = $class_to_omega->{$2};
      my $type = 'negative1';
      $type = 'positive1' if ($omg > 1);

      $sitewise_hash->{$site} = {
        aln_position => $site,
        omega_lower  => $omg,
        omega_upper  => $omg,
        omega        => $omg,
        ncod => $nongaps,
        type => $type
      };
    }
  }
  close(IN);

  unlink $aln_f;
  unlink $rate_f;
  unlink $ctrl_f;
  unlink $output_f;
  chdir($cwd);

  return {
    aln => $aln,
    omegas => $sitewise_hash
  };
}

sub get_submodel_string {
  my $self   = shift;
  my $params = shift;

  my @bins  = $self->get_equally_spaced_bins($params);
  my @probs = $self->get_equally_spaced_probs($params);
  my @final_bins;
  my @final_probs;
  for ( my $i = 0 ; $i < scalar(@probs) ; $i++ ) {
    my $prob = $probs[$i];
    next unless ( $prob > 0 );    # I think this works...
    push @final_bins,  $bins[$i];
    push @final_probs, $prob;
  }
  @probs = @final_probs;
  @bins  = @final_bins;

  print "probs: @probs\n";
  print "bins: @bins\n";

  pop @probs;  # Remove the last probability; Indelible only wants (n-1) probabilities in its input!

  my $class_to_omega;
  for ( my $i = 0 ; $i < scalar(@bins) ; $i++ ) {
    $class_to_omega->{$i} = $bins[$i];
  }
  $self->param( 'class_to_omega', $class_to_omega );

  my $kappa = $params->{phylosim_kappa};
  my $submodel_str = $kappa . "\n\t" . join( " ", @probs ) . "\n\t" . join( " ", @bins );
  return $submodel_str;
}

sub domains_to_indelible_setup {
  my $self           = shift;
  my $domains_object = shift;
  my $tree           = shift;
  my $params         = shift;

  my $newick = $tree->newick_format;

  # Add some padding to zero-length branches.
  $newick =~ s/(:0\.?0+)([;,()])/:0.0005$2/g;

  # Get rid of branch length on root.
  $newick =~ s/:\d\.?\d+;/;/g;

  my $defaults = {
    insertmodel => $params->{phylosim_insertmodel},
    deletemodel => $params->{phylosim_deletemodel},
    insertrate  => $params->{phylosim_insertrate},
    deleterate  => $params->{phylosim_deleterate}
  };

  my $str = "";
  $str .= "[TREE] tree $newick\n";

  # Print out the model definition.
  my $i = 0;
  foreach my $obj (@$domains_object) {
    my $model_name = "model" . $i++;
    $str .= "[MODEL] $model_name\n";
    foreach my $key (qw^insertmodel deletemodel insertrate deleterate^) {
      my $value = $defaults->{$key};
      $value = $obj->{$key} if ( defined $obj->{$key} );
      $str .= "  [$key]\t$value\n";
    }
    $str .= "  [submodel] " . $self->get_submodel_string($params) . "\n";
  }

  my $i = 0;
  $str .= "[PARTITIONS] partition\n";
  foreach my $obj (@$domains_object) {
    my $length     = $obj->{length};
    my $model_name = "model" . $i++;
    $str .= "  [tree $model_name $length]\n";
  }
  return $str;
}

sub simulate_alignment_phylosim {
  my $self   = shift;
  my $tree   = shift;
  my $params = shift;

  my $aln_out = $self->worker_temp_directory . "sim.fas";

  my $newick    = $tree->newick_format;
  my $omega_var = "omegaVarM7(seq,p,p0=0.5,p=1,q=2,mean=1.5,sd=0.5);";

  my $rcmd = qq^
library(R.oo);
library(ape);
source("/homes/greg/lib/phylosim/FullSource.R");

p<-NY98(kappa=2);
seq<-CodonSequence(length=100);
attachProcess(seq,p);
$omega_var

sampleStates(seq);

sim<-PhyloSim(
  phylo=read.tree(text="$newick"),
  root.seq=seq
);

Simulate(sim);

saveAlignment(
  sim,
  file="$aln_out"
);

^;

  Bio::Greg::EslrUtils->run_r($rcmd);

}

sub _store_sitewise_omegas {
  my $self         = shift;
  my $output_table = shift;
  my $sitewise_ref = shift;
  my $sa           = shift;
  my $params       = shift;

  my $tree = $self->get_tree;

  my @blocks           = ();
  my $data_id          = $self->data_id;
  my $node_id          = $self->node_id;
  my $parameter_set_id = $self->parameter_set_id;

  #  eval {

  my @insert_strings = ();
  foreach my $hr ( @{$sitewise_ref} ) {

    #print "hr: $hr\n";
    #printf "%d %d\n",$hr->{aln_position},$hr->{omega};

    my $type;
    if ( $hr->{omega} <= 1 ) {
      $type = "negative1";
    } else {
      $type = "positive1";
    }

    my $ncod = Bio::EnsEMBL::Compara::AlignUtils->get_nongaps_at_column( $sa, $hr->{aln_position} );

    push @insert_strings,
      sprintf(
      '(%d,%d,%d,%d,%.5f,%.5f,%.5f,"%s",%d)',
      $hr->{aln_position}, $node_id,     $data_id,
      $parameter_set_id,   $hr->{omega},           $hr->{omega}, $hr->{omega},
      $type,               $ncod
      );

    if ( scalar(@insert_strings) >= 1000 ) {
      my $insert = join( ",", @insert_strings );
      my $cmd =
        "INSERT INTO $output_table (aln_position,node_id,data_id,parameter_set_id,omega,omega_lower,omega_upper,type,ncod) values $insert;";
 #     print "$cmd\n";
      $self->dbc->do($cmd);
      @insert_strings = ();
    }
  }

  if ( scalar(@insert_strings) > 0 ) {
    my $insert = join( ",", @insert_strings );
    my $cmd =
      "INSERT INTO $output_table (aln_position,node_id,data_id,parameter_set_id,omega,omega_lower,omega_upper,type,ncod) values $insert;";
#    print "$cmd\n";
    $self->dbc->do($cmd);
  }

  #  };
}

sub get_m3_bins {
  my $self   = shift;
  my $params = shift;

  my $i = 0;
  my @bins;
  while ( $self->param( 'phylosim_w' . $i ) ) {
    push @bins, $self->param( 'phylosim_w' . $i );
    $i++;
  }
  return @bins;
}

sub get_m3_probs {
  my $self   = shift;
  my $params = shift;

  my $i = 0;
  my @probs;
  while ( $self->param( 'phylosim_p' . $i ) ) {
    push @probs, $self->param( 'phylosim_p' . $i );
    $i++;
  }
  return @probs;
}

sub get_equally_spaced_bins {
  my $self   = shift;
  my $params = shift;

  if ( $self->param('phylosim_omega_distribution') =~ m/M3/i ) {
    return $self->get_m3_bins($params);
  }

  my $k  = $self->param('phylosim_num_bins');
  my $lo = 0;
  my $hi = 3;

  my $r_cmd = qq^
   x = seq(from=$lo,to=$hi,by=($hi-$lo)/($k-1));
   print(sprintf("%.4f",x));
^;
  my @nums = $self->get_r_numbers($r_cmd);
  return @nums;
}

sub get_equally_spaced_probs {
  my $self   = shift;
  my $params = shift;

  my $function = $self->param('phylosim_omega_distribution');
  my $k        = $self->param('phylosim_num_bins');
  my $lo       = 0;
  my $hi       = 3;

  my $r_pre = qq^

ppaml_m8 = function(x,p0,p,q,w) {
  y = pbeta(x,shape1=p,shape2=q)*p0
  y[x>=w] = y[x>=w] + (1-p0)
  return(y)
}
pconst = function(x,w) {
  y = rep(0,length(x))
  y[x>=w] = 1
  return(y)
}
  ^;
  my $f = "";

  if ( $function =~ m/M3/i ) {
    return $self->get_m3_probs($params);
  } elsif ( $function =~ m/M8/i ) {
    my $p0 = $self->param('phylosim_p0');
    my $p  = $self->param('phylosim_p');
    my $q  = $self->param('phylosim_q');
    my $w  = $self->param('phylosim_w');
    $f = qq^paml_m8(x,p0=${p0},p=$p,q=$q,w=$w)^;
  } elsif ( $function eq "constant" ) {
    my $w = $self->param('phylosim_w');
    $f = qq^const(x,w=$w)^;
  } elsif ( $function eq "uniform" ) {
    my $min = $self->param('phylosim_min');
    my $max = $self->param('phylosim_max');
    $f = qq^unif(x,min=$min,max=$max)^;
  } elsif ( $function eq "gamma" ) {
    my $shape = $self->param('phylosim_shape');
    my $rate  = $self->param('phylosim_rate');
    $f = qq^gamma(x,shape=$shape,rate=$rate)^;
  } elsif ( $function eq "beta" ) {
    my $shape1 = $self->param('phylosim_shape1');
    my $shape2 = $self->param('phylosim_shape2');
    $f = qq^beta(x,shape1=$shape1,shape2=$shape2)^;
  } elsif ( $function eq "lognormal" ) {
    my $meanlog = $self->param('phylosim_meanlog');
    my $sdlog   = $self->param('phylosim_sdlog');
    $f = qq^lnorm(x,meanlog=$meanlog,sdlog=$sdlog)^;
  }

  my $r_cmd = qq^
   $r_pre
   
   x = seq(from=$lo,to=$hi,by=($hi-$lo)/($k-1))
   y = p$f
   z = diff(y)
   x = x[length(x)]
   z = c(z,1-p$f)
   print(sprintf("%.6f",z))
^;
  my @nums = $self->get_r_numbers($r_cmd);
  return @nums;
}

sub get_r_numbers {
  my $self  = shift;
  my $r_cmd = shift;

  my @lines = Bio::Greg::EslrUtils->get_r_values( $r_cmd, $self->worker_temp_directory );
  my $line = join( " ", @lines );

  #print "values:". $line."\n";
  #$line =~ s/\[.+?\]//g;
  $line =~ s/"//g;
  my @toks = split( " ", $line );
  return @toks;
}

# Calls R to extract the equiprobable bins for a given distribution.
# GJ 20010-01-15 this is currently NOT used!!!
sub get_equiprobable_bins {
  my $self   = shift;
  my $params = shift;

  my $function = $self->param('phylosim_omega_distribution');
  my $k        = $self->param('phylosim_num_bins');

  my $r_cmd = qq^
    dd = discretize_distribution(n=$k)
    #means = dd\$means
    print(sprintf("%.4f",dd\$means));
    ^;

  my $r_pre = "";
  if ( $function eq "uniform" ) {
    my $a = $self->param('phylosim_min');
    my $b = $self->param('phylosim_max');

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
  } elsif ( $function eq "gamma" ) {
    my $a = $self->param('phylosim_shape');
    my $b = $self->param('phylosim_rate');

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
  } elsif ( $function eq "beta" ) {
    my $p = $self->param('phylosim_shape1');
    my $q = $self->param('phylosim_shape2');

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
  } elsif ( $function eq "lognormal" ) {
    my $meanlog = $self->param('phylosim_meanlog');
    my $sdlog   = $self->param('phylosim_sdlog');
    $r_pre = qq^
      meanlog = $meanlog;
      sdlog = $sdlog;
      p1 = function(x) {
        return(plnorm(x,meanlog=meanlog,sdlog=sdlog));
      }
      d1 = function(x) {
        return(dlnorm(x,meanlog=meanlog,sdlog=sdlog));
      }
      q1 = function(x) {
        return(qlnorm(x,meanlog=meanlog,sdlog=sdlog));
      }
    ^;
  }

  my $r_prefix = qq^
    $r_pre
    discretize_distribution = function(n,p=p1,d=d1,q=q1) {
      num_categories = n;
	quantile_bounds = seq(from=0,to=1,length.out=num_categories+1)
	x_bounds = q(quantile_bounds)
	
	if(is.infinite(x_bounds[1])) {
	  x_bounds[1] = -100
	}
      if(is.infinite(x_bounds[num_categories+1])) {
	x_bounds[num_categories+1] = x_bounds[num_categories]*2;
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

  my @toks = $self->get_r_numbers($r_cmd);
  return @toks;
}

1;
