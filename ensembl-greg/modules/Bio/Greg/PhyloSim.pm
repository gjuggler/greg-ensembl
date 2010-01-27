package Bio::Greg::PhyloSim;

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
use Bio::Greg::ProcessUtils;

use Bio::AlignIO;

our @ISA = qw(Bio::EnsEMBL::Hive::Process Bio::Greg::ProcessUtils);

my $p_pre = "phylosim_";

my $dba;
my $pta;

my $params;

my $tree;
my @sitewise_omegas;
my $aln;

sub fetch_input {
  my $self = shift;

  # Load up the Compara DBAdaptor.
  if ($self->{dba}) {
    $dba = $self->dba;
  } else {
    $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-DBCONN=>$self->db->dbc);
  }
  $pta = $dba->get_ProteinTreeAdaptor;


  $params = {
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

    phylosim_num_bins => 20,
    phylosim_omega_distribution => 'M3',

    phylosim_seq_length => 500,
    phylosim_kappa => 4,
    phylosim_ins_rate => 0.05,
    phylosim_del_rate => 0.05,
    phylosim_simulation_program => 'indelible',
    phylosim_insertmodel => 'NB 0.3 2',
    phylosim_deletemodel => 'NB 0.3 2',
  };

  my $node_id;
  
  my $p_params = $self->get_params($self->parameters);
  my $i_params = $self->get_params($self->input_id);
  $node_id = $i_params->{'protein_tree_id'};
  $node_id = $i_params->{'node_id'} if (!defined $node_id);
  my $t_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_tree_tags($dba,$node_id);

  $params = $self->replace_params($params,$p_params,$i_params,$t_params);
  Bio::EnsEMBL::Compara::ComparaUtils->hash_print($params);

  # Load the tree.
  $pta->table_base($params->{'input_table_base'});
  $tree = $pta->fetch_node_by_node_id($node_id);

  print "Simulation tree: ".$tree->newick_format."\n";
}


sub run {
  my $self = shift;

  if (get('simulation_program') eq 'indelible') {
    $aln = $self->simulate_alignment_indelible($tree,$params);
    #print "ALN : $aln\n";
  } elsif (get('simulation_program') eq 'phylosim') {
    $aln = $self->simulate_alignment_phylosim($tree,$params);
  }
  
}

sub write_output {
  my $self = shift;

  my $final_cdna = $aln;
  Bio::EnsEMBL::Compara::AlignUtils->pretty_print($final_cdna,{length => 200});
  my $final_aa = Bio::EnsEMBL::Compara::AlignUtils->translate($final_cdna);
  
  my $out_table = "protein_tree_member";
  print "STORING ALIGNMENT\n";
  Bio::EnsEMBL::Compara::ComparaUtils->store_SimpleAlign_into_table($out_table,$tree,$final_aa,$final_cdna);
  print "STORING OMEGAS\n";
  $self->_store_sitewise_omegas("sitewise_omega",\@sitewise_omegas,$final_aa,$params);

  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
  my $aa_aln = $tree->get_SimpleAlign();
  my $cdna_aln = $tree->get_SimpleAlign(-cdna => 1);
  $cdna_aln = Bio::EnsEMBL::Compara::AlignUtils->sort_by_tree($cdna_aln,$treeI);
  my $alnout = Bio::AlignIO->new
    ('-format'      => 'phylip',
     '-file'          => $self->worker_temp_directory,
     '-interleaved' => 0,
     '-idlinebreak' => 1,
     '-idlength'    => $cdna_aln->maxdisplayname_length + 1);
  $alnout->write_aln($cdna_aln);
  $alnout->close();

  #Bio::EnsEMBL::Compara::AlignUtils->pretty_print($aa_aln,{length => 1500});
}

sub throw_error {
  my $self = shift;
  my $msg = shift;

  print "ERROR!!!\n";
  open(OUT,">/tmp/".$tree->node_id.".txt");
  print OUT $msg."\n";
  close(OUT);
}

sub simulate_alignment_indelible {
  my $self = shift;
  my $tree = shift;
  my $params = shift;

  #$tree = Bio::EnsEMBL::Compara::TreeUtils->scale($tree,1.95);

  my $newick = $tree->newick_format;
  # Add some padding to zero-length branches.
  $newick =~ s/(:0\.?0+)([;,()])/:0.0005$2/g;
  # Get rid of branch length on root.
  $newick =~ s/:\d\.?\d+;/;/g;
  print $newick."\n";

  my $length = get('seq_length');
  if ($length eq 'sample') {
    ($length) = $self->get_r_numbers("rlnorm(1,meanlog=6.03,sdlog=0.685);");
    $length = int($length);
  }
  my $ins_rate = get('ins_rate');
  my $del_rate = get('del_rate');

  print "Simulating sequence with $length codons...\n";
  
  my $tmp_dir = $self->worker_temp_directory . $tree->node_id . '/';
  print "Temp dir: $tmp_dir\n";
  #use File::Path qw(mkpath rmtree);
  mkpath([$tmp_dir]);

  my $output_f = $tmp_dir."sim.txt";
  my $ctrl_f = $tmp_dir."control.txt";

  my @bins = $self->get_equally_spaced_bins($params);
  my @probs = $self->get_equally_spaced_probs($params);
  my @final_bins;
  my @final_probs;
  for (my $i=0; $i < scalar(@probs); $i++) {
    my $prob = $probs[$i];
#    if ($prob > 0) {
      push @final_bins, $bins[$i];
      push @final_probs, $prob;
#    }
  }
  @probs = @final_probs;
  @bins = @final_bins;
  print "probs: @probs\n";
  print "bins: @bins\n";
  pop @probs; # Remove the last probability; Indelible only wants (n-1) probabilities in its input!
  
  my $kappa = get('kappa');
  my $submodel_str = $kappa."\n\t". join(" ",@probs) . "\n\t" . join(" ",@bins);

  my $class_to_omega;
  for (my $i=0; $i < scalar(@bins); $i++) {
    $class_to_omega->{$i} = $bins[$i];
  }

  my $node_id = $tree->node_id . time();

  my $ctrl_str = qq^
[TYPE] CODON 1
[SETTINGS]
  [output] FASTA
  [printrates] TRUE
  [randomseed] $node_id

[MODEL] model1
  [submodel] $submodel_str
  [insertmodel] NB 0.35 1
  [deletemodel] NB 0.35 1
  [insertrate] $ins_rate
  [deleterate] $del_rate

[TREE] tree1 $newick


[PARTITIONS] partition1
  [tree1 model1 $length]
  [EVOLVE] partition1 1 $output_f
  ^;
  print $ctrl_str."\n";
  open(OUT,">$ctrl_f");
  print OUT $ctrl_str;
  close(OUT);

  use Cwd;
  my $cwd = getcwd;
  chdir($tmp_dir);
  my $output = `indelible $ctrl_f\r\n`;

  my $aln_f = $output_f."_TRUE.fas";
#  if (!-e $aln_f) {
#    open(LOG,'LOG.txt');
#    while (<LOG>) {
#      chomp;
#      print "$_\n";
#    }
#    close(LOG);
#  }

  my $aln = Bio::EnsEMBL::Compara::AlignUtils->from_file($aln_f);

  print "$output\n";

  # Collect the omega values.
  my $rate_f = $output_f."_RATES.txt";
  open(IN,"$rate_f");
  while (<IN>) {
    if (/(\d+)\t(\d+)\t(\d+)/) {
      my $omg = $class_to_omega->{$2};
      push @sitewise_omegas, {
        aln_position => int($1), node_id => $tree->node_id,
        omega_lower => $omg, omega_upper => $omg, omega => $omg
      };
    }
  }
  close(IN);

  unlink $aln_f;
  unlink $rate_f;
  unlink $ctrl_f;
  unlink $output_f;
  chdir($cwd);

  return $aln;
}

sub simulate_alignment_phylosim {
  my $self = shift;
  my $tree = shift;
  my $params = shift;

  my $aln_out = $self->worker_temp_directory."sim.fas";

  my $newick = $tree->newick_format;
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
  my $self = shift;
  my $output_table = shift;
  my $sitewise_ref = shift;
  my $sa = shift;
  my $params = shift;

  my @blocks = ();
  my $node_id=$tree->node_id;
  my $parameter_set_id = 1;
  $parameter_set_id = $params->{'parameter_set_id'} if (defined $params->{'parameter_set_id'});

 
  eval {
 # Insert new omegas into the node.
#  my $sth = $self->dbc->prepare("INSERT IGNORE INTO $output_table (aln_position,aln_position_fraction,node_id,parameter_set_id,omega,omega_lower,omega_upper,type,ncod) values (?,?,?,?,?,?,?,?,?);");
  #print "sw: $sitewise_ref\n";
#  $self->dbc->do("LOCK TABLES $output_table WRITE");
  my @insert_strings = ();
  foreach my $hr (@{$sitewise_ref}) {
    #print "hr: $hr\n";
    #printf "%d %d\n",$hr->{aln_position},$hr->{omega};
    my $ncod = Bio::EnsEMBL::Compara::AlignUtils->get_nongaps_at_column($sa,$hr->{aln_position});

    my $type;
    if ($hr->{omega} <= 1) {
      $type = "negative1";
    } else {
      $type = "positive1";
    }

    my $aln_position_fraction = $hr->{aln_position} / $sa->length;

    push @insert_strings, sprintf('(%d,%.5f,%d,%d,%.5f,%.5f,%.5f,"%s",%d)',
                              $hr->{aln_position},
                              $aln_position_fraction,
                              $hr->{node_id},
                              $parameter_set_id,
                              $hr->{omega},
                              $hr->{omega_lower},
                              $hr->{omega_upper},
                              $type,
                              $ncod);
#    $sth->execute(
#$hr->{aln_position},
#                  $aln_position_fraction,
#		  $hr->{node_id},
#		  $parameter_set_id,
#		  $hr->{omega},
#		  $hr->{omega_lower},
#		  $hr->{omega_upper},
#                  $type,
#		  $ncod);
#    sleep(0.1);
  }
  
  my $insert = join(",",@insert_strings);
  $self->dbc->do("INSERT IGNORE INTO $output_table (aln_position,aln_position_fraction,node_id,parameter_set_id,omega,omega_lower,omega_upper,type,ncod) values $insert ;");
#  $sth->finish();
#  $self->dbc->do("UNLOCK TABLES;");
  };
  if ($@) {
    open(OUT,">/homes/greg/tmp/".$tree->node_id.".txt");
    print OUT $@."\n";
    close(OUT);
    die($@);
  }
}

sub get_m3_bins {
  my $self = shift;
  my $params = shift;
  
  my $i=0;
  my @bins;
  while (get('w'.$i)) {
    push @bins,get('w'.$i);
    $i++;
  }
  return @bins;
}

sub get_m3_probs {
  my $self = shift;
  my $params = shift;

  my $i=0;
  my @probs;
  while (get('p'.$i)) {
    push @probs,get('p'.$i);
    $i++;
  }
  return @probs;
}

sub get_equally_spaced_bins {
  my $self = shift;
  my $params = shift;

  if (get('omega_distribution') =~ m/M3/i) {
    return $self->get_m3_bins($params);
  }

  my $k = get('num_bins');
  my $lo = 0;
  my $hi = 3;

  my $r_cmd = qq^
   x = seq(from=$lo,to=$hi,by=($hi-$lo)/($k-1));
   print(sprintf("%.4f",x));
^;
  my @nums = $self->get_r_numbers($r_cmd);
  return @nums;
}

sub get {
  my $key = shift;
  return $params->{$p_pre.$key};
}

sub get_equally_spaced_probs {
  my $self = shift;
  my $params = shift;

  my $function = get('omega_distribution');
  my $k = get('num_bins');
  my $lo = 0;
  my $hi = 3;

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

  if ($function =~ m/M3/i) {
    return $self->get_m3_probs($params);
  } elsif ($function =~ m/M8/i) {
    my $p0 = get('p0');
    my $p = get('p');
    my $q = get('q');
    my $w = get('w');
    $f = qq^paml_m8(x,p0=${p0},p=$p,q=$q,w=$w)^;
  } elsif ($function eq "constant") {
    my $w = get('w');
    $f = qq^const(x,w=$w)^;
  } elsif ($function eq "uniform") {
    my $min = get('min');
    my $max = get('max');
    $f = qq^unif(x,min=$min,max=$max)^;
  } elsif ($function eq "gamma") {
    my $shape = get('shape');
    my $rate = get('rate');
    $f = qq^gamma(x,shape=$shape,rate=$rate)^;
  } elsif ($function eq "beta") {
    my $shape1 = get('shape1');
    my $shape2 = get('shape2');
    $f = qq^beta(x,shape1=$shape1,shape2=$shape2)^;
  } elsif ($function eq "lognormal") {
    my $meanlog = get('meanlog');
    my $sdlog = get('sdlog');
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
  my $self = shift;
  my $r_cmd = shift;

  my @lines = Bio::Greg::EslrUtils->get_r_values($r_cmd,$self->worker_temp_directory);
  my $line = join(" ",@lines);
#  print $line."\n";
  $line =~ s/\[.+?\]//g;
  $line =~ s/"//g;
  my @toks = split(" ",$line);
  return @toks;
}

# Calls R to extract the equiprobable bins for a given distribution.
# GJ 20010-01-15 this is currently NOT used!!!
sub get_equiprobable_bins {
  my $self = shift;
  my $params = shift;

  my $function = get('omega_distribution');
  my $k = get('num_bins');

  my $r_cmd = qq^
    dd = discretize_distribution(n=$k)
    #means = dd\$means
    print(sprintf("%.4f",dd\$means));
    ^;
  
  my $r_pre = "";
  if ($function eq "uniform") {
    my $a = get('min');
    my $b = get('max');
    
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
    my $a = get('shape');
    my $b = get('rate');

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
    my $p = get('shape1');
    my $q = get('shape2');

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
  } elsif ($function eq "lognormal") {
    my $meanlog = get('meanlog');
    my $sdlog = get('sdlog');
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
