use strict;
use File::Path;
use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::AlignUtils; 
use Bio::EnsEMBL::Compara::TreeUtils;
use Bio::Greg::EslrUtils;

my $cmd = qq^

#!/usr/bin/env Rscript
library(R.oo)
library(ape)
#system("make cat");
source("./FullSource.R");
####################################################
#
# Example script file showing how to set up a very simple simulation with codon sequences:
#

# Create the codon substituion process:
p<-NY98(kappa=4);

# Create a sequence of length 200:
seq<-CodonSequence(length=200); 

# Attach substitution process to root sequence:
attachProcess(seq,p);
# After attaching the NY98 process we have the M0 (equal omegas) model by default.

# Set an arbitrary omega pattern:
setOmegas(seq,p,value=c(0.5,1,2));
# Now this is not the M0 model any more.

# Copy a portion (for insert template):
ins.template<-copySubSequence(seq,1:10);

# Sample the states from the attached substitution process(es):
sampleStates(seq);

# Create the insertion process:
ins<-DiscreteInsertor(rate=0.05,sizes=c(1:5),probs=c(1,1,1,1,1));
# Create the deletion process:
del<-DiscreteDeletor(rate=0.05,sizes=c(1:5),probs=c(1,1,1,1,1));

# Set the template sequence for the insertion process:
ins.template\$processes<-list(list(p,ins,del));

# Unset write protection on insertion process:
ins\$writeProtected<-FALSE;
# Set insert template:
ins\$templateSeq<-ins.template;

# Set an insert hook for the insertion process:
ins\$insertHook<-function(seq){

  # For some reason have to specify the full method name here. To be fixed.
  setOmegas.CodonSequence(seq,p,value=c(0.5,0.6,0.7));
  return(seq);

}

# Attach the indel processes to root sequence:
seq\$processes<-list(list(p));

# Create the simulation object:
sim<-PhyloSim(
  phylo=read.tree(file="smalldemotree.nwk"),# tree as a phylo object.
  root.seq=seq# the root sequence.
  );

print(summary(sim))# Summary for the simulation object.

# Set the log file:
  sim\$logFile<-"phylosim.log";

# Set the log level to 0, only messages passed through Log() will be logged.
sim\$logLevel<-0;

# Run the simulation:
Simulate(sim);

# And finally save the resulting alignment:
saveAlignment(
  sim,  # the phylo object
  file="codon_example.fas"# filename for the saved alignment
  );


^;

Bio::Greg::EslrUtils->run_r($cmd);
