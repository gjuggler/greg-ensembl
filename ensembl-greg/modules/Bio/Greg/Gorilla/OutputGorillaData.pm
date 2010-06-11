package Bio::Greg::Gorilla::OutputGorillaData;

use strict;
use Time::HiRes qw(sleep);

use POSIX qw(strftime mktime);
use Cwd;
use File::Path;

use Bio::EnsEMBL::Hive::Process;

use Bio::Greg::EslrUtils;

use base ('Bio::Greg::Hive::Process');

sub run {
  my $self = shift;

  my $base = '/nfs/users/nfs_g/gj1/scratch/gorilla/output';
  $self->get_output_folder($base);

  $self->export_topologies;
  $self->export_likelihoods;

  #$self->export_counts;
  #$self->analyze_counts;
  #$self->get_enriched_genes;
}

sub export_topologies {
  my $self = shift;

  my $data_file = $self->get_output_folder . "/stats_topology.Rdata";
  my $folder = $self->get_output_folder;

  if (!-e $data_file) {
    my $cmd = qq^
dbname="gj1_gor_57"
source("../../scripts/collect_sitewise.R");

gcinfo(TRUE)
stats.topology <- get.vector(con,"SELECT * from stats_topology;")
save(stats.topology,file="${data_file}")
^;    
#    print "$cmd\n";
    my $params = {};
    Bio::Greg::EslrUtils->run_r($cmd,$params);
  }   
}

sub export_likelihoods {
  my $self = shift;

  my $data_file = $self->get_output_folder . "/stats_lnl.Rdata";
  my $folder = $self->get_output_folder;
  
  my $force = 0;
  if (!-e $data_file || $force) {
    my $cmd = qq^
dbname="gj1_gor_57"
source("../../scripts/collect_sitewise.R");

stats.lnl <- get.vector(con,"SELECT * from stats_lnl;")

# Likelihoods key:
# a: (H, G, others)
# b: (H#1, G, others)
# c: (H, G#1, others)
# d: (H#1, G#2, others)
# e: (H#1, G#1, others)
# Which omegas are which for each test?
# pval.1: fg=b_omega_1, bg=b_omega_0 # human
# pval.2: fg=c_omega_1, bg=c_omega_0 # gorilla
# pval.3: fg=d_omega_2, bg=d_omega_0 [and d_omega_1] # human
# pval.4: fg=d_omega_1, bg=d_omega_0 [and d_omega_2] # gorilla
# pval.5: fg=e_omega_1, bg=e_omega_0 # both

# Add the p-values
stats.lnl[,'pval.1'] = with(stats.lnl,1 - pchisq(2*(b_lnL-a_lnL),df=1))
stats.lnl[,'pval.2'] = with(stats.lnl,1 - pchisq(2*(c_lnL-a_lnL),df=1))
stats.lnl[,'pval.3'] = with(stats.lnl,1 - pchisq(2*(d_lnL-b_lnL),df=1))
stats.lnl[,'pval.4'] = with(stats.lnl,1 - pchisq(2*(d_lnL-c_lnL),df=1))
stats.lnl[,'pval.5'] = with(stats.lnl,1 - pchisq(2*(e_lnL-a_lnL),df=1))

# Save the data
save(stats.lnl,file="${data_file}")

^;    
#    print "$cmd\n";
    my $params = {};
    Bio::Greg::EslrUtils->run_r($cmd,$params);
    
  }   
}


sub export_counts {
  my $self = shift;

  my $counts_file = $self->get_output_folder . "/gorilla_counts.Rdata";
  my $test_file = $self->get_output_folder . "/gorilla_counts_test.Rdata";
  my $folder = $self->get_output_folder;

  if (!-e $counts_file || !-e $test_file) {
    my $cmd = qq^
dbname="gj1_gor_57"
source("../../scripts/collect_sitewise.R");
gcinfo(TRUE)
counts.sites <- get.vector(con,"SELECT * from outgroup_sites;")
counts.genes <- get.vector(con,"SELECT * from outgroup_genes;")
stats.genes <- get.vector(con,"SELECT * from stats_genes;")
save(counts.sites,counts.genes,stats.genes,file="${counts_file}")

counts.sites = counts.sites[1:100000,]
save(counts.sites,counts.genes,stats.genes,file="${test_file}")
^;
    print "$cmd\n";
    my $params = {};
    Bio::Greg::EslrUtils->run_r($cmd,$params);
  }

  my $genes_file = $self->get_output_folder . "/genes_merged.Rdata";
  if (!-e $genes_file) {
    my $cmd = qq^
dbname="gj1_gor_57"
source("../../scripts/collect_sitewise.R");
genes.merged = get.genes.merged()
save(genes.merged,file="${genes_file}")
^;
    print "$cmd\n";
    my $params = {};
    Bio::Greg::EslrUtils->run_r($cmd,$params);
  }

}

sub analyze_counts {
  my $self = shift;

  my $counts_file = $self->get_output_folder . "/gorilla_counts.Rdata";
  my $output_table = $self->get_output_folder . "/top_genes.csv";
  my $output_data = $self->get_output_folder . "/top_genes.Rdata";
  
  my $cmd = qq^
gcinfo(TRUE)
source("analyze_gorilla.R")
load("$counts_file")
require(plyr)
require(doBy)
summary.fn = function(data) {
  df = data
  df = subset(df,has_gap==0) # Remove sites with any gaps.
  df = subset(df,n_nucleotide_diffs<=1) # Remove sites with 2 or more nucleotide diffs.
  df = subset(df,mut_cpg==0) # Remove sites where a CpG dinucleotide is mutated.

  n.sites = nrow(df) # Get all sites, even those without mutations.
  n.g.muts = nrow(subset(df,species_a=="G")) # All sites with the one-off mutant is gorilla.
  n.h.muts = nrow(subset(df,species_a=="H")) # "" human.
  n.c.muts = nrow(subset(df,species_a=="C")) # "" chimp.
  return(
    data.frame(
      n.sites=n.sites,
      g.muts=n.g.muts,
      h.muts=n.h.muts,
      c.muts=n.c.muts
    )
  )
}

sites = counts.sites
mut.counts = ddply(sites,c("data_id"),summary.fn)
mut.counts = subset(mut.counts,n.sites > 50) # Filter out genes with tiny numbers of sites!

total.sites = sum(mut.counts[,'n.sites'])
total.c.muts = sum(mut.counts[,'c.muts'])
total.g.muts = sum(mut.counts[,'g.muts'])
total.h.muts = sum(mut.counts[,'h.muts'])

print(paste(total.c.muts/total.sites,total.g.muts/total.sites,total.h.muts/total.sites))

c.length = total.c.muts/total.sites
g.length = total.g.muts/total.sites
h.length = total.h.muts/total.sites

mut.counts[,'g.muts.corrected'] = mut.counts[,'g.muts'] / g.length / mut.counts[,'n.sites']
mut.counts[,'h.muts.corrected'] = mut.counts[,'h.muts'] / h.length / mut.counts[,'n.sites']
mut.counts[,'c.muts.corrected'] = mut.counts[,'c.muts'] / c.length / mut.counts[,'n.sites']

merged = merge(mut.counts,counts.genes,by=c('data_id'))
top.genes = orderBy(~-g.muts.corrected,data=merged)
save(top.genes,file="${output_data}")
write.csv(top.genes,file="${output_table}",row.names=F)
^;
  print "$cmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r($cmd,$params);

}


sub get_enriched_genes {
  my $self = shift;

  my $f = $self->get_output_folder;
  my $genes_csv = $f . "/top_genes.csv";
  my $genes_merged = $f ."/genes_merged.Rdata";
  my $scripts_dir = $self->base . "/scripts";
  my $enrich_file = $f . "/enriched.csv";

  my $cmd = qq^
dbname="gj1_gor_57"
source("${scripts_dir}/collect_sitewise.R")
source("${scripts_dir}/go_enrichments.R")

top.genes = read.csv("${genes_csv}")
load("${genes_merged}")
merged = merge(top.genes,genes.merged,by=c('node_id','human_gene'))

all.go = merged
print(nrow(all.go))
sub.go = subset(merged,g.muts.corrected > 5 & m_slr.dnds < 0.7)
print(nrow(sub.go))

enriched = get.enrich.df(sub.go[,'human_protein'],all.go[,'human_protein'],go.hs,go.field.name='stable_id')
write.csv(enriched,file="${enrich_file}")
^;
  print "$cmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r($cmd,$params);

}

1;
