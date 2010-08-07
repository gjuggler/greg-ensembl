package Bio::Greg::Gorilla::OutputGorillaData;

use strict;
use Time::HiRes qw(sleep);

use POSIX qw(strftime mktime);
use Cwd;
use File::Path;

use Bio::EnsEMBL::Hive::Process;

use Bio::Greg::EslrUtils;
use Bio::EnsEMBL::Compara::ComparaUtils;

use base ('Bio::Greg::Hive::Process');

sub run {
  my $self = shift;

  my $base = '/nfs/users/nfs_g/gj1/scratch/gorilla/output';
  $self->get_output_folder($base);

  $self->genomic_align_test;

  #$self->export_likelihoods('stats_branch','stats_branch_human.Rdata','human');
  #$self->export_likelihoods('stats_chimp_branch','stats_branch_chimp.Rdata','chimp');

  #$self->branch_enrichments('stats_branch','stats_branch_human.Rdata','human');
  #$self->branch_enrichments('stats_chimp_branch','stats_branch_chimp.Rdata','chimp');

  #$self->bryndis_enrichments;
  #$self->bryndis_dup_counts;

  #$self->export_counts;
  #$self->analyze_counts;
  #$self->get_enriched_genes;
}

sub genomic_align_test {
  my $self = shift;
  
  my $node_id = 1686351;

  my $tree = $self->pta->fetch_node_by_node_id($node_id);

  my ($gen_cdna,$gen_aa) = Bio::EnsEMBL::Compara::ComparaUtils->genomic_align_for_tree($tree,9606);

  my $params = $self->params;
  $params->{tree} = $tree;
  my $hom_cdna = $self->get_cdna_aln($params);
  my $hom_aa = $self->get_aln($params);
  ($hom_aa) = Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns($hom_aa);

  Bio::EnsEMBL::Compara::AlignUtils->pretty_print($gen_aa);
  Bio::EnsEMBL::Compara::AlignUtils->pretty_print($hom_aa);

}

sub bryndis_dup_counts {
  my $self = shift;

  my $folder = $self->get_output_folder;
  my $base = Bio::Greg::EslrUtils->baseDirectory;
  my $eslr_script = $self->base."/scripts/collect_sitewise.R";
  my $bryndis_output = "$folder/gene_copy_counts.csv";

  my $cmd = qq^
dbname = 'gj1_gor_58'
source("${eslr_script}",echo=F)
# Get gene copy count data and save as CSV.
stats.dups <- get.vector(con,"SELECT * from stats_dups;")
write.csv(stats.dups,file="$bryndis_output")
^;
  Bio::Greg::EslrUtils->run_r($cmd,{});
}


sub bryndis_enrichments {
  my $self = shift;

  my $folder = $self->get_output_folder;
  my $base = Bio::Greg::EslrUtils->baseDirectory;
  my $bryndis_list = "${base}/projects/gorilla/bryndis_genes.txt";
  my $go_script = Bio::Greg::EslrUtils->baseDirectory."/scripts/go_enrichments.R";
  my $eslr_script = Bio::Greg::EslrUtils->baseDirectory."/scripts/collect_sitewise.R";

  my $cmd = qq^
dbname = 'gj1_gor_58'
source("${eslr_script}",echo=F)

bryndis.genes <- read.csv("${bryndis_list}",header=T)
names(bryndis.genes) = c('gorilla_gene')
stats.dups <- get.vector(con,"SELECT * from stats_dups;")
stats.branch <- get.vector(con,"SELECT * from stats_branch;")

source("${go_script}")

stats.branch[,'stable_id'] = stats.branch[,'human_protein']
interesting.genes = merge(bryndis.genes,stats.branch,by='gorilla_gene')
go.df = go.hs

all.ids = stats.branch[,'stable_id']
interesting.ids = interesting.genes[,'stable_id']

get.df.subset = function(subset) {
  print(paste("size: ",nrow(subset)))
  return(get.enrich.by.subset(
    subset = subset,
    all = all.ids,
    go.df = go.df,
    go.field.name = 'stable_id',
    nodeSize = 3
  ))
}

get.go.data = function(subset,all) {
  geneSelectionFun = function(score){return(score >= 1)}

  go.vec = strsplit(go.df[,'go'],split=",",fixed=T)
  names(go.vec) = go.df[,'stable_id']

  scores = as.integer(all %in% subset)
  names(scores) = all

  GOdata <- new("topGOdata",
    ontology = 'BP',
    allGenes = scores,
    annot = annFUN.gene2GO,
    gene2GO = go.vec,
    geneSelectionFun = geneSelectionFun,
    nodeSize = 1,
    description = ''
  )
  return(GOdata)
}

go.data = get.go.data(interesting.ids,all.ids)
enrichment.table = get.df.subset(interesting.ids)

save(go.data,interesting.genes,enrichment.table,file="${folder}/bryndis_enrichments.Rdata")
write.csv(enrichment.table,file="${folder}/bryndis_enrichments.csv",row.names=F)

^;

    my $params = {};
    Bio::Greg::EslrUtils->run_r($cmd,$params);
}

sub export_likelihoods {
  my $self = shift;
  my $table = shift;
  my $output_filename = shift;
  my $short_prefix = shift;

  my $data_file = $self->get_output_folder . "/$output_filename";

  my $folder = $self->get_output_folder;
  mkpath([$folder]);
  print "FOLDER: $folder\n";

  my $collect_script = Bio::Greg::EslrUtils->baseDirectory."/scripts/collect_sitewise.R";
  my $lrt_script = Bio::Greg::EslrUtils->baseDirectory."/projects/gorilla/lrt_analysis.R";
  my $go_script = Bio::Greg::EslrUtils->baseDirectory."/scripts/go_enrichments.R";
  
  my $force = 1;
  if (!-e $data_file || $force) {
    my $cmd = qq^
dbname="gj1_gor_58"
source("${collect_script}");

stats.lnl <- get.vector(con,"SELECT * from $table;")
stats.dups <- get.vector(con,"SELECT * from stats_dups;")

stats.lnl = merge(stats.lnl,stats.dups[,c('data_id','name')],by='data_id')

# Likelihoods key:
# a: (H, G, others)
# b: (H#1, G, others)
# c: (H, G#1, others)
# d: (H#1, G#2, others)
# e: (H#1, G#1, others)
# Which omegas are which for each test?
# pval.1: fg=b_omega_1, bg=b_omega_0 # human
# pval.2: fg=c_omega_1, bg=c_omega_0 # gorilla
# pval.3: fg=d_omega_1, bg=d_omega_0 # gorilla
# pval.4: fg=d_omega_2, bg=d_omega_0 # human
# pval.5: fg=e_omega_1, bg=e_omega_0 # both

# Add the p-values
stats.lnl[,'pval.1'] = with(stats.lnl,1 - pchisq(2*(b_lnL-a_lnL),df=1))
stats.lnl[,'pval.2'] = with(stats.lnl,1 - pchisq(2*(c_lnL-a_lnL),df=1))
stats.lnl[,'pval.3'] = with(stats.lnl,1 - pchisq(2*(d_lnL-b_lnL),df=1))
stats.lnl[,'pval.4'] = with(stats.lnl,1 - pchisq(2*(d_lnL-c_lnL),df=1))
stats.lnl[,'pval.5'] = with(stats.lnl,1 - pchisq(2*(e_lnL-a_lnL),df=1))
stats.lnl[,'pval.6'] = with(stats.lnl,1 - pchisq(2*(d_lnL-e_lnL),df=1))
stats.lnl[,'pval.7'] = with(stats.lnl,1 - pchisq(2*(f_lnL-h_lnL),df=1))
stats.lnl[,'pval.8'] = with(stats.lnl,1 - pchisq(2*(g_lnL-h_lnL),df=1))

method = 'BH'
stats.lnl[,'pval.1.bh'] = with(stats.lnl,p.adjust(pval.1,method=method))
stats.lnl[,'pval.2.bh'] = with(stats.lnl,p.adjust(pval.2,method=method))
stats.lnl[,'pval.3.bh'] = with(stats.lnl,p.adjust(pval.3,method=method))
stats.lnl[,'pval.4.bh'] = with(stats.lnl,p.adjust(pval.4,method=method))
stats.lnl[,'pval.5.bh'] = with(stats.lnl,p.adjust(pval.5,method=method))
stats.lnl[,'pval.6.bh'] = with(stats.lnl,p.adjust(pval.6,method=method))
stats.lnl[,'pval.7.bh'] = with(stats.lnl,p.adjust(pval.7,method=method))
stats.lnl[,'pval.8.bh'] = with(stats.lnl,p.adjust(pval.8,method=method))

# Save the data
save(stats.lnl,file="${data_file}")
^;
#    print "$cmd\n";
    my $params = {};
    Bio::Greg::EslrUtils->run_r($cmd,$params);
  }
}

sub branch_enrichments {
  my $self = shift;
  my $table = shift;
  my $output_filename = shift;
  my $short_prefix = shift;

  my $data_file = $self->get_output_folder . "/$output_filename";

  my $folder = $self->get_output_folder . "/$short_prefix";
  mkpath([$folder]);
  print "FOLDER: $folder\n";

  my $collect_script = Bio::Greg::EslrUtils->baseDirectory."/scripts/collect_sitewise.R";
  my $lrt_script = Bio::Greg::EslrUtils->baseDirectory."/projects/gorilla/lrt_analysis.R";
  my $go_script = Bio::Greg::EslrUtils->baseDirectory."/scripts/go_enrichments.R";
  
  my $force = 1;
  if (!-e $data_file || $force) {
    my $cmd = qq^
dbname="gj1_gor_58"
source("${collect_script}");

stats.lnl <- get.vector(con,"SELECT * from $table;")
stats.dups <- get.vector(con,"SELECT * from stats_dups;")

# Likelihoods key:
# a: (H, G, others)
# b: (H#1, G, others)
# c: (H, G#1, others)
# d: (H#2, G#1, others)
# e: (H#1, G#1, others)
# f: (((H,C),G)#1, others)
# g: (((H,C),G)$1, others)
# h: (((H,C),G), others)
# Which omegas are which for each test?
# pval.1: fg=b_omega_1, bg=b_omega_0 # human
# pval.2: fg=c_omega_1, bg=c_omega_0 # gorilla
# pval.3: fg=d_omega_1, bg=d_omega_0 # gorilla
# pval.4: fg=d_omega_2, bg=d_omega_0 # human
# pval.5: fg=e_omega_1, bg=e_omega_0 # both
# pval.6: fg=either, bg=e_omega_0 # diff. b/t hum&gor
# pval.7: fg=f_omega_1, bg=f_omega_0 # great ape branch
# pval.8: fg=g_omega_1, bg=g_omega_0 # great ape clade

# Add the p-values
stats.lnl[,'pval.1'] = with(stats.lnl,1 - pchisq(2*(b_lnL-a_lnL),df=1))
stats.lnl[,'pval.2'] = with(stats.lnl,1 - pchisq(2*(c_lnL-a_lnL),df=1))
stats.lnl[,'pval.3'] = with(stats.lnl,1 - pchisq(2*(d_lnL-b_lnL),df=1))
stats.lnl[,'pval.4'] = with(stats.lnl,1 - pchisq(2*(d_lnL-c_lnL),df=1))
stats.lnl[,'pval.5'] = with(stats.lnl,1 - pchisq(2*(e_lnL-a_lnL),df=1))
stats.lnl[,'pval.6'] = with(stats.lnl,1 - pchisq(2*(d_lnL-e_lnL),df=1))
stats.lnl[,'pval.7'] = with(stats.lnl,1 - pchisq(2*(f_lnL-h_lnL),df=1))
stats.lnl[,'pval.8'] = with(stats.lnl,1 - pchisq(2*(g_lnL-h_lnL),df=1))

method = 'none'
stats.lnl[,'pval.1.bh'] = with(stats.lnl,p.adjust(pval.1,method=method))
stats.lnl[,'pval.2.bh'] = with(stats.lnl,p.adjust(pval.2,method=method))
stats.lnl[,'pval.3.bh'] = with(stats.lnl,p.adjust(pval.3,method=method))
stats.lnl[,'pval.4.bh'] = with(stats.lnl,p.adjust(pval.4,method=method))
stats.lnl[,'pval.5.bh'] = with(stats.lnl,p.adjust(pval.5,method=method))
stats.lnl[,'pval.6.bh'] = with(stats.lnl,p.adjust(pval.6,method=method))
stats.lnl[,'pval.7.bh'] = with(stats.lnl,p.adjust(pval.7,method=method))
stats.lnl[,'pval.8.bh'] = with(stats.lnl,p.adjust(pval.8,method=method))

# Save the data
save(stats.lnl,file="${data_file}")

source("${lrt_script}",echo=T)

q()
source("${go_script}")

t = 0.05 # pval / FDR threshold.

print(paste("Before: ",nrow(genes)))
print(paste("After: ",nrow(stats.lnl)))

stats.lnl[,'stable_id'] = stats.lnl[,'human_protein']
go.df = go.hs
all.ids = stats.lnl[,'stable_id']

get.df.subset = function(subset) {
  print(paste("size: ",nrow(subset)))
  return(get.enrich.by.subset(
    subset = subset[,'stable_id'],
    all = all.ids,
    go.df = go.df,
    go.field.name = 'stable_id',
    nodeSize = 3
  ))
}
get.df.scores = function(score.field.name,gene.universe) {
  scores = gene.universe[,score.field.name]
  names(scores) = gene.universe[,'stable_id']
  return(get.enrich.by.score(
    named.scores=scores,
    go.df = go.df,
    go.field.name = 'stable_id',
    nodeSize = 3
  ))
}

get.go.data = function() {
  scores = stats.lnl[,'pval.3.bh']
  names(scores) = stats.lnl[,'stable_id']
  geneSelectionFun = function(score){return(score < 0.1)}

  go.vec = strsplit(go.df[,'go'],split=",",fixed=T)
  names(go.vec) = go.df[,'stable_id']

  GOdata <- new("topGOdata",
    ontology = 'BP',
    allGenes = scores,
    annot = annFUN.gene2GO,
    gene2GO = go.vec,
    nodeSize = 3,
    description = '',
    geneSelectionFun = geneSelectionFun
  )
  return(GOdata)
}

GOdata = get.go.data()
save(GOdata,file="${folder}/go_data.Rdata")

# Useful things to do with the GOdata object:

#usedTerms <- usedGO(GOdata)
#term <- usedTerms[1]
#annotated.genes <- genesInTerm(GOdata,term)[[1]]
#sig.genes <- sigGenes(GOdata) # Get the significant genes according to the score threshold.

#q()

`1.up`   = subset(stats.lnl,pval.1.bh < t & b_omega_1 > b_omega_0)
`2.up`     = subset(stats.lnl,pval.2.bh < t & c_omega_1 > c_omega_0)
`3.up`   = subset(stats.lnl,pval.3.bh < t & d_omega_1 > d_omega_0)
`4.up` = subset(stats.lnl,pval.4.bh < t & d_omega_2 > d_omega_0)
`5.up` = subset(stats.lnl,pval.5.bh < t & e_omega_1 > e_omega_0)
`6.up` = subset(stats.lnl,pval.6.bh < t & d_omega_1 > d_omega_2)
`7.up` = subset(stats.lnl,pval.7.bh < t & f_omega_1 > f_omega_0)
`8.up` = subset(stats.lnl,pval.8.bh < t & g_omega_1 > g_omega_0)

`1.down`   = subset(stats.lnl,pval.1.bh < t & b_omega_1 < b_omega_0)
`2.down`     = subset(stats.lnl,pval.2.bh < t & c_omega_1 < c_omega_0)
`3.down` = subset(stats.lnl,pval.3.bh < t & d_omega_1 < d_omega_0)
`4.down`   = subset(stats.lnl,pval.4.bh < t & d_omega_2 < d_omega_0)
`5.down` = subset(stats.lnl,pval.5.bh < t & e_omega_1 < e_omega_0)
`6.down` = subset(stats.lnl,pval.6.bh < t & d_omega_1 < d_omega_2)
`7.down` = subset(stats.lnl,pval.7.bh < t & f_omega_1 < f_omega_0)
`8.down` = subset(stats.lnl,pval.8.bh < t & g_omega_1 < g_omega_0)

# Save the top N genes for each of the tests.
n = 100
all.genes = merge(stats.lnl,stats.dups,by=c('human_gene'))
print(nrow(all.genes))
top.genes = data.frame()
for (i in c(1:8)) {
  up.name = paste(i,'.up',sep="")
  down.name = paste(i,'.down',sep="")

  up = get(up.name)
  down = get(down.name)

  pval.name = paste('pval.',i,'.bh',sep='')
  up[,'lrt'] = up[,pval.name]
  down[,'lrt'] = down[,pval.name]

  up.genes = merge(up,stats.dups)
  down.genes = merge(down,stats.dups)

  up.top <- orderBy(~lrt,data=up.genes)
  down.top <- orderBy(~lrt,data=down.genes)

  up.top[,'label'] <- paste(pval.name,'up')
  top.genes = rbind(top.genes,up.top[1:n,])
  
  down.top[,'label'] <- paste(pval.name,'down')
  top.genes = rbind(top.genes,down.top[1:n,])
}
top.genes = subset(top.genes,!is.na(label))
write.csv(top.genes,file="${folder}/top.genes.csv",row.names=F)

#q()

`tbl.1.up` = get.df.subset(subset=`1.up`)
`tbl.2.up` = get.df.subset(subset=`2.up`)
`tbl.3.up` = get.df.subset(subset=`3.up`)
`tbl.4.up` = get.df.subset(subset=`4.up`)
`tbl.5.up` = get.df.subset(subset=`5.up`)
`tbl.6.up` = get.df.subset(subset=`6.up`)
`tbl.7.up` = get.df.subset(subset=`7.up`)
`tbl.8.up` = get.df.subset(subset=`8.up`)

`tbl.1.down` = get.df.subset(subset=`1.down`)
`tbl.2.down` = get.df.subset(subset=`2.down`)
`tbl.3.down` = get.df.subset(subset=`3.down`)
`tbl.4.down` = get.df.subset(subset=`4.down`)
`tbl.5.down` = get.df.subset(subset=`5.down`)
`tbl.6.down` = get.df.subset(subset=`6.down`)
`tbl.7.down` = get.df.subset(subset=`7.down`)
`tbl.8.down` = get.df.subset(subset=`8.down`)

for (i in c(1:8)) {
  up.name = paste('tbl.',i,'.up',sep="")
  down.name = paste('tbl.',i,'.down',sep="")
  
  up = get(up.name)
  down = get(down.name)

  write.csv(up,file=paste("${folder}/",up.name,".csv",sep=""),row.names=F)
  write.csv(down,file=paste("${folder}/",down.name,".csv",sep=""),row.names=F)
}

gor.up = subset(stats.lnl,d_omega_1 > d_omega_0)
tbl.gor.up = get.df.scores('pval.3.bh',gor.up)
gor.down = subset(stats.lnl,d_omega_1 < d_omega_0)
tbl.gor.down = get.df.scores('pval.3.bh',gor.down)

hum.up = subset(stats.lnl,d_omega_2 > d_omega_0)
tbl.hum.up = get.df.scores('pval.4.bh',hum.up)
hum.down = subset(stats.lnl,d_omega_2 < d_omega_0)
tbl.hum.down = get.df.scores('pval.4.bh',hum.down)

write.csv(tbl.gor.up,file="${folder}/gor.up.ks.csv",row.names=F)
write.csv(tbl.gor.down,file="${folder}/gor.down.ks.csv",row.names=F)
write.csv(tbl.hum.up,file="${folder}/hum.up.ks.csv",row.names=F)
write.csv(tbl.hum.down,file="${folder}/hum.down.ks.csv",row.names=F)

# Using uncorrected p-values here!!
gor.up   = subset(stats.lnl,pval.3 < t & d_omega_1 > d_omega_0)
hum.up = subset(stats.lnl,pval.4 < t & d_omega_2 > d_omega_0)
both.up = merge(hum.up,gor.up)
print(nrow(both.up))
tbl.both.up = get.df.subset(subset=both.up)

stats.lnl[,'pval.x'] = stats.lnl[,'pval.3'] * stats.lnl[,'pval.4']
tbl.both.up.scores = get.df.scores('pval.x',stats.lnl)

write.csv(tbl.both.up,file="${folder}/both.up.csv",row.names=F)
write.csv(tbl.both.up.scores,file="${folder}/both.up.ks.csv",row.names=F)

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
dbname="gj1_gor_58"
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
dbname="gj1_gor_58"
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
dbname="gj1_gor_58"
source("${scripts_dir}/collect_sitewise.R")
source("${scripts_dir}/go_enrichments.R")

top.genes = read.csv("${genes_csv}")
load("${genes_merged}")
merged = merge(top.genes,genes.merged,by=c('node_id','human_gene'))

all.go = merged
print(nrow(all.go))
sub.go = subset(merged,g.muts.corrected > 5 & m_slr_dnds < 0.7)
print(nrow(sub.go))

enriched = get.enrich.df(sub.go[,'human_protein'],all.go[,'human_protein'],go.hs,go.field.name='stable_id')
write.csv(enriched,file="${enrich_file}")
^;
  print "$cmd\n";
  my $params = {};
  Bio::Greg::EslrUtils->run_r($cmd,$params);

}

1;
