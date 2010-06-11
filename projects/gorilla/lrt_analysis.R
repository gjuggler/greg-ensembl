### This part of the script gets the data from my local DB and saves it to a file.
## I'm commenting it all out for now.
#dbname="gj1_gor_57"
#source("../../scripts/collect_sitewise.R");

#stats.lnl <- get.vector(con,"SELECT * from stats_lnl;")

## Likelihoods key:
## a: (H, G, others)
## b: (H#1, G, others)
## c: (H, G#1, others)
## d: (H#1, G#2, others)
## e: (H#1, G#1, others)
## Which omegas are which for each test?
## pval.1: fg=b_omega_1, bg=b_omega_0 # human
## pval.2: fg=c_omega_1, bg=c_omega_0 # gorilla
## pval.3: fg=d_omega_2, bg=d_omega_0 [and d_omega_1] # human
## pval.4: fg=d_omega_1, bg=d_omega_0 [and d_omega_2] # gorilla
## pval.5: fg=e_omega_1, bg=e_omega_0 # both

## Add the p-values
#stats.lnl[,'pval.1'] = with(stats.lnl,1 - pchisq(2*(b_lnL-a_lnL),df=1))
#stats.lnl[,'pval.2'] = with(stats.lnl,1 - pchisq(2*(c_lnL-a_lnL),df=1))
#stats.lnl[,'pval.3'] = with(stats.lnl,1 - pchisq(2*(d_lnL-b_lnL),df=1))
#stats.lnl[,'pval.4'] = with(stats.lnl,1 - pchisq(2*(d_lnL-c_lnL),df=1))
#stats.lnl[,'pval.5'] = with(stats.lnl,1 - pchisq(2*(e_lnL-a_lnL),df=1))

## Save the data
#save(stats.lnl,file="stats_lnl.Rdata")

### Quick example of how to use the data.

# Load the data
load("stats_lnl.Rdata")

ls()
str(stats.lnl)

t = 0.05 # Choose a threshold for p-values.

# Look for genes with significantly different w in the H or G terminal branch.
# [ Right now we select for relaxed or increased selection relative to the 
#   background, but this could be made stricter by requiring > 1 ]
sig.human.up   = subset(stats.lnl,pval.1 < t & b_omega_1 > b_omega_0)
sig.human.down = subset(stats.lnl,pval.1 < t & b_omega_1 < b_omega_0)
sig.gor.up     = subset(stats.lnl,pval.1 < t & c_omega_1 > c_omega_0)
sig.gor.down   = subset(stats.lnl,pval.1 < t & c_omega_1 < c_omega_0)

nrow(sig.human.up)
nrow(sig.human.down)
nrow(sig.gor.up)
nrow(sig.gor.down)

# Look for genes with significantly different w using the more nuanced test
independent.human.up = subset(stats.lnl,pval.3 < t & d_omega_2 > d_omega_0)
independent.gor.up   = subset(stats.lnl,pval.4 < t & d_omega_1 > d_omega_0)

nrow(independent.human.up)
nrow(independent.gor.up)

# Look at overlap between tests 1-3 and 2-4
human.overlap = merge(sig.human.up,independent.human.up)
gor.overlap   = merge(sig.gor.up,  independent.gor.up)

nrow(human.overlap)
nrow(gor.overlap)
# [it's strange how gorilla has much stronger overlap between the two tests,
#  even though human has many more results from the 'independent' method]

# Look at overlap between tests 1-2 and 3-4 (independent acceleration in H+G)
sig.hg.overlap         = merge(sig.human.up,sig.gor.up)
independent.hg.overlap = merge(independent.human.up,independent.gor.up)

nrow(sig.hg.overlap)
nrow(independent.hg.overlap)

# ...and so on. Have fun!
