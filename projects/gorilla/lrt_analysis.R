### This part of the script gets the data from my local DB and saves it to a file.
## I'm commenting it all out for now.
### Quick example of how to use the data.

# Load the data
if (!exists("stats.lnl")) {
  print("Need to load the stats.lnl dataset!")
  q()
}

#ls()
#str(stats.lnl)

t = 0.05 # Choose a threshold for p-values.

# Look for genes with significantly different w in the H or G terminal branch.
# [ Right now we select for relaxed or increased selection relative to the 
#   background, but this could be made stricter by requiring > 1 ]
sig.human.raw   = subset(stats.lnl,pval.1 < t)
sig.human   = subset(stats.lnl,pval.1.bh < t)
sig.human.up   = subset(stats.lnl,pval.1.bh < t & b_omega_1 > b_omega_0)
sig.human.down = subset(stats.lnl,pval.1.bh < t & b_omega_1 < b_omega_0)
sig.gor.raw   = subset(stats.lnl,pval.2 < t)
sig.gor   = subset(stats.lnl,pval.2.bh < t)
sig.gor.up     = subset(stats.lnl,pval.2.bh < t & c_omega_1 > c_omega_0)
sig.gor.down   = subset(stats.lnl,pval.2.bh < t & c_omega_1 < c_omega_0)

diff.gor = subset(stats.lnl,pval.6.bh < t)
diff.gor.up = subset(stats.lnl,pval.6.bh < t & d_omega_1 > d_omega_2)
diff.gor.down = subset(stats.lnl,pval.6.bh < t & d_omega_1 < d_omega_2)

print.options = function(df,omega_a,omega_b) {
  a = df[,omega_a]
  b = df[,omega_b]

  print(nrow(df))
  print(nrow(df[a > 1,]))
  print(nrow(df[a > b,]))
  print(nrow(df[b > a,]))
}

print.options(sig.human.raw,'b_omega_1','b_omega_0')
print.options(sig.human,'b_omega_1','b_omega_0')
print.options(sig.gor.raw,'c_omega_1','c_omega_0')
print.options(sig.gor,'c_omega_1','c_omega_0')

print.options(diff.gor,'d_omega_1','d_omega_2')

# Look for genes with significantly different w using the more nuanced test
independent.gor.raw = subset(stats.lnl,pval.3 < t)
independent.gor = subset(stats.lnl,pval.3.bh < t)
independent.human.raw = subset(stats.lnl,pval.4 < t)
independent.human = subset(stats.lnl,pval.4.bh < t)

print.options(independent.gor.raw,'d_omega_1','d_omega_0')
print.options(independent.gor,'d_omega_1','d_omega_0')
print.options(independent.human.raw,'d_omega_2','d_omega_0')
print.options(independent.human,'d_omega_2','d_omega_0')

combined.hg.raw = subset(stats.lnl,pval.5 < t)
print.options(combined.hg.raw,'e_omega_1','e_omega_0')
combined.hg = subset(stats.lnl,pval.5.bh < t)
print.options(combined.hg,'e_omega_1','e_omega_0')

# Overlaps.
nrow(merge(combined.hg.raw,sig.gor.raw))
nrow(merge(combined.hg.raw,sig.human.raw))


independent.gor.up   = subset(stats.lnl,pval.3 < t & d_omega_1 > d_omega_0)
independent.human.up = subset(stats.lnl,pval.4 < t & d_omega_2 > d_omega_0)
nrow(independent.human.up)
nrow(independent.gor.up)

library(doBy)
independent.gor.up <- orderBy(~pval.4,data=independent.gor.up)
#independent.gor.up[1:10,]

# Look at overlap between tests 1-3 and 2-4
human.overlap = merge(sig.human.up,independent.human.up)
gor.overlap   = merge(sig.gor.up,  independent.gor.up)
nrow(human.overlap)
nrow(gor.overlap)

# Look at overlap between tests 1-2 and 3-4 (independent acceleration in H+G)
sig.hg.overlap         = merge(sig.human.up,sig.gor.up)
independent.hg.overlap = merge(independent.human.up,independent.gor.up)
nrow(sig.hg.overlap)
nrow(independent.hg.overlap)
nrow(merge(sig.hg.overlap,independent.hg.overlap))

# Look at results from the single test for human and/or gorilla acceleration.
combined.hg.up = subset(stats.lnl,pval.5 < t & e_omega_1 > e_omega_0)
nrow(combined.hg.up)

# ...and so on. Have fun!
