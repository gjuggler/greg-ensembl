source("slr.tools.R")

node_id=21
base = sprintf("http://www.ebi.ac.uk/~greg/eslr/root_nodes/%s",node_id)
file = sprintf("%s/%s",base,"origSlr.txt")
plot.slr(file)