# Weighted TPM
# Depends on tpmwR.cpp; first compile as 'R CMD SHLIB tpmwR.cpp'
# p: array of p-values
# w: weights
# tau: truncation parameter [0-1] (1 for Fisher, i.e. no truncation)
# loops: number of Monte-Carlo loops
# seed: random seed 
Tpmw <- function(p, w, tau, loops, seed=0) {
    n <- length(p)
    if(seed == 0) seed <- as.integer(runif(1,0,1e9))
    ppv <- as.double(1)
    dyn.load("tpmwR.so")
    tpm <- .C("tpmwR", as.double(p), as.double(w), as.integer(n), 
	as.double(tau), as.integer(loops), as.double(seed), pv=ppv)
    ppv <- as.double(tpm$pv)
    dyn.unload("tpmwR.so")
    return(ppv)
}
