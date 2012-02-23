uname  <- Sys.getenv("USER")
lib <- ifelse(uname == 'gj1', 'src', 'lib')
source(sprintf("~/%s/greg-ensembl/projects/subclades/subclades_sim.R", lib))

args <- commandArgs(trailingOnly=T)
main <- function() {
  if (!is.na(args[1])) {
    fn.name <- args[1]
    print(paste("Calling", fn.name))

    jobindex <- Sys.getenv('LSB_JOBINDEX')
    if (jobindex != '' && as.integer(jobindex) > 0) {
      args <- c(args, as.integer(jobindex))
    }
    print(args)

    error.f <- function(e) {
      print("###ERROR###")
      print(as.character(e))
    }

    tryCatch(
      do.call(fn.name, as.list(args[2:length(args)])),
      error = error.f
    )
  }
}

ifebi <- function(a, b) {
  ifelse(Sys.getenv("USER") == 'greg', return(a), return(b))
}

bsub.function <- function(fn_name, queue='normal', mem=3, extra.args='', jobarray=NULL, jobarray_id=fn_name, drop.output=F) {
  array_s <- ''
  if (!is.null(jobarray)) {
    array_s <- paste('-J ', jobarray_id, '[1-', jobarray, ']', sep='')
  }
  output_s <- ifebi(
    sprintf('-o "/homes/greg/scratch/lsf_logs/%s_%%J_%%I.txt"', fn_name),
    sprintf('-o "/nfs/users/nfs_g/gj1/scratch/lsf_logs/%s_%%J_%%I.txt"', fn_name)
  )
  if (drop.output) {
    output_s <- '-o /dev/null'
  }

  queue <- ifebi(
    'research-rh6',
    queue
  )

  args_s <- paste(fn_name, ' ', extra.args, sep='')
  cmd <- ifebi(
    sprintf('bsub -q %s -R "select[mem>%d000] rusage[mem=%d000]" -M%d000 %s %s "R-2.14.0 --vanilla --args %s < lsf_job.R"',
      queue, mem, mem, mem, 
      array_s, output_s, args_s
    ),
    sprintf('bsub -q %s -R "select[mem>%d000] rusage[mem=%d000]" -M%d000000 %s %s "/software/R-2.13.0/bin/R --vanilla --args %s < lsf_job.R"',
      queue, mem, mem, mem,
      array_s, output_s, args_s
    )
  )
  print(cmd)
  system(cmd)
}

main()
