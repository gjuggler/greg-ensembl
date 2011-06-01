library(plyr)

get.orang.refseqs <- function() {
  # Read in the Kosiol list, split out the vega/ucsc/refseq IDs.
  tbl <- read.delim("geneset-masked-v2-ponAbe2.txt", sep=" ", stringsAsFactors=F, header=F)
  print(str(tbl))
  id.col <- tbl[,1]
  tokens <- strsplit(id.col, "\\.")
  #print(head(tokens))
  ids <- laply(tokens, function(x){ x[2] })
  #print(head(ids))
  
  # Read in an Ensembl Biomart export with ensembl / etc. IDs
  ens.tbl <- read.delim("ens_60_id_export.txt", sep="\t", stringsAsFactors=F, header=T)
  print(str(ens.tbl))

  ens.g <- ens.tbl[,'Ensembl.Gene.ID']
  vega <- ens.tbl[,'VEGA.transcript.ID.s...OTTT.']
  ref <- ens.tbl[,'RefSeq.DNA.ID']
  ucsc <- ens.tbl[,'UCSC.ID']

  # Remove trailing '.1' from ucsc labels.
  tokens <- strsplit(ucsc, "\\.")
  ucsc <- laply(tokens, function(x){ x[1] })  

  assign('ucsc',ucsc,envir=.GlobalEnv)
  assign('vega',vega,envir=.GlobalEnv)
  assign('ref',ref,envir=.GlobalEnv)

  ens.ids <- c()
  lost.ids <- c()

  ens.ids <- ens.g[ vega %in% ids | ucsc %in% ids | ref %in% ids ]
  lost.ids <- ids[ !( ids %in% vega | ids %in% ucsc | ids %in% ref ) ]

  assign('ids',ids,envir=.GlobalEnv)
  assign('ens.ids',ens.ids,envir=.GlobalEnv)
  assign('lost.ids',lost.ids,envir=.GlobalEnv)
  write.table(ens.ids, file="orang_ens_ids.txt")
  write.table(lost.ids, file="orang_lost_ids.txt")

  return()

  for (id in ids) {
    print(id)
    found.vega <- grep(id, vega)
    found.ucsc <- grep(id, ucsc)
    found.ref <- grep(id, ref)

    if (length(found.vega) > 0) {
      ens.ids <- c(ens.ids, ens.g[found.vega[1]])
    } else if (length(found.ucsc) > 0) {
      ens.ids <- c(ens.ids, ens.g[found.ucsc[1]])
    } else if (length(found.ref) > 0) {
      ens.ids <- c(ens.ids, ens.g[found.ref[1]])
    } else {
      lost.ids <- c(lost.ids, id)
    }
  }
  
}