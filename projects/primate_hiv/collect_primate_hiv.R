dbname <- 'gj1_hiv_58'
source("~/src/greg-ensembl/scripts/collect_sitewise.R",echo=T)

all.data <- data.frame()
table.name <- 'stats_windows'

cols <- 'parameter_set_id,aln_type,gene_name,n_leaves,n_pos_sites,n_sites,parameter_set_shortname,peptide_window_start,peptide_window_end,peptide_window_width,pval,stable_id_transcript'
query <- sprintf("SELECT * FROM %s",table.name)
cur.df <- get.vector(con,query,columns='all')
all.data <- rbind(all.data,cur.df)

print(nrow(all.data))

psets <- unique(all.data$parameter_set_shortname)
aln.types <- unique(all.data$aln_type)
print(psets)
print(aln.types)

if (length(psets) > 1 || length(aln.types) > 1) {
  for (pset in psets) {
    print(pset)
    for (aln.type in aln.types) {
      print(aln.type)
      data <- subset(all.data,parameter_set_shortname==pset & aln_type==aln.type)
  
      file.name <- paste(output_folder,"/windows_",pset,"_",aln.type,sep="")
      print(file.name)
      save(data,file=paste(file.name,".Rdata",sep=""))
      #write.csv(data,file=paste(file.name,".csv",sep=""),row.names=F)
    }
  }
}

data <- all.data
save(data,file=paste(output_folder,"/windows_all.Rdata",sep=""))