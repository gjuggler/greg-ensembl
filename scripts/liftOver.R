lift.over <- function(df, chainName) {
  temp.in <- tempfile()
  temp.out <- tempfile()
  temp.err <- tempfile()

  chain.file <- paste(chainName, ".over.chain", sep='')
  if (!file.exists(chain.file)) {
    chain.url <- paste(
      'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/', 
      chainName, 
      '.over.chain.gz',
      sep=''
    )
    system(paste('wget -nc', chain.url))
    system(paste('gunzip ', chainName, '.over.chain.gz', sep=''))
  }
  
  # Store an ID field to merge the results by.
  df$tmp_id <- 1:nrow(df)
  
  # Only keep non-NA chr_name rows.
  df <- subset(df, !is.na(chr_name) & !is.na(chr_start) & !is.na(chr_end))

  # Output the data frame as chr, start, end tuples
  sub.df <- subset(df, select=c('chr_name', 'chr_start', 'chr_end', 'tmp_id'))
  sub.df <- format(sub.df, scientific=F)
  write.table(sub.df, file=temp.in, row.names=F, quote=F, col.names=F)

  exec <- '/software/ensembl/compara/bin/liftOver'
  cmd <- paste(exec, temp.in, chain.file, temp.out, temp.err)
  system(cmd)
  
  in.df <- read.table(temp.out, header=F)
  names(in.df) <- c('chr_name', 'chr_start', 'chr_end', 'tmp_id')

  merged <- merge(in.df, df, by=c('tmp_id'), suffixes=c('', '.old'))
  merged$tmp_id <- NULL

  file.remove(temp.in, temp.out, temp.err)
  return(merged)
}