lift.over <- function(df, chainName, start.s='chr_start', end.s='chr_end', chr.s='chr_name') {
  temp.in <- tempfile()
  temp.out <- tempfile()
  temp.err <- tempfile()

  chain.file <- paste(chainName, ".over.chain", sep='')
  if (!file.exists(chain.file)) {
    if (chainName == 'hg18ToHg19') {
      base.url <- 'http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/'  
    } else {
      base.url <- 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/'
    }
    chain.url <- sprintf("%s%s.over.chain.gz", base.url, chainName)
    system(paste('wget -nc', chain.url))
    system(sprintf("gunzip %s.over.chain.gz", chainName))
  }
  
  # Store an ID field to merge the results by.
  df$tmp_id <- 1:nrow(df)
  
  # Only keep non-NA chr_name rows.
  starts <- df[, start.s]
  ends <- df[, end.s]
  chrs <- df[, chr.s]
  df <- df[!is.na(starts) & !is.na(ends) & !is.na(chrs),]
  rm(starts, ends, chrs) 

  # Output the data frame as chr, start, end tuples
  sub.df <- subset(df, select=c(chr.s, start.s, end.s, 'tmp_id'))
  sub.df <- format(sub.df, scientific=F)
  write.table(sub.df, file=temp.in, row.names=F, quote=F, col.names=F)

  exec <- 'liftOver'
  cmd <- paste(exec, temp.in, chain.file, temp.out, temp.err)
  system(cmd)
  
  in.df <- read.table(temp.out, header=F)
  names(in.df) <- c(chr.s, start.s, end.s, 'tmp_id')

  merged <- merge(in.df, df, by=c('tmp_id'), suffixes=c('', '.old'))
  merged$tmp_id <- NULL

  file.remove(temp.in, temp.out, temp.err)
  return(merged)
}