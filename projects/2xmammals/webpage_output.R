library(hwriter)

###
# The main HTML writing function for output.
#
tbl.count <- 0
write.rows <- function(page, rows) {
  # Create linked image thumbnails for any image columns
  img.fields <- grep('(plot)',colnames(rows),value=T)
  for (field in img.fields) {
    if(length(rows[,field]) > 0) {
      rows[,field] <- hwriteImage(rows[,field],table=F,width=300)
    }
  }

  # Store all file fields in the col.links list
  file.fields <- grep('(file)',colnames(rows),value=T)
  #print(file.fields)
  col.links <- list()
  for (field in file.fields) {
    col.links[[field]] <- rows[,field]
    rows[,field] <- field
  }
  col.links[['Alignment plot']] <- rows[,'aln_pdf']
  col.links[['SLR plot']] <- rows[,'sites_pdf']

  # Remove unnecessary fields from display.
  remove.fields <- c('aln_pdf','sites_pdf','node_id','parameter_set_id','data_id')
  for (field in remove.fields) {
    rows[,field] <- NULL
  }

  # Links to GeneCards and Ensembl
  rows[,'Gene'] <- hmakeTag('a',rows[,'Gene'],href=paste('http://www.genecards.org/cgi-bin/carddisp.pl?gene=',rows[,'Gene'],sep=''),target='_blank')
  rows[,'Ensembl ID'] <- hmakeTag('a',rows[,'Ensembl ID'],href=paste('http://www.ensembl.org/Homo_sapiens/protview?peptide=',rows[,'Ensembl ID'],sep=''),target='_blank')

  hwrite(rows,page,
    br=TRUE,
    table.class='width sortable',
    table.id=paste('table_',tbl.count,sep=''),
    col.link=col.links
  )
  tbl.count <- tbl.count + 1
}

write.rows.page <- function(page_base,rows,css.txt) {
  # Open the HTML file.
  file.name <- paste(page_base,'.html',sep='')
  file.name <- gsub('[- \\(\\)/\\\\]','_',file.name)

  print(paste("Writing page",file.name))
  page <- openPage(file.name,
    title=paste(page_base),
    link.javascript=c('sorttable.js','tablefilter.js'),
    css=css.txt
  )
  # Write rows.
  write.rows(page,rows)
  closePage(page)
  return(file.name)
}

bed.output <- function(rows,color.field,color.lo,color.hi) {
  # BED wants '+' or '-' strand.
  strand <- rows$chr_strand
  strand[strand == 1] <- '+'
  strand[strand == -1] <- '-'
  rows$chr_strand <- strand

  rows <- subset(rows,!is.na(Chr))
  rows[rows$Chr == 'chrMT',]$Chr <- 'chrM'

  color.values <- rows[,color.field]
  if (color.field == 'Pos-sel score') {
    color.values[color.values > 30] <- 30
  } else if (color.field == 'Overall dN/dS') {
    color.values[color.values > 1] <- 1
  } else if (color.field == 'Tree length') {
    color.values[color.values > 0.7] <- 0.7
  }
  color.values.no.na <- color.values[!is.na(color.values)] # Remove NAs.
  # Scale color values from 0 to 1.
  color.values[is.na(color.values)] <- 0
  color.range <- (max(color.values,na.rm=T)-min(color.values,na.rm=T))
  color.values <- (color.values - min(color.values,na.rm=T)) / color.range
  #print(summary(color.values))
  color.ramp <- colorRamp(c(color.lo,color.hi))
  colors <- color.ramp(color.values)
  colors.string <- apply(colors,1,function(x) {
    x <- as.integer(x)
    z <- paste(x[1],x[2],x[3],sep=',')
    return(z)
  })
  
  rows$color.tmp <- colors.string
  rows$score.tmp <- 0
  bed <- rows[,c('Chr','chr_start','chr_end','Gene','score.tmp','chr_strand','chr_start','chr_end','color.tmp','Ensembl ID',color.field)]

  bed[is.na(bed[,color.field]),color.field] <- 0
  
  # Generate the bigBed file.
  file.name <- paste(color.field,'.bb',sep='')
  file.name <- gsub('[- \\(\\)/\\\\]','_',file.name)
  full.file <- paste(folder,"/",file.name,sep='')
  if (!file.exists("chrom.sizes")) {
    system("fetchChromSizes hg19 > chrom.sizes")
  }
  print(paste("Outputting BED",file.name))
  write.table(bed,file="genes_tmp.bed",sep="\t",quote=F,row.names=F,col.names=F)  
  system("sort -k1,1 -k2,2n genes_tmp.bed > genes_tmp_sorted.bed")
  system(paste("bedToBigBed -unc -as=webpage_bed_format.as -bedFields=9 genes_tmp_sorted.bed chrom.sizes",full.file))

#  unlink("genes_tmp.bed")
#  unlink("genes_tmp_sorted.bed")

  rows$color.tmp <- NULL
  rows$score.tmp <- NULL
  return(file.name)
}

bigwig.output <- function(filename) {
  if (!exists("sites")) { 
    load("sites_1.Rdata")
    sites[,'signed_lrt'] = sites$lrt_stat
    sites[sites$omega < 1,]$signed_lrt = -sites[sites$omega < 1,]$signed_lrt
    sites <- subset(sites,!is.na(chr_name))
    sites[sites$chr_name == 'chrMT',]$chr_name <- 'chrM'
    assign("sites",sites,envir=.GlobalEnv)
  }

  sites$chr_end <- sites$chr_start+1
  unique.sites <- sites[!duplicated(sites[,c('chr_name','chr_start')]),]

  bedGraph <- unique.sites[,c('chr_name','chr_start','chr_end','signed_lrt')]
  bedGraph$chr_start <- sprintf("%d",bedGraph$chr_start)
  bedGraph$chr_end <- sprintf("%d",bedGraph$chr_end)
  write.table(bedGraph,file="sites.bedGraph",sep="\t",quote=F,row.names=F,col.names=F)  

  ### Generate the sites bigWig file.
  sites.filename <- "sitewise_dnds.bw"
  sites.full.file <- paste(folder,"/",sites.filename,sep='')
  print(paste("Outputting BigWig",sites.filename))
  system("sort -k1,1 -k2,2n sites.bedGraph > sites_sorted.bedGraph")
  system(paste("bedGraphToBigWig -unc sites_sorted.bedGraph chrom.sizes",sites.full.file))

#  unlink("sites.bedGraph")
#  unlink("sites_sorted.bedGraph")

  return(sites.filename)
}

# Load data from the table.
dbname <- 'gj1_2x_57'
source('~/src/greg-ensembl/scripts/collect_sitewise.R')
query <- 'SELECT * FROM web_data WHERE sites_csv IS NOT NULL;'
rows <- get.vector(con,query)


### Data clean-up
###
# Take neg-log of fisher p-values (easier sorting in web interface)
rows$pos_sel_score <- as.numeric(sprintf("%.3f",-log(rows$pval_fisher)))
rows$pval_fisher <- NULL
# Clean up gene names (remove dash at end)
rows$gene_name <- sub("-.*$","",rows$gene_name)
# Shorten some numbers.
for (field in c('gc_content_mean','slr_dnds','slr_kappa','seq_length_mean')) {
  rows[,field] <- as.numeric(sprintf("%.2f",rows[,field]))
}
for (field in c('f_neg','f_neut','f_pos')) {
  rows[,field] <- as.numeric(sprintf("%.4f",rows[,field]))
}
# Rename some columns.
rename.col <- function(rows,old,new) {
  rows[,new] <- rows[,old]
  rows[,old] <- NULL
  return(rows)
}
new.names <- list(
  c('stable_id','Ensembl ID'),
  c('gene_name','Gene'),
  c('chr_name','Chr'),
  c('aln_png','Alignment plot'),
  c('sites_png','SLR plot'),
  
  c('f_neg','Fraction negative'),
  c('f_pos','Fraction positive'),
  c('f_neut','Fraction neutral'),

  c('leaf_count','Sequence count'),
  c('tree_mean_path','Tree length'),
  c('slr_dnds','Overall dN/dS'),
  c('slr_kappa','Kappa'),
  c('gc_content_mean','GC'),
  c('seq_length_mean','Seq length'),
  c('pos_sel_score','Pos-sel score'),

  c('tree_file','Tree file'),
  c('aln_file','Alignment file'),
  c('data_file','Data file'),
  c('params_file','Debug file'),
  c('sites_csv','SLR file')
)
for (names in new.names) {
  rows <- rename.col(rows,names[1],names[2])
}
###
### /end Data clean-up


### BigBed / BigWig output
###
gene.links <- c()
color.fields <- c('Overall dN/dS','GC','Pos-sel score','Tree length')
color.lo <- c('black','black','black','black')
color.hi <- c('blue','red','green','orange')
for (i in 1:length(color.fields)) {
  cur.field <- color.fields[i]
  file.name <- bed.output(rows,cur.field,color.lo[i],color.hi[i])

  # Produce a link to the UCSC browser loading this track.
  hgct.custom.text <- paste('track',
                            'type=bigBed',
                            paste('name=',cur.field,sep=''),
                            'visibility=squish',
                            'itemRgb=on',
                            paste('bigDataUrl=','http://www.ebi.ac.uk/~greg/mammals/',file.name,sep=''),
                            sep=' ')
  hgct.custom.encoded <- URLencode(hgct.custom.text)
  hgct.url <- sprintf("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&position=chr22&hgct_customText=%s",hgct.custom.encoded)

  gene.links <- c(gene.links,hmakeTag('a', cur.field, href=hgct.url)  )
}

bigwig.file <- bigwig.output(rows)
hgct.custom.text <- paste('track',
                          'type=bigWig',
                          'name=Sitewise dnds',
                          'visibility=full',
                          paste('bigDataUrl=','http://www.ebi.ac.uk/~greg/mammals/',bigwig.file,sep=''),
                          sep=' ')
hgct.custom.encoded <- URLencode(hgct.custom.text)
hgct.url <- sprintf("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&position=chr22&hgct_customText=%s",hgct.custom.encoded)
sitewise.links <- hmakeTag('a','Sitewise dN/dS', href=hgct.url)
###
### /end UCSC track output


# Constants.
folder <- '~/scratch/mammals/output/2010-04-30/2010-04-30_01/web'
file <- 'tables.html'

# Chdir into the output folder.
orig.wd <- getwd()
#print(orig.wd)
setwd(folder)

# Sort and filter javascripts.
js <- 'http://www.kryogenix.org/code/browser/sorttable/sorttable.js'
system(paste('wget -nc',js))
js2 <- 'http://www.javascriptkit.com/script/script2/tablefilter.js'
system(paste('wget -nc',js2))

# CSS text.
css.txt <- "
body {
  font-size: 9px;
}
table {
  border-collapse: collapse;
  background-color: white;
  border-spacing: 2px;
}
.odd {
#  background-color: #def;
}
"

### Pages by chromosome.
chrs <- unique(rows[,'Chr'])
chrs <- chrs[!is.na(chrs)]
chrs <- sub('chr','',chrs)
numeric.chrs <- sort(as.numeric(chrs[!is.na(as.numeric(chrs))]))
chrs <- c(numeric.chrs,chrs[is.na(as.numeric(chrs))])
chrs <- paste('chr',chrs,sep='')
for (chr in c('chr22')) {
#for (chr in chrs) {
  if (is.na(chr)){next}
  sub.rows <- subset(rows, Chr==chr)
  write.rows.page(chr,sub.rows,css.txt)
}
print(chrs)
chr.links <- hmakeTag('a', chrs, href=paste('./',chrs,'.html',sep=''))

### Top / bottom 100 by field.
all.fields <- colnames(rows)
numeric.fields <- c()
top.files <- c()
bot.files <- c()
for (field in c('Fraction neutral')) {
#for (field in all.fields) {
  if (!is.numeric(rows[,field])){next}

  numeric.fields <- c(numeric.fields,field)

  x <- rows[order(-rows[,field]),]
  sub.rows <- head(x,n=200)
  f <- write.rows.page(paste(field,'top'),sub.rows,css.txt)
  top.files <- c(top.files,f)

  x <- rows[order(rows[,field]),]
  sub.rows <- head(x,n=200)
  f <- write.rows.page(paste(field,'bottom'),sub.rows,css.txt)
  bot.files <- c(bot.files,f)
}
top.links <- hmakeTag('a', numeric.fields, href=paste('./',top.files,sep=''))
bot.links <- hmakeTag('a', numeric.fields, href=paste('./',bot.files,sep=''))

###
### HTML template output.
###
html.template <- file("~/src/greg-ensembl/projects/2xmammals/webpage_template.html")
html.lines <- readLines(html.template)
close(html.template)
html.string <- paste(html.lines,collapse='\n',sep='')

html.string <- sub('CHROMOSOME_LINKS',paste(chr.links,collapse=' '),html.string)
html.string <- sub('TOP_LINKS',paste(top.links,collapse=' '),html.string)
html.string <- sub('BOTTOM_LINKS',paste(bot.links,collapse=' '),html.string)
html.string <- sub('SITEWISE_LINKS',paste(sitewise.links,collapse=' '),html.string)
html.string <- sub('GENE_LINKS',paste(gene.links,collapse=' '),html.string)

writeChar(html.string,paste(folder,"/webpage_output.html",sep=''))

# Back to original working dir.
setwd(orig.wd)

