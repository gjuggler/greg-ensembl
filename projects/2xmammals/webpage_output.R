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

  print(paste("Writing page",file.name
))
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

# Load data from the table.
dbname <- 'gj1_2x_57'
source('~/src/greg-ensembl/scripts/collect_sitewise.R')
query <- 'SELECT * FROM web_data WHERE sites_csv IS NOT NULL;'
rows <- get.vector(con,query)


### Data clean-up
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

### /end Data clean-up

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
chrs <- sort(unique(rows[,'Chr']))
for (chr in c('chr22')) {
#for (chr in chrs[1:2]) {
  if (is.na(chr)){next}
  sub.rows <- subset(rows, Chr==chr)
  write.rows.page(chr,sub.rows,css.txt)
}
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

print(bot.links)

###
### HTML template output.
###
html.template <- file("~/src/greg-ensembl/project/2xmammals/webpage_template.html")
html.lines <- readLines(html.template)
close(html.template)
html.string <- paste(html.lines,collapse='\n',sep='')



# Back to original working dir.
setwd(orig.wd)

