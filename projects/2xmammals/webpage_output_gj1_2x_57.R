library(hwriter)
library(plyr)

###
# Some important constants used throughout the script.
#
parameter_set_shortname <- '57'
folder <- paste('~/scratch/gj1_2x_57/2011-01-21_01/web_data_',parameter_set_shortname,sep='')

setwd(folder)

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
  col.links <- list()
  for (field in file.fields) {
    col.links[[field]] <- c(NA,rows[,field])
    rows[,field] <- field
  }
  col.links[['Alignment plot']] <- c(NA,rows[,'aln_pdf'])

  # Links to GeneCards, Ensembl, UCSC
  rows[,'Gene'] <- hmakeTag('a',rows[,'Gene'],name=rows[,'Ensembl ID'],href=paste('http://www.genecards.org/cgi-bin/carddisp.pl?gene=',rows[,'Gene'],sep=''),target='_blank')
  peptide_ids <- hmakeTag('a', rows[,'Ensembl ID'], href=paste('http://www.ensembl.org/Homo_sapiens/protview?peptide=', rows[, 'Ensembl ID'], sep=''), target='_blank')
  transcript_ids <- hmakeTag('a', rows[,'Ensembl transcript ID'], href=paste('http://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=', rows[, 'Ensembl transcript ID'], sep=''), target='_blank')
  gene_ids <- hmakeTag('a', rows[,'Ensembl gene ID'], href=paste('http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=', rows[, 'Ensembl gene ID'], sep=''), target='_blank')
  rows[,'Ensembl ID'] <- paste(peptide_ids, transcript_ids, gene_ids, sep="<br/>")
  rows[,'UCSC'] <- hmakeTag('a','browse',href=rows[,'UCSC'],target='_blank')

  # Remove unnecessary fields from display.
  remove.fields <- c(
    'aln_pdf',
    'sites_pdf',
    'sites_png',
    'node_id',
    'parameter_set_id',
    'data_id',
    'Data file',
    'Debug file',
    'gene_taxon_id',
    'job_id',
    'parameter_set_shortname',
    'Ensembl gene ID',
    'Ensembl transcript ID',
    'Chr',
    'Start',
    'End',
    'Strand',
    'Sequence count',
    'Number negative',
    'Number neutral',
    'GC CDS',
    'GC genomic',
    'pfam_domains'
  )
  for (field in remove.fields) {
    rows[,field] <- NULL
  }

  # Format output for some numeric rows.
  for (field in c('GC3', 'Overall dN/dS', 'Kappa')) {
    rows[,field] <- sprintf("%.2f",as.numeric(rows[,field]))
  }
  for (field in c('Fraction negative', 'Fraction neutral', 'Fraction positive')) {
    rows[,field] <- sprintf("%.3f",as.numeric(rows[,field]))
  }

  # Add help links for the column names.
  header <- paste(colnames(rows)," <a class='help' target='_blank' href='./index.html#",colnames(rows),"'>?</a>",sep='')
  rows <- rbind(header,rows) # Note: rows now contains the column names as the first row!

  hwrite("<em>Click row headers to sort. Click '?' for more information on columns.</em><br/>",page)
  hwrite(rows,page,
    br=TRUE,
    col.names =FALSE,
    row.names =FALSE,
    row.bgcolor=list('1'='#ddffdd'),
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
    link.css='webpage_output.css',
    css=css.txt
  )
  # Write rows.
  write.rows(page,rows)
  closePage(page)
  
  full.file <- paste(folder,'/',file.name,sep='')
  file.copy(file.name,full.file)
  return(file.name)
}

bed.output <- function(rows,color.field,color.lo,color.hi) {
  # BED wants '+' or '-' strand.
  strand <- rows$Strand
  strand[strand == 1] <- '+'
  strand[strand == -1] <- '-'
  rows$Strand <- strand

  rows <- subset(rows,!is.na(Chr))
  if (nrow(subset(rows, Chr == 'chrMT')) > 0) {
    rows[rows$Chr == 'chrMT',]$Chr <- 'chrM'
  }

  color.values <- rows[,color.field]
  if (color.field == 'Pos-sel score') {
    color.values[color.values > 10] <- 10
  } else if (color.field == 'Overall dN/dS') {
    color.values[color.values > 0.5] <- 0.5
  } else if (color.field == 'Tree length') {
    color.values[color.values > 0.5] <- 0.5
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
  bed <- rows[,c('Chr','Start','End','Gene','score.tmp','Strand','Start','End','color.tmp','Ensembl ID',color.field)]

  bed[is.na(bed[,color.field]),color.field] <- 0

  print(str(bed))
  
  # Generate the bigBed file.
  bed.name <- paste(color.field,'.bed',sep='')
  bed.name <- gsub('[- \\(\\)/\\\\]','_',bed.name)
  bb.name <- paste(color.field,'.bb',sep='')
  bb.name <- gsub('[- \\(\\)/\\\\]','_',bb.name)
  full.bed <- paste(folder,"/",bed.name,sep='')
  full.bb <- paste(folder,"/",bb.name,sep='')

  # Bail out if the file exists.
  if (file.exists(full.bb)) {
    print(paste("Skipping BED file [",full.bb,"] -- delete it if you want to regenerate the data."))
    return(full.bb)
  }
  if (!file.exists("chrom.sizes")) {
    system("fetchChromSizes hg18 > chrom.sizes")
  }
  bed.format.file <- paste(folder, '/', 'webpage_bed_format.as', sep='')
  if (!file.exists(bed.format.file)) {
    file.copy(from="~/src/greg-ensembl/projects/2xmammals/webpage_bed_format.as", to=bed.format.file)
  }
  
  #print(unique(rows[,'Chr']))

  print(paste("Outputting BED",bed.name))
  write.table(bed,file=full.bed,sep="\t",quote=F,row.names=F,col.names=F)  
  write.table(bed,file="genes_tmp.bed",sep="\t",quote=F,row.names=F,col.names=F)  
  system("sort -k1,1 -k2,2n genes_tmp.bed > genes_tmp_sorted.bed")
  system(paste("bedToBigBed -unc -as=webpage_bed_format.as -bedFields=9 genes_tmp_sorted.bed chrom.sizes",full.bb))

  file.remove("genes_tmp.bed", "genes_tmp_sorted.bed")

  rows$color.tmp <- NULL
  rows$score.tmp <- NULL
  return(file.name)
}


bigwig.output <- function() {
  if (!exists("sites")) { 
    print("  loading sites file")
    sites <- sites_m
#    sites.file <- paste("sites_", parameter_set_shortname, ".Rdata", sep='')
#    if (!file.exists(sites.file)) {
#      file.copy(paste("../", sites.file, sep=""), sites.file)
#    }
#    load(sites.file)
    print("  fixing up sites data")
    sites[,'signed_lrt'] = sites$lrt_stat
    sites[sites$omega < 1,]$signed_lrt = -sites[sites$omega < 1,]$signed_lrt

    sites <- subset(sites,!is.na(chr_name))

    if (nrow(subset(sites, chr_name=='chrM')) > 0) {
      # Change 'MT' to 'M', and don't forget to add a new factor level.
      levels(sites$chr_name) <- c(levels(sites$chr_name),'chrM')
      sites[sites$chr_name == 'chrMT',]$chr_name <- 'chrM'
    }

    assign("sites",sites,envir=.GlobalEnv)
  }

  output.filename <- 'sites.bw'
  full.path <- paste(folder,'/',output.filename,sep='')
  if (!file.exists(full.path)) {
    print(paste("Outputting BigWig to [",full.path,"]..."))

    if (!file.exists('sites_sorted.bedGraph')) {
      sites$chr_end <- sites$chr_start+3

      print("  sorting and un-overlapping")
      sites <- ddply(sites,.(chr_name),function(x) {
        x <- x[order(x[,'chr_start']),] # Sort sites by start pos.
        start.n <- nrow(x)
        x$chr_end_max <- pmax(x$chr_end)
        # Remove duplicate starting positions.
        x <- x[!duplicated(x[,'chr_start']),]
        # Remove sites starting before the previous end point.
        x <- x[x$chr_start >= c(0,x[-nrow(x),'chr_end_max']),]
        print(paste("  ",as.character(x[1,'chr_name']),start.n,nrow(x)))
        return(x)
      })

      print("  writing bedgraph")
      bedGraph <- sites[,c('chr_name','chr_start','chr_end','signed_lrt')]
      bedGraph$chr_start <- sprintf("%d",bedGraph$chr_start)
      bedGraph$chr_end <- sprintf("%d",bedGraph$chr_end)
      write.table(bedGraph,file="sites.bedGraph",sep="\t",quote=F,row.names=F,col.names=F)  

      ### Generate the sites bigWig file.
      print("  sorting")
      system("sort -k1,1 -k2,2n sites.bedGraph > sites_sorted.bedGraph")
    }
    print("  bg2bw")
    system(paste("bedGraphToBigWig -unc sites_sorted.bedGraph chrom.sizes",output.filename))
    file.copy(output.filename,full.path)
  } else {
    print(paste("Skipping BigWig file [",full.path,"] -- delete it if you want to re-generate."))
  }

  return(output.filename)
#  unlink("sites.bedGraph")
#  unlink("sites_sorted.bedGraph")
}

# Load data from the table.
dbname <- 'gj1_2x_57'

parameter_set_query_string <- ''

source('~/src/greg-ensembl/scripts/collect_sitewise.R')
query <- paste(
  'SELECT * FROM web_data', parameter_set_query_string, 
  ' WHERE sitewise_value_count > 0 ', 
#  ' AND gene_taxon_id=9606;',
  sep=''
)
rows <- get.vector(con,query)

# Just for 57.
rows <- genes_m
rows$protein_id <- rows$stable_id

fill.in.blank.fields <- c(
  'pfam_domains',
  'gene_description'
)

for (field in fill.in.blank.fields) {
  if (!field %in% colnames(rows)) {
    #print(paste("Filling in field",field,"with blanks!"))
    rows[, field] <- ''
  }
}


### Load Pfam annotations.
###
if (!file.exists("pfamA.txt")) {
  #pf.txt <- 'ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/database_files/pfamA.txt.gz'
  pf.txt <- 'ftp://ftp.sanger.ac.uk/pub/databases/Pfam/releases/Pfam24.0/database_files/pfamA.txt.gz'
  system(paste('wget -nc',pf.txt))
  system('gunzip pfamA.txt.gz')
}
if (!exists("pfam")) {
  pfam <- read.table('pfamA.txt', sep="", quote="'", stringsAsFactors=F)
  pfam <- pfam[,c(2,3,5)]
  colnames(pfam) <- c('id', 'name', 'description')
}
pf.domains <- rows[,'pfam_domains']
pf.domains <- aaply(pf.domains,1, function(x) {
  if (!is.na(x) && x == '') {return('')}
  domains <- strsplit(x, ',')
  links.list <- c()
  for (domain in domains) {
    # Find the row for this domain in the pfam data frame
    pf.row <- subset(pfam, id==domain)
    if (nrow(pf.row) > 0) {
      d.link <- hmakeTag('a', pf.row$name, href=paste('http://pfam.sanger.ac.uk/family/',pf.row$id,sep=''), title=pf.row$description, target='_blank')
      links.list <- c(links.list,d.link)
    }
  }
  return(paste(links.list, collapse=', '))
})
rows[,'pfam_domain_links'] <- pf.domains

### Data clean-up
###
# Take neg-log of fisher p-values (easier sorting in web interface)
rows$pos_sel_score <- as.numeric(sprintf("%.3f",-log(rows$pval_fisher)))
if (nrow(subset(rows, pos_sel_score > 999)) > 0) {
  big.rows <- rows$pos_sel_score > 999 & !is.na(rows$pos_sel_score)
  rows[big.rows,]$pos_sel_score <- 999
}
rows$pval_fisher <- NULL

# Clean up gene names (remove dash at end)
rows$gene_name <- sub("-.*$","",rows$gene_name)


# Remove bracketed text from description.
rows[,'gene_description'] <- sub('\\[.+\\]','', rows[,'gene_description'])


# Add UCSC browser links.
rows$ucsc_link <- NA
non.na <- rows[!is.na(rows$chr_name),]
rows[!is.na(rows$chr_name),'ucsc_link'] <- sprintf("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=%s:%s-%s",non.na$chr_name,non.na$chr_start,non.na$chr_end)

# Limit to rows with > 30 sitewise values.
rows <- subset(rows, sitewise_value_count > 30)

# Rename some columns.
rename.cols <- function(data, old.new.list) {
  for (old.new in old.new.list) {
    old <- old.new[1]
    new <- old.new[2]
    if (old %in% colnames(rows)) {
      rows[,new] <- rows[,old]
      rows[,old] <- NULL
    } else {
      warning(paste("Column",old,"not found when trying to rename to",new))
      rows[,old] <- ''
      rows[,new] <- rows[,old]
      rows[,old] <- NULL
    }
  }
  return(rows)
}
new.names <- list(
  c('protein_id','Ensembl ID'),
  c('gene_id', 'Ensembl gene ID'),
  c('transcript_id', 'Ensembl transcript ID'),

  c('gene_name','Gene'),
  c('gene_description','Gene description'),
  c('pfam_domain_links', 'Pfam domains'),
  c('ucsc_link','UCSC'),

  c('aln_png','Alignment plot'),
  
  c('slr_dnds','Overall dN/dS'),
  c('pos_sel_score','Pos-sel score'),

  c('f_neg','Fraction negative'),
  c('f_pos','Fraction positive'),
  c('f_neut','Fraction neutral'),
  c('n_neg','Number negative'),
  c('n_pos','Number positive'),
  c('n_neut','Number neutral'),

  c('sitewise_value_count','Sitewise value count'),
  c('leaf_count','Sequence count'),
  c('tree_mean_path','Tree length'),
  c('duplication_count', 'Duplication count'),
  c('slr_kappa','Kappa'),
  c('gc_3','GC3'),
  c('gc_cds','GC CDS'),
  c('gc_genomic','GC genomic'),
  c('aln_length', 'Alignment length'),
  c('seq_length','Seq length'),

  c('chr_name','Chr'),
  c('chr_start','Start'),
  c('chr_end','End'),
  c('chr_strand','Strand'),

  c('tree_file','Tree file'),
  c('aln_file','Alignment file'),
  c('pep_aln_file','Protein alignment file'),
  c('sites_file','SLR file'),
  c('data_file','Data file'),
  c('params_file','Debug file')
)
rows <- rename.cols(rows, new.names)

notes.list <- list(

  c('Ensembl ID',"The Ensembl ID for the reference peptide. The first
  human or mouse peptide was used as the reference if possible;
  otherwise, the first peptide ID was used. Hyperlinks to the Ensembl
  overview of the reference peptide are provided."),

  c('Gene',"The Ensembl gene name for the reference peptide (see
  above), if available. Hyperlinks to the GeneCards page for the
  reference gene are provided."),

  c('Gene description',"A description of the protein, from Ensembl."),

  c('Chr',"The chromosome or contig on which the given reference
  peptide resides. Note: this is relative to the <em>reference</em>
  peptide, so an alignment with a non-human reference (see above) will
  have non-human chromosome coordinates."),

  c('Start',"The start position of the coding sequence."),
  
  c('End',"The end position of the coding sequence."),

  c('Strand',"The strand of the coding sequence."),

  c('Alignment plot',"A plot of the alignment and mammalian sitewise
  dN/dS data. The alignment is below, with Ensembl IDs of all
  constituent peptides on the left and the aligned sequences on the
  right. The Taylor coloring scheme is used except that missing or
  filtered sequence is colored black. On top, the sitewise data is
  plotted as a series of bars, where the lower and upper ranges of the
  bar corresponds to the lower and upper bounds of the 95% confidence
  interval on dN/dS, plotted on a logarithmic scale. Bars are colored
  according to the sitewise SLR statistic, ranging from blue
  (confident purifying selection) to gray (neutral selection) and red
  (positive selection).<p>The thumbnail is linked to a PDF rendering
  of the same alignment, spread over a number of lines to make a
  closer analysis possible, and including sitewise data from all
  clades.<br/><em>Note: scale bars are not shown, but the complete SLR data
  for each gene can be accessed by clicking the link in the 'SLR data'
  column.</em>"),

  c('Fraction negative',"The fraction of sites, over all sites in the alignment which
  were analyzed, that show significant evidence for negative selection at a
  nominal 1% false positive rate (e.g., not corrected for multiple
  testing)."),
  
  c('Fraction positive',"The fraction of sites, over all sites in the alignment which
  were analyzed, that show significant evidence for positive selection at a nominal
  1% false positive rate (e.g., not corrected for multiple testing)."),

  c('Fraction neutral',"The fraction of sites, over all sites in the
  alignment which were analyzed, that show no significant evidence for
  either positive or negative selection at a nominal 1% false positive
  rate (e.g., not corrected for multiple testing)."),
  
  c('Sequence count',"The number of sequences in the alignment."),

  c('Tree length',"The mean root-to-tip path length of the tree."),

  c('Overall dN/dS',"The overall dN/dS parameter obtained from SLR,
  resulting from the optimization of the data under a M0 (one-rate)
  codon model of evolution."),

  c('Kappa',"The kappa parameter obtained from SLR, resulting from the
  optimization of the data under a M0 (one-rate) codon model of evolution."),

  c('GC',"The mean GC content of all coding sequences in the alignment."),

  c('Seq length',"The average length of all coding sequences in the alignment."),

  c('Pos-sel score',"A combined score for positive selection in the
  alignment. Fisher's method for combining independent p-values was
  employed: the nominal p-values for all sites with dN/dS > 0 were
  multiplied together and Bonferroni corrected, yielding a combined
  p-value for positive selection. This score is the negative logarithm
  of that p-value."),

  c('Tree file',"The Ensembl gene tree for the current alignment, in Newick format."),

  c('Alignment file',"The current alignment in Fasta format."),

  c('SLR file',"The sitewise values from SLR in tab-delimited format.")
)

notes.list.items <- c()
for (items in notes.list) {
  cur.str <- paste(
    "<li><a name='",items[1],"'></a>",
    "<b>",items[1],": </b>",
    items[2],
    "</li>",
    sep=''
  )
  notes.list.items <- c(notes.list.items,cur.str)
}
table.notes <- paste(notes.list.items,collapse=' ')

### /end Data clean-up


### BigBed / BigWig output
###
genes.file.links <- c()
color.fields <- c('Overall dN/dS','Pos-sel score','Tree length')
color.lo <- c('black','black','black','black')
color.hi <- c('red','red','red','red')

for (i in 1:length(color.fields)) {
  cur.field <- color.fields[i]
  print(cur.field)
  file.name <- bed.output(rows,cur.field,color.lo[i],color.hi[i])
  genes.file.links <- c(genes.file.links,hmakeTag('a', cur.field, href=paste('./',file.name,sep='')))
}

get.genes.url <- function(positions,fields) {
  output.txt <- c()
  for (i in 1:length(fields)) {
    cur.pos <- positions[i]
    cur.field <- fields[i]
    file.name <- paste(cur.field,'.bb',sep='')
    file.name <- gsub('[- \\(\\)/\\\\]','_',file.name)

    hgct.custom.text <- paste('track',
                              'type=bigBed',
                              'name=%s',
                              paste('description="','Genes colored by ',cur.field,'"',sep=''),
                              'visibility=pack',
                              'itemRgb=on',
                              'bigDataUrl=%s',
                              sep=' ')
    hgct.custom.text <- sprintf(hgct.custom.text,cur.field,paste('http://www.ebi.ac.uk/~greg/mammals/',file.name,sep=''))
    hgct.custom.encoded <- URLencode(hgct.custom.text)
    cur.url <- sprintf("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=%s&hgct_customText=%s",cur.pos,hgct.custom.encoded)
    output.txt <- c(output.txt,cur.url)
  }
  return(output.txt)
}

bigwig.file <- bigwig.output()
#bigwig.file <- 'asdf.bw'

hgct.custom.text <- paste('track',
                          'type=bigWig',
                          'name=Sitewise_SLR',
                          'description="Sitewise SLR statistic in mammals"',
                          'yLineMark=0',
                          'yLineOnOff=on',
                          'autoScale=off',
                          'viewLimits=-20:10',
                          'minLimit=-20',
                          'maxLimit=10',
                          'visibility=full',
                          paste('bigDataUrl=','http://www.ebi.ac.uk/~greg/mammals/',bigwig.file,sep=''),
                          sep=' ')
hgct.custom.encoded <- URLencode(hgct.custom.text)
get.sites.url <- function(position) {
  return(sprintf("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=%s&hgct_customText=%s",position,hgct.custom.encoded))
}
###
### /end UCSC track output


# Constants.

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
css.txt <- ""

### Page of interesting examples.
###

interesting.genes <- c('TRIM5','GSTM1','EIF2AK2','EIF2S1','BST2','CEP152','TLR5')
interesting.rows <- rows[rows$Gene %in% interesting.genes,]
#print(interesting.rows)

### Pages by chromosome.
###
chrs <- unique(rows[,'Chr'])
chrs <- chrs[!is.na(chrs)]
chrs <- sub('chr','',chrs)
numeric.chrs <- sort(as.numeric(chrs[!is.na(as.numeric(chrs))]))
chrs <- c(numeric.chrs,chrs[is.na(as.numeric(chrs))])
chrs <- paste('chr',chrs,sep='')
#for (chr in c('chr22')) {
for (chr in chrs) {
  if (is.na(chr)){next}
  sub.rows <- subset(rows, Chr==chr)
  sub.rows <- sub.rows[order(sub.rows[,'Start']),]
  write.rows.page(chr,sub.rows,css.txt)
}
chr.links <- hmakeTag('a', chrs, href=paste('./',chrs,'.html',sep=''), class='chr')

### Top / bottom 100 by field.
all.fields <- colnames(rows)
numeric.fields <- c()
top.files <- c()
bot.files <- c()
#for (field in c('Fraction neutral','Fraction positive','Fraction negative')) {
for (field in all.fields) {
  if (!is.numeric(rows[,field])){next}

  remove.numerics <- c(
    'parameter_set_id',
    'data_id',
    'node_id',
    'job_id',
    'gene_taxon_id',
    'Number negative',
    'Number neutral',
    'Sequence count'
  )
  if (field %in% remove.numerics) { next }

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
top.links <- hmakeTag('a', 'Top', href=paste('./',top.files,sep=''))
bottom.links <- hmakeTag('a', 'Bottom', href=paste('./',bot.files,sep=''))

top.bottom.links <- paste('<div class="top-bot">',numeric.fields,'<br/>',top.links,' | ',bottom.links,'</div>')

### Sitewise links for chromosomes.
sites.file.links <- hmakeTag('a', 'SLR statistic in mammals', href=paste('./',bigwig.file,sep=''))
sites.chr.links <- hmakeTag('a', chrs, class='chr', href=get.sites.url(chrs))

sites.genes <- list(
  c('RHD','[1]'),
  c('PTPRCAP','[2]'),
  c('PRM1','[3]'),
  c('BRCA1','[4]'),
  c('IL8','[5]'),
  c('IL6','[5]'),
  c('IL4','[5]'),
  c('ESX1','[6]'),
  c('EIF2AK2','(PKR) [7]'),
  c('BST2','[8]'),
  c('CEP152','[9]')
)

sites.genes.links <- c()
dnds.ecdf <- ecdf(rows[,'Overall dN/dS'])
pos.ecdf <- ecdf(rows[,'Pos-sel score'])
neg.ecdf <- ecdf(rows[,'Fraction negative'])
for (couplet in sites.genes) {
  row <- subset(rows,Gene==couplet[1])
  line.ensp <- row[,'Ensembl ID']
  line.name <- row[,'Gene']
  line.chr <- row[,'Chr']
  line.start <- row[,'Start']
  line.end <- row[,'End']
  line.stats <- sprintf("dnds %.2f / pos-sel %.2f / neg-sel %.2f", dnds.ecdf(row[,'Overall dN/dS']), pos.ecdf(row[, 'Pos-sel score']), neg.ecdf(row[, 'Fraction negative']))

  link.page <- hmakeTag('a','Summary',href=paste("./",line.chr,".html#",line.ensp,sep=''))
  line.pos <- paste(line.chr,':',line.start,'-',line.end,sep='')
  link.ucsc <- hmakeTag('a','UCSC',href=get.sites.url(line.pos))

  cite <- couplet[2]

  line.all <- paste(line.name,' ',cite,"<br/>",link.page,' | ',link.ucsc,sep='')
  line.all <- paste("<div class='top-bot'>",line.all,"</div>")
  sites.genes.links <- c(sites.genes.links,line.all)
}

# Tip -- find rows of a specific gene: rows[grep('^cep152',rows$Gene,ignore.case=T),]

### Genewise links for chromosomes.
genes.links.chrs <- c('chr22','chr19','chr10','chrX')
genes.links.fields <- c('GC','Overall dN/dS','Pos-sel score','Tree length')
genes.links.labels <- paste(genes.links.chrs,genes.links.fields,sep=' ')
genes.links <- hmakeTag('a', genes.links.labels, class='top-bot', href=get.genes.url(genes.links.chrs,genes.links.fields))

### Output a page with all genes by name.
print("Writing page all.html")
page <- openPage('all.html',
  title='All genes',
  link.css='webpage_output.css',
)
hwrite("<div style='width:800px;margin:auto;'>",page)
rows <- rows[order(rows[,'Gene']),] # Sort by gene name.
all.genes.links <- hmakeTag('a', rows[,'Gene'], class='top-bot', href=paste("./",rows[,'Chr'],".html#",rows[,'Ensembl ID'],sep=''))
all.genes.links.string <- paste(all.genes.links,collapse="\n")
hwrite(all.genes.links.string,page)
hwrite("</div>",page)
closePage(page)


###
### HTML template output.
###
html.template <- file("~/src/greg-ensembl/projects/2xmammals/webpage_template.html")
html.lines <- readLines(html.template)
css.file <- file("~/src/greg-ensembl/projects/2xmammals/webpage_template.css")
css.lines <- readLines(css.file)
close(css.file)
close(html.template)

css.string <- paste(css.lines,collapse='\n',sep='')
html.string <- paste(html.lines,collapse='\n',sep='')

html.string <- sub('CHROMOSOME_LINKS',paste(chr.links,collapse=' '),html.string)
html.string <- sub('TOP_BOTTOM_LINKS',paste(top.bottom.links,collapse=' '),html.string)
html.string <- sub('SITES_FILES',paste(sites.file.links,collapse=', '),html.string)
html.string <- sub('SITES_LINKS',paste(sites.chr.links,collapse=' '),html.string)
html.string <- sub('SITES_GENES_LINKS',paste(sites.genes.links,collapse=' '),html.string)
html.string <- sub('GENES_FILES',paste(genes.file.links,collapse=', '),html.string)
html.string <- sub('GENES_LINKS',paste(genes.links,collapse=' '),html.string)
html.string <- sub('TABLE_NOTES',table.notes,html.string)

writeChar(css.string,paste(folder,"/webpage_output.css",sep=''))
writeChar(html.string,paste(folder,"/index.html",sep=''))

# Back to original working dir.
setwd(orig.wd)
