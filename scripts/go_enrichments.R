filename <- NA
for (e in commandArgs(trailingOnly=T)) {
  ta = strsplit(e,"=",fixed=TRUE)
  f <- ta[[1]]
  v <- f[1]
  if (v=='file') {
    filename <- f[2]
  } else if (v=='dir') {
    direction <- f[2]
  }
}
#print(filename)

library(topGO)
library(doBy)
library(plyr)
library(Hmisc)
library(goseq)

getScoresForGenes <- function(GOdata, whichGenes) {
  all.scores <- geneScore(GOdata)
  all.names <- genes(GOdata)
  names(all.scores) <- all.names

  all.scores <- all.scores[names(all.scores) %in% whichGenes]
  return(all.scores)
}

getRowsForTerms <- function(GOdata,data,whichTerms) {
  genes <- getGenesForTerms(GOdata,whichTerms)

  rows <- data[data$Hsap_protein %in% genes,]
  return(rows)
}

getGenesForTerms <- function(GOdata,whichTerms) {
  whichTerms <- unlist(whichTerms)
  genes.within.terms <- c()
  for (term in whichTerms) {
    #print(term)
    cur.genes <- genesInTerm(GOdata,term)
    if (length(cur.genes) > 0) {
      cur.genes <- cur.genes[[1]]
      genes.within.terms <- c(genes.within.terms,cur.genes)
    }
  }
  genes.within.terms <- unique(genes.within.terms)
  return(genes.within.terms)
}

getScoresForTerms <- function(GOdata, whichTerms, useRank=FALSE, aboveZero=FALSE) {
  all.scores <- geneScore(GOdata)
  all.names <- genes(GOdata)
  names(all.scores) <- all.names
  if (aboveZero) {
    all.scores <- all.scores[all.scores >= 0]
  }
  if (useRank) {
    all.scores <- rank(all.scores,ties.method="random")
  }

  genes.within.terms <- getGenesForTerms(GOdata,whichTerms)

  scores.within.terms <- all.scores[names(all.scores) %in% genes.within.terms]
  return(scores.within.terms)
}

getTermsDefinition <- function(whichTerms, ontology, numChar = 20, multipLines = FALSE) {

  qTerms <- paste(paste("'", whichTerms, "'", sep = ""), collapse = ",")
  retVal <- dbGetQuery(GO_dbconn(), paste("SELECT term, go_id FROM go_term WHERE ontology IN",
                                          "('", ontology, "') AND go_id IN (", qTerms, ");", sep = ""))

  termsNames <- retVal$term
  names(termsNames) <- retVal$go_id

  if(!multipLines)
    shortNames <- paste(substr(termsNames, 1, numChar),
                        ifelse(nchar(termsNames) > numChar, '...', ''), sep = '')
  else
    shortNames <- sapply(termsNames,
                         function(x) {
                           a <- strwrap(x, numChar)
                           return(paste(a, sep = "", collapse = "\\\n"))
                         })

  names(shortNames) <- names(termsNames)
  return(shortNames[whichTerms])
}

## methodsSig contains a named vector of p-values for each run method
sigAllMethods <- function (methodsSig)
{
  names.index <- names(methodsSig[[1]])
  retval <- as.data.frame(lapply(methodsSig, function(x) x[names.index]))
  names(retval) <- names(methodsSig)
  return(retval)
}


get.go.table = function(scores,symbols.to.go,ontology,nodeSize=8,description='',scoresAreBinary=FALSE,p.adjust='none',scoreOrder="decreasing") {
  if (scoresAreBinary) {
    geneSelectionFun = function(score){return(score >= 1)}
  } else {
    geneSelectionFun = function(score){return(score < 0.05)}
  }

  GOdata <- new("topGOdata",
    ontology = ontology,
    allGenes = scores,
    annot = annFUN.gene2GO,
    gene2GO = symbols.to.go,
    geneSelectionFun = geneSelectionFun,
    nodeSize = nodeSize,
    description = description
  )

  print(paste("  # of genes", length(scores)))

  n.terms <- length(usedGO(GOdata))

  if (scoresAreBinary) {
    
    fis <- runTest(GOdata,algorithm="classic",statistic="Fisher")
    topgo <- runTest(GOdata,algorithm="weight01",statistic="Fisher")

    print("Generating table...")

    res.list <- list(
      pval.fis=fis,
      pval.topgo=topgo
    )
    res.list <- lapply(res.list, score)
    res.df <- sigAllMethods(res.list)
    whichTerms <- rownames(res.df)[1:n.terms]

    shortNames <- getTermsDefinition(whichTerms, ontology(GOdata), numChar = 50)
    annoStat <- termStat(GOdata, whichTerms)
    infoMat <- data.frame('GO ID' = whichTerms, 'Term' = shortNames, stringsAsFactors = FALSE)
    scores.df <- data.frame(annoStat,infoMat,res.df,check.names=F,stringsAsFactors=F)
    test.df <- scores.df
    test.df <- orderBy(~pval.topgo,data=test.df)

  } else {
    # scoreOrder = 'increasing' if we're using p-values, 'decreasing' if we're using a higher-is-better score.
    topgo <- runTest(GOdata,algorithm="weight01",statistic="KS",scoreOrder=scoreOrder)
    ks <- runTest(GOdata,algorithm="classic",statistic="KS",scoreOrder=scoreOrder)

    res.list <- list(
      pval.topgo=topgo,
      pval.ks=ks
    )
    res.list <- lapply(res.list, score)
    res.df <- sigAllMethods(res.list)
    whichTerms <- rownames(res.df)[1:n.terms]

    shortNames <- getTermsDefinition(whichTerms, ontology(GOdata), numChar = 50)
    annoStat <- termStat(GOdata, whichTerms)
    infoMat <- data.frame('GO ID' = whichTerms, 'Term' = shortNames, stringsAsFactors = FALSE)
    scores.df <- data.frame(annoStat,infoMat,res.df,check.names=F,stringsAsFactors=F)
    test.df <- scores.df

    test.df <- orderBy(~pval.ks,data=test.df)
  }

  assign("last.godata", GOdata, envir=.GlobalEnv)
  return(test.df)
}

get.enrich.by.score = function(named.scores,go.df=NULL,go.vec=NULL,go.field.name='node_id',nodeSize=8) {

  if (is.null(go.vec) && !is.null(go.df)) {
    print("Splitting df into vec!")
    go.vec = strsplit(go.df$go,split=",",fixed=T)
    names(go.vec) = go.df[,go.field.name]
  }

  bp = get.go.table(
    scores = named.scores,
    symbols.to.go = go.vec,
    ontology = 'BP',
    nodeSize = nodeSize,
    description = '',
    scoresAreBinary = FALSE,
    scoreOrder = "decreasing"
  )
  return(bp)
}

get.enrich.by.subset = function(subset,all,go.df=NULL,go.vec=NULL,go.field.name='node_id',nodeSize=8) {
  
  if (is.null(go.vec) && !is.null(go.df)) {
    print("Splitting df into vec!")
    go.vec = strsplit(go.df$go,split=",",fixed=T)
    names(go.vec) = go.df[,go.field.name]
  }

  scores = as.integer(all %in% subset)
  names(scores) = all

  bp = get.go.table(
    scores = scores,
    symbols.to.go = go.vec,
    ontology = 'BP',
    nodeSize = nodeSize,
    description = '',
    scoresAreBinary = TRUE
  )
  return(bp)
}

enrich.list <- function(subset.ids,all.ids=NULL,include.iea=TRUE) {
  if (is.null(all.ids)) {
    print("No background given -- using all Ensembl Protein IDs!")
    all.ids <- ens.genes$Ensembl.Protein.ID
  } 
  #print(paste("Number of foreground IDs:",length(unique(subset.ids))))
  #print(paste("Number of background IDs:",length(unique(all.ids))))

  print(paste("  enriching", length(unique(subset.ids)), "out of", length(unique(all.ids))))

  go.df <- NULL
  go.vec <- NULL
  if (include.iea) {
    go.vec <- go.ens
  } else {
    go.vec <- go.ens.excl.iea
  }
  #print(head(go.vec,n=5))

  tbl <- get.enrich.by.subset(
    subset=unique(subset.ids),
    all=unique(all.ids),
    go.vec=go.vec,
    go.df=NULL,
    go.field.name='stable_id'
  )
  return(tbl)
}

enrich.dir <- function(dir) {
  old.wd <- getwd()
  setwd(dir)
  files <- list.files(path=dir,pattern="*_sorted.csv$")
  for (cur.file in files) {
    enrich.file(cur.file,'up')
    enrich.file(cur.file,'down')
  }
  setwd(old.wd)
}

load.go.csv <- function(csv.f, out.f) {
  print("Loading GO data...")

  go.csv = read.csv(csv.f, header=T, stringsAsFactors=F)
  print(head(go.csv))
  go.ens = by(go.csv,go.csv$Ensembl.Protein.ID, function(df) {
    df$GO.Term.Accession
  })
  go.csv.excl.iea = subset(go.csv, GO.Term.Evidence.Code != 'IEA')
  go.ens.excl.iea = by(go.csv.excl.iea, go.csv.excl.iea$Ensembl.Protein.ID, function(df) {
      df$GO.Term.Accession
  })
  save(go.ens, go.ens.excl.iea, file=out.f)
}

enrich.file <- function(cur.file,dir) {
  if (!exists("go.ens")) {
    if (file.exists("go_60.Rdata")) {
      print("Loading R data...")
      load("go_60.Rdata",envir=.GlobalEnv)
    } else {
      print("Loading GO data...")
      # Download link:
      # http://www.ensembl.org/biomart/martview/11303929625047b175098c75bacc9456/11303929625047b175098c75bacc9456/11303929625047b175098c75bacc9456?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_peptide_id|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id|hsapiens_gene_ensembl.default.feature_page.go_biological__id|hsapiens_gene_ensembl.default.feature_page.go_biological_process_linkage_type&FILTERS=hsapiens_gene_ensembl.default.filters.biol_process_evidence_code."IC,IDA,IEA,IEP,IGI,IMP,IPI,ISS,NAS,ND,TAS,EXP"&VISIBLEPANEL=filterpanel
      go.csv = read.csv("ens_go_60.csv",header=T,stringsAsFactors=F)
      go.ens = by(go.csv,go.csv$Ensembl.Protein.ID,function(df){df$GO.Term.Accession..bp.})
      go.csv.excl.iea = subset(go.csv,GO.Term.Evidence.Code..bp. != 'IEA')
      go.ens.excl.iea = by(go.csv.excl.iea,go.csv.excl.iea$Ensembl.Protein.ID,function(df){df$GO.Term.Accession..bp.})
      assign('go.ens',go.ens,envir=.GlobalEnv)
      assign('go.ens.excl.iea',go.ens.excl.iea,envir=.GlobalEnv)

      print("Loading Ensembl gene list...")
      # http://www.ensembl.org/biomart/martview/11303929625047b175098c75bacc9456/11303929625047b175098c75bacc9456/11303929625047b175098c75bacc9456?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_peptide_id|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id|hsapiens_gene_ensembl.default.feature_page.clone_based_ensembl_gene_name|hsapiens_gene_ensembl.default.feature_page.external_gene_id&FILTERS=&VISIBLEPANEL=resultspanel
      ens.genes = read.csv("ens_genes_60.csv",header=T,stringsAsFactors=F)
      assign('ens.genes',ens.genes,envir=.GlobalEnv)
    }
  }

  rdata.file <- load.csv.save.rdata(cur.file,dir)
  load(rdata.file)

  # Do a KS enrichment by score.
  dir.lbl = dir
  print(dir.lbl)

  named.scores <- data[,'lrt.signed']
  names(named.scores) <- data[,'Hsap_protein']
  above.zero <- named.scores[named.scores > 0]

  print(str(named.scores))

  ks.tbl <- get.enrich.by.score(
    named.scores=above.zero,
    go.vec=go.ens,
    go.field.name='stable_id'
  )
  ks.tbl$category <- ''
  ks.tbl[ks.tbl$pval.ks < 0.05,]$category <- '05'
  ks.tbl[ks.tbl$pval.ks < 0.01,]$category <- '01'
  ks.tbl <- orderBy(~pval.topgo,data=ks.tbl)
  write.csv(ks.tbl,file=paste(cur.file,"_enrich_ks_",dir.lbl,".csv",sep=""))

  # Take the top 100, 250, 500 genes
  for (top.n in c(50,100,250,500)) {
    d.sub <- data[1:top.n,]
    d.sub <- subset(d.sub,lrt.signed > 0)
    print(nrow(d.sub))
    subset.ids <- d.sub[,c('Hsap_protein')]
    enrich.tbl <- enrich.list(subset.ids=subset.ids,all.ids=data[,c('Hsap_protein')])

    enrich.tbl$category <- ''
    enrich.tbl[enrich.tbl$pval.fis < 0.05,]$category <- '05'
    enrich.tbl[enrich.tbl$pval.fis < 0.01,]$category <- '01'
    
    write.csv(enrich.tbl,file=paste(cur.file,"_enrich_",dir.lbl,'_',top.n,".csv",sep=""))
  }
}

load.csv.save.rdata <- function(file,direction) {
  rdata.file <- paste(file,"_go_data_",direction,".Rdata",sep="")
  if (file.exists(rdata.file)) {
    return(rdata.file)
  }

  data <- read.csv(file,header=T)
  if (direction == 'down') {
    data$lrt.signed = -data$lrt.signed
  }

  data <- orderBy(~-lrt.signed,data=data)

  named.scores <- data[,'lrt.signed']
  names(named.scores) <- data[,'Hsap_protein']

  # Save a GOdata object.
  go.data <- new("topGOdata",
    ontology = 'BP',
    allGenes = named.scores,
    annot = annFUN.gene2GO,
    gene2GO = go.ens,
    nodeSize = 8,
    geneSelectionFun = function(score){return(score <= 0.05)}
  )
  save(data,go.data,file=rdata.file)
  return(rdata.file)
}

get.term.groups <- function() {
  ## Gorilla-specific acceleration
  male.aggression = c(
    'GO:0001820', # serotonin secretion
    'GO:0014063', # neg. reg. ser. secr.
    'GO:0006837', # ser. transport
    'GO:0042427', # Ser. biosyn. proc.
    'GO:0006587', # Ser. biosyn. from. tryptophan
    'GO:0051583', # dopamine uptake
    'GO:0015872' # dopamine transport
  )

  sexual.dimorphism = c(
    'GO:0046661', # Male sex diffrntn.
    'GO:0019102', # Male somatic sex diffrntn.
    'GO:0007530', # Sex differentiation
    'GO:0008209', # Androgen metabolic process
    'GO:0030518', # steroid hormone receptor signaling
    'GO:0035258', # steroid hormone receptor binding
    'GO:0048806', # genitalia devel.
    'GO:0005497'  # androgen binding
  )

  growth.hormone = c(
    'GO:0030252', # gr. hr. secre.
    'GO:0060416', # resp. to gr. hr. stimul.
    'GO:0060398', # regul. of gr. hr. secr. pthwy.
    'GO:0060123' # regul. of gr. hr. secr.
  )

  sperm.competition = c(
    'GO:0046546', # devel. of primary male sex characteristics
    'GO:0042713', # sperm ejaculation
    'GO:0007286', # spermatid development
    'GO:0030317', # sperm motility
    'GO:0048232', # male gamete generation
    'GO:0007283', # spermatogenesis
    'GO:0008584', # male gonad development
    'GO:0009566'  # fertilization
  )

  diet.change = c(
    'GO:0007586', # digestion
    'GO:0044241', # lipid digestion
    'GO:0022600', # digestive system process
    'GO:0055123'  # digestive system development
  )

  tooth.development = c(
    'GO:0042476', # odontogenesis
    'GO:0070166', # enamel mineralization
#    'GO:0046848', # hydroxyapatite binding
    'GO:0070172', # positive regulation of tooth mineralization
    'GO:0034505', # tooth mineralization
    'GO:0030345', # structural constituent of tooth enamel
#    'GO:0032967', # pos. regl. of collagen biosynthetic process
    'GO:0031214'  # biomineral tissue development
  )

  heart.development = c(
    'GO:0007507', # heart devel.
    'GO:0001570', # vasculogenesis
    'GO:0001525', # angiogen.
    'GO:0035050', # embryonic heart tube devel.
    'GO:0003007'  # heart morphogenesis
  )

  placental.development = c(
    'GO:0001890', # plac. devel.
    'GO:0001893', # mat. plac. devel.
    'GO:0060709', # glycog. cell devel. involved in plac. devel.
    'GO:0060716', # labyrinthine layer blood vessel devel.
    'GO:0060713'  # labyrinthine layer morphogenesis
  )

  brain.size = c(
    'GO:0007420', # brain development
    'GO:0007417', # central nervous system development
    'GO:0048854', # brain morphogen.
    'GO:0021987', # cerebral cortex devel.
    'GO:0021904', # dorsal / ventral neural tube patterning
    'GO:0021536', # diencephalon devel.
    'GO:0021871', # forebrain regionalization
    'GO:0021953', # CNS neuron differentiation
    'GO:0021872' # generation of neurons in forebrain
  )

  body.fluid = c(
    'GO:0050878' # regulation of body fluid levels
  )

  innate.immunity = c(
    'GO:0045087' # innate immune response
  )

  adaptive.immunity = c(
    'GO:0002250' # adaptive immune response
  )

  defense.response = c(
    'GO:0006952', # defense response
    'GO:0009615', # resp. to vir.
    'GO:0050832', # def. resp. to fungus
    'GO:0042742' # def. resp. to bact.
  )

  sound.perception = c(
    'GO:0007605', # sensory perception of sound
    'GO:0042472', # inner ear morphogenesis
    'GO:0048839' # innear ear development
  )


  hair.skin = c(
    'GO:0001942', # hair follicle devel.
    'GO:0030855', # epithelial cell differentiation
    'GO:0060429', # epithelium development
    'GO:0050678'  # regulation of epithelial cell differentiation
  )

  keratin = c(
    'GO:0030216' # keratinocyte development
  )

  fat = c(
    'GO:0045444' # fat cell devel.
  )

  ### Go through and get the subsets of genes
  list.of.term.lists <- list(
#    male.aggression = male.aggression,
#    sexual.dimorphism = sexual.dimorphism,
#    growth.hormone = growth.hormone,
#    placental.development = placental.development,
#    diet.change = diet.change,
    sperm.competition = sperm.competition,
#    tooth.development = tooth.development,
#    heart.development = heart.development,
#    body.fluid = body.fluid,
    hair.skin = hair.skin,
    brain.size = brain.size,
    innate.immunity = innate.immunity,
    adaptive.immunity = adaptive.immunity,
    defense.response = defense.response,
    sound.perception = sound.perception,
    keratin = keratin,
    fat = fat
  )
  return(list.of.term.lists)
}

group.presence.absence <- function(name,go.terms) {
  list.of.term.lists <- get.term.groups()

  pres.abs <- llply(list.of.term.lists, function(x) {
    ret.val <- as.numeric(any(x %in% go.terms))
    if (is.na(ret.val)) {
      ret.val <- 0
    }
    return(ret.val)
  })
  names(pres.abs) <- paste('GO',names(pres.abs),sep='.')
  return(pres.abs) 
}


subset.densities <- function(files,direction='up') {

  list.of.term.lists <- get.term.groups()
  list.of.term.lists[['all']] <- c('')
#  list.of.term.lists[['clark']] <- c('')
#  print(list.of.term.lists)

  all.data <- data.frame()
  for (file in files) {
    load(file)
#    data <- read.csv(file,header=T)
    data <- data[,c('Hsap_protein','name','lrt.signed')]
    data$filename <- file
    if (direction == 'none') {
      data$rank <- rank(abs(data$lrt.signed))
    } else {
      data$rank <- rank(data$lrt.signed)
    }
    all.data <- rbind(all.data,data)
  }
  file <- 'subset.densities'
  clark.genes <- read.table('clark 2005 genes.txt',sep="\t",stringsAsFactors=F,header=T)

  each.df <- ddply(all.data, .(filename), function(data) {
    ks <- data.frame(p.value = 1)
    ws <- data.frame(p.value = 1)

    term.names <- names(list.of.term.lists)
    terms.df <- ldply(term.names, function(x) {
      ret.df <- data.frame()
      if (x == 'all') {
        rows <- data
      } else if (x == 'clark') {
        rows <- subset(data, Hsap_gene %in% clark.genes[,'ens_gene'])
        print(nrow(rows)) 
      } else {
        term.list <- list.of.term.lists[x]
        rows <- getRowsForTerms(go.data,data,term.list)
        print(paste(data[1,]$filename,x,nrow(rows)))
      }
      if (nrow(rows) > 0) {
        if (direction == 'none') {
          ks <- ks.test(rows$rank,data$rank)
          ws <- wilcox.test(rows$rank,data$rank)
        } else if (direction == 'up') {
          rows <- subset(rows,lrt.signed > 0)
          if (nrow(rows) > 0) {
            ks <- ks.test(rows$rank,data[data$lrt.signed > 0,]$rank,alternative='l')
            ws <- wilcox.test(rows$rank,data[data$lrt.signed > 0,]$rank,alternative='g')
          }
        } else {
          rows <- subset(rows,lrt.signed < 0)
          if (nrow(rows) > 0) {
            ks <- ks.test(rows$rank,data[data$lrt.signed < 0,]$rank,alternative='g')
            ws <- wilcox.test(rows$rank,data[data$lrt.signed < 0,]$rank,alternative='l')
          }
        }
        
        group.lbl <- paste(x,"\nN = ",nrow(rows),sep="")
        group.lbl <- paste(group.lbl,"\nKS = ",format.pval(ks$p.value,digits=3),sep="")
        group.lbl <- paste(group.lbl,"\nMW = ",format.pval(ws$p.value,digits=3),sep="")
  
#        print(group.lbl)
        if (nrow(rows) > 0) {
          ret.df <- data.frame(rows,group=x,filename=data[1,]$filename,lbl=group.lbl,sig=ks$p.value,sig.ks=ks$p.value,sig.mw=ws$p.value,n=nrow(rows))
        } else {
          ret.df <- data.frame(group=x,filename=data[1,]$filename,lbl=group.lbl,sig=ks$p.value,sig.ks=ks$p.value,sig.mw=ws$p.value,n=nrow(rows))
        }
      }
      return(ret.df)
    })
    return(terms.df)
  })

  library(ggplot2)
  library(plyr)

  file.diffs <- data.frame()
  files <- unique(each.df$filename)
  print(files)
  groups <- unique(each.df$group)
  for (i in files) {
    for (j in files) {
      for (k in groups) {
      if (i == '1' || j == '1') {
        next
      }
#        print(i)
         print(paste(i,j,k))
        rows.i <- subset(each.df,filename==i & group==k)
        rows.j <- subset(each.df,filename==j & group==k)
        if (nrow(rows.i) > 0 && nrow(rows.j) > 0) {
          ks <- ks.test(rows.i$rank,rows.j$rank)
          file.diffs <- rbind(file.diffs,data.frame(
            file.i = i,
            file.j = j,
            group = k,
            p = ks$p.value
          ))
#         print(format.pval(ks$p.value,digits=3))
        }
      }
    }
  }
#  print(head(file.diffs))
  file.diffs$p.bh <- p.adjust(file.diffs$p,'BH')
  write.csv(file.diffs,file=paste(file,'_',direction,'_diffs.csv',sep=''),row.names=F,quote=F)

  noted.genes <- ddply(each.df,c('group','filename'),function(x) {
    x <- orderBy(~-rank,data=x)
    if (direction == 'up') {
      row <- x[1:min(10,nrow(x)),]
    } else {
      lo <- max(0,nrow(x)-10)
      hi <- nrow(x)
      row <- x[lo:hi,]
    }
    return(data.frame(group=row$group,filename=row$filename,rank=row$rank,id=row$name))
  })
  write.csv(noted.genes,file=paste(file,'_',direction,'_genes.csv',sep=''),row.names=F,quote=F)

  n.tests <- nrow(unique(each.df[,c('filename','group')]))
  n.tests <- n.tests - length(files)

  print(n.tests)
  summary.df <- ddply(each.df,c('filename','group'), function(x) {
    row <- x[1,]
    ret <- data.frame(
      filename = row$filename,
      group = row$group,
      n = row$n,
      ks = row$sig.ks,
      mw = row$sig.mw,
      ks.bonf = pmin(1,row$sig.ks * n.tests),
      mw.bonf = pmin(1,row$sig.mw * n.tests)
    )
    return(ret)
  })

  summary.df$ks.bonf2 = p.adjust(summary.df$ks,'bonf')
  summary.df$ks.bh = p.adjust(summary.df$ks,'BH')

  write.csv(summary.df,file=paste(file,'_',direction,'_pvals.csv',sep=''),row.names=F,quote=F)

  #w <- 10
  #h <- length(list.of.term.lists) * w / 3
  w <- 20
  h <- 5

  file <- 'subset.densities'
  pdf.file <- paste(file,'_',direction,'.pdf',sep='')
  pdf(pdf.file,height=h,width=w,pointsize=8)
  #pdf.file <- paste(file,'_',direction,'.png',sep='')
  #png(pdf.file,height=2400,width=3200)
  p <- ggplot(each.df,aes(x=rank))
  p <- p + scale_fill_gradient2(low="red",mid="grey",high="grey",midpoint=-1,trans="log10")
  p <- p + stat_density(aes(
    ymax = ..scaled.., ymin = -..scaled.., fill=sig),
    colour = "grey50",
    geom = "ribbon", position = "identity")
  p <- p + geom_point(aes(y=0),colour='black',size=1,alpha=0.25,position=position_jitter(height=0))
#  p <- p + geom_text(data=noted.genes,aes(label=id),y=0,fill='black',size=4,angle=-45,hjust=1)
  p <- p + geom_text(aes(label=lbl),x=-4,y=-1,hjust=0,vjust=0,size=2)
  p <- p + facet_grid(filename ~ group)
  p <- p + opts(strip.text.y = theme_text())
  p <- p + opts(title=pdf.file)
#  p <- p + theme_gray(base_size=36)
  print(p)
  dev.off()  
}


do.a.bunch.of.enrichments <- function() {
  file.names <- c()
  for (species in c('Hsap','Ptro')) {
    for (test in c(1,2,3,4,5,6,7,8)) {
      file <- paste('lnl_m_',species,'_',test,'_sorted.csv',sep='')
      file.names <- c(file.names,file)
    }
  }

  file.names <- c(file.names,
#  file.names <- c(
    'chimp_sorted.csv',
    'human_sorted.csv',
    'gorilla_sorted.csv',
    'gc_parallel_sorted.csv',
    'hc_parallel_sorted.csv',
    'gh_parallel_sorted.csv'
  ) 

  for (file in file.names) {
    for (direction in c('down','up')) {
      if (file.exists(file)) {
        cmd <- paste(
         "bsub -o /dev/null 'R-2.12.0 --slave --args file=",
         file,
         " dir=",direction,
         " < ~/lib/greg-ensembl/scripts/go_enrichments.R '",
         sep="")
         print(cmd)
         system(cmd)
       }
    }
  }
}

test.sig.diffs <- function() {
  files <- c(
#    'chimp_sorted.csv',
#    'human_sorted.csv',
#    'gorilla_sorted.csv',
#    'gc_parallel_sorted.csv',
#    'hc_parallel_sorted.csv',
#    'gh_parallel_sorted.csv',
    "lnl_m_Hsap_3_sorted.csv",
    "lnl_m_Hsap_4_sorted.csv",
    "lnl_m_Ptro_4_sorted.csv",
    "lnl_m_Hsap_7_sorted.csv",
    "lnl_m_Hsap_8_sorted.csv"
  )

  files <- paste(files,'_go_data_up.Rdata',sep='')

  subset.densities(files,'up')
  subset.densities(files,'down')
#  subset.densities(files,'none')
}

cor.enrichments <- function(file.a,file.b) {
  data.a <- read.csv(file.a,header=T)
  data.b <- read.csv(file.b,header=T)

  data.a$pval.a <- data.a$pval.fis
  data.b$pval.b <- data.b$pval.fis

  merged <- merge(data.a,data.b,by='GO.ID')

#  merged <- subset(merged,pval.a != 1 || pval.b != 1)
  ab.cor <- cor(merged$pval.a,merged$pval.b,use='complete.obs',method='spearman')
  return(ab.cor)
}

compare.all.tests <- function() {
#  files <- list.files(path='./',pattern=".*_enrich_ks_up.*.csv$")
  files <- list.files(path='./',pattern=".*(gorilla|Hsap_3).*(up|down).*.csv$")
  res.matrix <- matrix(nrow=length(files),ncol=length(files))

  # Get all enrichments into merged data frame for corrgram.
  all.enrichments <- data.frame()
  for (file in files) {
    cur.file <- read.csv(file,header=T)
    if ('pval.ks' %in% names(cur.file)) {
      cur.file <- cur.file[,c('GO.ID','pval.ks')]
    } else if ('pval.fis' %in% names(cur.file)) {
      cur.file <- cur.file[,c('GO.ID','pval.fis')]    
    }
    file.name <- file
    file.name <- paste(substr(file,1,20),substr(file,nchar(file)-20,nchar(file)),sep="\n")
    names(cur.file) <- c('GO.ID',file.name)
    print(file)
    if (nrow(all.enrichments) > 0) {
      all.enrichments <- merge(all.enrichments,cur.file,by='GO.ID',all.x=T,all.y=F)
    } else {
      all.enrichments <- cur.file
    }
  }

  library(corrgram)
  pdf(file="enrichment_corrgram.pdf")
  print("Plotting corrgram...")
  data <- all.enrichments
  nums <- unlist(lapply(data,is.numeric))
  data <- data[,nums]
  corrgram(data,order=TRUE,upper.panel=panel.pie,cex.labels=0.5,label.pos=c(0.5,0.5))
  dev.off()

#  for (i in 1:length(files)) {
#    for (j in 1:length(files)) {
#      cur.cor <- cor.enrichments(files[i],files[j])
#      res.matrix[i,j] <- cur.cor
#      print(cur.cor)
#    }
#  }

#  write.csv(res.matrix,file='enrichment_correlation_matrix.csv',quote=F,row.names=F) 
}

do.enrichments <- function(x, lbl, go.ens) {
  # x must have the following:
  #  binary.score = a binarized score.
  #  score = the actual score
  #  length = the "gene length"
  #  id = the ID of some sort
  #  name = a more meaningful name (e.g., id=ENSP000xyz, name=BRCA2)

  geneSelFn <- function(score){return(score == 1)}
  godata.f <- scratch.f(sprintf("godata_%s.Rdata", lbl))
  if (!file.exists(godata.f)) {
    print("Godata file doesn't exist -- calculating stuff...")
    scores <- x$binary.score
    names(scores) <- x$id
    root.godata <- new("topGOdata",
      ontology = "BP",
      allGenes = scores,
      annot = annFUN.gene2GO,
      gene2GO = go.ens,
      geneSelectionFun = geneSelFn,
      nodeSize = 8,
      description = ''
    )
    go.terms <- usedGO(root.godata)
    go.pairs <- ldply(go.terms, function(x) {
      cur.term <- x
      cur.genes <- as.character(unlist(genesInTerm(root.godata, cur.term)))
      data.frame(
        gene.id = cur.genes,
        go.id = cur.term
      )
    })
    print(head(go.pairs))
    save(root.godata, go.pairs, file=godata.f)
  }
  load(godata.f)

  #go.terms <- data.frame(go_id=c('GO:0000003', 'GO:0000018', 'GO:0000070', 'GO:0006458'), stringsAsFactors=F)
  go.terms <- data.frame(go_id=usedGO(root.godata), stringsAsFactors=F)

  scores <- x$binary.score
  names(scores) <- x$id
  geneSelFn <- function(score){return(score == 1)}

  print(paste("N sig:", sum(geneSelFn(scores))))

  go.data <- root.godata
  total.genes <- length(genes(go.data))
  sub.x <- x[x$id %in% genes(go.data), ]

  pwf.scores <- as.integer(x$binary.score)
  names(pwf.scores) <- x$id
  pwf <- nullp(pwf.scores, bias.data=x$length)

  print(head(pwf))

  print("Thinning go.pairs...")
  print(nrow(go.pairs))
  go.pairs <- subset(go.pairs, gene.id %in% genes(go.data))
  go.pairs <- subset(go.pairs, go.id %in% usedGO(go.data))
  print(nrow(go.pairs))

  print(head(pwf))
  go.corr <- goseq(pwf, gene2cat=go.pairs)
  print(head(go.corr))
  go.uncorr <- goseq(pwf, gene2cat=go.pairs, method='Hypergeometric')
  print(head(go.uncorr))

  topgo.res <- runTest(go.data, algorithm="weight01", statistic="Fisher")
  topgo.pvals <- score(topgo.res)
  topgo.ids <- names(topgo.pvals)
  topgo.df <- data.frame(id=topgo.ids, pval=as.numeric(topgo.pvals), stringsAsFactors=F)

  topgo.fis.res <- runTest(go.data, algorithm="classic", statistic="Fisher")
  topgo.fis.pvals <- score(topgo.fis.res)
  topgo.fis.ids <- names(topgo.fis.pvals)
  topgo.fis.df <- data.frame(id=topgo.fis.ids, pval=as.numeric(topgo.fis.pvals), stringsAsFactors=F)

  sig.indices <- geneSelFn(sub.x$binary.score)

  go.res <- ddply(go.terms, 'go_id', function(y) {
    cur.id <- y[1, 'go_id']
    cur.id <- as.character(cur.id)
    cur.genes <- unlist(genesInTerm(go.data, cur.id))
    n.genes <- length(cur.genes)

    row.indices <- sub.x$id %in% cur.genes
    cur.rows <- sub.x[row.indices, ]
    sig.rows <- cur.rows[geneSelFn(cur.rows$binary.score), ]
    sig.rows <- sig.rows[order(sig.rows$score, decreasing=T),]

    f.in.term <- sum(row.indices) / total.genes
    f.sig <- sum(sig.indices) / total.genes
    total.n.sig <- sum(sig.indices)      
    n.expected <- f.in.term * f.sig * total.genes
    n.sig <- sum(geneSelFn(cur.rows$binary.score))

    min.n <- 1

    fisher.pval <- 1
    # Only run FET for terms with > N significant genes.
    if (n.sig >= min.n) {
      a.factor <- factor(row.indices)
      b.factor <- factor(sig.indices)
      if (length(levels(a.factor)) < 2 || length(levels(b.factor)) < 2) {
        fisher.pval <- 1
      } else {
        fisher.res <- fisher.test(a.factor, b.factor, alternative='g')
        fisher.pval <- fisher.res$p.value
      }
    }

    hyper.pval <- 1
    if (n.sig >= min.n) {
      n.white <- sum(row.indices)
      n.black <- length(row.indices) - sum(row.indices)
      n.balls.drawn <- sum(sig.indices)
      n.white.drawn <- sum(sig.indices & row.indices)
      hyper.pval <- phyper(n.white.drawn, n.white, n.black, n.balls.drawn, lower.tail=F)
    }

    mwu.pval <- 1
    ks.pval <- 1
    contains.all <- all(row.indices, TRUE)
    if (!contains.all) {
      cur.scores <- sub.x[row.indices, 'score']
      other.scores <- sub.x[!(row.indices), 'score']
      if (length(cur.scores) > 0) {
        ks.res <- ks.test(cur.scores, other.scores, alternative='l')
        ks.pval <- ks.res$p.value
      }
      if (length(cur.scores) > 0) {
        mwu.res <- wilcox.test(cur.scores, other.scores, alternative='g')
        mwu.pval <- mwu.res$p.value
      }
    }

    score.lo <- NA
    score.mid <- NA
    score.hi <- NA
    cur.scores <- sub.x[row.indices, 'score']
    if (length(cur.scores) > 0) {
      qs <- quantile(cur.scores, c(0.25, 0.5, 0.75))
      
      score.lo <- qs[1]
      score.mid <- qs[2]
      score.hi <- qs[3]
    }

    t.pval <- 1
    t.dir <- NA
    has.t.values <- any(!is.na(cur.rows$t_test_a))
    if (has.t.values) {
      t.res <- t.test(cur.rows$t_test_a, cur.rows$t_test_b, paired=T)
      t.pval <- t.res$p.value
      t.dir <- t.res$estimate
    }

    mean.length <- mean(cur.rows$length)
    goseq.corr.pval <- go.corr[go.corr$category == cur.id, ]$over
    goseq.uncorr.pval <- go.uncorr[go.uncorr$category == cur.id, ]$over
    topgo.pval <- topgo.df[topgo.df$id == cur.id, ]$pval
    topgo.fis.pval <- topgo.fis.df[topgo.fis.df$id == cur.id, ]$pval
    descr <- getTermsDefinition(cur.id, 'BP', numChar=40)
    genes.str <- paste(sig.rows$name, collapse=' ')

    cur.df <- data.frame(
      'go_id'=cur.id,
      'mean_length' = mean.length,
      'total_n' = total.genes,
      'total_sig' = total.n.sig,
      'n' = sum(row.indices),
      'n_sig' = n.sig,
      'n_exp' = n.expected,
      'fis_pval' = fisher.pval,
      'hyper_pval' = hyper.pval,
      'ks_pval' = ks.pval,
      'mwu_pval' = mwu.pval,
      'goseq_wall' = goseq.corr.pval,
      'goseq_fis' = goseq.uncorr.pval,
      'topgo_pval' = topgo.pval,
      'topgo_fis' = topgo.fis.pval,
      t_pval = t.pval,
      t_dir = t.dir,
      score_lo = score.lo,
      score_mid = score.mid,
      score_hi = score.hi,      
      descr=descr,
      genes=genes.str,
      stringsAsFactors=F
    )
    #print(sprintf("%s", cur.df$descr))
    cur.df
  })

  adj.f <- function(x, fld) {
    adj.s <- paste(fld, '_adj', sep='')
    large.terms <- x[, 'n'] >= 20
    x[, adj.s] <- 1
    x[, adj.s] <- p.adjust(x[, fld], method='BH')
    x
  }

  go.res <- adj.f(go.res, 'fis_pval')
  go.res <- adj.f(go.res, 'hyper_pval')
  go.res <- adj.f(go.res, 'ks_pval')
  go.res <- adj.f(go.res, 'mwu_pval')
  go.res <- adj.f(go.res, 'goseq_wall')
  go.res <- adj.f(go.res, 'topgo_fis')

  go.res$label <- lbl
  go.res
}
