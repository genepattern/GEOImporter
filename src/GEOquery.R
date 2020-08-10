.packageName <- "GEOquery"
"GDS2MA" <-
  function(GDS,do.log2=FALSE,GPL=NULL) {
    require(limma)
    if(is.null(GPL)) {
      GPL <- getGEO(Meta(GDS)$platform)
    }
    ord.table <- match(Table(GDS)[,1],Table(GPL)[,1])
                                        # exclude non-numeric columns
    inc.columns <- grep('GSM',colnames(Table(GDS)))
    mat <- suppressWarnings(as.matrix(apply(Table(GDS)[,inc.columns],2,function(x) {as.numeric(as.character(x))})))
    if(do.log2) {
      M <- log2(mat)
    } else {
      M <- mat
    }
    MA <- new('MAList',list(M=M,
                            A=NULL,
                            targets=Columns(GDS),
                            genes=Table(GPL)[ord.table,],
                            notes=Meta(GDS)
                            ))
    return(MA)
  }

"GDS2eSet" <-
  function(GDS,do.log2=FALSE,GPL=NULL) {
    require(Biobase)
                                        # exclude non-numeric columns
    if(is.null(GPL)) {
      GPL <- getGEO(Meta(GDS)$platform)
    }
    ord.table <- match(Table(GDS)[,1],Table(GPL)[,1])
    inc.columns <- grep('GSM',colnames(Table(GDS)))
    mat <- suppressWarnings(as.matrix(apply(Table(GDS)[,inc.columns],2,
                                            function(x) {as.numeric(as.character(x))})))
    if(do.log2) {
      expr <- log2(mat)
    } else {
      expr <- mat
    }
    rownames(expr) <- as.character(Table(GDS)$ID_REF)
    tmp <- Columns(GDS)
    rownames(tmp) <- as.character(tmp$sample)
    pheno <- new("AnnotatedDataFrame",data=tmp)
    mabstract=ifelse(is.null(Meta(GDS)$description),"",Meta(GDS)$description)
    mpubmedids=ifelse(is.null(Meta(GDS)$pubmed_id),"",Meta(GDS)$pubmed_id)
    mtitle=ifelse(is.null(Meta(GDS)$title),"",Meta(GDS)$title)
    dt <- Table(GPL)
    rownames(dt) <- as.character(dt$ID)
    featuredata <- new('AnnotatedDataFrame',data=dt[ord.table,],
                       varMetadata=data.frame(Column=Columns(GPL)[,1],
                         labelDescription=Columns(GPL)[,2]))
    eset <- new('ExpressionSet',exprs=expr,phenoData=pheno,
                featureData=featuredata,
                experimentData=new("MIAME",
                  abstract=mabstract,
                  title=mtitle,
                  pubMedIds=mpubmedids,
                  other=Meta(GDS)))
    return(eset)
  }
### $Id: classes.R 20111 2006-09-25 20:40:09Z sdavis2@mail.nih.gov $

### Generic GEO classes:

setClass("GEODataTable",
         representation(
                        columns="data.frame",
                        table="data.frame"
                        ),
         prototype=list(
           columns=data.frame(matrix(nr=0,nc=0)),
           table=data.frame(matrix(nr=0,nc=0))
           )
         )


setClass("GEOData",
         representation(
                        header="list"
                        ),
         prototype=list(
           header=list()
         )
         )

setClass("GSE",
         representation(
                        header="list",
                        gsms="list",
                        gpls="list"
                        ),
         prototype=list(
           header=list(),
           gsms=list(),
           gpls=list()
         )
         )

setClass("GPL",
         representation(
                        dataTable='GEODataTable'
                        ),
         prototype=list(
           dataTable=new("GEODataTable")
           ),
         contains="GEOData")
setClass("GSM",
         representation(
                        dataTable='GEODataTable'
                        ),
         prototype=list(
           dataTable=new("GEODataTable")
           ),
         contains="GEOData")
setClass("GDS",
         representation(
                        gpl="GPL",
                        dataTable='GEODataTable'
                        ),
         prototype=list(
           gpl=new("GPL"),
           dataTable=new("GEODataTable")
           ),
         contains="GEOData")

printHead <- function(x)
#  Print leading 5 elements or rows of atomic object
#  From limma and Gordon Smyth
{
        if(is.atomic(x)) {
                d <- dim(x)
                if(length(d)<2) which <- "OneD"
                if(length(d)==2) which <- "TwoD"
                if(length(d)>2) which <- "Array"
        } else {
                if(inherits(x,"data.frame")) {
                        d <- dim(x)
                        which <- "TwoD"
                } else
                        which <- "Recursive"
        }
        switch(which,
        OneD={
                n <- length(x)
                if(n > 20) {
                        print(x[1:5])
                        cat(n-5,"more elements ...\n")
                } else
                        print(x)
        },
        TwoD={
                n <- d[1]
                if(n > 10) {
                        print(x[1:5,])
                        cat(n-5,"more rows ...\n")
                } else
                        print(x)
        },
        Array={
                n <- d[1]
                if(n > 10) {
                        dn <- dimnames(x)
                        dim(x) <- c(d[1],prod(d[-1]))
                        x <- x[1:5,]
                        dim(x) <- c(5,d[-1])
                        if(!is.null(dn[[1]])) dn[[1]] <- dn[[1]][1:5]
                        dimnames(x) <- dn
                        print(x)
                        cat(n-5,"more rows ...\n")
                } else
                        print(x)
        },
        Recursive=print(x)
        )
}

chopColumns <- function(x) {
  apply(x,2,substring,first=1,last=35)
}

setGeneric("Meta",
           function(object) standardGeneric("Meta"))
setGeneric("Accession",
           function(object) standardGeneric("Accession"))
setGeneric("dataTable",
           function(object) standardGeneric("dataTable"))
setGeneric("Columns",
           function(object) standardGeneric("Columns"))
setGeneric("GPL",
           function(object) standardGeneric("GPL"))
setGeneric("GSM",
           function(object) standardGeneric("GSM"))
setGeneric("GPLList",
           function(object) standardGeneric("GPLList"))
setGeneric("GSMList",
           function(object) standardGeneric("GSMList"))
setGeneric("Table",
           function(object) standardGeneric("Table"))


setMethod("Meta","GEOData",
          function(object) {
            return(object@header)
          }
          )
setMethod("Meta","GSE",
          function(object) {
            return(object@header)
          }
          )
          
setMethod("Accession","GEOData",
		 function(object) {
		   return(Meta(object)$geo_accession)
		 }
		 )

setMethod("dataTable","GEOData",
          function(object) {
            return(object@dataTable)
          }
          )

setMethod("Table","GEODataTable",
          function(object) {
            return(object@table)
          }
          )
setMethod("Columns","GEODataTable",
          function(object) {
            return(object@columns)
          }
          )

setMethod("Table","GEOData",
          function(object) {
            return(Table(dataTable(object)))
          }
          )

setMethod("Columns","GEOData",
          function(object) {
            return(Columns(dataTable(object)))
          }
          )

setMethod("GPLList","GSE",
          function(object) {
            return(object@gpls)
          }
          )
setMethod("GPL","GDS",
          function(object) {
            return(object@gpl)
          }
          )
setMethod("GSMList","GSE",
          function(object) {
            return(object@gsms)
          }
          )
          

setMethod("show","GEODataTable",
          function(object) {
            cat("An object of class \"",class(object),"\"\n",sep="")
            cat("****** Column Descriptions ******\n")
            print(Columns(object))
                  cat("****** Data Table ******\n")
            printHead(Table(object))
          })

setMethod("show","GEOData",
          function(object) {
            cat("An object of class \"",class(object),"\"\n",sep="")
            for (i in names(Meta(object))) {
              cat(i,"\n")
              print(Meta(object)[[i]])
            }
            print(dataTable(object))
          })
getGEO <- function(GEO=NULL,
                   filename=NULL,
                   destdir=tempdir(),
                   GSElimits=NULL) {
  filename <- filename
  if(!is.null(GSElimits)) {
    if(length(GSElimits)!=2) {
      stop('GSElimits should be an integer vector of length 2, like (1,10) to include GSMs 1 through 10')
    }
  }
  if(is.null(GEO) & is.null(filename)) {
    stop("You must supply either a filename of a GEO file or a GEO accession")
  }
  if(is.null(filename)) {
    GEO <- toupper(GEO)
    filename <- getGEOfile(GEO,destdir=destdir)
  }
  if(length(grep('\\.gz$',filename,perl=TRUE))>0) {
    con <- gzfile(filename,'r')
  } else {
    con <- file(filename,'r')
  }
  ret <- parseGEO(con,GSElimits)
  close(con)
  return(ret)
}
                   
                   
getGEOfile <- function(GEO,destdir=tempdir(),
                       amount=c('full','brief','quick','data'))
  {
    amount <- match.arg(amount)
    geotype <- toupper(substr(GEO,1,3))
    mode <- 'wb'
    if (geotype == 'GDS') {
      gdsurl <- 'ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/GDS/'
      myurl <- paste(gdsurl,GEO,'.soft.gz',sep="")
      destfile <- file.path(destdir,paste(GEO,'.soft.gz',sep=""))
    }
    if (geotype == 'GSE' & amount=='full') {
      gseurl <- 'ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_series/'
      myurl <- paste(gseurl,GEO,'/',GEO,'_family.soft.gz',sep="")
      destfile <- file.path(destdir,paste(GEO,'.soft.gz',sep=""))
    }
    if (geotype == 'GSE' & amount!='full' & amount!='table') {
      gseurl <- "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
      myurl <- paste(gseurl,'?targ=self&acc=',GEO,'&form=text&view=',amount,sep='')
      destfile <- file.path(destdir,paste(GEO,'.soft',sep=""))
      mode <- 'w'
    }
    if (geotype == 'GPL') {
      gseurl <- "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
      myurl <- paste(gseurl,'?targ=self&acc=',GEO,'&form=text&view=',amount,sep='')
      destfile <- file.path(destdir,paste(GEO,'.soft',sep=""))
      mode <- 'w'
    }
    if (geotype == 'GSM') {
      gseurl <- "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
      myurl <- paste(gseurl,'?targ=self&acc=',GEO,'&form=text&view=',amount,sep='')
      destfile <- file.path(destdir,paste(GEO,'.soft',sep=""))
      mode <- 'w'
    }
    download.file(myurl,destfile,mode=mode, quiet=T)
   # writeLines('File stored at: ')
   # writeLines(destfile)
    invisible(destfile)
  }

getGEORaw <- function(GEO,destdir=tempdir()) {
  geotype <- toupper(substr(GEO,1,3))
  if(geotype=='GSE') {
    GEOurl <- 'ftp://ftp.ncbi.nih.gov/pub/geo/data/geo/raw_data/series/'
    myurl <- paste(GEOurl,GEO,'/',GEO,'_RAW.tar',sep="")
    destfile <- file.path(destdir,paste(GEO,'_RAW.tar',sep=""))
    download.file(myurl,destfile)
    # writeLines('File stored at: ')
    # writeLines(destfile)
    invisible(destfile)
  } else {
    stop('Fetching raw data supported for GSE only....')
  }
}
parseGeoData <- function(txt) {
	# Get the stuff between table_begin and table_end
	# Arguments:
	#    txt: text associated with a single entity
	# Returns:
	#    data.frame containing just the table and associated
	#    column names
		
	# Get start and end of data table
	# Thanks to GEO folks for adding table_begin and end
	tbl.start <- grep('!\\w+_table_begin',txt,perl=TRUE)+1
	tbl.end <- grep('!\\w+_table_end',txt,perl=TRUE)-1
        if(length(tbl.end)==0) {
          tbl.end <- length(txt)
        }
#       This is slower, but correctly deals with column types
        txtcon <- textConnection(txt[tbl.start:tbl.end])
        ret <- read.delim(txtcon,header=TRUE,sep="\t",na.strings="NULL")
        close(txtcon)
        return(ret)
	# Grab first line after table_begin for column names
#	tbl.colnames <- strsplit(txt[tbl.start],"\t")[[1]]
#	tbl.tmp <- do.call('rbind',strsplit(txt[(tbl.start+1):tbl.end],"\t"))
#	colnames(tbl.tmp) <- tbl.colnames
#	return(data.frame(tbl.tmp))
}

parseGeoMeta <- function(txt) {
		
	leader <- strsplit(grep('!\\w*_',txt,perl=TRUE,value=TRUE)[1],'_')[[1]][1]
     # pull out only lines that are in the header
  	tmp <- txt[grep(leader,txt)]
  	tmp <- gsub(paste(leader,'_',sep=""),'',tmp)
  	first.eq <- regexpr(' = ',tmp)
  	tmp <- cbind(substring(tmp,first=1,last=first.eq-1),substring			    (tmp,first=first.eq+3))
     # remove blank lines
     #tmp <- tmp[-grep('^\\s?$',tmp[,2],perl=TRUE),]
	header <- split(tmp[,2],tmp[,1])
	return(header)
}

parseGDSSubsets <- function(txt) {
  # takes GDS text as input
  # returns a data frame suitable for inclusion as a "Column" slot
  #   in a GDS GeoDataTable object
  numSubsets <- length(grep('^\\^subset',txt,ignore.case=TRUE))
  subset.lists <- list()
  if (numSubsets>0) {
    subset.types <-
      do.call('rbind',strsplit(txt[grep("!subset_type",txt)],' = '))[,2]
    subset.descriptions <-
      do.call('rbind',strsplit(txt[grep("!subset_description",txt)],' = '))[,2]
    subset.samples <-
      strsplit(do.call('rbind',strsplit(txt[grep("!subset_sample_id",txt)],' = '))[,2],',')
    for(i in 1:length(subset.types)) {
      for(j in 1:length(subset.samples[[i]])) {
        subset.lists[[subset.types[i]]][[subset.samples[[i]][j]]] <- subset.descriptions[i]
      }
    }
  }
  sample.descriptions <-
    do.call('rbind',strsplit(txt[grep('#GSM',txt)],' = '))
  sample.descriptions[,1] <- gsub('#','',sample.descriptions[,1])
  samp.table <- data.frame(sample=as.character(sample.descriptions[,1]))
  if(length(subset.lists)>0) {
    for(i in 1:length(subset.lists)) {
      samp.table[match(names(subset.lists[[i]]),samp.table[,1]),
                 names(subset.lists)[i]] <- as.factor(unlist(subset.lists[[i]]))
    }
    colnames(samp.table) <- c('sample',gsub(' ','.',names(subset.lists)))
  } else {
    colnames(samp.table) <- 'sample'
  }
  samp.table[match(sample.descriptions[,1],samp.table[,1]),'description'] <-
    sample.descriptions[,2]
  return(as.data.frame(samp.table))
}


splitOnFirst <- function(x,pattern) {
  patlen <- nchar(pattern)
    matches <- regexpr(pattern,x)
    leftside <- substr(x,start=1,stop=matches-1)
    rightside <- substr(x,start=matches+patlen,stop=10000000)
    return(data.frame(leftside,rightside))
  }


parseGeoColumns <- function(txt) {
	cols <- as.data.frame(splitOnFirst(txt[grep('^#',txt,perl=TRUE)],' = '))
    	cols[,1] <- sub('#','',as.character(cols[,1]))
    	colnames(cols) <- c('Column','Description')
    	return(cols)
}

parseGeoDataTable <- function(txt) {
	cols <- parseGeoColumns(txt)
	tab <- parseGeoData(txt)
	return(new('GEODataTable',columns=cols,table=tab))
}

parseGSM <- function(txt) {
  geoDataTable <- parseGeoDataTable(txt)
  meta <- parseGeoMeta(txt)
  gsm <- new('GSM',
             header=meta,
             dataTable = geoDataTable)
  return(gsm)
}

parseGPL <- function(txt) {
	geoDataTable <- parseGeoDataTable(txt)
	meta <- parseGeoMeta(txt)
	gsm <- new('GPL',
			  header=meta,
			  dataTable = geoDataTable)
		
	return(gsm)
}

parseGSE <- function(con,GSElimits) {
  gsmlist <- list()
  gpllist <- list()
  GSMcount <- 0
 # writeLines('Parsing....')
  ## This gets the header information for the GSE
  lines <- 1
  a <- vector()
  finished <- FALSE
  nextEntity <- ""
  while(!finished) {
    line <- readLines(con,1)
    if(length(line)==0) finished <- TRUE
    a[lines] <- line
    lines <- lines+1
    b <- grep('^\\^(SAMPLE|PLATFORM)',line,value=TRUE,perl=TRUE)
    if(length(b)>0) {
      nextEntity <- b
      # writeLines(b)
      finished <- TRUE
      lines <- 1
      header=parseGeoMeta(a)
    }
  }
  finished <- FALSE
  while(!finished) {
    line <- readLines(con,1)
    if(length(line)==0) {
      finished <- TRUE
    } else {
      a[lines] <- line
      lines <- lines+1
      b <- grep('^\\^(SAMPLE|PLATFORM)',line,value=TRUE,perl=TRUE)
    }
    if(length(b)>0 | finished) {
      lines <- 1
      #new SAMPLE
      if(length(grep('SAMPLE',nextEntity))>0) {
        accession <- strsplit(nextEntity,' = ')[[1]][2]
        GSMcount <- GSMcount+1
        # Look to see if limits should be used, otherwise, proceed
        if(is.null(GSElimits)) {
          gsmlist[[accession]] <- parseGSM(a)
          # writeLines(b)
        } else {
          if((GSMcount>=GSElimits[1]) &
             (GSMcount<=GSElimits[2])) {
            gsmlist[[accession]] <- parseGSM(a)
            # writeLines(b)
          } else {
            cat('Skipping sample',GSMcount,': Accession',accession,'at user request\n')
          }
        }
      }
      if(length(grep('PLATFORM',nextEntity))>0) {
        accession <- strsplit(nextEntity,' = ')[[1]][2]
        gpllist[[accession]] <- parseGPL(a)
      }
      if(!is.null(GSElimits)) {
        if(GSMcount+1>GSElimits[2]) {
                                        # end if beyond GSElimits[2]
          cat('Stopping here at user request\n')
          finished <- TRUE
        }
      }
      if(!finished) {
        nextEntity <- b
      }
      a <- vector()
    }
  }
  gse <- new("GSE",
             header=header,
             gsms = gsmlist,
             gpls = gpllist
             )
  return(gse)
}


parseGDS <- function(txt) {
  # writeLines("parsing geodata")
  tab <- parseGeoData(txt)
  # writeLines("parsing subsets")
  cols <- parseGDSSubsets(txt)
  
  # writeLines("ready to return")
  return(new('GDS',
             header=parseGeoMeta(txt),
             dataTable=new('GEODataTable',
               table=tab,
               columns=cols
               )
             ))
	}

findFirstEntity <- function(con) {
  while(TRUE) {
    line <- readLines(con,1)
    if(length(line)==0) return(0)
    entity.line <- grep('^\\^(DATASET|SAMPLE|SERIES|PLATFORM)',
                        line,ignore.case=TRUE,value=TRUE,perl=TRUE)
    if(length(entity.line)>0) {
      ret <- c(tolower(sub('\\^','',strsplit(entity.line,' = ')[[1]][1])),
               strsplit(entity.line,' = ')[[1]][2])
      return(ret)
    }
  }
}

parseGEO <- function(con,GSElimits) {
  first.entity <- findFirstEntity(con)
  ret <- switch(as.character(first.entity[1]),
                sample= {
                  txt <- readLines(con)
                  parseGSM(txt)
                },
                series= parseGSE(con,GSElimits),
                dataset= {
                  txt <- readLines(con)
                  parseGDS(txt)
                },
                platform= {
                  txt <- readLines(con)
                  parseGPL(txt)
                },
                )
  return(ret)
}
	
"readUrl" <-
function (url) 
{
    options(show.error.messages = FALSE)
    con <- try(url(url, open = "r"))
    options(show.error.messages = TRUE)
    if (inherits(con, "try-error")) {
        stop(paste("Can't connect to url", url))
    }
    temp <- readLines(con)
    close(con)
    return(temp)
}
.onLoad <- function(lib,pkg) require(methods)
