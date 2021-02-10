message <- function (..., domain = NULL, appendLF = TRUE) {

}

# GSExxx, Series
GseToGct <- function(gse=NULL, data.column.name='VALUE',
					gct.output.filename=NULL, omit.na=FALSE) {
	t <- Table(GPLList(gse)[[1]])
	ids <- t$ID # ids in GPL
	
	desc <- t$DESCRIPTION
	if(is.null(desc)) {
		desc <- t$GENE_NAME
	} 
	if(is.null((desc))) {
		desc <- t$Name
	}
	if(is.null((desc))) {
		desc <- t$GB_ACC
	}
	
	data.matrix <- do.call("cbind", lapply(GSMList(gse), function(x) {
	tab <- Table(x)
		mymatch <- match(ids, tab$ID_REF)
		return(tab[,data.column.name][mymatch])
	}))

	if(omit.na==TRUE) {
		filter <- which(rowSums(is.na(data.matrix)) != ncol(data.matrix))
		data.matrix <- data.matrix[filter,]
		desc <- desc[filter]
		ids<- ids[filter]
	}

	row.names(data.matrix) <- ids
	gct <- list(data=data.matrix, row.descriptions=desc)
	write.gct(gct, gct.output.filename)
}



# GDSxxx, e.g. GDS1, GDS2577
GdsToGct <- function(gds=NULL, gct.output.filename,
					omit.na=FALSE) {
	eset <- GDS2eSet(gds, do.log2 = FALSE)
	f <- eset@featureData
	annotations <- row.names(varMetadata(f))
	
	row.descriptions=NULL
	
	if(!is.na(annotations['CLONE_ID'])) {
		row.descriptions <- f$CLONE_ID
	}
	
	if(is.null(row.descriptions) && !is.na(annotations['Gene Title'])) {
		row.descriptions <- f[['Gene Title']]
	}
	
	if(is.null(row.descriptions) && !is.na(annotations['GENE'])) {
		row.descriptions <- f$GENE
	}
	row.descriptions <- as.vector(row.descriptions)
	
	gct <- list(data=exprs(eset), row.descriptions=row.descriptions)
	write.gct(gct, gct.output.filename)
}
	install.required.packages <- function(libdir) {
	if(!is.package.installed(libdir, "BiocGenerics")) {
		info("installing BiocGenerics")
		install.package(libdir, "BiocGenerics_0.36.0.tar.gz")
	}        
	if(!is.package.installed(libdir, "Biobase")) {
		info("installing Biobase")
		suppressMessages(suppressWarnings(install.package(libdir, "Biobase_2.50.0.tar.gz")
		))
	}        
}

run <- function(libdir, args) {
	suppressMessages(.run(libdir, args))
}

.run <- function(libdir, args) {
	
	library(methods)
	library(tools)
	
	source(paste(libdir, "common.R", sep=''))
	
	if(libdir!='') {
		setLibPath(libdir)
		install.required.packages(libdir)
	}
	
	library(Biobase)
	
	source(paste(libdir, "GEOquery.R", sep=''))
	DEBUG <<- F
	geo.id <- NULL
	filename <- NULL
	data.column.name <- 'VALUE'
	gct.output.filename <- NULL
	
	ftp.proxy.server <- NULL
	ftp.proxy.user <- NULL
	ftp.proxy.password <- NULL
	http.proxy.server <- NULL
	http.proxy.user <- NULL
	http.proxy.password <- NULL
	
	 
	for(i in 1:length(args)) {
		flag <- substring(args[[i]], 0, 2)
		value <- trim(substring(args[[i]], 3, nchar(args[[i]])))
		if(value=='') {
			next
		}
		if(flag=='-a') {
			geo.id <- value
		} else if(flag=='-f') {
			filename <- value
		} else if(flag=='-n') {
			omit.na <- as.logical(value)
		} else if(flag=='-d') {
			data.column.name <- value
		} else if(flag=='-o') {
			gct.output.filename <- value
		} else if(flag=='-1') {
			ftp.proxy.server <- value
		} else if(flag=='-2') {
			ftp.proxy.user <- value
		} else if(flag=='-3') {
			ftp.proxy.password <- value
		} else if(flag=='-4') {
			http.proxy.server <- value
		} else if(flag=='-5') {
			http.proxy.user <- value
		} else if(flag=='-6') {
			http.proxy.password <- value
		}
	}

	if(!is.null(ftp.proxy.server)) {
		info("set FTP proxy")
		Sys.setenv("ftp_proxy", ftp.proxy.server)
	}
	if(!is.null(ftp.proxy.user)) {
		info("set FTP user")
		Sys.setenv("ftp_proxy_user", ftp.proxy.user)
	}
	if(!is.null(ftp.proxy.password)) {
		info("set FTP password")
		Sys.setenv("ftp_proxy_password", ftp.proxy.password)
	}
	if(!is.null(http.proxy.server)) {
		info("set HTTP proxy")
		Sys.setenv("http_proxy", http.proxy.server)
	}
	if(!is.null(http.proxy.user)) {
		s <- http.proxy.user
		if(!is.null(http.proxy.password)) {
			s <- paste(http.proxy.user, ":", http.proxy.password, sep='')
		} 
		info("set HTTP username, password")
		Sys.setenv("http_proxy_user", s)
	}

	if(is.null(geo.id) && is.null(filename)) {
		exit("Either a GEO accession or a GEO SOFT file must be specified.")
	}
	
	if(is.null(filename)) {
		info("download GEO id ", geo.id, "...")
		geo.query <- getGEO(geo.id)
	} else {
		info("loading GEO file ", filename, "...")
		geo.query <- getGEO(filename=filename)
	}
	info("loaded GEO")
	if (class(geo.query) == "GSE") {
		info("converting GSE...") 
		GseToGct(gse=geo.query, data.column.name=data.column.name,
			gct.output.filename=gct.output.filename, omit.na=omit.na)
	} else if (class(geo.query) == "GDS") {
		info("converting GDS...")
		GdsToGct(gds=geo.query, gct.output.filename=gct.output.filename)
	} else {
		exit("Can only retrieve data for GEO Datasets and GEO Series")
	}
}

init = commandArgs(trailingOnly = TRUE)
run(libdir = init[1], args = as.list(init))
