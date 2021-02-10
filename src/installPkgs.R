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
}

install.required.packages <- function(libdir) {
	if(!is.package.installed(libdir, "BiocGenerics")) {
		info("installing BiocGenerics")
		install.package(libdir, "BiocGenerics_0.36.0.tar.gz")
	}
	if(!is.package.installed(libdir, "Biobase")) {
		info("installing Biobase")
		install.package(libdir, "Biobase_2.50.0.tar.gz")
	}
}

is.package.installed <- function(libdir, pkg) {
	f <- paste(libdir, pkg, sep='')
	return(file.exists(f) && file.info(f)[["isdir"]])
}

install.package <- function(dir, other) {
	f <- paste(dir, other, sep="")
	.install.unix(f)
}

.install.unix <- function(pkg) {
	if(DEBUG) {
		info("Installing package ", pkg)
	}
    lib <- .libPaths()[1]
   # cmd <- paste(file.path(R.home(), "bin", "R"), "CMD INSTALL --with-package-versions")
	cmd <- paste(file.path(R.home(), "bin", "R"), "CMD INSTALL")
    cmd <- paste(cmd, "-l", lib)
    cmd <- paste(cmd, " '", pkg, "'", sep = "")
    status <- system(cmd)
    if (status != 0)
    	cat("\tpackage installation failed\n")
}

  init = commandArgs(trailingOnly = TRUE)
run(libdir = init[1] , args = as.list(init))