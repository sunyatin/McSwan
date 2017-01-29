

#' @title Initialize the working session
#' @param tempDir path/name of a temporary directory (will be recursively created if it does not exist)
#' @param pythonPath path to the \code{python} executable, or set "python" if python is callable as an environment variable
#' @param javaPath path to the \code{java} executable, or set "java" if python is callable as an environment variable
#' @details If python is not set as an environment variable, please see <https://docs.python.org/2/using/windows.html>, alternatively you can specify the python executable ("java.exe" for instance on Windows). Java is to be specified only if you simulate pseudo-observed datasets for the validation procedure, otherwise you can leave its default value.
#' @examples Please refer to the vignette.
#' @export
set_session <- function(tempDir, pythonPath = "python", javaPath = "java") {

  tempDir <<- tempDir
  if (!file.exists(tempDir)) {
	#stop(paste(tempDir,"does not exist"))
	dir.create(tempDir, recursive=TRUE)
	cat("Temporary directory ",normalizePath(tempDir)," does not exist and has been created.\n", sep="")
  } else {
    cat("Please note that the temporary directory \"",normalizePath(tempDir),"\" already exists.\n", sep="")
  }
  cat("The temporary directory has been set (",normalizePath(tempDir),")\n", sep="")
  
  pythonPath <<- pythonPath
  if ( grepl(".", pythonPath, fixed=T) && !file.exists(pythonPath) ) stop(paste(pythonPath,"does not exist"))
  #if ( !grepl(".", pythonPath, fixed=T) ) {
    #WPM <- tryCatch(system("python --help", intern=T), error=function(e) e, warning=function(w) w)
    WPM <- tryCatch(system(paste(pythonPath,"--help"), intern=T), error=function(e) e, warning=function(w) w)
    if ("simpleError" %in% attributes(WPM)$class) stop("Python does not seem to be set as an environment variable, please see https://docs.python.org/2/using/windows.html or set the path to a python executable")
  #}
  cat("Python has been correctly attached.\n")
  
  javaPath <<- javaPath
  if ( grepl(".", javaPath, fixed=T) && !file.exists(javaPath) ) stop(paste(javaPath,"does not exist"))
  #if (!file.exists(javaPath)&&javaPath!="java") stop(paste(javaPath,"does not exist"))
  
  pyPath <<- system.file("data", "multisfs_v4.py", package = "McSwan")
  #if (!file.exists(pyPath)&&pyPath!="python") stop(paste(pyPath,"does not exist"))
  
  msmsPath <<- system.file("data", "msms3.2rc-b163.jar", package = "McSwan")
  if (!file.exists(msmsPath)) stop(paste(msmsPath,"does not exist"))
  
  pyMSMSPath <<- system.file("data", "pos_sas.py", package = "McSwan")
  if (!file.exists(pyMSMSPath)) stop(paste(pyMSMSPath,"does not exist"))
  
  pyCLEANMSPath <<- system.file("data", "clean_msOutput.py", package = "McSwan")
  if (!file.exists(pyCLEANMSPath)) stop(paste(pyCLEANMSPath,"does not exist"))
  
  pyVCF2PACPath <<- system.file("data", "vcf2pac.py", package = "McSwan")
  if (!file.exists(pyVCF2PACPath)) stop(paste(pyVCF2PACPath,"does not exist"))
  
  cat("All paths have been set.\n")
}



#' @title Sample n classes
#' @keywords internal
discr <- function(n, classes) {
  if (!is.numeric(classes)) classes <- as.numeric(as.character(classes))
  return(sample(classes, n, replace=TRUE))
}


#' @title Format a manually defined prior distribution as it is.
#' @keywords internal
exact <- function(n, classes) {
  if (!is.numeric(classes)) classes <- as.numeric(as.character(classes))
  if (length(classes)!=n) stop("the number of elements in your exact array must match nSim!")
  return(classes)
}


#' @title Generates log10-uniform distribution
#' @description This function samples values out of a log10-uniform distribution. This distribution is useful to uniformly sample wide-ranged values across all the orders of magnitudes (e.g. ranges of values checking log10(upperBound) - log10(lowerBound) > 3).
#' @param n number of values to sample
#' @param a lower limit of the distribution
#' @param b upper limit of the distribution
#' @details \code{Y} follows a log10-uniform distribution <=> \code{Y=10^(x)} where \code{x} follows a uniform distribution between \code{log10(a)} and \code{log10(b)}.
#' @return A vector of \emph{n} values drawn from a log10-uniform distribution.
#' @export
rlogunif <- function(n, a, b) {
  e <- runif(n, log10(a), log10(b))
  return(10**e)
}

#' @title Calculate the mode of a univariate distribution
#' @param x a vector of numeric values
#' @return A numeric value corresponding to the mode of the interpolated empirical distribution.
#' @keywords internal
getmode <- function(x, na.rm = TRUE) {
  if (all(is.na(x))) return(NA)
  if (length(!is.na(x))==1) return(x)
  z <- density(x, na.rm = na.rm)
  return(z$x[which.max(z$y)])
}


#' @title Print a summary of a \code{referenceTable} object content
#' @param X a \code{referenceTable} object
#' @seealso \code{\link{generate_priors}}
#' @export
print.referenceTable <- function(X) {
  if (class(X)!="referenceTable") stop("X must be a valid referenceTable object")
    
    cat(paste("\n>>> GENERAL\n\n"))
    cat(paste("Demograpy:\t",X$GENERAL$msDemography,"\n\n"))
    cat(paste("Diploid effective size of the reference population (No):            ",X$GENERAL$No,"\n"))
    cat(paste("Number of simulations per model:                                    ",X$GENERAL$nSimul,"\n"))
    cat(paste("Window size (in base pairs):                                        ",X$GENERAL$windowSize,"\n"))
    cat(paste("Population sample sizes (number of gametes):                        ",paste(X$GENERAL$islandSizes, collapse=" "),"\n"))
    cat("\n")
                                                                                                                           
    cat("\n>>> PRIORS\n\n")
	cat("Selective models: ",names(X$PRIORS),"\n")
	cat("Summary of the sweep age prior distributions (in generations before present):\n")
    print(t(simplify2array(lapply(X$PRIORS, function(x) summary(x$sweepAge)))))
    cat("\n")
    
    cat(paste("\n>>> SFS (joint multidimensional site frequency spectra):\n"))
	if (is.null(X$SFS)) {
		cat("EMPTY.\n")
	} else {
		cat("Format: ",ifelse(class(X$SFS[[1]])=="dgCMatrix", "sparseMatrix", "matrix"),"\n")
		cat("Models: ",names(X$SFS),"\n")
		cat("SFS properties:",ifelse(X$GENERAL$fold, "folded (= minor alleles frequencies)", "unfolded (= derived alleles frequencies)"), "and",ifelse(X$DIMREDUC$GENERAL$relativeSFS, "relative", "absolute"),"\n")
		cat("Number of SFS bins: ",ncol(X$SFS[[1]]),"\n")
	}
    cat("\n")
    
    cat(paste("\n>>> DIMENSION REDUCTION\n"))
	if (is.null(X$DIMREDUC)) {
		cat("EMPTY.\n")
	} else {
		cat(paste("\t Removed collinear cells:",X$DIMREDUC$GENERAL$removedCollinearColumns,"\n"))
		cat(paste("\t Standardized for PLS:",X$DIMREDUC$GENERAL$PLS_normalize,"\n"))
		cat("\n")
		cat("\t_LDA_\t",ifelse(class(X$DIMREDUC$LDA$model)=="lda","MASS::lda","HiDimDA::Dlda"),"\n")
		cat(paste("\t\t |__ Number of retained components:",ncol(X$DIMREDUC$LDA$scores),"\n"))
		cat("\t_PLS_\t","pls::plsr","\n")
		cat(paste("\t\t |__ Number of retained components:",paste(simplify2array(lapply(X$DIMREDUC$PLS, function(x) ncol(x$scores))), collapse=" & "),"\n"))
    }
	cat("\n")
    
}


#' @title Print a summary of a \code{validationTable} object content
#' @param X a \code{validationTable} object
#' @seealso \code{\link{generate_pseudoobs}}
#' @export
print.validationTable <- function(X) {
  if (class(X)!="validationTable") stop("X must be a valid validationTable object")
    
    cat(paste("\n>>> GENERAL\n\n"))
    cat(paste("Demograpy:\t",X$GENERAL$msDemography,"\n\n"))
    cat(paste("Diploid effective size of the reference population (No):            ",X$GENERAL$No,"\n"))
    cat(paste("Number of simulations per model:                                    ",X$GENERAL$nSimul,"\n"))
    cat(paste("Window size (in base pairs):                                        ",X$GENERAL$windowSize,"\n"))
    cat(paste("Population sample sizes (number of gametes):                        ",paste(X$GENERAL$islandSizes, collapse=" "),"\n"))
    cat("\n")
                                                                                                                           
    cat("\n>>> PRIORS\n\n")
	cat("Selective models: ",names(X$PRIORS),"\n")
	cat("Sweep position: ",X$PRIORS[[1]]$sweepPos[1],"\n")
	cat("Summary of the sweep age distributions:\n")
	nama <- names(X$PRIORS); nama <- nama[nama!="i0"]
    print(t(simplify2array(lapply(X$PRIORS[nama], function(x) summary(x$sweepAge)))))
	cat("Summary of the recombination rate distributions:\n")
	print(t(simplify2array(lapply(X$PRIORS[nama], function(x) summary(x$recRate)))))
    cat("\n")
    
    cat(paste("\n>>> SFS (joint multidimensional site frequency spectra):\n"))
	if (is.null(X$SFS)) {
		cat("EMPTY.\n")
	} else {
		cat("Models: ",names(X$SFS),"\n")
		cat("SFS properties:",ifelse(X$GENERAL$fold, "folded (= minor alleles frequencies)", "unfolded (= derived alleles frequencies)"),"\n")
		cat("Number of SFS bins: ",length(X$SFS[[1]]),"\n")
	}
    cat("\n")
    
    cat(paste("\n>>> ANALYSIS\n"))
	if (!"ANALYSIS"%in%names(X)) {
		cat("EMPTY.\n")
	} else {
		cat("DONE!\n")
    }
	cat("\n")
    
}


#' @title Print a summary of a \code{observedDataset} object content
#' @param X a \code{observedDataset} object
#' @seealso \code{\link{get_SFS}}
#' @export
print.observedDataset <- function(X) {
  if (class(X)!="observedDataset") stop("X must be a valid observedDataset object")
    
    cat("Number of valid SNPs: ", nrow(X$obsData),"\n")
	cat("Allele type (minor/derived): ",ifelse(X$folded,"minor","derived"),"\n")
	cat("Expected number of bins in the multidimensional joint site frequency spectra: ",length(X$template),"\n")
    
}



#' @title Genomic SNP density in the observed dataset
#' @description Given a window length, compute the distribution of within-window SNP counts in the observed dataset. This function is useful to get an idea about the SNP density in the observed dataset to fix (i) an appropriate window size and (ii) the minimum number of SNPs to retain a window for downstream analyses.
#' @param obs an \code{observedDataset} object, generated by \code{\link{get_SFS}}
#' @param windowSize (integer) a window size (in base pairs)
#' @return A summary of the distribution of within-window SNP counts.
#' @seealso \code{\link{generate_priors}}
#' @examples load(system.file("data", "UNFOLDED_POLARIZED_2dSFS_CEU-LWK_20-20_chr2", package = "McSwan"))
#' summary(obs, windowSize = 10000)
#' @export
summary.observedDataset <- function(obs, windowSize) {

  L <- windowSize
  p <- obs$obsData$POS
  
  d <- c(); start <- p[1]; td <- 0
  for (i in seq_along(p)) {
    end <- start + L
    if (p[i]>=end) {
      start <- end
      d <- c(d, td)
      td <- 0
    }
    td <- td + 1
  }
  d <- c(d, td)
  
  hist(d, main="Histogram of the within-window SNP amount", xlab="Amount of SNPs within windows", ylab="Number of windows")
  cat("Summary of the within-window SNP density:\n")
  return(summary(d))
}



