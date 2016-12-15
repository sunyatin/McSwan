

#' @title Set mandatory paths to external softwares
#' @param tempDir path to a temporary directory
#' @param pyPath path to your \code{python} executable, or the OS environment variable
#' @param javaPath path to your \code{java} executable, or the OS environment variable
#' @param msmsPath path to the \code{msms} executable, or the OS environment variable
#' @export
set_session <- function(tempDir, pythonPath="python", javaPath) {

  tempDir <<- tempDir
  if (!file.exists(tempDir)) {
	#stop(paste(tempDir,"does not exist"))
	dir.create(tempDir)
	cat(paste("Temporary directory ",normalizePath(tempDir)," does not exist and has been created.\n",sep=""))
  } else {
    cat("Please note that temporary directory ",normalizePath(tempDir)," already exists.\n")
  }
  cat("Temporary directory has been set as: ",normalizePath(tempDir),"\n")
  
  pythonPath <<- pythonPath
  if ( grepl(".", pythonPath, fixed=T) && !file.exists(pythonPath) ) stop(paste(pythonPath,"does not exist"))
  #if ( !grepl(".", pythonPath, fixed=T) ) {
    #WPM <- tryCatch(system("python --help", intern=T), error=function(e) e, warning=function(w) w)
    WPM <- tryCatch(system(paste(pythonPath,"--help"), intern=T), error=function(e) e, warning=function(w) w)
    if ("simpleError" %in% attributes(WPM)$class) stop("Python does not seem to be set as an environment variable, please see https://docs.python.org/2/using/windows.html or set the path to a python executable")
  #}
  cat("Python has been correctly attached.\n")
  
  javaPath <<- javaPath
  #if (!file.exists(javaPath)&&javaPath!="java") stop(paste(javaPath,"does not exist"))
  
  pyPath <<- system.file("data", "multisfs_v4.py", package = "McSwan")
  #if (!file.exists(pyPath)&&pyPath!="python") stop(paste(pyPath,"does not exist"))
  
  msmsPath <<- system.file("data", "msms3.2rc-b163.jar", package = "McSwan")
  if (!file.exists(msmsPath)) stop(paste(msmsPath,"does not exist"))
  
  pyMSMSPath <<- system.file("data", "pos_sas.py", package = "McSwan")
  if (!file.exists(pyMSMSPath)) stop(paste(pyMSMSPath,"does not exist"))
  
  pyCLEANMSPath <<- system.file("data", "clean_msOutput.py", package = "McSwan")
  if (!file.exists(pyCLEANMSPath)) stop(paste(pyCLEANMSPath,"does not exist"))
  
  cat("All paths have been set.\n")
}



#' @title Sampling
discr <- function(n, classes) {
  if (!is.numeric(classes)) classes <- as.numeric(as.character(classes))
  return(sample(classes, n, replace=TRUE))
}


#' @title Exact prior distrib
exact <- function(n, classes) {
  if (!is.numeric(classes)) classes <- as.numeric(as.character(classes))
  if (length(classes)!=n) stop("the number of elements in your exact array must match nSim!")
  return(classes)
}


#' @title Random log-uniform distribution
#' @export
rlogunif <- function(n, a, b) {
  e <- runif(n, log10(a), log10(b))
  return(10**e)
}

#' @title Get mode
#' @export
getmode <- function(x, na.rm = TRUE) {
  if (all(is.na(x))) return(NA)
  if (length(!is.na(x))==1) return(x)
  z <- density(x, na.rm = na.rm)
  return(z$x[which.max(z$y)])
}

#' @title Get confusion matrix
#' @export
confusion_matrix <- function(cvd) {
  # first col: true ; second col: predicted
  # cvd is a 2-col matrix of strings representing demes under sel
  
  cvd <- apply(cvd, 2, as.character)
  
  var <- unique(c(cvd))
  cm <- matrix(0, nrow=length(var), ncol=length(var))
  rownames(cm) <- colnames(cm) <- var
  cm <- as.data.frame(cm)
  
  for (i in 1:nrow(cvd)) {
    cm[cvd[i,1],cvd[i,2]] <- cm[cvd[i,1],cvd[i,2]] + 1
  }
  
  return(cm)
}


#' @title explain a \code{referenceTable} object
#' @param X a \code{referenceTable} object
#' @return A summary of the content of the reference table.
#' @seealso \code{\link{generate_priors}}
#' @export
explain <- function(X) {
  if (class(X)=="referenceTable") {
    
    cat(paste("\nGENERAL\n"))
    cat(paste("\t ms-formatted demography:",X$GENERAL$msDemography,"\n"))
    cat(paste("\t Reference No:",X$GENERAL$No,"\n"))
    cat(paste("\t Number of simulations per model:",X$GENERAL$nSimul,"\n"))
    cat(paste("\t Windows size:",X$GENERAL$windowSize,"\n"))
    cat(paste("\t Island sizes:",paste(X$GENERAL$islandSizes, collapse=" "),"\n"))
    cat("\n")
    
    cat(paste("\nPRIORS\n"))
    cat(str(X$PRIORS))
    cat("\n")
    
    cat(paste("\nSFS\n"))
    cat(str(X$SFS))
    cat("\n")
    
    cat(paste("\nDIMENSION REDUCTION TECHNIQUES\n"))
    cat(paste("\t Relative SFS:",X$DIMREDUC$GENERAL$relativeSFS,"\n"))
    cat(paste("\t Removed collinear cells of the SFS:",X$DIMREDUC$GENERAL$removedCollinearColumns,"\n"))
    cat(paste("\t Normalizing for PLS:",X$DIMREDUC$GENERAL$PLS_normalize,"\n"))
    cat("\n\n")
    cat(paste("\t LDA:","LDA"%in%names(X$DIMREDUC),"\n"))
    cat(paste("\t\t |__ Number of retained components:",ncol(X$DIMREDUC$LDA$scores),"\n"))
    cat(paste("\t PLS:","PLS"%in%names(X$DIMREDUC),"\n"))
    cat(paste("\t\t |__ Number of retained components:",ncol(X$DIMREDUC$PLS[[1]]$scores),"\n"))
    cat("\n")
    
  }
}


# statistics on observed dataset
#' @title Statistics on observed dataset (e.g. SNP density)
#' @export
obs_stats <- function(obs, reftb) {
  L <- reftb$GENERAL$windowSize
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
  
  hist(d, main="Histogram of the SNP density along the genomic region", xlab="SNP counts within windows")
  return(summary(d))
}



