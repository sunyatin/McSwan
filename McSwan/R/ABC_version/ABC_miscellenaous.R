

#' @title Print a summary of a \code{referenceTable} object content
#' @param X a \code{referenceTable} object
#' @seealso \code{\link{generate_priors}}
#' @export
ABC_print.referenceTable <- function(X) {
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


#' @title Genomic SNP density in the observed dataset
#' @description Given a window length, compute the distribution of within-window SNP counts in the observed dataset. This function is useful to get an idea about the SNP density in the observed dataset to fix (i) an appropriate window size and (ii) the minimum number of SNPs to retain a window for downstream analyses.
#' @param obs an \code{observedDataset} object, generated by \code{\link{get_SFS}}
#' @param windowSize (integer) a window size (in base pairs)
#' @return A summary of the distribution of within-window SNP counts.
#' @seealso \code{\link{generate_priors}}
#' @examples load(system.file("data", "UNFOLDED_POLARIZED_2dSFS_CEU-LWK_20-20_chr2", package = "McSwan"))
#' summary(obs, windowSize = 10000)
#' @export
ABC_summary.observedDataset <- function(obs, windowSize) {

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



