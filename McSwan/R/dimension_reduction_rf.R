

#' @title LDA & RF supervised machine learning of the selection signals
#' @description Performs supervised learning algorithms (LDA and random forests) on the summary statistics (joint multidimensional allele frequency spectra) prior to genome scan. This function is to be run after \code{\link{coalesce}}.
#' @param x a \code{referenceTable} object with non-empty \code{PRIORS} and \code{SFS} slots

#' @param removeCollinearCols (logical) whether to remove collinear columns from the site frequency spectra prior to LDA & PLS dimension reduction (TRUE is recommended)

#' @param LDA_xPC (float between 0 and 1) number of principal components to retain as a fraction of the original number of cells in the SFS; LDA will be trained on the PCA-projected simulations
#' @param LDA_sizePerModel (integer) number of simulations per model to consider for training the LDA; \emph{NULL} (by default) forces to use all data

#' @param PLS_normalize (logical) whether to normalize the spectra prior to PLS analysis (TRUE is recommended)
#' @param PLS_ncomp (integer) if NULL, automatically select the fittest number of PLS components, else an integer to manually set the number of retained components
#' @param PLS_maxncomp (integer) if automatic selection of PLS features, the maximum number of PLS components to consider before downsampling to the fittest set
#' @param PLS_ncomp_method (string) the method to automatically select the best number of features: either "\emph{elbow}", "\emph{onesigma}" or "\emph{randomization}" (cf. \code{selectNcomp} in package \code{pls})
#' @param plot_PLScv (logical) if TRUE, plots the PLS cross-validation graphs
#' @param RF_ntree (integer) number of trees to grow for the random forest (has a quadratic effect on computation time)
#' @param RF_pcomp (numeric, between 0 and 1) proportion of features to retain for the random forest modelling
#' @param ... other arguments passed to \code{randomForest}
#' @return An object of class \code{referenceTable} with filled-in \code{DIMREDUC} slot.
#' @seealso \code{\link{coalesce}}, \code{\link{plsr}}, \code{\link{lda}}
#' @examples Please refer to the vignette.
#' @import Matrix
#' @import randomForest
#' @export
dim_reduction_rf <- function(x,
						# REFTB PRE-PROCESSING
                          removeCollinearCols = TRUE,
						  compress_SFS = TRUE,
						#LDA
						  LDA_xPC = .1,
						  LDA_sizePerModel = NULL,
						# dating method
						dating_method = "rf",
						##> "pls" == PLS
                          PLS_normalize = TRUE,
                          PLS_ncomp = NULL,
						  PLS_ncomp_method = "elbow",
                          PLS_maxncomp = 100,
                          PLS_plotCV = TRUE,
						##> "rf" == random forests
						  RF_ntree = 500,
						  RF_pcomp = 1, # proportion of retained features (with best importance score)
						  ...) {
						  
  warning("dim_reduction_rf is in development")
  
  if (class(x)!="referenceTable"||is.null(x$SFS)) stop("x is not a referenceTable object or its SFS slot is empty")
  
  # internally set options
  relativeSFS = TRUE # TRUE is now MANDATORY for McSwan 1.1 BECAUSE VARYING WINDOW SIZES WOULD INTRODUCE BIASES IF NOT SCALING SFS BY TOTAL SNP COUNT
  QR_multicollinearity = TRUE # whether multicollinearity should be detected by QR decomposition or not (TRUE is recommended)
  doPLSrecRate = FALSE # whether to regress (in addition to the sweep ages) onto the recombination parameter in the PLS analyses (since the recombination rate is considered a nuisance parameter, it can be sometimes useful to discard it from the PLS analysis)
  genetic_sumstats = FALSE # whether to calculate additional summary statistics. NOTA BENE: this is still in infancy
  
  # McSwan 1.1 :
  # following options are now obsolete, have no impact (only valid for previous version McSwan v1.0), but for record:
  LDA_minVariance = 0
  LDA_minVarPerModel = FALSE
  LDA_NonNullP = .05
  LDA_tol = 1e-9
  LDA_fastDiag = FALSE
  GOF_n = 2000 # because no goodness of fit tests anymore (only relevant if ABC)
  GOF_abcTol = .01
  
  # fun
  
  x$DIMREDUC <- list()
  x$DIMREDUC$GENERAL <- list(relativeSFS = relativeSFS, 
                             removedCollinearColumns = removeCollinearCols,
                             PLS_normalize = PLS_normalize)

  if ( class(x$SFS[[1]]) == "dgCMatrix" ) {
    x$SFS <- lapply(x$SFS, as.matrix)
  }
  
  # relativize multiSFS, if requested
  if (relativeSFS) {
    cat("\nRELATIVE SFS. The joint site frequency spectra will be scaled by their SNP content.\n")
    sfs <- lapply(x$SFS, function(z) {
      t(apply(z, 1, function(y) {
        sm <- sum(y, na.rm = TRUE)
        if (sm==0) return(y)
        return(y/sm)
      }))
    })
  } else {
    cat("\nABSOLUTE SFS. The joint site frequency spectra will _not_ be scaled by their SNP content.\n")
    sfs <- x$SFS
  }

  
  ###############################
  ### ROBUST LDA (by PCA-LDA) ###
  ###############################
  
  cat("\n\n>>> Robust PCA-LDA dimension reduction.\n\n") 

    cat("\n\nMatrix preparation.\n")
    lsfs <- linearize(sfs, LDA_sizePerModel = LDA_sizePerModel)
	modelIndices <- lsfs$model_indices; lsfs <- lsfs$sfs
    invisible(gc(F))
  
    nPCs = as.integer(LDA_xPC * ncol(lsfs))

    PCA.means = colMeans(lsfs)
    lsfs = sweep(lsfs, 2, PCA.means)

    PCA.euCols = apply(lsfs, 2, function(zk) length(unique(zk))>1)
    lsfs = lsfs[,PCA.euCols,drop=F]

	cat("PCA computation.\n")
	# for PCA, only centered var (not scaled because SFS has same unit across vars)
    PCA.model = prcomp(lsfs, center = FALSE, scale = FALSE, retx = FALSE)$rotation[,1:nPCs,drop=F]

    lsfs = lsfs %*% PCA.model
	
	cat("LDA computation.\n")
	LDA = lda(lsfs, modelIndices)
	
	x$DIMREDUC$LDA <- list("PCA.model" = PCA.model,
							"PCA.means" = PCA.means,
							"euCols" = PCA.euCols,
							"model" = LDA,
							"modelIndices" = factor(modelIndices, levels=names(reftb$SFS)))
							
	rm(list=c("PCA.means", "lsfs", "LDA", "PCA.euCols", "PCA.model")); invisible(gc(F))
  
  
  if (dating_method == "pls") {
  
      ###########
      ### PLS ###
      ###########
	  
	  cat("\n>>> PLS\n")
	  
	  colz <- c("sweepAge", "recRate")
	  if (doPLSrecRate==FALSE) colz <- c("sweepAge")

	  PLS <- lapply(as.list(names(x$PRIORS)), function(z) {
		get_pls(param = x$PRIORS[[z]][,colz,drop=FALSE],
				ss = sfs[[z]],
				PLS_normalize = PLS_normalize,
				removeCollinearCols = removeCollinearCols,
				PLS_ncomp = PLS_ncomp,
				PLS_ncomp_method = PLS_ncomp_method,
				plot_PLScv = PLS_plotCV,
				maxPLSncomp = PLS_maxncomp,
				deme = z)
	  })
	  names(PLS) <- names(x$PRIORS)
	  
	  # return all
	  x$DIMREDUC$PLS <- PLS
	  rm("PLS"); invisible(gc(F))
	  
	} else {
	
      ###########
      ###  RF ###
      ###########
	  
	  cat("\n>>> random Forests\n")
	  
	  colz <- c("sweepAge", "recRate")
	  if (doPLSrecRate==FALSE) colz <- c("sweepAge")
	  
	  PLS <- lapply(as.list(names(x$PRIORS)), function(z) {
	    # remove non variable SFS categories (intra pop)
		sds <- apply(sfs[[z]], 2, sd, na.rm=T)
		euCols <- sds != 0
		cat(paste("\t - Removed",sum(!euCols),"columns with zero variance.\n"))
		
		if (RF_pcomp != 1) {
		stop("in development, for now, RF_pcomp must be set to 1")
		# fast RF (with approximations)
			cat("Computing feature importance\n")
			print(table(euCols))
			ori_name <- names(euCols)
			rF_1 <- randomForest::randomForest(x = sfs[[z]][,euCols,drop=FALSE], 
									   y = as.vector(x$PRIORS[[z]][,colz]),
									   ntree = 100,
									   importance = TRUE,
									   proximity = FALSE)
			rF_1 <- rF_1$importance[,1]
			rF_1 <- order(rF_1, decreasing=TRUE)
			rF_1 <- rF_1[1:floor(RF_pcomp*ncol(sfs[[z]]))]
			names(euCols) <- paste0("i",seq_along(euCols))
			rF_1_name <- names(euCols)[euCols]
			rF_1_name <- rF_1_name[rF_1]
			euCols[!names(euCols) %in% rF_1_name] <- FALSE
			names(euCols) <- ori_name
			print(table(euCols))
			## random Forest
			rF <- randomForest::randomForest(x = sfs[[z]][,euCols,drop=FALSE], 
									   y = as.vector(x$PRIORS[[z]][,colz]),
									   ntree = RF_ntree,
									   importance = FALSE,
									   proximity = FALSE,
									   do.trace = TRUE,
									   ...)
		} else {
		# most reliable RF modelling
			rF <- randomForest::randomForest(x = sfs[[z]][,euCols,drop=FALSE], 
									   y = as.vector(x$PRIORS[[z]][,colz]),
									   ntree = RF_ntree,
									   importance = FALSE,
									   proximity = FALSE,
									   do.trace = TRUE,
									   ...)
		}
		
        invisible(gc(FALSE))
        return(list("model" = rF, "euCols" = euCols, "mean" = rep(0, ncol(sfs[[z]])), "sd" = rep(1, ncol(sfs[[z]]))))
	  })
	  names(PLS) <- names(x$PRIORS)
	  
	  # return all
	  x$DIMREDUC$PLS <- PLS
	  rm("PLS"); invisible(gc(FALSE))
	
	}
  
  if (compress_SFS) {
    cat("The SFS slot was compressed into a sparseMatrix object. You can access it like a regular matrix.\n")
    x$SFS <- lapply(x$SFS, function(z) Matrix::Matrix(z, sparse = TRUE))
    invisible(gc(F))
  }
  
  # 4 test
  if (genetic_sumstats) {
	cat("Additional genetic summary statistics to be computed.\n")
    x$DIMREDUC$genetic_sumstats <- TRUE
	TAJ <- tajima(sfs$i0, internal = TRUE)
	for (i in seq_along(x$DIMREDUC$PLS)) {
		nama <- names(x$DIMREDUC$PLS)[i]
		taj <- tajima(sfs[[nama]], internal = TRUE, do_H = !x$GENERAL$folded)
		x$DIMREDUC$PLS[[i]]$scores <- cbind(x$DIMREDUC$PLS[[i]]$scores, taj)
		TAJ <- rbind(TAJ, taj)
	}
	x$DIMREDUC$LDA$scores <- cbind(x$DIMREDUC$LDA$scores, TAJ)
  } else {
    #x$DIMREDUC$genetic_sumstats <- FALSE
  }
  
  x$GENERAL$call.dim_reduction <- match.call()
  return(x)
}




