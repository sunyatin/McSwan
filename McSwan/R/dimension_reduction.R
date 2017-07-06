

#' @title Merge model-specific multiSFSs into one big matrix
#' @param sfs a list of model-specific multiSFS matrices
#' @return A list of two elements: a vector of model indices; a big matrix containing all merged multiSFSs.
#' @keywords internal
linearize <- function(sfs, LDA_sizePerModel = NULL) {
  model_indices <- gsfs <- c()
  if (is.null(LDA_sizePerModel)) {
	for (i in seq_along(sfs)) {
		model_indices <- c(model_indices, rep(names(sfs)[i], nrow(sfs[[i]])))
		#gsfs <- rbind(gsfs, sfs[[i]])
	}
	gsfs = do.call("rbind", sfs)
  } else {
	for (i in seq_along(sfs)) {
		nr <- nrow(sfs[[i]])
		if (nr < LDA_sizePerModel) stop("LDA_sizePerModel is greater than the number of simulations in your model(s)")
		if (nr == LDA_sizePerModel) {
			model_indices <- c(model_indices, rep(names(sfs)[i], nrow(sfs[[i]])))
			gsfs <- rbind(gsfs, sfs[[i]])
		} else {
			model_indices <- c(model_indices, rep(names(sfs)[i], LDA_sizePerModel))
			gsfs <- rbind(gsfs, sfs[[i]][sample(1:nrow(sfs[[i]]), LDA_sizePerModel, replace=FALSE),])
		}
	}
  }
  return(list(model_indices=model_indices, sfs=gsfs))
}



#' @title Get model-specific PLS components
#' @param PLS_ncomp if NULL, automatically select the fittest number of PLS components, else an integer to set the number of components to a specific value
#' @keywords internal
#' @import pls
#' @import MASS
get_pls <- function(param, ss, PLS_normalize, removeCollinearCols, PLS_ncomp, PLS_ncomp_method = "elbow", plot_PLScv = FALSE,
                    maxPLSncomp = 100, deme) {
	# as of Feb 19th 2016, 3 possible methods to select the number of components: "elbow" / "onesigma" / "randomization"
  
  #if (!is.null(dim(param)) && ncol(param)>1) stop("PLS is computed against a single parameter andyou provided more than one parameter")
  param <- as.matrix(param)
  cat("\n\n\t Model-specific PLS.\n")

  # normalize
  means <- colMeans(ss)
  sds <- apply(ss, 2, sd, na.rm=T)
  euCols <- sds != 0
  cat(paste("\t - Removed",sum(!euCols),"columns with zero variance.\n"))
  if (PLS_normalize) {
    cat("\t - Matrix standardization.\n")
    ss <- t(apply(ss, 1, function(z) (z-means)/sds))
  }
  
  if (removeCollinearCols) {
    
    # perform QR decomposition
    colnames <- names(euCols)[euCols]
    QR <- qr(x = ss[,euCols,drop=F], tol = 1e-7, LAPACK = FALSE)
    rank <- QR$rank
    cat(paste("\t - Removed",sum(euCols)-rank,"collinear columns, based on QR decomposition.\n"))
    euQRcols <- QR$pivot[1:rank]
    euQRcols <- colnames[euQRcols]
    euCols[!names(euCols)%in%euQRcols] <- FALSE
    
  } #else {
  
    # detect collinear columns, if requested
    #if (removeCollinearCols) {
    #  cat("\t - Removing collinear columns using caret::findCorrelation algorithm.\n")
    #  nonCollin <- rep(TRUE, ncol(ss)); names(nonCollin) <- colnames(ss)
    #  tmp <- caret::findCorrelation(cor(ss[,euCols]), cutoff=.90, names=TRUE)
    #  nonCollin[names(nonCollin)%in%tmp] <- FALSE
    #  euCols <- euCols & nonCollin
    #  rm("tmp"); invisible(gc(F))
    #}
    
  #}
  
  ss <- ss[,euCols,drop=F]

  ##############################
  ##############################  
  ##### K SELECTION
  
  if (is.null(PLS_ncomp)) {
    cat("\t - Automatic selection of the fittest number of PLS components.\n")
    
    if ( maxPLSncomp >= ncol(ss)-1 ) maxPLSncomp <- ncol(ss)-1
    plsa <- pls::plsr(param ~ ss, scale = FALSE, validation = "CV", ncomp = maxPLSncomp)
	
	if (PLS_ncomp_method == "elbow") {
    
		r <- unlist(pls::RMSEP(plsa, intercept=F, estimate="CV", se=F)$val)[1,,]
		
		if (which.min(r)<=1) {
		  PLS_ncomp <- 1
		} else if (which.min(r)==2) {
		  PLS_ncomp <- 2
		} else {
		
		  # line1 is the line crossing the L-shaped curve
		  slope <- (min(r) - r[1])/(which.min(r) - 1)
		  intercept <- r[1] - slope*1
		  lin <- function(z) slope * z + intercept
		  
		  # find slope & intercept of line2, orthogonal to line1 (itself defined by its slope and intercept),
		  # and passing through a given point (defined by xcross and ycross)
		  find_perp <- function(xcross, ycross, slope, intercept) {
			pslope <- -1 / slope
			pintercept <- ycross - pslope * xcross
			return(c("slope" = pslope, "intercept" = pintercept))
		  }
		  
		  # K is selected as the one maximizing the distance between a point on the L-shaped curve and
		  # its projection onto line1
		  d <- c()
		  for (x in 1:which.min(r)) {
			p <- find_perp(x, r[x], slope, intercept)
			xcross <- 1 / slope * ( p[1] * x + p[2] - intercept)
			ycross <- p[1] * xcross + p[2]
			d <- c(d, sqrt( (x-xcross)**2 + (lin(x)-ycross)**2 ))
		  }
		  
		  # plots
		  if (plot_PLScv) {
			par(mfrow=c(1,2))
			  plot(r, xlim=c(0,which.min(r)), type="l", ylim=c(min(r),max(r[1:which.min(r)])), main=paste("Deme:",deme))
			  lines(1:which.min(r), lin(1:which.min(r)), col="red")
			  plot(d, type="b", ylab="Elbowity")
			par(mfrow=c(1,1))
		  }

		  ##### end:  K SELECTION
		  ##############################
		  ##############################  

		  PLS_ncomp <- which.max(d)

		}
	} else if (PLS_ncomp_method == "onesigma") {
		PLS_ncomp = selectNcomp(plsa, method = "onesigma", plot = plot_PLScv)
	} else {
		PLS_ncomp = selectNcomp(plsa, method = "randomization", alpha = .05, plot = plot_PLScv)
	}
	
	rm("plsa"); invisible(gc(FALSE))
	
    cat(paste("\t - Number of retained PLS components:",PLS_ncomp,"(based on the ",PLS_ncomp_method," method).\n"))
  }

  # PLS
  cat("\t - PLS computation.\n")
  pls <- pls::plsr(param ~ ss, scale = FALSE, validation = "none", ncomp = PLS_ncomp, model = F, x = F, y = F)
  #pls = plsRglm::plsRglm(dataY = c(param), dataX = ss, scaleX = F, scaleY = F, nt = PLS_ncomp, modele = "pls-glm-logistic")
  
  #pls.scores <- pls$scores
  #class(pls.scores) = NULL; rownames(pls.scores) = NULL

  # remove unused
  pls$scores = NULL
  pls$Yscores = NULL
  pls$fitted.values = NULL
  pls$residuals = NULL
  
  # output
  PLS <- list("model" = pls,
              "euCols" = euCols,
              #"scores" = pls.scores, 
              "mean" = means, 
              "sd" = sds)
  rm(list=c("pls", "means", "sds", "param")); invisible(gc(FALSE))
  
  cat("\n\n")
  return(PLS)    
}


#' @title Supervised machine learning of the selection signals
#' @description Performs supervised learning algorithms (LDA and PLSR) on the summary statistics (joint multidimensional allele frequency spectra) prior to genome scan. This function is to be run after \code{\link{coalesce}}.
#' @param x a \code{referenceTable} object with non-empty \code{PRIORS} and \code{SFS} slots

#' @param removeCollinearCols (logical) whether to remove collinear columns from the site frequency spectra prior to LDA & PLS dimension reduction (TRUE is recommended)

#' @param LDA_xPC (float between 0 and 1) number of principal components to retain as a fraction of the original number of cells in the SFS; LDA will be trained on the PCA-projected simulations
#' @param LDA_sizePerModel (integer) number of simulations per model to consider for training the LDA; \emph{NULL} (by default) forces to use all data

#' @param PLS_normalize (logical) whether to normalize the spectra prior to PLS analysis (TRUE is recommended)
#' @param PLS_ncomp (integer) if NULL, automatically select the fittest number of PLS components, else an integer to manually set the number of retained components
#' @param PLS_maxncomp (integer) if automatic selection of PLS features, the maximum number of PLS components to consider before downsampling to the fittest set
#' @param PLS_ncomp_method (string) the method to automatically select the best number of features: either "\emph{elbow}", "\emph{onesigma}" or "\emph{randomization}" (cf. \code{selectNcomp} in package \code{pls})
#' @param plot_PLScv (logical) if TRUE, plots the PLS cross-validation graphs

#' @return An object of class \code{referenceTable} with filled-in \code{DIMREDUC} slot.
#' @seealso \code{\link{coalesce}}, \code{\link{plsr}}, \code{\link{lda}}
#' @examples Please refer to the vignette.
#' @import Matrix
#' @export
dim_reduction <- function(x,
						# REFTB PRE-PROCESSING
                          removeCollinearCols = TRUE,
						  compress_SFS = TRUE,
						#LDA
						  LDA_xPC = .1,
						  LDA_sizePerModel = NULL,
						#PLS
                          PLS_normalize = TRUE,
                          PLS_ncomp = NULL,
						  PLS_ncomp_method = "elbow",
                          PLS_maxncomp = 100,
                          PLS_plotCV = TRUE) {
  
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




