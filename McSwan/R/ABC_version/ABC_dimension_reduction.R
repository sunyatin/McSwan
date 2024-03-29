

#' @title Get the null distribution of distances within a model
#' @param x a model-specific multiSFS matrix
#' @param n_dist_bstp_per_model number of distances to calculate per model to build the null distribution
#' @return A vector of distances, giving an empirical null distribution.
#' @keywords internal
ABC_get_dists <- function(x, n_dist_bstp_per_model, abcTolerance) {
  targets <- sample(1:nrow(x), n_dist_bstp_per_model, replace=T)
  d <- rep(NA, n_dist_bstp_per_model)
  for (i in 1:length(targets)) {
    tp <- suppressWarnings(abc::abc(x[targets[i],], 
                               rep(1, nrow(x)-1), 
                               x[-c(targets[i]),], 
                               tol = abcTolerance, 
                               corr = T,
                               method = "rejection"))
    d[i] <- mean(tp$dist)
  }
  return(d)
}



#' @title Get model-specific PLS components
#' @param PLS_ncomp if NULL, automatically select the fittest number of PLS components, else an integer to set the number of components to a specific value
#' @keywords internal
#' @import pls
#' @import MASS
ABC_get_pls <- function(param, ss, PLS_normalize, removeCollinearCols, PLS_ncomp, plot_PLScv,
                    maxPLSncomp = 100, deme, QR_multicollinearity = T) {
  
  #if (!is.null(dim(param)) && ncol(param)>1) stop("PLS is computed against a single parameter andyou provided more than one parameter")
  param <- as.matrix(param)
  cat("\n\n\t Model-specific PLS.\n")

  # normalize
  means <- apply(ss, 2, mean, na.rm=T)
  sds <- apply(ss, 2, sd, na.rm=T)
  euCols <- sds != 0
  cat(paste("\t - Removed",sum(!euCols),"columns with zero variance.\n"))
  if (PLS_normalize) {
    cat("\t - Matrix standardization.\n")
    ss <- t(apply(ss, 1, function(z) (z-means)/sds))
  }
  
  if (QR_multicollinearity) {
    
    # perform QR decomposition
    colnames <- names(euCols)[euCols]
    QR <- qr(x = ss[,euCols], tol = 1e-7, LAPACK = FALSE)
    rank <- QR$rank
    cat(paste("\t - Removed",sum(euCols)-rank,"collinear columns, based on QR decomposition.\n"))
    euQRcols <- QR$pivot[1:rank]
    euQRcols <- colnames[euQRcols]
    euCols[!names(euCols)%in%euQRcols] <- FALSE
    
  } else {
  
    # detect collinear columns, if requested
    if (removeCollinearCols) {
      cat("\t - Removing collinear columns using caret::findCorrelation algorithm.\n")
      nonCollin <- rep(TRUE, ncol(ss)); names(nonCollin) <- colnames(ss)
      tmp <- caret::findCorrelation(cor(ss[,euCols]), cutoff=.90, names=TRUE)
      nonCollin[names(nonCollin)%in%tmp] <- FALSE
      euCols <- euCols & nonCollin
      rm("tmp"); invisible(gc(F))
    }
    
  }
  
  ss <- ss[,euCols]

  ##############################
  ##############################  
  ##### K SELECTION
  
  if (is.null(PLS_ncomp)) {
    cat("\t - Automatic selection of the fittest number of PLS components.\n")
    
    if ( maxPLSncomp >= ncol(ss)-1 ) maxPLSncomp <- ncol(ss)-1
    plsa <- pls::plsr(param ~ ss, scale = FALSE, validation = "CV", ncomp = maxPLSncomp)
    
    r <- unlist(pls::RMSEP(plsa, intercept=F, estimate="CV", se=F)$val)[1,,]
    
    rm("plsa"); invisible(gc(FALSE))
    
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
    cat(paste("\t - Number of retained PLS components:",PLS_ncomp,"\n"))
  }

  # PLS
  cat("\t - PLS computation.\n")
  pls <- pls::plsr(param ~ ss, scale = FALSE, validation = "none", ncomp = PLS_ncomp, model = F, x = F, y = F)
  #pls = plsRglm::plsRglm(dataY = c(param), dataX = ss, scaleX = F, scaleY = F, nt = PLS_ncomp, modele = "pls-glm-logistic")
  pls$scores <- NULL
  pls.scores <- predict(pls, ss, type="scores")
  
  # output
  PLS <- list("model" = pls,
              "euCols" = euCols,
              "scores" = pls.scores, 
              "mean" = means, 
              "sd" = sds)
  rm(list=c("pls", "pls.scores", "means", "sds", "param")); invisible(gc(FALSE))
  
  cat("\n\n")
  return(PLS)    
}


#' @title Dimension reduction prior to approximate Bayesian computation
#' @description Reduces the dimension of the summary statistics (joint multidimensional allele frequency spectra) prior to ABC model selection (through linear discriminant analysis) and ABC parameter estimation (through model-specific partial least square regressions). This function is to be run after \code{\link{coalesce}}.
#' @param x a \code{referenceTable} object with non-empty \code{PRIORS} and \code{SFS} slots

#' @param removeCollinearCols (logical) whether to remove collinear columns of the site frequency spectra prior to LDA & PLS dimension reduction (TRUE is recommended)
#' @param relativeSFS (logical) whether to rescale the site frequency spectra by their SNP content (see Details) 

#' @param LDA_minVariance (numeric) minimal variance across models' simulations for a given SFS bin (if the variance is inferior or equal to this value, the SFS bin will be discarded from LDA analysis)
#' @param LDA_tol (numeric) a tolerance to decide if a matrix is singular, see \code{?lda}

#' @param PLS_normalize (logical) whether to normalize the spectra prior to PLS analysis (TRUE is recommended)
#' @param PLS_ncomp (integer) if NULL, automatically select the fittest number of PLS components, else an integer to manually set the number of retained components
#' @param PLS_maxncomp (integer) maximum number of PLS components to calculate before selecting the best set of features
#' @param plot_PLScv (logical) if TRUE, will render the cross-validation graphs of the PLS

#' @param GOF_n (integer) estimate size of the null distribution for later window goodness-of-fit tests
#' @param GOF_abcTol (numeric between 0 and 1) ABC tolerance ratio for the calculation of null distance distributions, ideally should equal the ABC tolerance for model selection (see \code{\link{abc}})

#' @return An object of class \code{referenceTable} with filled-in \code{DIMREDUC} slot.
#' @seealso \code{\emph{coalesce}}, \code{\emph{plsr}}, \code{\emph{lda}}, \code{\emph{abc}}
#' @examples Please refer to the vignette.
#' @import Matrix
#' @export
ABC_dim_reduction <- function(x,
						# REFTB PRE-PROCESSING
                          relativeSFS,
						  genetic_sumstats = TRUE,
                          removeCollinearCols = TRUE,
						  compress_SFS = TRUE,
						#LDA
                          LDA_minVariance = 0,
						  LDA_minVarPerModel = FALSE,
						  LDA_NonNullP = .05,
						  LDA_tol = 1e-9,
						  LDA_fastDiag = TRUE,
						  LDA_sizePerModel = 10000,
						#PLS
                          PLS_normalize = TRUE,
                          PLS_ncomp = NULL,
                          PLS_maxncomp = 100,
                          PLS_plotCV = FALSE,
						# NULL DISTRIB OF DISTANCES
                          GOF_n = 2000,
                          GOF_abcTol = .01 ) {
  
  if (class(x)!="referenceTable"||is.null(x$SFS)) stop("x is not a referenceTable object or its SFS slot is empty")
  
  # internally set options
  QR_multicollinearity = TRUE # whether multicollinearity should be detected by QR decomposition or not (TRUE is recommended)
  doPLSrecRate = FALSE # whether to regress (in addition to the sweep ages) onto the recombination parameter in the PLS analyses (since the recombination rate is considered a nuisance parameter, it can be sometimes useful to discard it from the PLS analysis)
  
  
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

  
  ###########
  ### LDA ###
  ###########
  
  cat("\n\n>>> LDA\n\nMatrix preparation.\n")
  lsfs <- linearize(sfs, LDA_sizePerModel = LDA_sizePerModel); modelIndices <- lsfs$model_indices; lsfs <- lsfs$sfs
  invisible(gc(F))
  
  if (QR_multicollinearity) {
    
	if (LDA_minVarPerModel==FALSE) {
		# remove 0-varianced columns
		sds <- apply(lsfs, 2, function(z) sqrt(var(z, na.rm=T)))
		euCols <- sds > LDA_minVariance
		cat(paste("\t- Removed",sum(!euCols),"columns with too low variance __across__ all models.\n"))
		cat(paste("\t- Kept",sum(euCols),"columns.\n"))
	} else {
		sds <- simplify2array(by(lsfs, modelIndices, function(z) apply(z, 2, function(zz) sqrt(var(zz)) > LDA_minVariance)))
		euCols <- apply(sds, 1, function(z) all(z==TRUE))
		cat(paste("\t- Removed",sum(!euCols),"columns with too low variance __in at least__ one model.\n"))
		cat(paste("\t- Kept",sum(euCols),"columns.\n"))
	}
	
	if (!is.null(LDA_NonNullP)) {
		nnp <- apply(lsfs, 2, function(x) sum(x!=0)/length(x) >= LDA_NonNullP)
		euCols = as.logical(euCols * nnp)
		cat(paste0("\t- Found ",sum(!nnp)," columns having less than ",LDA_NonNullP,"% non-null values.\n"))
	}
    
    # perform QR decomposition
    colnames <- names(euCols)[euCols]
    QR <- qr(x = lsfs[,euCols], tol = 1e-7, LAPACK = FALSE)
    rank <- QR$rank
    cat(paste("\t- Removed",sum(euCols)-rank,"collinear columns, based on QR decomposition.\n"))
    euQRcols <- QR$pivot[1:rank]
    euQRcols <- colnames[euQRcols]
    euCols[!names(euCols)%in%euQRcols] <- FALSE
    
  } else {
	warning("LDA_minVarPerModel not taken into acount when QR_multicollinearity=FALSE")
	
    # detect zero-varianced columns
    sds <- apply(lsfs, 2, sd, na.rm=T)
    euCols <- sds > LDA_minVariance
    cat(paste("\t- Removed",sum(!euCols),"columns with near zero variance.\n"))
    
    if (F) {
    print(table(euCols))
    TMP <- rep(TRUE, length(euCols))
    md <- unique(modelIndices)
    for (i in seq_along(md)) {
      print(i)
      TMP <- TMP * apply(lsfs[modelIndices==md[i],], 2, function(z) sd(z, na.rm=T) > LDA_minVariance)
    }
    euCols <- as.logical(euCols * TMP)
    print(table(euCols))
    }
    
    
    # detect collinear columns, if requested
    if (removeCollinearCols) {
      cat("\t - Removing collinear columns using caret::findCorrelation algorithm.\n")
      nonCollin <- rep(TRUE, ncol(lsfs)); names(nonCollin) <- colnames(lsfs)
      tmp <- caret::findCorrelation(cor(lsfs[,euCols]), cutoff=.9, names=TRUE)
      nonCollin[names(nonCollin)%in%tmp] <- FALSE
      euCols <- euCols & nonCollin
      rm("tmp"); invisible(gc(F))
    }
    
  }
  
  # LDA
  cat('LDA computation.\n')
  if (!LDA_fastDiag) {
    cat("\t- Using MASS::lda algorithm.\n")
    lda <- MASS::lda(x = lsfs[,euCols], grouping = modelIndices, tol = LDA_tol)
    ldaScores <- predict(lda, lsfs[,euCols])$x
  } else {
    cat("\t- Using HiDimDA::Dlda algorithm based on a diagonal covariance matrix estimator.\n")
    lda <- HiDimDA::Dlda(data = lsfs[,euCols], grouping = as.factor(modelIndices))
    ldaScores <- predict(lda, lsfs[,euCols])$Z
  }
  if (F) {
  cat('PCA computation.\n')
  lda <- prcomp(lsfs[,euCols], retx=F, center=F, scale.=F, tol=.1)
  ldaScores <- predict(lda, lsfs[,euCols])
  cat(dim(ldaScores))
  }
  
  # compute null distrib of distances for each model
  cat("\nComputation of per-model null distance distributions, for future goodness-of-fit tests.\n")
  
  dists <- by(ldaScores, modelIndices, function(z) get_dists(z, GOF_n, GOF_abcTol))  
  
  # store results
  x$DIMREDUC$LDA <- list("model" = lda, 
                          "euCols" = euCols,
                          "scores" = ldaScores,
                          "modelIndices" = modelIndices,
                          "nullDists" = dists)
  
  rm(list=c("lsfs", "lda", "ldaScores", "sds", "dists")); invisible(gc(F))

  
  ###########
  ### PLS ###
  ###########
  
  cat("\n>>> PLS\n")
  
  colz <- c("sweepAge", "recRate")
  if (doPLSrecRate==FALSE) colz <- c("sweepAge")

  PLS <- lapply(as.list(names(x$PRIORS)), function(z) {
    get_pls(x$PRIORS[[z]][,colz,drop=FALSE],
            sfs[[z]],
            PLS_normalize,
            removeCollinearCols,
            PLS_ncomp,
            PLS_plotCV,
            maxPLSncomp = PLS_maxncomp,
            deme = z,
            QR_multicollinearity = QR_multicollinearity)
  })
  names(PLS) <- names(x$PRIORS)
  
  # return all
  x$DIMREDUC$PLS <- PLS
  rm("PLS"); invisible(gc(F))
  
  if (compress_SFS) {
    cat("The SFS slot was compressed into a sparseMatrix object. You can access it like a common matrix.\n")
    x$SFS <- lapply(x$SFS, function(z) Matrix::Matrix(z, sparse = TRUE))
    invisible(gc(F))
  }
  
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
    x$DIMREDUC$genetic_sumstats <- FALSE
  }
  
  return(x)
}


#' @title Project an observed multiSFS onto LDA/PLS axes, stored in a referenceTable object
#' @param reftb a reference table (class \emph{referenceTable}) with non-NULL DIMREDUC element
#' @param target the observed multiSFS
#' @param method (string) either "LDA" (model selection) or "PLS" (parameter estimation)
#' @param focalIsland index of the focal island
#' @export
#! @keywords internal
ABC_project_target <- function(reftb, target, method, focalIsland = NULL) {
  
  if (is.null(reftb$DIMREDUC)) return(NULL)
  
  target <- c(unlist(target))
  
  if (reftb$DIMREDUC$GENERAL$relativeSFS) {
    sm <- sum(target, na.rm = TRUE) 
    if (sm != 0) target <- target / sm
  }
  
  if (method=="LDA") {
    
    euc <- reftb$DIMREDUC$LDA$euCols
    obs <- matrix(target[euc], nrow=1)
    colnames(obs) <- names(euc)[euc]
    if (class(reftb$DIMREDUC$LDA$model) == "lda") {
        obs <- predict(reftb$DIMREDUC$LDA$model, obs)$x
    } else {
        obs <- predict(reftb$DIMREDUC$LDA$model, obs)$Z
    }

  } else if (method=="PLS") {
    
    if (is.null(focalIsland)) stop("When switching method to PLS, you must provide the index of the focal island.")
    obs <- target
    if (reftb$DIMREDUC$GENERAL$PLS_normalize) {
      obs <- (target - reftb$DIMREDUC$PLS[[focalIsland]]$mean) / reftb$DIMREDUC$PLS[[focalIsland]]$sd
    }
    obs <- matrix(obs[reftb$DIMREDUC$PLS[[focalIsland]]$euCols], nrow=1)
    obs <- predict(object = reftb$DIMREDUC$PLS[[focalIsland]]$model, 
                   newdata = obs, 
                   type = "scores")
    
  } else {
    
    return(NULL)
    
  }

  if (reftb$DIMREDUC$genetic_sumstats==TRUE) obs <- cbind(obs, tajima(target, internal = TRUE, do_H = !x$GENERAL$folded))
  return(obs)
}




