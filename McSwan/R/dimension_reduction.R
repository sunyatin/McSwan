

#' @title Merge model-specific multiSFSs into one big matrix
#' @param sfs a list of model-specific multiSFS matrices
#' @return A list of two elements: a vector of model indices; a big matrix containing all merged multiSFSs.
#' @keywords internal
linearize <- function(sfs) {
  model_indices <- gsfs <- c()
  for (i in seq_along(sfs)) {
    model_indices <- c(model_indices, rep(names(sfs)[i], nrow(sfs[[i]])))
    gsfs <- rbind(gsfs, sfs[[i]])
  }
  return(list(model_indices=model_indices, sfs=gsfs))
}


#' @title Get the null distribution of distances to some pseudo-observed datastes
#' @param x a model-specific multiSFS matrix
#' @param n_dist_btsp_per_model number of distances to perform per model for computation of the null distribution
#' @return A vector of null distances.
#' @keywords internal
get_dists <- function(x, n_dist_btsp_per_model, abcTolerance) {
  targets <- sample(1:nrow(x), n_dist_btsp_per_model, replace=T)
  d <- rep(NA, n_dist_btsp_per_model)
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
get_pls <- function(param, ss, PLS_normalize, removeCollinearCols, PLS_ncomp, plot_PLScv,
                    maxPLSncomp = 100, deme, QR_multicollinearity = T) {
  
  #if (!is.null(dim(param)) && ncol(param)>1) stop("PLS is computed against a single parameter andyou provided more than one parameter")
  param <- as.matrix(param)

  # normalize
  means <- apply(ss, 2, mean, na.rm=T)
  sds <- apply(ss, 2, sd, na.rm=T)
  euCols <- sds != 0
  cat(paste("removed",sum(!euCols),"columns with zero variance\n"))
  if (PLS_normalize) {
    cat("\nmodel-specific standardization of multiSFS for PLS analysis\n")
    ss <- t(apply(ss, 1, function(x) (x-means)/sds))
  }
  
  if (QR_multicollinearity) {
    
    # perform QR decomposition
    colnames <- names(euCols)[euCols]
    QR <- qr(x = ss[,euCols], tol = 1e-7, LAPACK = FALSE)
    rank <- QR$rank
    cat(paste("removed",sum(euCols)-rank,"collinear columns based on QR decomposition\n"))
    euQRcols <- QR$pivot[1:rank]
    euQRcols <- colnames[euQRcols]
    euCols[!names(euCols)%in%euQRcols] <- FALSE
    
  } else {
  
    # detect collinear columns, if requested
    if (removeCollinearCols) {
      cat("removing collinear columns in multiSFS\n")
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
    cat("\nAutomatic selection of the fittest number of PLS components\n")
    
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
      lin <- function(x) slope * x + intercept
      
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
    cat(paste("\nFittest K for PLS:",PLS_ncomp,"\n"))
  }

  # PLS
  cat("\nPLS\n")
  pls <- pls::plsr(param ~ ss, scale = FALSE, validation = "none", ncomp = PLS_ncomp, model = F, x = F, y = F)
  pls$scores <- NULL
  pls.scores <- predict(pls, ss, type="scores")
  
  # output
  PLS <- list("model" = pls,
              "euCols" = euCols,
              "scores" = pls.scores, 
              "mean" = means, 
              "sd" = sds)
  rm(list=c("pls", "pls.scores", "means", "sds", "param")); invisible(gc(FALSE))
  
  return(PLS)    
}


#' @title Dimension reduction for multidimensional SFSs
#' @description Performs a series of dimensional reductions: LDA (for later model selection) and model-specific PLSs (for later parameter estimation).
#' @param x (referenceTable) a reference table with \code{SFS} elements
#' @param PLS_ncomp (integer) if NULL, automatically select the fittest number of PLS components, else an integer to manually set the number of components to retain
#' @param PLS_normalize (logical) whether to normalize the multiSFSs prior to PLS analysis (TRUE is recommended)
#' @param removeCollinearCols (logical) whether to remove collinear columns across the multiSFSs (TRUE is recommended for LDA analysis)
#' @param doPLSrecRate (logical) whether to regress (in addition to the sweep ages) onto the recombination parameter in the PLS analyses (since the recombination rate is considered a nuisance parameter, it can be sometimes useful to discard it from the PLS analysis)
#' @param n_dist_bstp_per_model (integer) number of jacknifes to perform per model to compute the null distribution of distances for later goodness-of-fit tests
#' @param relativeSFS (logical) if TRUE, will use relative SFS, else will use absolute SFS (with absolute SNP counts per frequency class) (FALSE is recommended)
#' @param nonvarTolerance (numeric) variance threshold across the multiSFSs for a given frequency class, if the variance is below this threshold, the frequency class will be discarded for LDA analysis
#' @param abcTolerance (numeric < 1) ABC tolerance ratio for the distance-based goodness-of-fit test (\emph{cf. n_dist_btsp_per_model} and see \code{\link{abc}})
#' @param plot_PLScv (logical) if TRUE, will plot the cross-validation graphs of the PLS
#' @return An object of class \code{referenceTable} with filled-in \code{DIMREDUC} element.
#' @export
dim_reduction <- function(x,
                          removeCollinearCols = TRUE,
                          tolLDA = 1e-4,
                          minSDlda = 0,
                          
                          PLS_normalize = TRUE,
                          PLS_ncomp = NULL,
                          doPLSrecRate = TRUE,
                          
                          n_dist_bstp_per_model = 2e3,
                          relativeSFS = FALSE,
                          nonvarTolerance = 1.0e-4,
                          abcTolerance = .01,
                          
                          QR_multicollinearity = TRUE,
                          
                          plot_PLScv = FALSE) {
  
  if (class(x)!="referenceTable"||is.null(x$SFS)) stop("x is not a referenceTable object or its SFS element is NULL")
  
  x$DIMREDUC <- list()
  x$DIMREDUC$GENERAL <- list(relativeSFS = relativeSFS, 
                             removedCollinearColumns = removeCollinearCols,
                             PLS_normalize = PLS_normalize)
  
  # relativize multiSFS, if requested
  if (relativeSFS) {
    cat("\nLower-dimensionned multiSFS is computed based on relative fractions of SNPs")
    sfs <- lapply(x$SFS, function(z) {
      t(apply(z, 1, function(y) {
        sm <- sum(y, na.rm = TRUE)
        if (sm==0) return(y)
        return(y/sm)
      }))
    })
  } else {
    sfs <- x$SFS
  }

  
  ###########
  ### LDA ###
  ###########
  
  cat("\nAcross-model LDA\n")
  lsfs <- linearize(sfs); modelIndices <- lsfs$model_indices; lsfs <- lsfs$sfs
  invisible(gc(F))
  
  if (QR_multicollinearity) {
    
    # remove 0-varianced columns
    sds <- apply(lsfs, 2, function(x) sqrt(var(x, na.rm=T)))
    euCols <- sds > minSDlda
    cat(paste("removed",sum(!euCols),"columns with zero variance\n"))
    
    # perform QR decomposition
    colnames <- names(euCols)[euCols]
    QR <- qr(x = lsfs[,euCols], tol = 1e-7, LAPACK = FALSE)
    rank <- QR$rank
    cat(paste("removed",sum(euCols)-rank,"collinear columns based on QR decomposition\n"))
    euQRcols <- QR$pivot[1:rank]
    euQRcols <- colnames[euQRcols]
    euCols[!names(euCols)%in%euQRcols] <- FALSE
    
  } else {
  
    # detect zero-varianced columns
    sds <- apply(lsfs, 2, sd, na.rm=T)
    euCols <- sds >= nonvarTolerance
    cat(paste("removed",sum(!euCols),"columns with near zero variance\n"))
    
    if (F) {
    print(table(euCols))
    TMP <- rep(TRUE, length(euCols))
    md <- unique(modelIndices)
    for (i in seq_along(md)) {
      print(i)
      TMP <- TMP * apply(lsfs[modelIndices==md[i],], 2, function(x) sd(x, na.rm=T) > nonvarTolerance)
    }
    euCols <- as.logical(euCols * TMP)
    print(table(euCols))
    }
    
    
    # detect collinear columns, if requested
    if (removeCollinearCols) {
      cat("removing collinear columns in multiSFS\n")
      nonCollin <- rep(TRUE, ncol(lsfs)); names(nonCollin) <- colnames(lsfs)
      tmp <- caret::findCorrelation(cor(lsfs[,euCols]), cutoff=.9, names=TRUE)
      nonCollin[names(nonCollin)%in%tmp] <- FALSE
      euCols <- euCols & nonCollin
      rm("tmp"); invisible(gc(F))
    }
    
  }
  
  # LDA
  print('lda')
  lda <- MASS::lda(x=lsfs[,euCols], grouping=modelIndices, tol=tolLDA)
  ldaScores <- predict(lda, lsfs[,euCols])$x
  if (F) {
  print('pca')
  lda <- prcomp(lsfs[,euCols], retx=F, center=F, scale.=F, tol=.1)
  print('predict pca x')
  ldaScores <- predict(lda, lsfs[,euCols])
  print(dim(ldaScores))
  }
  
  # compute null distrib of distances for each model
  cat("computing null distribution of distances for future goodness of fit tests\n")
  
  dists <- by(ldaScores, modelIndices, function(x) get_dists(x, n_dist_bstp_per_model, abcTolerance))  
  
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
  
  cat("\nModel-specific PLS's\n")
  
  colz <- c("sweepAge", "recRate"); if (!doPLSrecRate) colz <- c("sweepAge")

  PLS <- lapply(as.list(names(x$PRIORS)), function(z) {
    get_pls(x$PRIORS[[z]][,colz,drop=FALSE],
            sfs[[z]],
            PLS_normalize,
            removeCollinearCols,
            PLS_ncomp,
            plot_PLScv,
            deme = z,
            QR_multicollinearity = QR_multicollinearity)
  })
  names(PLS) <- names(x$PRIORS)
  
  # return all
  x$DIMREDUC$PLS <- PLS
  rm("PLS"); invisible(gc(F))
  
  return(x)
}


#' @title Project an observed multiSFS onto LDA/PLS components stored in a referenceTable object
#' @param reftb a reference table (class \emph{referenceTable}) with non-NULL DIMREDUC element
#' @param target the observed multiSFS
#' @param method (string) either "LDA" (model selection) or "PLS" (parameter estimation)
#' @param focalIsland index of the focal island
#' @keywords internal
project_target <- function(reftb, target, method, focalIsland = NULL) {
  
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
    obs <- predict(reftb$DIMREDUC$LDA$model, obs)$x
    
  } else if (method=="PLS") {
    
    if (is.null(focalIsland)) stop("We switching method to PLS, you must provide the index of the focal island.")
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

  return(obs)
}




