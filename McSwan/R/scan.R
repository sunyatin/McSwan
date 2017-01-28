
#' @title Assign a model based on the Bayes Factors
#' @details postproba: a vector of model posterior probabilities, the first element must be the neutral model
#' @details this vector must be NAMED by the model names ("i0", "i1", "i2", ... by default)
#' @keywords internal
impute <- function(postproba, cutoff) {
  sBF <- postproba / postproba[1]
  i <- which(sBF >= cutoff)
  if (length(i)==0) return(list("i0"))
  if (length(i)>1) return("failedDISCRIMINABILITY")
  bf <- ifelse(sBF[i]==Inf, 300, sBF[i])
  
bf <- ifelse(sBF[i]>300, 300, sBF[i])
  
  #bf <- postproba[i]
  return(list(sweepingIsland=names(postproba)[i], BayesFactor=bf))
}


#' @title Sweep detection & sweep age estimation for a single genomic region
#' @param target multidimensional joint site-frequency-spectra of the target genomic window
#' @param reftb a \code{referenceTable} object, with non-empty \code{PRIORS}, \code{SFS} and \code{DIMREDUC} elements
#' @param plot_simCloud (logical) plot the observed dataset amid the simulated ones, in the space of the retained LDA components
#' @param verbose (logical) verbose mode
#' @param minSNP (integer) minimum number of within-window SNPs; if the number of within-window SNPs \code{< minSNP}, the function will return \code{NA}
#' @param tolABC (numeric or an array of 2 elements, values between 0 and 1) if a single numeric value, the tolerance ratio for both sweep detection and sweep estimation (see \code{\link{postpr}} and \code{\link{abc}}); if an array of 2 numeric elements, the first one will be the tolerance ratio for sweep detection (see \code{\link{postpr}}) and the second one for sweep estimation (see \code{\link{abc}})
#' @param tolGFIT (numeric between 0 and 1) significance level for the goodness-of-fit test, ie. if \code{p-value(GoF) < tolGFIT}, we reject Ho: "The observed dataset falls within the null distribution of distances"; i.e. the observed dataset is not explainable by a given evolutionary model; if you want to disable the GoF test, set \code{tolGFIT = 0}
#' @param cutoff (numeric) the Bayes Factor cut-off value above which a window is assigned to the fittest model
#' @return A list of 4 elements:
#' \itemize{
#' \item \code{sweepingIsland} the population which has experienced a recent positive sweep
#' \item \code{BayesFactor} the Bayes Factor in favour of the model
#' \item \code{params} point estimate for (\code{sweepAge})
#' \item \code{ic} a vector giving the lower and upper limits of the 95\% credible intervals for (\code{sweepAge}
#' }
#' @seealso \code{\link{gscan}} for iterating this function by sliding windows along the genome
#' @keywords internal
analyze <- function(target, 
                    reftb,
                    minSNP = 2, 
                    tolABC = .01, 
                    tolGFIT = .05,
                    cutoff = 5,
                    plot_simCloud = FALSE,
                    verbose = F) {

  abc_summary_stat <- "median"
  doESTIMATION <- TRUE
  #plot_simCloud <- FALSE

  # target conditions
  if (is.null(reftb$DIMREDUC)) stop("You have not performed the dimension reduction!")
  if (is.data.frame(target)) target <- unlist(target)
  target <- c(target)
  if (length(target) != ncol(reftb$SFS$"i0")) stop("target SFS has not the same number of cells as simulated SFS")
  if (any(is.na(target))) stop("target SFS contains NA values")
  if (sum(target)==1 && length(unique(target))>2) stop("target SFS must contain absolute SNP counts")
  if (sum(target) < minSNP) return("failedMINSNP")
  if (length(tolABC)==1) tolABC <- rep(tolABC, 2)
  
  # get vars
  ss <- reftb$DIMREDUC$LDA$scores
  modelIndices <- reftb$DIMREDUC$LDA$modelIndices
  
  # target treatment
  tproj <- project_target(reftb, target, method="LDA")
  colnames(ss) <- colnames(tproj) <- paste("L",1:ncol(ss),sep="")
  
  # simulation cloud
  if (plot_simCloud) {
    if (ncol(ss)==1) { yy <- ss[,1]*0; yyo <- 0 } else { yy <- ss[,2]; yyo <- tproj[2] }
    dtc <- data.frame(group=reftb$DIMREDUC$LDA$modelIndices, x=ss[,1], y=yy)
    gz=ggplot2::ggplot(data=dtc) + geom_point(aes(x=x, y=y, col=factor(group)), alpha=.5) +
      geom_point(x=tproj[1], y=yyo, cex=6, shape="+", col="purple") +
      xlab("LDA 1") + ylab("LDA 2")
    print(gz)
    Sys.sleep(1)
  }

  # goodness of fit test
  if (tolGFIT!=0) {
    models <- unique(modelIndices)
    gfitOK <- rep(FALSE, length(models))
    for ( m in seq_along(models) ) {
      pp <- ss[modelIndices==models[m],, drop=F]
      tp <- suppressWarnings(abc::abc(tproj, 
                                 rep(1,nrow(pp)),
                                 pp, 
                                 tol = tolABC, 
                                 corr = TRUE, 
                                 method = "rejection"))
      d <- mean(tp$dist)
      dists <- reftb$DIMREDUC$LDA$nullDists[[models[m]]]
      gfitOK[m] <- (sum(dists > d)/length(dists)) > tolGFIT
    }
    gfitOK <- ifelse(any(gfitOK==TRUE), TRUE, FALSE)
    if (gfitOK==FALSE && verbose==TRUE) cat("\n>> Out of all model clouds!\n")
  } else {
    gfitOK <- TRUE
  }
  
  # detect if sweep
  if (!gfitOK) return("failedGFIT")

  pp <- try(suppressWarnings(summary(suppressWarnings(abc::postpr(target = tproj,
                                index = modelIndices,
                                sumstat = ss,
                                tol = tolABC[1],
                                method = "mnlogistic",
                                corr = TRUE,
                                kernel = "epanechnikov")), print=F)), silent = TRUE)
  
  att <- attributes(pp)$"class"
  if (!is.null(att) && att=="try-error") return("failedABC")
  
  if ("mnlogistic"%in%names(pp)) {
    pp <- pp$mnlogistic$Prob
  } else {
    #print("failed POSTPR zeroVariance")
    pp <- pp$Prob
  }
  
  imp <- impute(pp, cutoff)
  if (!is.list(imp) && imp=="failedDISCRIMINABILITY") return("failedDISCRIMINABILITY")
  sweepIsl <- imp[[1]]
  if (sweepIsl=="i0") return("i0")
  
  if (doESTIMATION) {
  
    # estimate if not "i0"
    lb <- t(apply(reftb$PRIORS[[sweepIsl]], 2, range))
    ss <- reftb$DIMREDUC$PLS[[sweepIsl]]$scores
    tproj <- project_target(reftb, target, method="PLS", focalIsland = sweepIsl)
    colnames(ss) <- colnames(tproj) <- 1:ncol(ss)
    
    par <- reftb$PRIORS[[sweepIsl]]
    #par$recRate <- NULL
  
    invisible(capture.output(aa <- try(suppressWarnings(abc::abc(target = tproj,
                               param = par,
                               sumstat = ss,
                               tol = tolABC[2],
                               method = "loclinear",
                               hcorr = TRUE,
                               transf = "logit",
                               logit.bounds = lb )), silent = TRUE)))
    
    att <- attributes(aa)$"class"
    if (!is.null(att) && att=="try-error") {
      return("failedABC")
    } else {
      if ("adj.values"%in%names(aa)) {
        prm <- apply(aa$adj.values, 2, get(abc_summary_stat), na.rm=TRUE)
        adjvals <- aa$adj.values[,1]
      } else {
  #print("failed ABC zeroVariance")
  #prm <- apply(aa$unadj.values, 2, mean, na.rm=TRUE)
  prm <- c(NA, NA)
  adjvals <- NA
      }
      prm <- list(prm)
      adjvals <- list(adjvals)
     
      ic <- apply(aa$adj.values, 2, quantile, c(.025, .975), na.rm=T)
      IC <- c(ic)
      names(IC) <- paste(rep(colnames(aa$adj.values), each=2), rep(c("low","up"), ncol(aa$adj.values)), sep=".")
      
      return(c(imp, params=prm, ic=list(IC), ageAdjValues=adjvals))
    }
  } else {
    return(c(imp, params=list(c(NA,NA)), ic=list(c(NA,NA,NA,NA)), ageAdjValues=list(c(NA))))
  }
}


#' @title Genome scan to detect selective sweeps and estimate sweep ages
#' @description This function scans a genome to detect selective sweeps and estimate their age using sliding (possibly overlapping) windows of constant length. For each window, it estimates the Bayes factor of the most likely evolutionary model. If the window is successfully assigned to a model of selective sweep, the function estimates the age of the sweep.
#' @param X either an object of class \code{observedDataset} produced by \code{\link{get_SFS}} or of class \code{validationTable} produced by \code{\link{generate_pseudoobs}}
#' @param reftb a \code{referenceTable} object, with non-empty \code{PRIORS}, \code{SFS} and \code{DIMREDUC} elements
#' @param minSNP (integer) minimum number of within-window SNPs; if the number of within-window SNPs \code{< minSNP}, the function will return \code{NA}
#' @param tolABC (numeric or an array of 2 elements, values between 0 and 1) if a single numeric value, the tolerance ratio for both sweep detection and sweep estimation (see \code{\link{postpr}} and \code{\link{abc}}); if an array of 2 numeric elements, the first one will be the tolerance ratio for sweep detection (see \code{\link{postpr}}) and the second one for sweep estimation (see \code{\link{abc}})
#' @param tolGFIT (numeric between 0 and 1) significance level for the goodness-of-fit test, ie. if \code{p-value(GoF) < tolGFIT}, we reject Ho: "The observed dataset falls within the null distribution of distances"; i.e. the observed dataset is not explainable by a given evolutionary model; if you want to disable the GoF test, set \code{tolGFIT = 0}
#' @param cutoff (numeric) the Bayes Factor cut-off value above which a window is assigned to the fittest model
#' @param startPos (integer) the start position of the genome scan (if \code{NULL} will start at the first SNP of your observed dataset)
#' @param lastPos (integer) the stop position of the genome scan (if \code{NULL} will stop at the last SNP of your observed dataset)
#' @param windowSlide (integer) the sliding gap, if \code{NULL} then \emph{windowSlide} will be set equal to the window size (windows will thus slide without overlap)
#' @return A list of two elements: \code{res} and \code{adjVals}. The \code{adjVals} element contains a list of adjusted posterior values for the sweep ages, this information is used internally by the function \code{\link{thin}}. The \code{res} element contains a dataframe with genomic windows in rows. In columns: \itemize{
#' \item \code{start.pos} the first position of the window
#' \item \code{end.pos} the last position of the window
#' \item \code{deme.under.sel} the population which has been detected to have experienced a recent positive sweep 
#' \item \code{BayesFactor} the Bayes factor in favour of the selective model for \code{deme.under.sel}
#' \item \code{sweepAge} the posterior estimate of the sweep age
#' \item \code{sweepAge.IC.low} the lower limit of the 95\% credible interval for the sweep age estimation
#' \item \code{sweepAge.IC.up} the upper limit of the 95\% credible interval for the sweep age estimation
#' \item \code{failure} if the sweep detection has failed, it can be (i) because of an insufficient minimum number of SNPs (\code{"minSNP"}); (ii) because the null hypothesis of the goodness-of-fit test was rejected (\code{"goodnessOfFit"}); (iii) because the ABC failed due to zero variance in a LDA/PLS component (\code{"ABC"}); (iv) because various models were not distinguishable based on their Bayes factor (i.e. case of equally likely models) (\code{"modelDiscrimination"}).
#' }
#' @seealso \code{\link{thin}}, \code{\link{get_SFS}}, \code{\link{generate_pseudoobs}}, \code{\link{abc}}, \code{\link{postpr}}
#' @export
gscan <- function(X, reftb, minSNP, startPos = NULL, lastPos = NULL, windowSlide = NULL, tolABC = .01, tolGFIT = .05, cutoff = 5, printProgress = T, plot_simCloud = FALSE) {
  
  # preciser et corriger si windowSize observed != windowSize reftb
  
  # internally set options
  #plot_simCloud <- FALSE
  
  if (class(X)=="validationTable") {
    
    # initialize linearized df for trues
    TR_PRIORS <- do.call(rbind, X$PRIORS)
    TR_PRIORS <- data.frame(sweepingIsland=sapply(rownames(TR_PRIORS), function(x) {
      unlist(strsplit(x, "\\."))[1]
    }), TR_PRIORS, stringsAsFactors=F)
    TR_SFS <- do.call(rbind, X$SFS)
    
    # initialize df for estimates
    ES_PRIORS <- as.data.frame(matrix(NA, ncol=ncol(TR_PRIORS), nrow=nrow(TR_PRIORS)))
    names(ES_PRIORS) <- names(TR_PRIORS)
    rownames(ES_PRIORS) <- rownames(TR_PRIORS)
    ES_PRIORS$BayesFactor <- NA
    
    FAILURES <- as.data.frame(matrix(FALSE, ncol=4, nrow=nrow(TR_PRIORS)))
    rownames(FAILURES) <- rownames(TR_PRIORS)
    names(FAILURES) <- c("failedMINSNP", "failedGFIT", "failedABC", "failedDISCRIMINABILITY")
    
    for (i in 1:nrow(TR_SFS)) {
      cat(i, ".")

      A <- analyze(TR_SFS[i,], reftb, minSNP = minSNP, tolABC = tolABC, tolGFIT = tolGFIT, cutoff = cutoff, plot_simCloud = plot_simCloud, verbose = FALSE)
     
      if (is.null(A)) stop("err!")
      if (!is.list(A) && A=="failedMINSNP") {
        FAILURES[i,"failedMINSNP"] <- TRUE
      } else if (!is.list(A) && A=="failedGFIT") {
        FAILURES[i,"failedGFIT"] <- TRUE
      } else if (!is.list(A) && A=="failedABC") {
        FAILURES[i,"failedABC"] <- TRUE
      } else if (!is.list(A) && A=="failedDISCRIMINABILITY") {
        FAILURES[i,"failedDISCRIMINABILITY"] <- TRUE
      } else {
        if (!is.list(A) && A=="i0") {
          ES_PRIORS[i,"sweepingIsland"] <- "i0"
        } else {
          ES_PRIORS[i,"sweepingIsland"] <- A$sweepingIsl
          ES_PRIORS[i,"BayesFactor"] <- A$BayesFactor
          ES_PRIORS[i,names(A$params)] <- A$params
        }
      }
      
    }
    
    names(FAILURES) <- c("failMINSNP","failGFIT","failABC","failDISCRIM")
    FAILURES <- apply(FAILURES, 2, as.integer)
if (any(FAILURES)!=0) print(FAILURES)
    ES_PRIORS <- data.frame(ES_PRIORS, FAILURES, stringsAsFactors=F)
    out <- list(trues=TR_PRIORS, estimates=ES_PRIORS, failures=FAILURES)
    class(out) <- "validationResult"
    return(out)
  }
  
  
  if (class(X) == "observedDataset") {
    
    L <- reftb$GENERAL$windowSize
    template <- X$template
    PAC <- X$obsData
    #POS <- X$obsData$POS
    #NSNP <- X$obsData$nSNP
    
    if (X$folded != reftb$GENERAL$folded) stop("folding is not shared between the reftb and the observed dataset, please make sure you have folded or unfolded both")
    
    if (is.null(windowSlide)) windowSlide <- L
    if (is.null(startPos)) startPos <- min(PAC$POS)
    if (is.null(lastPos)) lastPos <- max(PAC$POS)
    
    res <- c(); adjVals <- list()
    while (TRUE) {
      if (printProgress) cat(prettyNum(startPos, big.mark=",", scientific=FALSE), " >")
      
      endPos <- startPos + L
      
      # get local observed SFS
      oSFS <- template
      #tmp <- table(PAC[which(POS >= startPos & POS <= endPos)])
      #oSFS[names(tmp)] <- tmp
      localPAC <- subset(PAC, POS >= startPos & POS < endPos)
      tmp <- table(subset(localPAC, nSNP == 1)$PAC)
      oSFS[names(tmp)] <- tmp
      tmp <- table(subset(localPAC, nSNP == .5)$PAC)
      oSFS[names(tmp)] <- oSFS[names(tmp)] + tmp / 2
#oSFS <<- oSFS
      
      #A <- analyze(oSFS, reftb, minSNP = minSNP, ...)
	  A <- analyze(oSFS, reftb, minSNP = minSNP, tolABC = tolABC, tolGFIT = tolGFIT, cutoff = cutoff, plot_simCloud = plot_simCloud, verbose = FALSE)


      d <- NA; bf <- NA; est <- c(NA); ic <- rep(NA, 2); fail <- NA; agevals <- NA
      if (is.null(A)) stop("err!")
      if (!is.list(A) && A=="failedMINSNP") {
        d <- NA
        fail <- "minSNP"
      } else if (!is.list(A) && A=="failedGFIT") {
        d <- NA
        fail <- "goodnessOfFit"
      } else if (!is.list(A) && A=="failedABC") {
        d <- NA
        fail <- "ABC"
      } else if (!is.list(A) && A=="failedDISCRIMINABILITY") {
        d <- NA
        fail <- "modelDiscrimination"
      } else {
        if (!is.list(A) && A=="i0") {
          d <- "i0"
        } else {
          d <- A$sweepingIsl
          bf <- A$BayesFactor
          est <- A$params
          ic <- A$ic
          agevals <- A$ageAdjValues
        }
      }
      
      res <- rbind(res, c(startPos, endPos, d, bf, est, ic, fail))
      adjVals <- c(adjVals, list(agevals))
      
      # move window
      startPos <- startPos + windowSlide
      if (startPos >= lastPos) break()
      
    }
    
    res <- as.data.frame(res, stringsAsFactors=F)
    #names(res) <- c("start.pos", "end.pos", "deme.under.sel", "BayesFactor", "sweepAge", "recRate",
    #                "sweepAge.IC.low", "sweepAge.IC.up", "recRate.IC.low", "recRate.IC.up", "failure")
	names(res) <- c("start.pos", "end.pos", "deme.under.sel", "BayesFactor", "sweepAge", "sweepAge.IC.low", "sweepAge.IC.up", "failure")
    for (j in c(1,2,4,5,6,7)) res[,j] <- as.numeric(res[,j])
    
    return(list(res=res, adjVals=adjVals))
  }
  
}




