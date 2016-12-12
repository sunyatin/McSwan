
#' @title impute model based on Bayes Factors
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


#' @title For a single genomic window, detect sweep and estimate sweep age, if appropriates
#' @param target (array or single-row matrix or dataframe) the multiSFS of the target genomic window
#' @param reftb a \emph{referenceTable} object, with activated \code{PRIORS}, \code{SFS} and \code{DIMREDUC} elements
#' @param plot_simCloud (logical) whether to plot the observed dataset among the cloud of simulations datasets in the space of the retained LDA components
#' @param minSNP (integer) minimum number of SNPs required in the window to detect sweep; if the number of SNPs in the window is \eqn{n < minSNP}, the function will return \code{NA}
#' @param tolABC (numeric or an array of 2 elements) (values must be between 0 and 1) if a single numeric value, the tolerance ratio for both sweep detection and sweep estimation (see \code{\link{postpr}} and \code{\link{abc}}); if an array of 2 numeric elements, the first one will be the tolerance ratio for sweep detection (see \code{\link{postpr}}) and the second one for sweep estimation (see \code{\link{abc}})
#' @param tolGFIT (numeric between 0 and 1) the tolerance ratio for goodness-of-fit test (should be equal to the one specified in \code{\link{dim_reduction}})
#' @param verbose (logical) whether to print all messages (TRUE) or not (FALSE)
#' @param cutoff (numeric) the Bayes Factor cutoff value to decide that the sweep model under a specific deme is significantly more likely than the other competing models
#' @param doESTIMATION (logical) whether to do the sweep estimation in addition to the sweep detection
#' @return A list with 4 elements:
#' \itemize{
#' \item \code{sweepingIsland} the deme which has likely undergone a recent positive sweep
#' \item \code{BayesFactir} the Bayes Factor in favour of this corresponding model
#' \item \code{params} a vector of point estimates for the sweep age (\code{sweepAge}) and recombination rate (\code{recRate})
#' \item \code{ic} a vector of the lower and upper limits of the 95\% credible intervals for each point estimate (\code{sweepAge} and \code{recRate})
#' }
#' @seealso \code{\link{gscan}} for iterating this function along the genome
#' @export
# sweep_detection (postpr)
analyze <- function(target, 
                    reftb,
                    plot_simCloud = FALSE, 
                    minSNP = 2, 
                    tolABC = .01, 
                    tolGFIT = .1,
                    cutoff = 5,
                    doESTIMATION = T,
                    verbose = F,
                    forceESTIMATION = F) {
  
  abc_summary_stat <- "median"

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

#' @title Performs a genome scan, by detecting sweeps along sliding windows and estimating sweep ages
#' @description This function performs a genome scan along sliding windows, for each window it estimates the Bayes Factor of the most likely sweep model in a specific deme \emph{vs.} the neutral model. If a sweep is detected, it carries on estimating the age of the sweep.
#' @param X either an object with class \emph{observedDataset} produced by \code{get_SFS()} or with class \emph{validationTable}
#' @param reftb a \emph{referenceTable} object, with activated \code{PRIORS}, \code{SFS} and \code{DIMREDUC} elements
#' @param startPos (integer) the starting position of the genome scan
#' @param lastPos (integer) the ending position of the genome scan
#' @param windowSlide (integer) the sliding gap, if \code{NULL} then \emph{windowSlide} will be automatically set as equal the window size (therefore windows will not overlap)
#' @param ... other arguments passed to \code{\link{analyze}}
#' @return A dataframe with each row corresponding to one window and with the following columns: \itemize{
#' \item \code{start.pos} the starting position of the window
#' \item \code{end.pos} the ending position of the window
#' \item \code{deme.under.sel} the deme which has been detected to have undergone recent positive sweep 
#' \item \code{BayesFactor} the Bayes Factor for the sweep model of the sweeping deme (ie. \code{deme.under.sel})
#' \item \code{sweepAge} the point estimate of the sweep age
#' \item \code{recRate} the point estimate of the recombination rate between adjacent bases
#' \item \code{sweepAge.IC.low} the lower bound of the 95\% credible interval for the sweep age estimation
#' \item \code{sweepAge.IC.up} the upper bound of the 95\% credible interval for the sweep age estimation
#' \item \code{recRate.IC.low} the lower bound of the 95\% credible interval for the recombination rate estimation
#' \item \code{recRate.IC.up} the upper bound of the 95\% credible interval for the recombination rate estimation
#' \item \code{failure} if the sweep detection has failed, specifies the reason: either because the window failed the requirement for a minimum number of SNPs (\code{minSNP}); failed the goodness-of-fit test (\code{goodnessOfFit}); failed the ABC estimation because of zero variance in the LDA/PLS components for instance (\code{ABC}); failed the discriminability between equally likely models (\code{modelDiscrimination}).
#' }
#' @seealso \code{\link{thin}}, \code{\link{get_SFS}}
#' @export
gscan <- function(X, reftb, minSNP, startPos = NULL, lastPos = NULL, windowSlide = NULL, printProgress = T, ...) {
  
  # preciser et corriger si windowSize observed != windowSize reftb
  
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

      A <- analyze(TR_SFS[i,], reftb, minSNP = minSNP, ...)
     
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
      oSFS <<- oSFS
      
      A <- analyze(oSFS, reftb, minSNP = minSNP, ...)
      d <- NA; bf <- NA; est <- c(NA, NA); ic <- rep(NA, 4); fail <- NA; agevals <- NA
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
    names(res) <- c("start.pos", "end.pos", "deme.under.sel", "BayesFactor", "sweepAge", "recRate",
                    "sweepAge.IC.low", "sweepAge.IC.up", "recRate.IC.low", "recRate.IC.up", "failure")
    for (j in c(1,2,4,5,6,7,8,9,10)) res[,j] <- as.numeric(res[,j])
    
    return(list(res=res, adjVals=adjVals))
  }
  
}

if (F) {
  


#!' @title Timming the sweep detections to get a confident point estimate for sweep ages
#!' @export
thin_density <- function(x, reftb, stat = "mean", bwadjust = 1, minWindows = 1, epsilon = 1e-15) {
  
  dn <- paste("i",deme,sep="")
  
  y <- subset(x, !is.na(deme.under.sel) & deme.under.sel == dn)
  gr <- c(min(x$start.pos), max(x$end.pos))
  ymid <- (y$start.pos+y$end.pos)/2
  
  if (is.null(nrow(y))||nrow(y)<minWindows) return(NA)
  
  if (is.null(bwadjust)) {
    D <- density( ymid,
               from = gr[1],
               to = gr[2])
  } else {
    D <- density( ymid,
                  bw = reftb$GENERAL$windowSize * bwadjust,
                  from = gr[1],
                  to = gr[2])
  }
  
  plot(D, main="Gaussian spline of sweeped regions")
  points(ymid, rep(0, length(ymid)), pch="|")
  
  D$y[D$y<epsilon] <- 0
  
  # derivate
  dD <- list(x=D$x[-1], y=diff(D$y) / diff(D$x))
  
  # find near-zero intersection
  dDy <- dD$y
  zI <- sapply(2:(length(dDy)-1), function(j) dDy[j-1]>0&&dDy[j]>=0&&dDy[j+1]<0 ) # local maximum
  zI <- c(FALSE, zI, FALSE)
  
  xleft <- dD$x[zI]
  xright <- dD$x[c(FALSE, zI[-length(zI)])]
  xP <- ( xleft + xright ) / 2
  
  # discardOrphanWindow? => still missing in the code
  
  idx <- sapply(xP, function(p) {
    which.min( abs(ymid - p) )
  })
  
  dy <- D$y
  lims <- sapply(which(c(F,zI)), function(p) {
    lm <- which(dy[1:p]<=epsilon)
    lm <- lm[length(lm)]
    
    rm <- which(dy[p:length(dy)]<=epsilon) + (p-1)
    rm <- rm[1]
    
    L <- which.min( abs(ymid - D$x[lm]) )
    L <- ifelse(length(L)==0, min(y$start.pos), y$start.pos[L])
    R <- which.min( abs(ymid - D$x[rm]) )
    R <- ifelse(length(R)==0, max(y$end.pos), y$end.pos[R])
    
    return(c(L, R))
  })
  lims <- t(lims)
  
  ret <- data.frame(sweep.center=xP, 
                    sweep.lbound = lims[,1],
                    sweep.rbound = lims[,2],
                    deme=deme, sweepAge=y$sweepAge[idx], 
                    sweepAge.IC.low=y$sweepAge.IC.low[idx],
                    sweepAge.IC.up=y$sweepAge.IC.up[idx])
  return(list("density"=cbind("x"=D$x, "y"=D$y), "estimation"=ret))
}


#!' @title Thinning the sweep detections to get a confident point estimate for sweep ages
#!' @export
thin_contig <- function(x, reftb, stat = "mean", bwadjust = 1, minWindows = 1, epsilon = 1e-15) {
    
  dn <- paste("i",deme,sep="")
    
  dn <- paste("i",deme,sep="")
  y <- subset(x, !is.na(deme.under.sel) & deme.under.sel == dn)
  L <- reftb$GENERAL$windowSize
  
  if (is.null(nrow(y))||nrow(y)<minWindows) return(NA)

  store <- TRUE
  startp <- y[1,"start.pos"]
  age <- age.low <- age.up <- bf <- c()
  ret <- c()
  for (i in 1:nrow(y)) {
   
    if (i>1) store <- ifelse( y[i,"start.pos"] - y[i-1,"start.pos"] > bwadjust * L, FALSE, TRUE)
   
    if (!store) {
      if (stat=="mean") {
        pe.Age <- mean(age)
        pe.lAge <- mean(age.low)
        pe.rAge <- mean(age.up)
        pe.BF <- mean(bf)
      } else {
        ix <- which.min( abs(y$start.pos+y$end.pos)/2 - (startp+endp)/2 )
        pe.Age <- age[ix]
        pe.lAge <- age.low[ix]
        pe.rAge <- age.up[ix]
        pe.BF <- bf[ix]
      }
      add <- data.frame(sweep.center = (startp+endp)/2, 
                        sweep.lbound = startp,
                        sweep.rbound = endp,
                        BF = pe.BF,
                        deme=deme, 
                        sweepAge = pe.Age, 
                        sweepAge.IC.low = pe.lAge,
                        sweepAge.IC.up = pe.rAge)
      ret <- rbind(ret, add)
      startp <- y[i,"start.pos"]
      age <- age.low <- age.up <- bf <- c()
    }
    
    endp <- y[i,"end.pos"]
    age <- c(age, y[i,"sweepAge"])
    age.low <- c(age.low, y[i,"sweepAge.IC.low"])
    age.up <- c(age.up, y[i,"sweepAge.IC.up"])
    bf <- c(bf, y[i,"BayesFactor"])
    
  }
  
  if (store) {
    if (stat=="mean") {
      pe.Age <- mean(age)
      pe.lAge <- mean(age.low)
      pe.rAge <- mean(age.up)
      pe.BF <- mean(bf)
    } else {
      ix <- which.min( abs(y$start.pos+y$end.pos)/2 - (startp+endp)/2 )
      pe.Age <- age[ix]
      pe.lAge <- age.low[ix]
      pe.rAge <- age.up[ix]
      pe.BF <- bf[ix]
    }
    add <- data.frame(sweep.center = (startp+endp)/2, 
                      sweep.lbound = startp,
                      sweep.rbound = endp,
                      BF = pe.BF,
                      deme=deme, 
                      sweepAge = pe.Age, 
                      sweepAge.IC.low = pe.lAge,
                      sweepAge.IC.up = pe.rAge)
    ret <- rbind(ret, add)
  }
  
  rgg <- c(x$start.pos[1], x$end.pos[nrow(x)])
  ymid <- (y$start.pos+y$end.pos)/2
  plot(ymid, rep(0, length(ymid)), pch="|", type="p", ylim=c(0,1), xlim=rgg)
  segments(x0=ret$sweep.center, y0=0, y1=1, col="gray")
  abline(v=diff(rgg)*.8, col="red")
  
  return(list("density"=cbind("x"=ret$sweep.center, "y"=ret$BF), "estimation"=ret))
}

#!' @title Thinning the sweep detections to get a confident point estimate for sweep ages
#!' @export
thin <- function(x, reftb, method, stat = "mean", ...) {

  if (method=="density") {
    y <- thin_density(x, reftb, stat = stat, ...)
  } else if (method=="contig") {
    y <- thin_contig(x, reftb, stat = stat, ...)
  } else {
    stop("no method provided")
  }
  return(y)
}

}



