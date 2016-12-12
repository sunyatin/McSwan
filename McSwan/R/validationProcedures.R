# validation procedures

#' @title Validation procedure on pseudo-chromosomes
#' @export
sliding_validation <- function(valtb, 
                               reftb,
                               method,
                               stat = "mean",
                               minSNP = 10, 
                               cutoff = 10, 
                               tolABC = .01, 
                               windowSlide = NULL, 
                               tolGFIT = .1,
                               plot_simCloud = F,
                               plot_thinning = F,
                               bwadjust = 1,
                               minWindows = 1,
                               summary_stat = "mean",
                               ...) {
  
  valtb$ANALYSIS <- list()[rep(1,length(valtb$SFS))]
  names(valtb$ANALYSIS) <- names(valtb$SFS)
  
  for (i in 1:length(valtb$SFS)) {
    print(paste("> ",names(valtb$SFS)[i]))
    mSFS <- valtb$SFS[[i]]
    #OO <- list(density=c(), estimation=c())
    for (j in 1:length(mSFS)) {
      obs <- mSFS[[j]]
      print(paste("===",j))
      O <- suppressWarnings(gscan(obs, reftb, 
                                  minSNP = minSNP,
                                  cutoff = cutoff, 
                                  tolABC = tolABC,
                                  windowSlide = windowSlide, 
                                  plot_simCloud = plot_simCloud, 
                                  tolGFIT = tolGFIT,
                                  printProgress = F, ...))
      O <- thin(O, reftb, method=method, stat = stat, bwadjust = bwadjust, minWindows = minWindows, plot_thinning = plot_thinning, summary_stat = summary_stat, obs = obs) #####!!!!!! 06102016
      #OO$density <- rbind(OO$density, O$density)
      #OO$estimation <- rbind(OO$estimation, O$estimation)
      valtb$ANALYSIS[[i]] <- c(valtb$ANALYSIS[[i]], list(O))
      names(valtb$ANALYSIS[[i]])[length(valtb$ANALYSIS[[i]])] <- j
    }

  }
  #ret <- list("validationTable"=ALL, "trueParams"=trues)
  #attr(ret, "deme") <- paste("i",pseudoobs$deme,sep="")
  #attr(ret, "windowSize") <- pseudoobs$windowSize
  return(valtb)
}


#' @title Validation trimming procedure for sweep age point estimation
#' @export
sliding_validation_trim <- function(valtb, reftb, deme, method, which=NULL, ...) {
  
  if (is.null(which)) {
    indices <- seq_along(valtb$validationTable)
  } else {
    indices <- which
  }
  ALL <- valtb$validationTable
  trues <- valtb$trueParams
  
  ests <- trues * NA
  estPos <- rep(NA, nrow(trues))
  
  nConflicts <- 0
  for (i in indices) {
    #cat(i, "=")
    tmp <- trim(ALL[[i]], reftb, deme, method, ...)
    
    if ("estimation"%in%names(tmp)) {
      if (!is.null(tmp$estimation)) {
        if (nrow(tmp$estimation)>1) {
          nConflicts <- nConflicts + 1
        } else {
          ests[i,] <- tmp$estimation[,c("sweepAge")]
          if (method=="contig") {
            estPos[i] <- tmp$density[,"x"]
          } else {
            estPos[i] <- tmp$density[which.max(tmp$density[,"y"]),"x"]
          }
        }
      }
    }
  }
  
  if (nConflicts>0) warning(paste("Discared",nConflicts,"simulations due to conflicting multiple sweeps detected. Try to increase the bandwidth."))
  
  return(list("trues"=trues, "ests"=ests, "estPos"=estPos))
}


#' @title Validation trimming procedure for sweep age point estimation
#' @export
sliding_validation_trim_MODIF <- function(valtb, which=NULL, method, bwadjust = NULL, discardSingleWindow = FALSE, POD = NULL, minWindows = 1) {
  
  if (is.null(which)) {
    indices <- seq_along(valtb$validationTable)
  } else {
    indices <- which
  }
  ALL <- valtb$validationTable
  trues <- valtb$trueParams
  
  ests <- trues * NA
  sweep.center <- c()
  for ( i in indices ) {
    cat(i, ".")
    
    io <- ALL[[i]]
    
    if (sum(io$deme.under.sel==attr(valtb, "deme"), na.rm=T)==0) next()

    i1 <- subset(io, deme.under.sel==attr(valtb, "deme"))
    
    # rec
    if (!is.null(nrow(i1)) && nrow(i1)>=1) {
      if (grepl("^r", method, perl=T)) {
        #ior <- i1[which.min(i1$recRate),]
        if (nrow(i1)>=minWindows) {
          ior <- i1
          ests[i,] <- c(mean(ior$sweepAge), mean(ior$recRate))
          
          #ior <- i1[which.max(i1$BayesFactor),]
          #ests[i,] <- c(mean(ior$sweepAge), mean(ior$recRate))
        }
        
        if (nrow(i1)>1) {
          yy <- scale(i1$recRate) * scale(i1$sweepAge)
        } else {
          yy <- i1$recRate * i1$sweepAge
        }
        plot((i1$start.pos+i1$end.pos)/2, log10(i1$BayesFactor), pch="+", xlim=c(0,8e5), main="BF"); abline(v=.6*8e5, col="red")
        #plot((i1$start.pos+i1$end.pos)/2, log10(i1$recRate), pch="+", xlim=c(0,8e5), main="recRate"); abline(v=.6*8e5, col="red"); segments(x0=(i1$start.pos+i1$end.pos)/2, y0=log10(i1$recRate.IC.low), y1=log10(i1$recRate.IC.up))
        plot((i1$start.pos+i1$end.pos)/2, log10(i1$sweepAge.IC.up-i1$sweepAge.IC.low), pch="+", xlim=c(0,8e5), main="recRate"); abline(v=.6*8e5, col="red")
        
      } else if (grepl("^d", method, perl=T)) {
        if (nrow(i1)<minWindows) {
          #if (!discardSingleWindow) ests[i,] <- c(ior$sweepAge, ior$recRate)
        } else {
          if (!is.null(bwadjust)) {
            D <- with(i1, density( (start.pos+end.pos)/2, bw=as.numeric(attr(valtb, "windowSize"))*bwadjust))
            plot(D, main="Gaussian spline of sweeped regions", xlim=c(0,8e5))
            abline(v=.6*8e5, col="red")
            ymid <- with(i1, (start.pos+end.pos)/2)
            points(ymid, rep(0, length(ymid)), pch="|")
          } else {
            D <- with(i1, density((start.pos+end.pos)/2)) #start.pos ONLY before
          }
          x <- D$x[which.max(D$y)]
          sweep.center <- c(sweep.center, x); abline(v=sweep.center, col="gray")
          x <- which.min( abs( ( (io$start.pos+io$end.pos)/2  - x) ) )
          ior <- io[x,]
          ests[i,] <- c(mean(ior$sweepAge), mean(ior$recRate))
        }
      } else if (method == "ld") {
        pd <- POD; pd[[1]] <- pd[[1]][i]; pd[[2]] <- pd[[2]][i,]
        ff <- pd$SFS[[1]]$obsData
        ff$Q <- as.numeric(unlist(lapply(strsplit(as.character(ff$PAC), "\\."), function(x) x[[2]])))

        res <- tmp <- c()
        for ( u in 1:nrow(i1) ) {
          st <- i1$start.pos[u]
          ed <- i1$end.pos[u]
          tp <- subset(ff, POS >= st & POS <= ed)$Q
          #tp[tp>20] <- 40 - tp[tp>20]
          tmp <- c(tmp, mean(tp))
        }
        #cent <- i1[which.min(tmp),]
        cent <- i1[which.max(tmp),]
        cent <- (cent$start.pos+cent$end.pos)/2
        
        sweep.center <- c(sweep.center, cent)
        
        plot((i1$start.pos+i1$end.pos)/2, tmp, pch="+", xlim=c(0,8e5))
        ymid <- with(i1, (start.pos+end.pos)/2)
        points(ymid, rep(0, length(ymid)), pch="|"); abline(v=.6*8e5, col="red")
        abline(v=cent, col="gray")
        
          if (F) {
            ad <- ff$Q[i]
            # if folding
            if (ad>20) ad <- 40-ad
            
            if (ff$POS[i]-start > 20000) {
              tmp <- c(tmp, ad)
              res <- rbind(res, c(start, mean(tmp)))
              tmp <- c()
              start <- ff$POS[i+1]
            } else {
              tmp <- c(tmp, ad)
            }
            # if folding
            cent <- res[which.min(res[,2]),1]
            
            sweep.center <- c(sweep.center, cent)
    
            plot(res, pch="+")
            ymid <- with(i1, (start.pos+end.pos)/2)
            points(ymid, rep(0, length(ymid)), pch="|"); abline(v=.6*8e5, col="red")
            abline(v=cent, col="gray")
          }
      
      } else if (grepl("^c", method, perl=T)) {
        
      } else {
        stop("missing method")  
      }
    }
    
  }
  
  
  return(list("x"=sweep.center, "trues"=trues, "ests"=ests))
}


#' @title Validation trimming procedure for sweep age point estimation
#' @export
sliding_validation_trim_REF <- function(valtb, method, bwadjust = NULL, discardSingleWindow = FALSE) {
  
  ALL <- valtb$validationTable
  trues <- valtb$trueParams
  
  ests <- trues * NA
  sweep.center <- c()
  for ( i in seq_along(ALL) ) {
    cat(i, ".")
    
    io <- ALL[[i]]
    
    if (sum(io$deme.under.sel==attr(valtb, "deme"), na.rm=T)==0) next()
    
    i1 <- subset(io, deme.under.sel==attr(valtb, "deme"))

    # rec
    if (!is.null(nrow(i1)) && nrow(i1)>=1) {
      if (grepl("^r", method, perl=T)) {
        ior <- i1[which.min(i1$recRate),]
        ests[i,] <- c(mean(ior$sweepAge), mean(ior$recRate))
      } else if (grepl("^d", method, perl=T)) {
        if (nrow(i1)==1) {
          if (!discardSingleWindow) ests[i,] <- c(ior$sweepAge, ior$recRate)
        } else {
          if (!is.null(bwadjust)) {
            D <- with(i1, density( (start.pos+end.pos)/2, bw=as.numeric(attr(valtb, "windowSize"))*bwadjust))
            plot(D, main="Gaussian spline of sweeped regions")
            abline(v=.6*8e5, col="red")
            ymid <- with(i1, (start.pos+end.pos)/2)
            points(ymid, rep(0, length(ymid)), pch="|")
          } else {
            D <- with(i1, density((start.pos+end.pos)/2)) #start.pos ONLY before
          }
          x <- D$x[which.max(D$y)]
          sweep.center <- c(sweep.center, x); abline(v=sweep.center, col="gray")
          x <- which.min( abs( ( (io$start.pos+io$end.pos)/2  - x) ) )
          ior <- io[x,]
          ests[i,] <- c(mean(ior$sweepAge), mean(ior$recRate))
        }
      } else {
        stop("missing method")
      }
    }
      
  }
  
  
  return(list("x"=sweep.center, "trues"=trues, "ests"=ests))
}

#' @title Computes the normalized root-mean-square error/deviation
#' @export
NRMSE <- function(true, est) {
  # NA are omitted
  if (length(true)!=length(est)) stop("true and est must be vectors of same length")
  d <- na.omit(data.frame(true=true, est=est))
  nrmse <- with(d, 1/mean(true) * sqrt( sum((true-est)**2) /length(true)) )
  return(nrmse)
}

#' @title Computes performance rate statistics: True Positive Rate (TPR=sensitivity), True Negative Rate (TNR=specificity) and False Discovery Rate (FDR=1-precision)
#' @export
perf_rates <- function(CM) {
  # CM is a square matrix
  # rownames are TRUE models
  # colnames are ESTIMATED models
  # for a positive outcome i_x, negative outcomes will include (i_y, for all y!=x)

  models <- gsub("true.", "", rownames(CM))
  positive_outcomes <- models[models!="i0"]
  
  RR <- c()
  for (i in seq_along(positive_outcomes)) {
    pos <- positive_outcomes[i]
    neg <- models[models!=pos]
    
    TP <- sum(CM[paste("true.",pos,sep=""),paste("est.",pos,sep="")])
    FN <- sum(CM[paste("true.",pos,sep=""),colnames(CM)%in%paste("est.",neg,sep="")])
    FP <- sum(CM[rownames(CM)%in%paste("true.",neg,sep=""),paste("est.",pos,sep="")])
    TN <- sum(CM[rownames(CM)%in%paste("true.",neg,sep=""),colnames(CM)%in%paste("est.",neg,sep="")])
    
    TPR <- TP/(TP+FN)
    TNR <- TN/(FP+TN)
    FDR <- FP/(FP+TP)
    
    RR <- rbind(RR, c(pos, paste(neg, collapse=" || "), TPR, TNR, FDR))
  }
  
  # only versus "i0"
  for (i in seq_along(positive_outcomes)) {
    pos <- positive_outcomes[i]
    neg <- "i0"
    
    TP <- sum(CM[paste("true.",pos,sep=""),paste("est.",pos,sep="")])
    FN <- sum(CM[paste("true.",pos,sep=""),colnames(CM)%in%paste("est.",neg,sep="")])
    FP <- sum(CM[rownames(CM)%in%paste("true.",neg,sep=""),paste("est.",pos,sep="")])
    TN <- sum(CM[rownames(CM)%in%paste("true.",neg,sep=""),colnames(CM)%in%paste("est.",neg,sep="")])
    
    TPR <- TP/(TP+FN)
    TNR <- TN/(FP+TN)
    FDR <- FP/(FP+TP)
    
    RR <- rbind(RR, c(pos, paste(neg, collapse="-"), TPR, TNR, FDR))
  }
  
  RR <- as.data.frame(RR, stringsAsFactors=F)
  names(RR) <- c("positive.outcome", "negative.outcomes", "TPR", "TNR", "FDR")
  RR$TPR <- as.numeric(RR$TPR)
  RR$TNR <- as.numeric(RR$TNR)
  RR$FDR <- as.numeric(RR$FDR)
  
  names(RR) <- c("positive.outcome", "negative.outcomes", "TPR.sensitivity", "TNR.specificity", "FDR")
  return(RR)
}

#' @title Detect sweeps and estimate sweep age
#' @param X either the output from simulate_obs || a VCF filepath || a multiSFS
#' @export
summary_cv <- function(X, file) {
  
  # (class(X) == "validationResult") {
  
  pdf(file)
  
  Y <- X$ANALYSIS
  R <- list()
  demes <- names(Y)
  
  D <- data.frame("true.model" = character(),
                  "simID" = integer(),
                  "est.model" = character(),
                  "est.sweepAge" = numeric(),
                  "est.sweepAge.IClow" = numeric(),
                  "est.sweepAge.ICup" = numeric(),
                  "est.sweepPos" = numeric(),
                  "true.sweepPos" = numeric(),
                  "est.sweepStart" = numeric(),
                  "est.sweepEnd" = numeric(),
                  "true.sweepAge" = numeric(),
                  "true.recRate" = numeric(), stringsAsFactors=FALSE)
  for (i in seq_along(Y)) {
    deme <- names(Y)[i]
    print(deme)
    for (j in seq_along(Y[[i]])) {
      
      #if (is.null(Y[[i]][[j]]$density)) next()
      #if (all(is.na(Y[[i]][[j]]$estimation))) next()
      if (is.null(Y[[i]][[j]]$estimation)||nrow(Y[[i]][[j]]$estimation)==0) {
        add <- data.frame("true.model" = deme,
                          "simID" = names(Y[[deme]])[j],
                          "est.model" = "i0",
                          "est.sweepAge" = NA,
                          "est.sweepAge.IClow" = NA,
                          "est.sweepAge.ICup" = NA,
                          "est.sweepPos" = NA,
                          "true.sweepPos" = X$GENERAL$sweepPos,
                          "est.sweepStart" = NA,
                          "est.sweepEnd" = NA,
                          "true.sweepAge" = ifelse(deme=="i0", NA, X$PRIORS[[deme]]$sweepAge[j]),
                          "true.recRate" = ifelse(deme=="i0", NA, X$PRIORS[[deme]]$recRate[j]), stringsAsFactors=F)        
      } else {
        add <- data.frame("true.model" = deme,
                        "simID" = names(Y[[deme]])[j],
                        "est.model" = as.character(Y[[deme]][[j]]$estimation$deme),
                        "est.sweepAge" = Y[[deme]][[j]]$estimation$sweepAge,
                        "est.sweepAge.IClow" = Y[[deme]][[j]]$estimation$sweepAge.IC.low,
                        "est.sweepAge.ICup" = Y[[deme]][[j]]$estimation$sweepAge.IC.up,
                        "est.sweepPos" = Y[[deme]][[j]]$estimation$sweep.center,
                        "true.sweepPos" = X$GENERAL$sweepPos,
                        "est.sweepStart" = Y[[deme]][[j]]$estimation$sweep.lbound,
                        "est.sweepEnd" = Y[[deme]][[j]]$estimation$sweep.rbound,
                        "true.sweepAge" = ifelse(deme=="i0", NA, X$PRIORS[[deme]]$sweepAge[j]),
                        "true.recRate" = ifelse(deme=="i0", NA, X$PRIORS[[deme]]$recRate[j]), stringsAsFactors=F)
      }
      D <- rbind(D, add)
    }
  }
  
  ## remove i0 from estimations for i(!=0), otherwise absurdly inflates the FPR
  #D <- D[-c(which(D$true.model!="i0"&D$est.model=="i0")),]
  
  ##################
  ##################
  # DETECTION

  # absolute confusion matrix
  # take only unique(estimated demes)
  UNIQUE_DEMES <- D[,c("true.model","simID","est.model")]
  UNIQUE_DEMES <- UNIQUE_DEMES[!duplicated(UNIQUE_DEMES),]
  CM <- confusion_matrix(na.omit(data.frame(UNIQUE_DEMES$true.model, UNIQUE_DEMES$est.model)))
  suppressWarnings(barplot(t(CM), 
                           col = brewer.pal(nrow(CM), "Paired"), 
                           xlab = "Simulated models",
                           main = "Absolute confusion matrices (black line = #sims/model)",
                           ylab = "Absolute counts of estimated models",
                           legend = rownames(CM)))
  abline(h = X$GENERAL$nSimul, col = "black", lwd = 2)
  rn <- rownames(CM)
  rownames(CM) <- paste("true.", rn, sep="")
  colnames(CM) <- paste("est.", rn, sep="")
  R$absolute.confusion.matrix <- CM
  
  # performance rates
  R$performance.rates <- McSwan::perf_rates(CM)
  
  # relative confusion matrix
  CM <- confusion_matrix(na.omit(data.frame(UNIQUE_DEMES$true.model, UNIQUE_DEMES$est.model)))
  CM <- t(apply(CM, 1, function(x) x/sum(x)))
  suppressWarnings(barplot(t(CM), 
                           col = brewer.pal(nrow(CM), "Paired"), 
                           xlab = "Simulated models",
                           main = "Relative confusion matrices",
                           ylab = "Fraction of estimated models",
                           legend = rownames(CM)))
  rn <- rownames(CM)
  rownames(CM) <- paste("true.", rn, sep="")
  colnames(CM) <- paste("est.", rn, sep="")
  R$relative.confusion.matrix <- CM
  
  # number of sweeps
  cv <- list()
  demes_no_i0 <- demes[demes!="i0"]
  for (i in seq_along(demes_no_i0)) {
    d <- subset(D, true.model==demes_no_i0[i] & est.model==demes_no_i0[i])
    d <- table(d$simID)
    cv <- c(cv, list(d))
  }
  cvv <- matrix(rep(NA, length(demes_no_i0)*length(0:max(unlist(cv)))), nrow=length(demes_no_i0)); colnames(cvv) <- 0:max(unlist(cv)); rownames(cvv) <- demes_no_i0
  for (i in seq_along(cv)) {
    tmp <- table(cv[[i]])
    cvv[i,names(tmp)] <- tmp/sum(tmp)*100
  }
  suppressWarnings(barplot(cvv, 
                           col = brewer.pal(nrow(cvv), "Paired"), 
                           xlab = "Number of discrete sweeps detected",
                           main = "Discrete sweeps detected when a single sweep was simulated",
                           ylab = "Per-model percentage of sweeps detected",
                           legend = rownames(cvv), beside = F))
  R$number.of.sweeps <- cvv
  
  # distance to true sweep
  cv <- subset(D, true.model!="i0" & est.model==true.model)
  g=ggplot2::ggplot(data=cv, aes(x=abs(est.sweepPos-true.sweepPos), y=..count.., group=true.model, fill=true.model, col=true.model)) + geom_density(alpha=.05) + theme_bw() + xlab("Distance (bp) (log10 scale)") + ylab("Density (counts)") + ggtitle("Density of estimated distances to the beneficial mutation position") + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + annotation_logticks(sides="b", size=.6)
  print(g)
  
  # inclusion stats
  cvv <- simplify2array( by(cv, cv$true.model, function(y) {
    tt <- apply(y[,-c(1:3)], 1, function(x) x['est.sweepEnd'] >= x['true.sweepPos'] && x['est.sweepStart'] <= x['true.sweepPos'])
    return(c("FALSE"=sum(!tt), "TRUE"=sum(tt))) }) )
  cvv <- apply(cvv, 2, function(x) x/sum(x)*100)
  barplot(cvv, xlab = "Simulated models", ylab = "% of sweeps which include the beneficial mutation position", legend = rownames(cvv))
  R$sweep.inclusiveness <- cvv
  
  ##################
  ##################
  # ESTIMATION
  
  tmp <- subset(D, true.model!="i0" & est.model==true.model)
  
  l <- range(c(tmp$true.sweepAge, tmp$est.sweepAge), na.rm=T)
  g <- ggplot2::ggplot(data=tmp, aes(x = true.sweepAge, 
                            y = est.sweepAge, 
                            col = log10(true.recRate))) +
    geom_abline(intercept=0, slope=1, col="gray", size=1, alpha=.6) +
    geom_point(size=1, alpha=.8) +
    theme_bw() + #theme(legend.position = "top") +
    facet_wrap(~true.model) +
    scale_colour_distiller("log10(rho)", palette="Spectral") +
    coord_equal(xlim=c(l[1], l[2]), ylim=c(l[1], l[2])) +
    ggtitle("Parameter estimation performance") +
    xlab("True sweep ages (generations BP, gBP)") +
    ylab("Estimated sweep ages (gBP)")
  print(g)
  
  R$param.estimation <- tmp
  
  # with confidence intervals
  tmp <- subset(D, true.model!="i0" & est.model==true.model)
  
  l <- range(c(tmp$true.sweepAge, tmp$est.sweepAge), na.rm=T)
  g <- ggplot2::ggplot(data=tmp, aes(x = true.sweepAge, 
                            y = est.sweepAge, 
                            col = log10(true.recRate))) +
    geom_segment(aes(x=true.sweepAge, xend=true.sweepAge, y=est.sweepAge.IClow, yend=est.sweepAge.ICup, col=log10(true.recRate)), alpha=.7) +
    geom_abline(intercept=0, slope=1, col="gray", size=1, alpha=.6) +
    geom_point(size=1, alpha=.8) +
    theme_bw() + #theme(legend.position = "top") +
    facet_wrap(~true.model) +
    scale_colour_distiller("log10(rho)", palette="Spectral") +
    coord_equal(xlim=c(l[1], l[2]), ylim=c(l[1], l[2])) +
    ggtitle("Parameter estimation performance with 95% confidence intervals") +
    xlab("True sweep ages (generations BP, gBP)") +
    ylab("Estimated sweep ages (gBP)")
  print(g)
  
  # NRMSE
  prmD <- R$param.estimation
  R$NRMSE <- as.matrix(by(prmD, prmD$true.model, function(x) McSwan::NRMSE(x$true.sweepAge, x$est.sweepAge)))

  dev.off()
  return(R)
}







