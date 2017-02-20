
#' @title Computes the mean-normalized (or standardized) root mean square error (or deviation)
#' @param true a vector of true numeric values
#' @param est a vector of estimated numeric values
#' @return A numeric value corresponding to the NRMSE.
#' @references Walther, B.A. & Moore, J.L. (2005). The concepts of bias, precision and accuracy, and their use in testing the performance of species richness estimators, with a literature review of estimator performance. Ecography, 28, 815-829.
#' @keywords internal
NRMSE <- function(true, est) {
  # NA are omitted
  if (length(true)!=length(est)) stop("true and est must be vectors of same length")
  d <- na.omit(data.frame(true=true, est=est))
  nrmse <- with(d, 1/mean(true) * sqrt( sum((true-est)**2) /length(true)) )
  return(nrmse)
}


#' @title Computes performance statistics: True Positive Rate (TPR=sensitivity), True Negative Rate (TNR=specificity) and False Discovery Rate (FDR=1-precision)
#' @param CM a confusion matrix
#' @return A dataframe of performance statistics.
#' @export
#' @keywords internal
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


#' @title Compute a confusion matrix
#' @param cvd a matrix of strings with two columns: (1) the true model; (2) the predicted model
#' @return A confusion matrix.
#' @keywords internal
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


#' @title Summary of the sliding validation
#' @param X a \code{validationTable} object with non-empty \code{ANALYSIS} slot (ie. obtained after executing the \code{sliding_validation} function)
#' @param file the name of a file where validation graphics are plotted
#' @import ggplot2
#' @return A list of 6 elements: \itemize{ 
#' \item \code{absolute.confusion.matrix} confusion matrix with absolute counts
#' \item \code{absolute.relative.matrix} confusion matrix with absolute counts
#' \item \code{performance rates} specificity (TNR), sensitivity (TPR) and false discovery rate (FDR)
#' \item \code{number.of.sweeps} number of contiguous sweeps detected in each simulated genomic fragment
#' \item \code{sweep.inclusiveness} number of times the detected sweeps encompass the beneficial mutation
#' \item \code{NRMSE} normalized/standardized root mean square error
#' }
#' @export
#' @seealso \code{\link{sliding_validation}}
#' @examples Please refer to the vignette.
summary.validationResult <- function(X, valtb, file) {
  
  # (class(X) == "validationResult") {
  if (class(valtb)!="validationTable") stop("valtb is not a valid validationTable object.")
  #if (!"ANALYSIS"%in%names(X)) stop("summary() can only be applied to validationTable with existing ANALYSIS slot, ie after executing the sliding_validation() funciton")
  
  pdf(file)
  
  Y <- X #X$ANALYSIS
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
    cat(deme,"\n")
    for (j in seq_along(Y[[i]])) {
      
      #if (is.null(Y[[i]][[j]]$density)) next()
      #if (all(is.na(Y[[i]][[j]]$estimation))) next()
      if (is.null(Y[[i]][[j]])||nrow(Y[[i]][[j]])==0) {
        add <- data.frame("true.model" = deme,
                          "simID" = names(Y[[deme]])[j],
                          "est.model" = "i0",
                          "est.sweepAge" = NA,
                          "est.sweepAge.IClow" = NA,
                          "est.sweepAge.ICup" = NA,
                          "est.sweepPos" = NA,
                          "true.sweepPos" = valtb$GENERAL$sweepPos,
                          "est.sweepStart" = NA,
                          "est.sweepEnd" = NA,
                          "true.sweepAge" = ifelse(deme=="i0", NA, valtb$PRIORS[[deme]]$sweepAge[j]),
                          "true.recRate" = ifelse(deme=="i0", NA, valtb$PRIORS[[deme]]$recRate[j]), stringsAsFactors=F)        
      } else {
        add <- data.frame("true.model" = deme,
                        "simID" = names(Y[[deme]])[j],
                        "est.model" = as.character(Y[[deme]][[j]]$deme),
                        "est.sweepAge" = Y[[deme]][[j]]$sweepAge,
                        "est.sweepAge.IClow" = Y[[deme]][[j]]$sweepAge.IC.low,
                        "est.sweepAge.ICup" = Y[[deme]][[j]]$sweepAge.IC.up,
                        "est.sweepPos" = Y[[deme]][[j]]$sweep.center,
                        "true.sweepPos" = valtb$GENERAL$sweepPos,
                        "est.sweepStart" = Y[[deme]][[j]]$sweep.lbound,
                        "est.sweepEnd" = Y[[deme]][[j]]$sweep.rbound,
                        "true.sweepAge" = ifelse(deme=="i0", NA, valtb$PRIORS[[deme]]$sweepAge[j]),
                        "true.recRate" = ifelse(deme=="i0", NA, valtb$PRIORS[[deme]]$recRate[j]), stringsAsFactors=F)
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
                           col = RColorBrewer::brewer.pal(nrow(CM), "Paired"), 
                           xlab = "Simulated models",
                           main = "Absolute confusion matrices (black line = #sims/model)",
                           ylab = "Absolute counts of estimated models",
                           legend = rownames(CM)))
  abline(h = valtb$GENERAL$nSimul, col = "black", lwd = 2)
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
                           col = RColorBrewer::brewer.pal(nrow(CM), "Paired"), 
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
  maxCV <- unlist(cv)
  emptyCV <- ifelse(length(maxCV)==0, TRUE, FALSE)
  maxCV <- ifelse(emptyCV, 0, max(maxCV))
  cvv <- matrix(rep(NA, length(demes_no_i0)*length(0:maxCV)), nrow=length(demes_no_i0))
  colnames(cvv) <- 0:maxCV
  rownames(cvv) <- demes_no_i0
  for (i in seq_along(cv)) {
    tmp <- table(cv[[i]])
    cvv[i,names(tmp)] <- tmp/sum(tmp)*100
  }
  suppressWarnings(barplot(cvv, 
                           col = RColorBrewer::brewer.pal(nrow(cvv), "Paired"), 
                           xlab = "Number of discrete sweeps detected",
                           main = "Discrete sweeps detected when a single sweep was simulated",
                           ylab = "Per-model percentage of sweeps detected",
                           legend = rownames(cvv), beside = F))
  R$number.of.sweeps <- cvv
  
  # distance to true sweep
  cv <- subset(D, true.model!="i0" & est.model==true.model)
  g=ggplot2::ggplot(data=cv, aes(x=abs(est.sweepPos-true.sweepPos), y=..count.., group=true.model, fill=true.model, col=true.model)) + 
  geom_density(alpha=.05) + theme_bw() + xlab("Distance (bp) (log10 scale)") + ylab("Density (counts)") + 
  ggtitle("Density of estimated distances to the beneficial mutation position") + 
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + 
  scale_x_log10() +
  annotation_logticks(sides="b", size=.6)
  print(g)
  
  # inclusion stats
  if (emptyCV) { 
	cvv <- NA 
  } else {
	  cvv <- simplify2array( by(cv, cv$true.model, function(y) {
		tt <- apply(y[,-c(1:3)], 1, function(x) x['est.sweepEnd'] >= x['true.sweepPos'] && x['est.sweepStart'] <= x['true.sweepPos'])
		return(c("FALSE"=sum(!tt), "TRUE"=sum(tt))) }) )
	  cvv <- apply(cvv, 2, function(x) x/sum(x)*100)
	  barplot(cvv, xlab = "Simulated models", ylab = "% of sweeps which include the beneficial mutation position", legend = rownames(cvv))
  }
  R$sweep.inclusiveness <- cvv
  
  ##################
  ##################
  # ESTIMATION
  
  tmp <- subset(D, true.model!="i0" & est.model==true.model)
  
  if (nrow(tmp)>0) {
  
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
		R$NRMSE <- as.matrix(by(prmD, prmD$true.model, function(x) NRMSE(x$true.sweepAge, x$est.sweepAge)))

  } else {
	R$NRMSE <- NA
  }
  
  
  dev.off()
  return(R)
}




