
#' @title Merge contiguous sweeps
#' @description This function merges contiguous regions detected by \code{\link{gscan}} as having experienced a selective sweep.
#' @param x a dataframe returned by \code{\link{gscan}}
#' @param reftb an initialized \code{referenceTable} object
#' @param bwadjust (numeric) bandwidth coefficient, the algorithm will extend contiguity as far as \code{bwadjust * windowSize} from any focal window
#' @param minWindows (integer) the minimum number of contiguous windows required to perform contiguity search
#' @param summary_stat (string) the central statistic to average the parameter estimates of merged contiguous windows ("mean" recommended)
#' @param plot_thinning (logical) plot a representation of the genomic fragment with the contiguous sweep regions
#' @return Returns a dataframe of the merged contiguous sweeps. One sweep per line, and in columns:
#' \itemize{ 
#' \item \code{\bold{sweep.center}} center of the contiguous sweep region
#' \item \code{\bold{sweep.lbound}} first position of the contiguous sweep region
#' \item \code{\bold{sweep.rbound}} last position of the contiguous sweep region
#' \item \code{\bold{BF}} aggregated Bayes Factor of the contiguous sweep region
#' \item \code{\bold{deme}} population detected to have experienced a selective sweep (given as "i" & index of the population in the \emph{MS} command)
#' \item \code{\bold{sweepAge}} posterior estimate for the sweep age
#' \item \code{\bold{sweepAge.IC.low}} lower boundary of the 95% credible interval for the sweep age estimation
#' \item \code{\bold{sweepAge.IC.up}} upper boundary of the 95% credible interval for the sweep age estimation
#' }
#' @seealso \code{\link{gscan}}
#' @export
#' @examples Please refer to the vignette.
thin <- function(scanResult, reftb, X, discard_extraRange = TRUE, max_L = 1e6, signif.threshold = .5, maxIter = 1000, stat = mean, weighted_metaEstimate = TRUE) {
	if (weighted_metaEstimate) cat("Meta-estimates are weighted by the SNP-wise posterior probability of sweep model.\n")

	if (class(X)=="observedDataset") {
		if (any(sapply(scanResult[[1]], is.list))) stop("scanResult must be a result obtained from an observedDataset object.")
		R <- c()
		demes4contig = levels(reftb$DIMREDUC$LDA$modelIndices); demes4contig = demes4contig[demes4contig!="i0"]
		for (dd in seq_along(demes4contig)) {
			ddeme <- demes4contig[dd]
			cat("\n==========================",ddeme,"==========================\n")
			II <- findContigs(scanResult = scanResult, reftb = reftb, POS = X$obsData$POS, PAC = X$obsData$PAC, wSNP = X$obsData$nSNP, deme = ddeme, max_L = max_L, signif.threshold = signif.threshold, maxIter = maxIter, stat = stat, discard_extraRange = discard_extraRange, weighted_metaEstimate = weighted_metaEstimate)
			R = rbind(R, II)
		}
		if (nrow(R)>=1) rownames(R) <- 1:nrow(R)
		return(R)
	} else if (class(X)=="validationTable") {
		if (!any(sapply(scanResult[[1]], is.list))) stop("scanResult must be a result obtained from a validationTable object.")
		ANALYSIS = scanResult
		for (d in seq_along(scanResult)) {
			deme <- names(scanResult)[d]
			cat("\n\n==========================\n",deme,"\n==========================\n")
			for (i in seq_along(scanResult[[deme]])) { cat("\n>")
				cat("\n\n>>> ",deme," ",i,"\n")
				I = scanResult[[deme]][[i]]
				demes4contig = names(scanResult); demes4contig = demes4contig[demes4contig!="i0"]
				R <- c()
				for (dd in seq_along(demes4contig)) {
					ddeme <- demes4contig[dd]
					II <- findContigs(scanResult = I, reftb = reftb, POS = X$SFS[[deme]][[i]]$obsData$POS, PAC = X$SFS[[deme]][[i]]$obsData$PAC, wSNP = X$SFS[[deme]][[i]]$obsData$nSNP, deme = ddeme, max_L = max_L, signif.threshold = signif.threshold, maxIter = maxIter, stat = stat, discard_extraRange = discard_extraRange, weighted_metaEstimate = weighted_metaEstimate)
					#true.age = if (deme=="i0") NA else X$PRIORS[[deme]]$sweepAge[i]
					#if (!is.na(II[1,1])) R = rbind(R, data.frame(true.model = deme, run = i, est.model = II$deme, true.age = true.age, II))
					R = rbind(R, II)
				}
				if (nrow(R)>=1) rownames(R) <- 1:nrow(R)
				ANALYSIS[[deme]][[i]] = R
			}
		}
		#R$deme = NULL
		#R$OK = R$true.model==R$est.model
		class(ANALYSIS) <- "validationResult"
		return(ANALYSIS)
	} else {
		stop("Object x is not an observedDataset or validationTable object.")
	}

}

findContigs = function(scanResult, reftb, POS, PAC, wSNP, deme, max_L, signif.threshold, maxIter = 1000, stat = mean, discard_extraRange = TRUE, weighted_metaEstimate = TRUE) {
	if (!is.character(deme)) stop("deme must be a string")

    v = if (deme != "i0") data.frame(pos = scanResult$pos, 
	n = scanResult$stability[,deme], 
	a = scanResult$param[,deme],
	w = scanResult$postpr[,deme]) else data.frame(pos = scanResult$pos, n = scanResult$stability[,deme])

	# note that n could be scanResult$postpr[deme]
	
	RES <- data.frame(sweep.center = integer(0),
					  sweep.lbound = integer(0),
					  sweep.rbound = integer(0),
					  BF = numeric(0),
					  deme = character(0),
					  sweepAge = numeric(0),
					  sweepAgeSynth = numeric(0),
					  sweepAge.IC.low = numeric(0),
					  sweepAge.IC.up = numeric(0),
					  nBaseSNPs = integer(0))

    Kontigs <- c()
    while (1) { cat(".")
        line <- which.max(v$n)

        # left
        iter <- 0; pp <- line
        while (iter < maxIter) {
            iter <- iter + 1
            poz <- v[pp,"pos"]
            inHorizon <- subset(v, pos>=(poz-max_L) & pos<poz & n>signif.threshold)
            if (nrow(inHorizon)==0) break
            pp <- as.integer(rownames(inHorizon)[which.min(inHorizon$pos)])
        }
        ppLeft <- pp

        # right
        iter <- 0; pp <- line
        while (iter < maxIter) {
            iter <- iter + 1
            poz <- v[pp,"pos"]
            inHorizon <- subset(v, pos<=(poz+max_L) & pos>poz & n>signif.threshold)
            if (nrow(inHorizon)==0) break
            pp <- as.integer(rownames(inHorizon)[which.max(inHorizon$pos)])
        }
        ppRight <- pp

        if (length(ppLeft)==0 || abs(ppLeft-ppRight) == 0) break

        # append
		do_NArm <- TRUE # important, sometimes weight-scaling can lead to NAs // added Feb 21st 2017
		if (weighted_metaEstimate) {
			subv <- na.omit(v[c(ppLeft:ppRight), c("a","w"), drop=F])
			subv$w <- subv$w / sum(subv$w)
			nsnp <- nrow(subv)
			central_estimate <- stat(subv$a, na.rm = do_NArm, weights = subv$w)
			IC_estimate <- c(coef(rq(subv$a~1, tau = c(.025, .975), weights = subv$w)))
		} else {
			nsnp <- sum(!is.na(v[c(ppLeft:ppRight),"a"]))
			central_estimate <- stat(v[c(ppLeft:ppRight),"a"], na.rm = do_NArm)
			IC_estimate <- quantile(v[c(ppLeft:ppRight),"a"], c(.025, .975), na.rm = do_NArm)
		}
        ad <- c(v[c(ppLeft, ppRight),"pos"], central_estimate, IC_estimate, nsnp)
        Kontigs <- rbind(Kontigs, ad)

        # substract
        v[c(ppLeft:ppRight),"n"] <- NA
    }

    if (is.null(Kontigs) || length(Kontigs)==0) return( RES )

    Kontigs <- as.data.frame(Kontigs)
    names(Kontigs) <- c("pos_left", "pos_right", "SNP_estim", "IC_low", "IC_up", "nBaseSNPs")
    Kontigs$pos_middle = with(Kontigs, 1/2*(pos_left+pos_right))
    Kontigs$length = with(Kontigs, pos_right-pos_left)
    Kontigs = Kontigs[order(Kontigs$pos_left),]

	if (FALSE) print(Kontigs)
	
	# contig-wise age inference
	Kontigs$meta_estim = NA; Kontigs$pp = NA
	cat("\n")
	if (deme != "i0") {
		for (i in 1:nrow(Kontigs)) {
			ad = scan_core(reftb, POS, PAC, wSNP, firstPos = Kontigs[i,"pos_left"], lastPos = Kontigs[i,"pos_right"], minSNP = 0, windowSizes = NULL, nSteps = 1, infer.age = TRUE, discard_extraRange = discard_extraRange)
			# meta_estim
			vl <- unique((ad$param)[,deme])
			if (length(vl)!=1) stop("error in meta-estimation: non-unique value (probably due to the presence of a NA)")
			Kontigs[i,"meta_estim"] = vl
			# pp
			vl <- unique((ad$postpr)[,deme])
			if (length(vl)!=1) stop("error in meta-postpr estimate: non-unique (value) (probably due to the presence of a NA)")
			Kontigs[i,"pp"] = vl
		}
	}
	
	RES <- rbind(RES, data.frame(
						sweep.center = Kontigs$pos_middle,
						sweep.lbound = Kontigs$pos_left,
						sweep.rbound = Kontigs$pos_right,
						BF = Kontigs$pp,
						deme = deme,
						sweepAge = Kontigs$SNP_est,
						sweepAgeSynth = Kontigs$meta_estim,
						sweepAge.IC.low = Kontigs$IC_low,
						sweepAge.IC.up = Kontigs$IC_up,
						nBaseSNPs = Kontigs$nBaseSNPs))
	rownames(RES) <- 1:nrow(RES)

    return(RES)
}
  
