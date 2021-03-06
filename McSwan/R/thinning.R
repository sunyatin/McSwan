
#' @title Merge contiguous swept loci (tiling path inference)
#' @description This function merges contiguous regions detected as having experienced a selective sweep (i.e. following a \code{\link{gscan}}).
#' @param scanResult a list returned by \code{\link{gscan}}
#' @param reftb an initialized \code{referenceTable} object
#' @param X the observed or pseudo-observed dataset used in the \code{\link{gscan}} call
#' @param max_L (integer) the "horizon distance" (in base pairs), ie. the maximum distance from any position above which loci detected as under selection cannot be assumed to be related to the same sweep region
#' @param signif.threshold (single or vector of numeric values) (between 0 and 1) a selection score cut-off above which we decide that a SNP is positive selected; if unique value, this cutoff will be used for all populations; alternatively, you can give a vector of cutoff values which will be specifically used for each population (e.g. signif.threshold=c(0.5, 0.2))
#' @param maxIter (integer) maximum number of iterations used for the tiling path inference
#' @param stat (function) a statistic used to center the SNP-wise sweep age estimates (\code{getmode} was shown to be the least biased but you may also try \code{median} or \code{mean})
#' @return Returns a dataframe of the estimated sweep contigs. One sweep region per line, and in columns:
#' \itemize{ 
#' \item \code{\bold{sweep.center}} the center of the sweep region contig
#' \item \code{\bold{sweep.lbound}} the first position of the sweep region contig
#' \item \code{\bold{sweep.rbound}} the last position of the sweep region contig
#' \item \code{\bold{RPP}} ratio of LDA-derived posterior probabilities: highest ratio of pp(is)/pp(i0) among all the SNPs within the sweep region, where pp(is) is the integrated posterior probability of the selective model is at a focal SNP and pp(i0) the integrated posterior probability of the neutral model at the same SNP
#' \item \code{\bold{deme}} the population detected to have experienced a selective sweep (given as "\emph{i}" and the index of the population in the \emph{MS} command)
#' \item \code{\bold{sweepAge}} the point estimate for the sweep age
#' \item \code{\bold{sweepAge.IC.low}} lower boundary of the 95\% confidence interval for the sweep age estimation
#' \item \code{\bold{sweepAge.IC.up}} upper boundary of the 95\% confidence interval for the sweep age estimation
#' \item \code{\bold{sweepAge.postDistrib}} posterior distribution of the sweep age estimation
#' }
#' @seealso \code{\link{gscan}}, \code{\link{summary}}, \code{\link{summary.validationResult()}}
#' @export
#' @examples Please refer to the vignette.
thin <- function(scanResult, reftb, X, max_L = 1e6, signif.threshold = .5, maxIter = 1000, stat = getmode) {
	output_region_wise_estimates = FALSE
	discard_extraRange = TRUE
	weighted_metaEstimate = TRUE

	if (weighted_metaEstimate) cat("Meta-estimates are weighted by the SNP-wise posterior probability of sweep model.\n")
	if (length(signif.threshold)>1) {
		names(signif.threshold) = paste0("i",seq_along(signif.threshold))
		cat("Multiple significance thresholds:",signif.threshold,"\n")
	}
	if (length(signif.threshold)==1) cat("Single significance threshold.\n")

	if (class(X)=="observedDataset") {
		if (any(sapply(scanResult[[1]], is.list))) stop("scanResult must be a result obtained from an observedDataset object.")
		R <- c()
		demes4contig = levels(reftb$DIMREDUC$LDA$modelIndices); demes4contig = demes4contig[demes4contig!="i0"]
		for (dd in seq_along(demes4contig)) {
			ddeme <- demes4contig[dd]
			cat("\n==========================",ddeme,"==========================\n")
			sth = ifelse(length(signif.threshold)==1, signif.threshold, signif.threshold[ddeme])
			### >>>
			II <- findContigs(scanResult = scanResult, reftb = reftb, POS = X$obsData$POS, PAC = X$obsData$PAC, wSNP = X$obsData$nSNP, deme = ddeme, max_L = max_L, signif.threshold = sth, maxIter = maxIter, stat = stat, discard_extraRange = discard_extraRange, weighted_metaEstimate = weighted_metaEstimate)
			R = rbind(R, II)
		}
		if (nrow(R)>=1) rownames(R) <- 1:nrow(R)
		if (output_region_wise_estimates==FALSE) {
			R$BFSynth = NULL
			R$sweepAgeSynth = NULL
			names(R)[which(names(R)=="BF")] <- "RPP"
		}
		return(R)
	} else if (class(X)=="validationTable") {
		if (!any(sapply(scanResult[[1]], is.list))) stop("scanResult must be a result obtained from a validationTable object.")
		ANALYSIS = scanResult
		for (d in seq_along(scanResult)) {
			deme <- names(scanResult)[d]
			cat("\n\n==========================\n",deme,"\n==========================\n")
			for (i in seq_along(scanResult[[deme]])) {
				#cat("\n\n>>> ",deme," ",i,"\n")
				I = scanResult[[deme]][[i]]
				demes4contig = names(scanResult); demes4contig = demes4contig[demes4contig!="i0"]
				R <- c()
				for (dd in seq_along(demes4contig)) {
					ddeme <- demes4contig[dd]
					sth = ifelse(length(signif.threshold)==1, signif.threshold, signif.threshold[ddeme])
					### >>>
					 capture.output({ II = findContigs(scanResult = I, reftb = reftb, POS = X$SFS[[deme]][[i]]$obsData$POS, PAC = X$SFS[[deme]][[i]]$obsData$PAC, wSNP = X$SFS[[deme]][[i]]$obsData$nSNP, deme = ddeme, max_L = max_L, signif.threshold = sth, maxIter = maxIter, stat = stat, discard_extraRange = discard_extraRange, weighted_metaEstimate = weighted_metaEstimate) })
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
		if (output_region_wise_estimates==FALSE) {
			ANALYSIS$BFSynth = NULL
			ANALYSIS$sweepAgeSynth = NULL
			names(R)[which(names(R)=="BF")] <- "RPP"
		}
		return(ANALYSIS)
	} else {
		stop("Object x is not an observedDataset or validationTable object.")
	}

}

findContigs = function(scanResult, reftb, POS, PAC, wSNP, deme, max_L, signif.threshold, maxIter = 1000, stat = getmode, discard_extraRange = TRUE, weighted_metaEstimate = TRUE) {
	if (!is.character(deme)) stop("deme must be a string")

	stat_BF <- max # max better represents the pr of the very core (ie at the very proximity of the advantageous mutation)

    v = if (deme != "i0") data.frame(pos = scanResult$pos, 
	n = scanResult$stability[,deme], 
	a = scanResult$param[,deme],
	w0 = scanResult$postpr[,"i0"],
	w = scanResult$postpr[,deme]) else data.frame(pos = scanResult$pos, n = scanResult$stability[,deme])

	# note that n could be scanResult$postpr[deme]
	
	RES <- data.frame(sweep.center = integer(0),
					  sweep.lbound = integer(0),
					  sweep.rbound = integer(0),
					  BF = numeric(0),
					  BFSynth = numeric(0),
					  deme = character(0),
					  sweepAge = numeric(0),
					  sweepAgeSynth = numeric(0),
					  sweepAge.IC.low = numeric(0),
					  sweepAge.IC.up = numeric(0),
					  nBaseSNPs = integer(0))

    Kontigs <- c()
	#postDistrib <- list()
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
			subv <- v[c(ppLeft:ppRight),,drop=F]
			if (discard_extraRange) {
				Range = unlist(reftb$GENERAL$sweepAgeDistrib[[deme]][2:3])
				subv[!is.na(subv[,"a"]) & (subv[,"a"]<Range[1] | subv[,"a"]>Range[2]),"a"] <- NA
			}
			subv <- na.omit(subv[, c("a","w"), drop=F])
			subv$w <- subv$w / sum(subv$w)
			nsnp <- nrow(subv)
			if (nsnp==0) {
				central_estimate = NA
				IC_estimate = c(NA, NA)
				#postDistrib = c(postDistrib, list(NA))
			} else {
				central_estimate <- stat(subv$a, na.rm = do_NArm, weights = subv$w)
				IC_estimate <- c(coef(quantreg::rq(subv$a~1, tau = c(.025, .975), weights = subv$w)))
				#postDistrib <- c(postDistrib, list(subv$a)) # should be adjusted, ideally
			}
		} else {
			nsnp <- sum(!is.na(v[c(ppLeft:ppRight),"a"]))
			central_estimate <- stat(v[c(ppLeft:ppRight),"a"], na.rm = do_NArm)
			IC_estimate <- quantile(v[c(ppLeft:ppRight),"a"], c(.025, .975), na.rm = do_NArm)
			#postDistrib <- c(postDistrib, list(v[c(ppLeft:ppRight),"a"]))
		}
		BF <- stat_BF( v[c(ppLeft:ppRight),"w"] / v[c(ppLeft:ppRight),"w0"], na.rm = TRUE )
        ad <- c(v[c(ppLeft, ppRight),"pos"], central_estimate, IC_estimate, nsnp, BF)
        Kontigs <- rbind(Kontigs, ad)
        # hide
        #v[c(ppLeft:ppRight),"n"] <- NA
		v[c(ppLeft:ppRight),c("n","w","w0","a")] <- NA
    }

    if (is.null(Kontigs) || length(Kontigs)==0) { 
	#RES$sweepAge.postDistrib <- list(); 
	return( RES ) }

    Kontigs <- as.data.frame(Kontigs)
    names(Kontigs) <- c("pos_left", "pos_right", "SNP_estim", "IC_low", "IC_up", "nBaseSNPs", "BF")
    Kontigs$pos_middle = with(Kontigs, 1/2*(pos_left+pos_right))
    Kontigs$length = with(Kontigs, pos_right-pos_left)
    Kontigs = Kontigs[order(Kontigs$pos_left),]

	if (FALSE) print(Kontigs)
	
	# contig-wise age inference
	Kontigs$meta_estim = NA; Kontigs$BFSynth = NA
	cat("\n")
	if (deme != "i0") {
		for (i in 1:nrow(Kontigs)) {
			ad = scan_core(reftb, POS, PAC, wSNP, firstPos = Kontigs[i,"pos_left"], lastPos = Kontigs[i,"pos_right"], minSNP = 0, windowSizes = NULL, nSteps = 1, infer.age = TRUE, discard_extraRange = discard_extraRange)
			# meta_estim
			vl <- unique((ad$param)[,deme])
			if (length(vl)!=1) stop("error in meta-estimation: non-unique value (probably due to the presence of a NA)")
			Kontigs[i,"meta_estim"] = vl
			# BF
			vl <- unique((ad$postpr)[,deme])
			vl_i0 <- unique((ad$postpr)[,"i0"])
			if (length(vl)!=1) stop("error in meta-postpr estimate: non-unique vl (probably due to the presence of a NA)")
			if (length(vl_i0)!=1) stop("error in meta-postpr estimate: non-unique vl_i0 (probably due to the presence of a NA)")
			Kontigs[i,"BFSynth"] = vl / vl_i0 # proxy of Bayes Factor, as postpr(sweep_model)/postpr(neutral_model)
		}
	}
	
	RES <- rbind(RES, data.frame(
						sweep.center = Kontigs$pos_middle,
						sweep.lbound = Kontigs$pos_left,
						sweep.rbound = Kontigs$pos_right,
						BF = Kontigs$BF,
						BFSynth = Kontigs$BFSynth,
						deme = deme,
						sweepAge = Kontigs$SNP_est,
						sweepAgeSynth = Kontigs$meta_estim,
						sweepAge.IC.low = Kontigs$IC_low,
						sweepAge.IC.up = Kontigs$IC_up,
						nBaseSNPs = Kontigs$nBaseSNPs))
	#RES$sweepAge.postDistrib <- postDistrib
	rownames(RES) <- 1:nrow(RES)

    return(RES)
}
  
