

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
#' @export
gscan = function(X, reftb, firstPos = NULL, lastPos = NULL, minSNP = 10, windowSizes = seq(1e4, 2e5, length.out = 20), nSteps = 20, n = NULL) {
	if (class(reftb)!="referenceTable") stop("reftb is not a valid referenceTable")
	if (is.null(reftb$DIMREDUC)) stop("You have not performed the dimension reduction.")
	
	infer.age = TRUE

	# by class
	if (class(X)=="observedDataset") { cat("Scanning an observedDataset object.\n")
		if (reftb$GENERAL$folded != X$folded) stop("Error: folding is not homogeneous across X and reftb.\n")
		Y = scan_core(reftb, POS = X$obsData$POS, PAC = X$obsData$PAC, wSNP = X$obsData$nSNP, firstPos = firstPos, lastPos = lastPos, minSNP = minSNP, windowSizes = windowSizes, nSteps = nSteps, infer.age = infer.age)
	} else if (class(X)=="validationTable") { cat("Scanning a validationTable object.\n")
		if (reftb$GENERAL$folded != X$GENERAL$folded) stop("Error: folding is not homogeneous across X and reftb.\n")
		Y <- lapply(X$SFS, function(x) NA)
		for (d in seq_along(X$SFS)) {
			deme <- names(X$SFS)[d]
			if (is.null(n)) n = length(X$SFS[[deme]]) else if (n>length(X$SFS[[deme]])) stop("n is superior to the actual number of simulations in the validationTable")
			cat("\n\n==========================\n",deme,"\n==========================\n")
			ad = list()
			for (i in 1:n) { cat("\n\n>>> ",deme," ",i,"\n")
				ad = c(ad, list(scan_core(reftb, POS = X$SFS[[deme]][[i]]$obsData$POS, PAC = X$SFS[[deme]][[i]]$obsData$PAC, wSNP = X$SFS[[deme]][[i]]$obsData$nSNP, firstPos = NULL, lastPos = NULL, minSNP = minSNP, windowSizes = windowSizes, nSteps = nSteps, infer.age = infer.age)))
			}
			names(ad) = 1:n
			Y[[deme]] = ad
		}
	} else {
		stop("X is not a valid observedDataset or referenceTable object.")
	}
	
	return(Y)
}

#' @title Core function for the genome scan
#' @keywords internal
scan_core = function(reftb, POS, PAC, wSNP, firstPos, lastPos, minSNP, windowSizes, nSteps, infer.age = TRUE, progressBar = TRUE) {
	if (class(reftb)!="referenceTable") stop("reftb is not a valid referenceTable")
 
	# if requested, contract the scan range
	if (!is.null(firstPos)) {
		inRange = POS >= firstPos
		POS = POS[inRange]
		PAC = PAC[inRange]
		wSNP = wSNP[inRange]
	}
	if (!is.null(lastPos)) {
		inRange = POS <= lastPos
		POS = POS[inRange]
		PAC = PAC[inRange]
		wSNP = wSNP[inRange]
	}
	minPos = min(POS)
	maxPos = max(POS)
	
	# initialization
    STABILITY = matrix(0, nrow=length(POS), ncol=nlevels(reftb$DIMREDUC$LDA$modelIndices))
    colnames(STABILITY) = levels(reftb$DIMREDUC$LDA$modelIndices)
    if (infer.age) {
        EST.AGE = matrix(0, nrow=length(POS), ncol=length(reftb$DIMREDUC$PLS))
        colnames(EST.AGE) <- names(reftb$DIMREDUC$PLS)
    }
    WEIGHTS = STABILITY
	N.TESTS = rep(0, length(POS))
	isFolded = (reftb$GENERAL$folded)

	# force the levels of PAC to the full range of the sfs, in the same order as the reftb sfs
	PAC <- factor(as.character(PAC), levels = colnames(reftb$SFS[[1]]))
	n = 0
	if (is.null(windowSizes)) windowSizes = maxPos - minPos + 100
	if (progressBar) { nBarElmts <- 50; nMax <- length(windowSizes) * nSteps; nLast <- 1; bar <- rep(" ", nBarElmts) }

    for (i in seq_along(windowSizes)) { if (!progressBar) cat("\n",i," ")
        starts = seq(minPos, minPos + windowSizes[i] - 1, length.out = nSteps)
        for (s in seq_along(starts)) { if (!progressBar) cat(".")
			n = n + 1

            # partition into windows
            breaks = seq(starts[s], maxPos, by=windowSizes[i])
            WIN = findInterval(POS, breaks)
            # if NA-full windows, labelling irrelevantly incremented, so we have to remove them:
            WIN = as.numeric(factor(WIN, levels=sort(unique(WIN))))

			if (isFolded) {
				local.sfs = do.call("rbind", by(PAC[wSNP==1], WIN, table)) + (1/2 * do.call("rbind", by(PAC[wSNP==0.5], WIN, table)))
			} else {
				local.sfs <- do.call("rbind", by(PAC, WIN, table))
			}
			nSNP = rowSums(local.sfs)
			local.sfs = sweep(local.sfs, MARGIN = 1, FUN = "/", STATS = nSNP)

            # LDA prediction
            if ("PCA.model" %in% names(reftb$DIMREDUC$LDA)) {
                if (!progressBar) cat(">")
                lslda = sweep(local.sfs, MARGIN = 2, FUN = "-", STATS = reftb$DIMREDUC$LDA$PCA.means)
                lslda = lslda[,reftb$DIMREDUC$LDA$euCols,drop=F]
                lslda = lslda %*% reftb$DIMREDUC$LDA$PCA.model
                LDA = predict(reftb$DIMREDUC$LDA$model, lslda)
            } else {
                LDA <- predict(reftb$DIMREDUC$LDA$model, local.sfs[,reftb$DIMREDUC$LDA$euCols,drop=F])
            }

            lvl = levels(LDA$class)
			classes = LDA$class[WIN]
			inferrableSNP = nSNP[WIN] >= minSNP
			N.TESTS[inferrableSNP] = N.TESTS[inferrableSNP] + 1

			if (infer.age && !progressBar) cat("!")
            for (l in 1:length(lvl)) {
				focus = (classes == lvl[l]) & inferrableSNP
				STABILITY[focus, l] = STABILITY[focus, l] + 1
				tmp = LDA$posterior[WIN, l]
                WEIGHTS[inferrableSNP, l] = WEIGHTS[inferrableSNP, l] + tmp[inferrableSNP]
				# estimate
				if (infer.age && lvl[l] != "i0") { 
					focal.sfs = local.sfs
                    focal.sfs = sweep( sweep(focal.sfs, 2, FUN="-", reftb$DIMREDUC$PLS[[lvl[l]]]$mean), 2, FUN="/", reftb$DIMREDUC$PLS[[lvl[l]]]$sd)
                    focal.sfs = focal.sfs[,reftb$DIMREDUC$PLS[[lvl[l]]]$euCols,drop=F]
                    pls = c(predict(reftb$DIMREDUC$PLS[[lvl[l]]]$model, focal.sfs, type="response", comps=1:reftb$DIMREDUC$PLS[[lvl[l]]]$model$ncomp))
					w = LDA$posterior[, lvl[l]]
                    pls = pls * w
                    pls = pls[WIN]
                    EST.AGE[inferrableSNP, lvl[l]] = EST.AGE[inferrableSNP, lvl[l]] + pls[inferrableSNP]
				}
            }

			if (progressBar) {
				if (floor(n/nMax*length(bar)) >= nLast) {
					nLast <- floor(n/nMax*length(bar))
					bar[1:nLast] <- "="
					if (n/nMax!=1) bar[nLast+1] <- ">"
					cat("[",paste0(bar, collapse=""),"]","\n", sep="")
				}
			}
        } #s
    } #i

    if (infer.age) EST.AGE = EST.AGE / WEIGHTS[,-1] # must be done before averaging the WEIGHTS
    STABILITY = sweep(STABILITY, MARGIN = 1, FUN = "/", STATS = N.TESTS)
    WEIGHTS = sweep(WEIGHTS, MARGIN = 1, FUN = "/", STATS = N.TESTS)

	return(list("pos" = POS,
	"stability" = STABILITY,
	"postpr" = WEIGHTS,
	"param" = EST.AGE,
	"inferentiality" = N.TESTS / n))
}


