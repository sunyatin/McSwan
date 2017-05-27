

#' @title Ensemble genome scan for sweep detection & sweep age estimation
#' @param X an observed dataset obtained with \code{get_SFS()} or a pseudo-observed dataset obtained with \code{generate_pseudoobs()}
#' @param reftb a \code{referenceTable} object, with non-empty \code{PRIORS}, \code{SFS} and \code{DIMREDUC} slots
#' @param firstPos (integer) the first position of the genome scan (NULL to automatically start at the first known position)
#' @param lastPos (integer) the last position of the genome scan (NULL to automatically stop at the last known position)
#' @param windowSizes (array of integers) a vector of window lengths over which genome scans will be performed iteratively, you may consider using the R function seq(minimum_length, maximum_length, length.out = number_of_values) to sample a set of equally-spaced lengths
#' @param nSteps (integer) number of overlapping shifts (the higher the number, the more overlapping windows will be introduced, increasing the scan resolution)
#' @param discard_extraRange (logical) if TRUE, will discard sweep age estimations going beyond the prior range (TRUE is recommended)
#' @param minSNP (integer) minimum number of SNPs in the window to consider it valid for inference
#' @return A list to be further processed with the \code{thin()} function.
#' @seealso \code{\link{get_SFS}}, \code{\link{generate_pseudoobs}}, \code{\link{thin}}
#' @export
gscan = function(X, reftb, firstPos = NULL, lastPos = NULL, minSNP = 10, windowSizes = seq(1e4, 2e5, length.out = 20), nSteps = 20, discard_extraRange = TRUE) {
	n = NULL # if n = integer, only 1:n POD simulations per model will be gscanned

	if (class(reftb)!="referenceTable") stop("reftb is not a valid referenceTable")
	if (is.null(reftb$DIMREDUC)) stop("You have not performed the dimension reduction.")
	if (minSNP < 1) stop("minSNP must be >= 1")
	if (!is.null(firstPos) && !is.null(lastPos) && lastPos <= firstPos) stop("lastPos must be strictly greater than firstPos")
	if (discard_extraRange) cat("Out-of-range parameter estimates will be discarded.\n")
 
	cat("\n")

	infer.age = TRUE
	windowSizes = sapply(windowSizes, as.integer)

	# by class
	if (class(X)=="observedDataset") { cat("Scanning an observedDataset object.\n")
		if (reftb$GENERAL$folded != X$folded) stop("Error: folding is not homogeneous across X and reftb.\n")
		if (!all.equal(colnames(reftb$SFS[[1]]), names(X$template))) stop("SFS bins are not the same between X and reftb")

		# go
		Y = scan_core(reftb, POS = X$obsData$POS, PAC = X$obsData$PAC, wSNP = X$obsData$nSNP, firstPos = firstPos, lastPos = lastPos, minSNP = minSNP, windowSizes = windowSizes, nSteps = nSteps, infer.age = infer.age, discard_extraRange = discard_extraRange)
	} else if (class(X)=="validationTable") { cat("Scanning a validationTable object.\n")
		if (reftb$GENERAL$folded != X$GENERAL$folded) stop("Error: folding is not homogeneous across X and reftb.\n")
		if (!all.equal(colnames(reftb$SFS[[1]]), names(X$SFS[[1]][[1]]$template))) stop("SFS bins are not the same between X and reftb")

		# go
		Y <- lapply(X$SFS, function(x) NA)
		for (d in seq_along(X$SFS)) {
			deme <- names(X$SFS)[d]
			if (is.null(n)) n = length(X$SFS[[deme]]) else if (n>length(X$SFS[[deme]])) stop("n is superior to the actual number of simulations in the validationTable")
			cat("\n\n==========================\n",deme,"\n==========================\n")
			ad = list()
			for (i in 1:n) { cat("\n\n>>> ",deme," ",i,"\n")
				ad = c(ad, list(scan_core(reftb, POS = X$SFS[[deme]][[i]]$obsData$POS, PAC = X$SFS[[deme]][[i]]$obsData$PAC, wSNP = X$SFS[[deme]][[i]]$obsData$nSNP, firstPos = NULL, lastPos = NULL, minSNP = minSNP, windowSizes = windowSizes, nSteps = nSteps, infer.age = infer.age, discard_extraRange = discard_extraRange)))
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
scan_core = function(reftb, POS, PAC, wSNP, firstPos, lastPos, minSNP, windowSizes, nSteps, infer.age = TRUE, discard_extraRange = TRUE, progressBar = TRUE) {
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
	if (discard_extraRange == TRUE) PARAM.WEIGHTS = WEIGHTS[,-1]
	N.TESTS = rep(0, length(POS))
	isFolded = (reftb$GENERAL$folded)

	# force the levels of PAC to the full range of the sfs, in the same order as the reftb sfs bins
	PAC <- factor(as.character(PAC), levels = colnames(reftb$SFS[[1]]))
	n = 0
	if (is.null(windowSizes)) windowSizes = maxPos - minPos + 100
	if (progressBar) { nBarElmts <- 100; nMax <- length(windowSizes) * nSteps; nLast <- 1; bar <- rep(" ", nBarElmts); cat("1% ") }

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
			local.sfs = sweep(local.sfs, MARGIN = 1, FUN = "/", STATS = replace(nSNP, nSNP==0, 1)) # replace() to avoid dividing by 0

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
				STABILITY[focus, lvl[l]] = STABILITY[focus, lvl[l]] + 1
				wgts = LDA$posterior[WIN, lvl[l]]
                WEIGHTS[inferrableSNP, lvl[l]] = WEIGHTS[inferrableSNP, lvl[l]] + wgts[inferrableSNP]
				# estimate
				if (infer.age && lvl[l] != "i0") { 
					focal.sfs = local.sfs
                    focal.sfs = sweep( sweep(focal.sfs, 2, FUN="-", reftb$DIMREDUC$PLS[[lvl[l]]]$mean), 2, FUN="/", reftb$DIMREDUC$PLS[[lvl[l]]]$sd)
                    focal.sfs = focal.sfs[,reftb$DIMREDUC$PLS[[lvl[l]]]$euCols,drop=F]
                    pls = c(predict(reftb$DIMREDUC$PLS[[lvl[l]]]$model, focal.sfs, type="response", comps=1:reftb$DIMREDUC$PLS[[lvl[l]]]$model$ncomp))
					w = LDA$posterior[, lvl[l]]
					if (discard_extraRange == TRUE) {
						Range = unlist(reftb$GENERAL$sweepAgeDistrib[[lvl[l]]][2:3])
						est_is_ok = (pls >= Range[1]) & (pls <= Range[2])
						pls[!est_is_ok] = NA
						pls = pls * w
						pls = pls[WIN]
						est_is_ok = est_is_ok[WIN]
						# append
						PARAM.WEIGHTS[inferrableSNP&est_is_ok, lvl[l]] = PARAM.WEIGHTS[inferrableSNP&est_is_ok, lvl[l]] + wgts[inferrableSNP&est_is_ok]
						EST.AGE[inferrableSNP&est_is_ok, lvl[l]] = EST.AGE[inferrableSNP&est_is_ok, lvl[l]] + pls[inferrableSNP&est_is_ok]
					} else {
						pls = pls * w
						pls = pls[WIN]
						# append
						EST.AGE[inferrableSNP, lvl[l]] = EST.AGE[inferrableSNP, lvl[l]] + pls[inferrableSNP]
					}
                    
				}
            }

			if (progressBar) {
				if (floor(n/nMax*length(bar)) > nLast) {
					nLast <- floor(n/nMax*length(bar))
					if (FALSE) { # bar
						bar[1:nLast] <- "="
						if (n/nMax!=1) bar[nLast+1] <- ">"
						cat("[",paste0(bar, collapse=""),"]","\n", sep="")
					} else { # only print %
						cat(floor(n/nMax*100),"%"," ",sep="")
					}
				}
			}
        } #s
    } #i

    if (infer.age) {
		if (discard_extraRange) {
			EST.AGE = EST.AGE / PARAM.WEIGHTS
		} else {
			EST.AGE = EST.AGE / WEIGHTS[,-1] # must be done before averaging the WEIGHTS
		}
	}
    STABILITY = sweep(STABILITY, MARGIN = 1, FUN = "/", STATS = N.TESTS)
    WEIGHTS = sweep(WEIGHTS, MARGIN = 1, FUN = "/", STATS = N.TESTS)
	
	# in case where the scaling by N.TESTS leads to Inf, set NA (because Inf could be interpreted as overconfidence, what Inf are NOT)
	STABILITY[!is.finite(STABILITY)] <- NA
	WEIGHTS[!is.finite(WEIGHTS)] <- NA
	EST.AGE[!is.finite(EST.AGE)] <- NA

	return(list("pos" = POS,
	"stability" = STABILITY,
	"postpr" = WEIGHTS,
	"param" = EST.AGE,
	"inferentiality" = N.TESTS / n))
}


