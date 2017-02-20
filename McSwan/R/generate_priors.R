
#' @title get_max_sweep_ages
#' @description Automatically determines the most ancient sweep age possible for local adaptation when >=2 populations.
#' @details Note that in the situation of a Kingman approximation to the real multicoalescent, we propose to underestimate the most ancient sweep age by \emph{dT(i)} (i.e. the multicoalescent delay approximation, by default \eqn{dT(i) = nSamples(i) - 1}) in order to avoid temporal conflicts with other neutral events specified in the \code{ms}-formatted demographic history.
#' @param ms (string) the \code{ms}-formatted neutral demographic history of the metapopulation
#' @param No (int) effective size of the reference population (ie. number of diploid individuals)
#' @seealso \code{\link{generate_priors}}
#' @return An array giving the most ancient sweep age (scaled in \eqn{4No}) for local adaptation in each deme.
#' @keywords internal
get_max_sweep_ages <- function(ms, No) {
  # if single deme, the user must specify max sweep age
  # note that -ej (forward admixture) should not be considered to put a bound to the focal deme-specific sweep events
  
  msarr <- unlist(strsplit(ms, " "))
  I <- which(msarr=="-I")
  n <- ifelse(length(I)==0, 1L, as.integer(msarr[I+1]))
  islandSizes <- as.integer(msarr[(I+2):(I+as.integer(msarr[I+1])+1)])
  names(islandSizes) <- seq_along(islandSizes)
  
  auto_max_sweep_age <- NULL
  if ( n > 1 ) {
    eji <- which(msarr=="-ej")
    ej <- matrix(sapply(eji, function(x) {
      return(as.numeric(c(msarr[x+1], msarr[x+2], msarr[x+3])))
    }), ncol=3, byrow=T)
    auto_max_sweep_age <- rep(NA, n)
    names(auto_max_sweep_age) <- 1:n
    # max sweep age is the min time to any divergence
    for ( p in 1:n ) {
      auto_max_sweep_age[p] <- min(ej[ej[,2]==p|ej[,3]==p,1])
    }
  }
  
  # underestimate to avoid conflicts
  auto_max_sweep_age <- auto_max_sweep_age - islandSizes/(4*No)
  
  return(auto_max_sweep_age)
}



#' @title Initialize the reference table
#' @description Semi-automatic initialization of the Bayesian prior distributions for all model parameters.
#' @param nSimul (int) number of simulations to perform for each model
#' @param msDemography (string) the neutral demographic history of the metapopulation, formatted according to Hudson's \emph{MS} conventions, where \eqn{\theta = 4N_oµ} with \eqn{N_o} the diploid effective size; times are scaled in units of \eqn{4N_o} and population size rescaling in units of \eqn{N_o}; note that the demographic history can involved unsampled (ghost) populations
#' @param No (integer) the \bold{diploid} effective size of the reference population
#' @param windowSize (integer) window length for the genome scan (number of base pairs)
#' @param restrictDemesTo (array of integers) by default, McSwan will detect population-specific sweeps across all the populations specified in the \code{MS} command; however, if you want to restrict the analysis to a particular subset of demes, provide their \bold{indices} in a vector; note that population indices correspond to their position under the \code{-I} switch of the \code{MS} command, and note that the index of the first population is \bold{1}
#' @param sweepAge (special list) prior distribution for the sweep ages (scaled in generations before present); if \code{NULL} and the number of populations is superior to 2, the prior distribution ranges will be automatically determined, otherwise it is mandatory to specify the distribution manually, see \code{Details})
#' @param fold (logical) whether allele counts are to be defined in reference to the \bold{minor} alleles (\code{fold = TRUE}, producing \bold{folded} site frequency spectra) or to the \bold{derived} alleles (\code{fold = FALSE), producing \bold{unfolded} site frequency spectra}. In the derived allele case, you must ensure that the SNPs in your VCF dataset are polarized relatively to some ancestral genomes.

#' @details Prior distributions must be specified using the following syntax: \code{list("P", arg1, arg2)} with \emph{P} the name of the distribution function (e.g. \code{\link{runif}} for the uniform distribution, \code{\link{rlogunif}} for the log-uniform; please make sure you have quoted the function name and removed the argument brackets); \emph{arg1} and \emph{arg2} respectively the first and second arguments of the function (e.g. for \code{runif} will be the lower and upper limits of the distribution). Note that if the upper and lower limits of the distribution are distant by more than 3 units in the \emph{log10} scale, it is recommended to use the \code{rlogunif} distribution.
#' 
#' If you manually provide the prior distributions for the \bold{sweep ages}, you have two options: \itemize{
#' \item either specify a single \code{list("P", arg1, arg2)} which will be the distribution set \bold{across} the demes,
#' \item or specify a list of distribution-lists, each distribution-list corresponding to the sweep age distribution for a specific deme (indexed as they appear in the \code{ms} command); for instance, for 2 demes, one would specify: \code{list(list("rlogunif", T_1, TT_1), list("runif", T_2, TT_2))}. Please note that you must specify distributions for every deme, even if you have restricted the analysis to specific demes \emph{via} the \code{restrictDemesTo} option.
#' }
#' @seealso \code{\link{coalesce}}
#' @return An object of class \emph{referenceTable} with initialized \code{PRIOR} slot.
#' @examples Please refer to the vignette.
#' @references Hudson, R.R. (2007) \code{ms} - a program for generating samples under neutral models \url{http://home.uchicago.edu/rhudson1/source/mksamples/msdir/msdoc.pdf}
#' @examples Please refer to the vignette.
#' @export
generate_priors <- function(msDemography,
                            No,
                            fold,
                            windowSize = 20000,
                            nSimul = 10000,
                            restrictDemesTo = NULL,
                            sweepAge = NULL) {
  
  # internally set options
  FUNnullSweepAge <- "rlogunif"
  recRate <- list("runif", 0, 0) #list("rlogunif", 1e-12, 1e-3)
  
  msDemography <- gsub("\\s+", " ", msDemography)
  
  # get island sizes
  msarr <- unlist(strsplit(msDemography, " "))
  
  # if single deme, force its appearance in ms command
  if (!"-I"%in%msarr) {
    cat("\nDetected you specified only one single deme.\n")
    th <- which(msarr=="-t")
    msarr <- c(msarr[1:(th+1)], "-I", 1, msarr[1], msarr[(th+2):length(msarr)])
    print(paste(msarr, collapse=" "))
  }
  
  I <- which(msarr=="-I")
  islands <- as.integer(msarr[which(msarr=="-t")-2])
  if ( length(I) > 0 ) islands <- as.integer(msarr[(I+2):(I+as.integer(msarr[I+1])+1)])
  names(islands) <- seq_along(islands)
  
  if (length(islands)==1 && is.null(sweepAge)) stop("if one single deme, you must specify the prior range of sweep ages")
  if (!is.null(sweepAge) && !is.list(sweepAge[[1]])) sweepAge <- lapply(seq_along(islands), function(x) sweepAge)
  if (!is.logical(fold)) stop("fold argument must be TRUE or FALSE")
  
  # island-specific priors
  priors <- list()
  saDistrib <- list()
  for ( i in seq_along(islands) ) {
    
    isl <- paste("i",names(islands)[i],sep="")
    priors[[isl]] <- data.frame()
    
    # sweepAge
    # default: automatically bounded by soonest divergence (LOCAL adaptation)
    # else: c("runif", absolute_time_lower, absolute_time_upper)
    if (is.null(sweepAge)) {
      maxT <- 4 * No * get_max_sweep_ages(msDemography, No)[i]
      attr(maxT, "names") <- NULL
      sa <- do.call(FUNnullSweepAge, 
                             list(nSimul, 
                             1, 
                             maxT))
      saDistrib <- c(saDistrib, list(list(FUNnullSweepAge, 1, maxT)))
    } else {
      if (length(sweepAge)!=1 && length(sweepAge)!=length(islands)) stop("unappropriate number of distribution-lists for the sweep age prior distribution")
      sa <- do.call(sweepAge[[i]][[1]], list(nSimul, sweepAge[[i]][[2]], sweepAge[[i]][[3]]))
      saDistrib <- sweepAge
    }
    names(saDistrib) <- paste("i",seq_along(saDistrib),sep="")
    
    # recombination
    #r <- do.call(recRate[[1]], c(nSimul, recRate[-1]))

    #priors[[isl]] <- data.frame(sweepAge=sa, recRate=r)
	priors[[isl]] <- data.frame(sweepAge=sa)
    
  }
  
  names(islands) <- paste("i",names(islands),sep="")
  if (!is.null(restrictDemesTo)) priors <- priors[restrictDemesTo]
  reftb <- list(GENERAL = list(msDemography = paste(msarr, collapse=" "),
                               No = as.integer(No),
                               nSimul = nSimul,
                               windowSize = as.integer(windowSize),
                               islandSizes = islands,
                               sweepAgeDistrib = saDistrib,
                               #recRateDistrib = recRate,
                               folded = fold,
							   call.generate_priors = match.call(),
							   creationDate = Sys.time()),
                PRIORS = priors,
                SFS = NULL,
                DIMREDUC = NULL)
  
  class(reftb) <- "referenceTable"
  return(reftb)
}

