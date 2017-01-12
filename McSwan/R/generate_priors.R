
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
#' @description Semi-automatic initialization of the Bayesian prior distributions of parameters across all possible demographic models: neutral and with local positive adaptation in each specific deme.
#' @param nSimul (int) number of simulations to perform for each model
#' @param msDemography (string) the neutral demographic history of the metapopulation, formatted according to Hudson's \code{\link{ms}} conventions, where \emph{theta} equals \eqn{4*No*mutRate}
#' @param No (int) effective size of the reference population (i.e. number of diploid individuals)
#' @param windowSize (int) window size for the genome scan (number of base pairs)
#' @param restrictDemesTo (array of int) by default, McSwan will detect positive sweeps across the recent history of \bold{every} deme specified in the \code{ms} command; however, if you want to restrict the analysis to specific demes, provide their \bold{indices} here in a vector (NB. the deme index corresponds to its position under the \code{-I} option of the \code{ms} command, note that indexation is 1-based)
#' @param sweepAge (special list) prior distribution for the sweep ages (scaled in generations before present) (if \code{NULL} and the number of demes is superior to 2, the prior distribution will be automatically set, otherwise it is mandatory to specify the distribution manually, see \link{Details})
#' @param recRate (special list) prior distribution for the recombination rates (see \link{Details})
#' @details Prior distributions must be specified using the following syntax: \code{list("P", arg1, arg2)} with \emph{P} the name of the distribution function (e.g. \code{\link{runif}} for the uniform distribution, \code{\link{rlogunif}} for the log-uniform; please make sure you have quoted the function name and removed the argument brackets); \emph{arg1} and \emph{arg2} respectively the first and second arguments of the function (e.g. for \code{runif} will be the lower and upper limits of the distribution). Note that if the upper and lower limits of the distribution are distant by more than 3 orders of magnitude in \emph{log10}, we advise to use the \code{rlogunif} distribution.
#' 
#' If you manually provide the prior distributions for the \bold{sweep ages}, you have two options: \itemize{
#' \item either specify uniquely \code{list("P", arg1, arg2)} which will be the distribution set \bold{across} the demes,
#' \item or specify a list of distribution-lists, each distribution-list corresponding to the sweep age distribution for a specific deme (indexed as they appear in the \code{ms} command); for instance, for 2 demes, one would specify: \code{list(list("rlogunif", T_1, TT_1), list("runif", T_2, TT_2))}. Please note that you must specify distributions for every deme, even if you have restricted the analysis to specific demes \emph{via} the \code{restrictDemesTo} option.
#' }
#' @seealso \code{\link{coalesce}}, \code{\link{get_max_sweep_ages}}
#' @return An object of class \emph{referenceTable} with initialized prior fields.
#' @examples mut <- 2.5e-8
#'No <- 2053
#'reftb <- generate_priors(msDemography = paste("80 1 -t",4*No*mut,
#'                                              "-I 2 40 40",
#'                                              "-en 0 2",29153/No,
#'                                              "-en",70/(4*No),2,11869/No,
#'                                              "-ej",2000/(4*No),"1 2",
#'                                              "-en",2683/(4*No),2,7083/No),
#'                         No = No, windowSize = 20000, nSimul = 1e3)
#' @references Hudson, R.R. (2007) \code{ms} - a program for generating samples under neutral models \url{http://home.uchicago.edu/rhudson1/source/mksamples/msdir/msdoc.pdf}
#' @export
generate_priors <- function(msDemography,
                            No,
                            fold,
                            windowSize = 20000,
                            nSimul = 10000,
                            restrictDemesTo = NULL,
                            sweepAge = NULL, 
                            recRate = list("rlogunif", 1e-12, 1e-3)) {
  
  FUNnullSweepAge <- "rlogunif"
  
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
    r <- do.call(recRate[[1]], 
                  c(nSimul, recRate[-1]))

    priors[[isl]] <- data.frame(sweepAge=sa, recRate=r)
    
  }
  
  names(islands) <- paste("i",names(islands),sep="")
  if (!is.null(restrictDemesTo)) priors <- priors[restrictDemesTo]
  reftb <- list(GENERAL = list(msDemography = paste(msarr, collapse=" "),
                               No = as.integer(No),
                               nSimul = nSimul,
                               windowSize = as.integer(windowSize),
                               islandSizes = islands,
                               sweepAgeDistrib = saDistrib,
                               recRateDistrib = recRate,
                               folded = fold),
                PRIORS = priors,
                SFS = NULL,
                DIMREDUC = NULL)
  
  class(reftb) <- "referenceTable"
  return(reftb)
}

