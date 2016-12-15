
#' @title Generates pseudo-observed genomic fragments under arbitrary demographic histories with and without positive sweeps
#' @description Simulation based on Ewing's \emph{et al.} \code{msms}. Assumes very strong selection coefficients.
#' @param reftb an initialized \emph{referenceTable} object
#' @param nSimul (integer) number of independent genomic fragments to simulate
#' @param L (integer) genomic fragment length (in base pairs); if \code{NULL}, the function will perform a non-sliding window validation; if not-\code{NULL}, the function will perform sliding window validation (closer to the genomic scan approach used to analyze observed dataset)
#' @param sweepAge (special list) prior distribution for the sweep ages (scaled in generations before present) (if \code{NULL} and the number of demes is superior to 2, the prior distribution will be automatically set, otherwise it is mandatory to specify the distribution manually, see \link{Details})
#' @param recRate (special list) prior distribution for the recombination rates (see \link{Details})
#' @param sweepPos (numeric) relative position of the beneficial mutation (eg. for \eqn{sweepPos=0.5}, the beneficial mutation will be placed at the center of the genomic region)
#' @param nReps (integer) number of simulations per parameter combination (by default 1, but you can increase it to reduce the stochasticity of the coalescent process, the function will automatically average the coalescent-generated multiSFSs)
#' @param verbose (logical) whether to print all messages (TRUE) or not (FALSE)
#' @param doSFS (logical) whether to compute the multiSFS (TRUE) or not (FALSE)
#' @param Smu forward mutation rate to the beneficial allele
#' @return An object with special class (see below) containing in the \code{SFS} element: \itemize{
#' \item multidimensional SFSs if \eqn{L = NULL} (class \code{validationTable})
#' \item positional allele frequencies if \eqn{L != NULL} (class \code{slidingValidationTable})
#' }
#'  
#' @details Prior distributions must be specified using the following syntax: \code{list("P", arg1, arg2)} with \emph{P} the name of the distribution function (e.g. \code{\link{runif}} for the uniform distribution, \code{\link{rlogunif}} for the log-uniform; please make sure you have quoted the function name and removed the argument brackets); \emph{arg1} and \emph{arg2} respectively the first and second arguments of the function (e.g. for \code{runif} will be the lower and upper limits of the distribution). Note that if the upper and lower limits of the distribution are distant by more than 3 orders of magnitude in \emph{log10}, we advise to use the \code{rlogunif} distribution.
#' 
#' If you manually provide the prior distributions for the \bold{sweep ages}, you have to respect the following format: \itemize{
#' \item specify a list of distribution-sublists, each distribution-sublist corresponding to the sweep age distribution for the specific deme(s) you provided in the \code{sweepingIsl} argument (indexed as they appear in the \code{ms} command) (if \eqn{sweepingIsl = NULL} you will need to provide the distribution-sublists for \bold{all} demes); for instance, for \eqn{sweepingIsl = c(1,2)}, one would specify: \code{list(list("rlogunif", T_1, TT_1), list("runif", T_2, TT_2))}.
#' }
#' 
#' @references Ewing et Hermisson (2010) MSMS: a coalescent simulation program including recombination, demographic structure and selection at a single locus. \emph{Bioinformatics}.
#' @export
generate_pseudoobs <- function(reftb, 
                               nSimul = 50, 
                               sweepingIsl = NULL,
                               recRate, 
                               sweepAge = NULL,
                               L = NULL,
                               sweepPos = .5,
                               Smu = .01, 
                               nReps = 1, 
                               verbose = FALSE, 
                               doSFS = TRUE,
                               default_sweepAge_prior_func = "rlogunif",
                               save_each_file = FALSE) {
  
  if (is.null(L)) { print("non-sliding validation") } else { print("sliding validation") }
  
  if (!is.list(recRate)) stop("recRate must be a single list with 3 elements, cf. documentation")
  if (is.list(recRate[[1]])) stop("recRate must be a single list with 3 elements, cf. documentation")
  
  # assumes a selection strength of 1 (complete sweep)
  # assumes half selection pressure on heterygote
  # assumes beneficial allele at freq of 0 before selection occurs
  
  # fill valtb
  valtb <- list(GENERAL = reftb$GENERAL,
                PRIORS = list(),
                SFS = list())
  class(valtb) <- ifelse(is.null(L), "validationTable", "slidingValidationTable")
  valtb$GENERAL$nSimul <- nSimul
  
  # permanent parameters
  sAA <- 1; sAa <- .5 * sAA; saa <- 0
  #sAa <- 1*sAA ####!!!!
  
  # wp
  No <- reftb$GENERAL$No
  windowSize <- ifelse(is.null(L), reftb$GENERAL$windowSize, L)
  islandSizes <- reftb$GENERAL$islandSizes
  nIsl <- length(islandSizes)
  
  # sweepAge format
  if (is.null(sweepAge)) {
    sweepAge <- reftb$GENERAL$sweepAgeDistrib
    sweepAge <- lapply(sweepAge, function(x) { x[[1]] <- default_sweepAge_prior_func; return(x) })
  } else {
    if (!is.list(sweepAge)) stop("sweepAge must be a list of sub-lists, cf. documentation")
    if (!is.list(sweepAge[[1]])) stop("sweepAge must be a list of sub-lists, cf. documentation")
    if (is.null(sweepingIsl)) {
      if (length(sweepAge)!=length(reftb$PRIORS)) stop("sweepAge must be a list of 'n' sub-lists with 'n' the number of islands, cf. documentation")
      names(sweepAge) <- names(reftb$PRIORS)
    } else {
      if (length(sweepAge)!=length(sweepingIsl[sweepingIsl!=0])) stop("sweepAge must be a list of 'n' sub-lists with 'n' the number of sweepingIsl, cf. documentation")
      names(sweepAge) <- paste("i",sweepingIsl,sep="")
    }
  }
  
  sweepAge <- lapply(sweepAge, function(x) { if (x[[2]]<=1) {
    x[[2]] <- 5 }
    return(x) })
  
  # definition of sweeping islands
  if (is.null(sweepingIsl)) {
    isl4sim <- names(reftb$PRIORS)
    if (!"i0"%in%isl4sim) isl4sim <- c("i0", isl4sim)
  } else {
    isl4sim <- paste("i",sweepingIsl,sep="")
    if (!"i0"%in%isl4sim) isl4sim <- c("i0", isl4sim)
    if (!all(isl4sim%in%c(names(reftb$PRIORS),"i0"))) stop("could not find this island index in the referenceTable")
  }

  # get priors
  ii <- isl4sim[isl4sim!="i0"]
  P <- lapply(as.list(ii), function(x) data.frame(sweepAge = do.call(sweepAge[[x]][[1]],
                                                                          list(nSimul,
                                                                          sweepAge[[x]][[2]],
                                                                          sweepAge[[x]][[3]])),
                                                       recRate = do.call(recRate[[1]], list(nSimul, recRate[[2]], recRate[[3]])),
                                                       sweepPos = sweepPos))
  names(P) <- ii
  P$i0 <- P$i0 * NA
  valtb$PRIORS <- P
  
  # set number of repetitions in MS command
  ms <- reftb$GENERAL$msDemography
  msarr <- unlist(strsplit(ms, " "))
  msarr[2] <- sprintf('%i', nReps)
  
  # altering THETA by incorporating windowSize specification
  msarr[4] <- sprintf('%f', as.numeric(msarr[4]) * windowSize)
  ms <- paste(msarr, collapse=" ")
  
  #################
  for ( i in 1:length(isl4sim) ) {
    
    isl <- isl4sim[i]
    cat(paste("\n\n\n>>>",isl,"\n\n\n"))
    
    SFS <- c()
    for ( s in 1:nSimul ) {
      cat(s, ".")
      
      if (isl=="i0") {
        msms <- paste("-N",No,"-ms",ms)
      } else {
        # initial frequency of beneficial allele in each deme
        #beneFreq <- paste(rep(1/No, nIsl), collapse=" ")
        #beneFreq <- paste(c(1/2053, 1/29153), collapse=" ")
        beneFreq <- sapply(seq_along(islandSizes), function(j) {
          1 / (2*No * get_size(ms, j, P[[isl]]$sweepAge[s]/(4*No)))
        })
        beneFreq <- paste(beneFreq, collapse=" ")
        
        # sweep start
        SI <- paste(P[[isl]]$sweepAge[s]/(4*No), nIsl, beneFreq, collapse=" ")
        
        # recombination
        rho <- paste(4 * No * P[[isl]]$recRate[s] * (windowSize-1), sprintf('%i', windowSize), sep=" ")
        msTmp <- paste(ms, "-r", rho, sep=" ")
        
        # deme-specific selection coefficient
        ## for deme under sweep
        Sc <- paste("-Sc", 0, gsub("i", "", isl), sAA*2*No, sAa*2*No, saa*2*No)
        ## for other demes (neutral evol)
        ScUnsweeped <- c()
        for (o in seq_along(islandSizes)) {
          if (paste("i",o,sep="")==isl) next()
          add <- paste("-Sc", 0, o, 0, 0, 0, sep=" ")
          ScUnsweeped <- c(ScUnsweeped, add)
        }
        Sc <- paste(Sc, ScUnsweeped, collapse=" ")
        
        # command
        msms <- paste("-N",No,"-ms",msTmp,Sc,"-SI",SI,"-Sp",P[[isl]]$sweepPos[s])
        if (!is.null(Smu)) msms <- paste(msms,"-Smu",Smu)
        
        ###!
        msms <- paste(msms,"-oOC","-SForceKeep")
        #msms <- paste(msms,"-oTrace")
      }
      
      msms <- gsub("\\s+", " ", msms)
      if (verbose) cat(paste("\n",msms,"\n\n"))
      
      # execution
      file.create(paste(tempDir,"/segsites.txt",sep=""), showWarnings=FALSE)
      
      if (grepl("windows", Sys.info()["sysname"], ignore.case = TRUE)) {
        cmd <- paste("\"",javaPath,"\" -jar ",msmsPath," ",msms, sep="")
        tmp <- system(cmd, intern = TRUE)
        write(tmp, file=paste(tempDir,"/segsites.txt",sep=""))
      } else {
        cmd <- paste("\"",javaPath,"\" -jar ",msmsPath," ",msms," > ",tempDir,"/segsites.txt", sep="")
        system(cmd, intern = FALSE)
      }

      if (save_each_file==TRUE) file.copy(paste(tempDir,"/segsites.txt",sep=""), paste(tempDir,"/segsites_",isl4sim[i],"_",s,".txt",sep=""), copy.date = TRUE, overwrite = TRUE)

      # SINGLE WINDOW
      if (is.null(L)) {
        if (doSFS) {
          cmd <- paste(pythonPath," ",pyPath," -i ",tempDir,"/segsites.txt -o ",tempDir,"/sfs.txt", sep="")
          if (valtb$GENERAL$folded==TRUE) cmd <- paste(cmd,"--fold")
          system(cmd, intern=F)
          SFS <- rbind(SFS, apply(read.table(paste(tempDir,"/sfs.txt",sep=""), header=T, sep="\t"), 2, mean))
        }
      # SLIDING VALIDATION
      } else {
        if (doSFS) {
          
          if (grepl("windows", Sys.info()["sysname"], ignore.case = TRUE)) {
            # cleaning
            cmd <- paste(pythonPath," ",pyCLEANMSPath," -i ",tempDir,"/segsites.txt -o ",tempDir,"/segsitesCLEAN.txt",sep="")
            system(cmd, intern=F)
          } else {
            file.rename(paste(tempDir,"/segsites.txt",sep=""), paste(tempDir,"/segsitesCLEAN.txt",sep=""))
          }
          
          # importing
          write("", file=paste(tempDir,"/MSMS.txt", sep=""))
          cmd <- paste(pythonPath," ",pyMSMSPath," -i ",tempDir,"/segsitesCLEAN.txt -o ",tempDir,"/MSMS.txt -L ",sprintf('%i',L), sep="")
          if (valtb$GENERAL$folded==TRUE) cmd <- paste(cmd,"--fold")
          system(cmd, intern=F)

          obs <- get_SFS(paste(tempDir,"/MSMS.txt", sep=""), ms)
          SFS <- c(SFS, list(obs))
          names(SFS)[[length(SFS)]] <- s
        }
      }
      
    }
    valtb$SFS[[isl]] <- SFS
  }
  
  valtb$GENERAL$sweepPos <- sweepPos * L
  return(valtb)
}









