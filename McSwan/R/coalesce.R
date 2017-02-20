

#' @title Extract the deme effective sizes just before an event time-scaled age
#' @param ms the ms-formatted demographic history
#' @param islandIndex the index of the focal island
#' @param eventAge the age of the event, scaled in units of 4No
#' @return the resizing coefficient, in reference to No
#' @keywords internal
get_size <- function(ms, islandIndex, eventAge) {
# IMPORTANT NOTE: McSwan will NOT re-establish Ne conditioned by population growth before (backward in time) the sweep

  msarr <- unlist(strsplit(ms, " "))

  df <- c()

  z <- which(msarr=="-n")
  if (length(z)>0) z <- z[sapply(z, function(x) msarr[x+1]==islandIndex)]
  if (length(z)>0) df <- rbind(df, cbind("n", z, 0, msarr[z+2]))

  z <- which(msarr=="-N")
  if (length(z)>0) df <- rbind(df, cbind("N", z, 0, msarr[z+1]))

  z <- which(msarr=="-en")
  if (length(z)>0) z <- z[sapply(z, function(x) as.numeric(msarr[x+1])<=eventAge &&
                                   as.numeric(msarr[x+2])==islandIndex)]
  if (length(z)>0) df <- rbind(df, cbind("en", z, msarr[z+1], msarr[z+3]))

  z <- which(msarr=="-eN")
  if (length(z)>0) z <- z[sapply(z, function(x) as.numeric(msarr[x+1])<=eventAge)]
  if (length(z)>0) df <- rbind(df, cbind("eN", z, msarr[z+1], msarr[z+2]))

  if (is.null(nrow(df)) || nrow(df)==0) return(1)

  df <- as.data.frame(df, stringsAsFactors=F)
  colnames(df) <- c("type", "index", "age", "k")
  df$age <- as.numeric(df$age)

  return(as.numeric(df[which.max(df$age),"k"]))
}


#' @title Extract all demographic events involving a given island, before an event time-scaled age
#' @description does not handle -es -eM -eG -eN -N -M -ema / this internal function is only called if method=="partition" or method=="nMigrants" at \code{coalesce()}
#' @param ms the ms-formatted demographic history
#' @param islandIndex the index of the focal island
#' @param eventAge the age of the event, scaled in units of 4No
#' @return a string of ms-formatted history events affecting the twinned/duplicated focal island
#' @keywords internal
get_events <- function(ms, islandIndex, eventAge) {
# ___/!\___ this internal function is only called if method=="partition" or method=="nMigrants" at \code{coalesce()}

  msarr <- unlist(strsplit(ms, " "))
  nOriIsl <- as.integer(msarr[which(msarr=="-I")+1])

  add <- c()

  # -n | -g
  ng <- which(msarr=="-n"|msarr=="-g")
  if (length(ng)>0) {
    ngI <- ng[sapply(ng, function(x) msarr[x+1]==islandIndex)]
    for (i in seq_along(ngI)) {
      d <- msarr[ngI[i]:(ngI[i]+2)]
      d[2] <- nOriIsl + 1
      add <- c(add, paste(d, collapse=" "))
    }
  }

  # -m
  m <- which(msarr=="-m")
  if (length(m)>0) {

    mI1 <- m[sapply(m, function(x) msarr[x+1]==islandIndex)]
    for (i in seq_along(mI1)) {
      d <- msarr[mI1[i]:(mI1[i]+3)]
      d[2] <- nOriIsl + 1
      add <- c(add, paste(d, collapse=" "))
    }

    mI2 <- m[sapply(m, function(x) msarr[x+2]==islandIndex)]
    for (i in seq_along(mI2)) {
      d <- msarr[mI2[i]:(mI2[i]+3)]
      d[3] <- nOriIsl + 1
      add <- c(add, paste(d, collapse=" "))
    }

  }

  # -en | -eg
  eng <- which(msarr=="-en"|msarr=="-eg")
  if (length(eng)>0) {
    engI <- eng[sapply(eng, function(x) as.numeric(msarr[x+1])<=eventAge & msarr[x+2]==islandIndex)]
    for (i in seq_along(engI)) {
      d <- msarr[engI[i]:(engI[i]+3)]
      d[3] <- nOriIsl + 1
      add <- c(add, paste(d, collapse=" "))
    }
  }

  # -em
  em <- which(msarr=="-em")
  if (length(em)>0) {

    emI1 <- em[sapply(em, function(x) as.numeric(msarr[x+1])<=eventAge & msarr[x+2]==islandIndex)]
    for (i in seq_along(emI1)) {
      d <- msarr[emI1[i]:(emI1[i]+4)]
      d[3] <- nOriIsl + 1
      add <- c(add, paste(d, collapse=" "))
    }

    emI2 <- em[sapply(em, function(x) as.numeric(msarr[x+1])<=eventAge & msarr[x+3]==islandIndex)]
    for (i in seq_along(emI2)) {
      d <- msarr[emI2[i]:(emI2[i]+4)]
      d[4] <- nOriIsl + 1
      add <- c(add, paste(d, collapse=" "))
    }

  }

  return(paste(add, collapse=" "))
}


#' @title Coalescent simulations , given a set of parameter priors
#' @description Given a \code{referenceTable} object containing prior parameter values and a neutral demographic history, this function generates the command lines of neutral and selective models and execute them with Hudson's \emph{MS} coalescent simulator.
#' @param x an initialized \code{referenceTable} object
#' @param execute (logical) whether to execute the simulations
#' @param verbose (logical) verbose mode
#' @return A reference table with all simulated site frequency spectra stored in the \code{SFS} slot.
#' @details Past population size rescalings are always set in reference to \code{No}. If you want to see the demographic command lines for each model before simulating, set \code{execute = FALSE} and \code{verbose = TRUE}.
#' @seealso \code{\link{generate_priors}}, \code{\link{dim_reduction}}
#' @references Hudson, R.R. (2002). Generating samples under a Wright-Fisher neutral model of genetic variation. Bioinformatics, 18, 337-338.
#' @examples Please refer to the vignette.
#' @export
coalesce <- function(x, execute = TRUE, verbose = FALSE) {
  # method: partition || nMigrants || NULL (NULL is recommended, RT 14012017)
  
  # { internally set options
  method = "NULL" # NO SUBPARTITION (unlike Przeworski) and therefore also NO MIGRATION
  nRep = 1 # nRep (optional) number of repetitions to approximate the sweep (default 100)
  SAA = 1
  phyclust = TRUE
  ## }

  G <- x$GENERAL
  P <- x$PRIORS
  SFS <- list()

  nIsl <- length(G$islandSizes)

  msarr <- unlist(strsplit(G$msDemography, " "))

  # altering theta by incorporating windowSize specification
  msarr[4] <- sprintf('%f', as.numeric(msarr[4]) * G$windowSize)

  ##############################
  ##############################
  # NEUTRAL KINGMAN

  if (execute) file.create(paste(tempDir,"/segsites.txt",sep=""), showWarnings=FALSE)
  if (execute) file.create(paste(tempDir,"/sfs",sep=""), showWarnings=FALSE)

  cat("\n> Neutral model\n")
  if (execute) {
    if (phyclust) {
      if (verbose) cat(paste(c(msarr[-c(1:2)],"\n"), collapse=" "))
      phyclust::ms(nsam = as.integer(msarr[1]),
                   nreps = as.integer(G$nSimul),
                   opts = paste(msarr[-c(1:2)], collapse=" "),
                   temp.file = paste(tempDir,"/segsites.txt",sep=""),
                   tbs.matrix = NULL)
    } else {
      msarr[2] <- sprintf('%i', G$nSimul)
      cmd <- paste(msPath,paste(msarr,collapse=" "),">",tempDir,"/segsites.txt")
      system(cmd, intern = F)
    }
  }
  cmd <- paste(pythonPath," ",pyPath," -i ",tempDir,"/segsites.txt -o ",tempDir,"/sfs.txt", sep="")
  if (G$folded==TRUE) cmd <- paste(cmd,"--fold")
  if (execute) system(cmd, intern=F)
  if (execute) file.create(paste(tempDir,"/segsites.txt",sep=""), showWarnings=FALSE)
  if (execute) SFS[["i0"]] <- as.matrix(read.table(paste(tempDir,"/sfs.txt",sep=""), header=T, sep="\t"))

  ##############################
  ##############################
  # MULTICOALESCENT
  cat("\n> Selective models\n\n")

  # bottleneck intensity to mimick multicoalescence
  Ib <- (1/2) /(G$No)

  for ( i in seq_along(G$islandSizes) ) {

    if (execute) file.create(paste(tempDir,"/segsites.txt",sep=""), showWarnings=FALSE)
    if (execute) file.create(paste(tempDir,"/sfs",sep=""), showWarnings=FALSE)

	if (!verbose) cat(">>> Population index:    ",i,"\n")
    for ( j in 1:G$nSimul ) {
      if (verbose) cat(paste(">>>     Population Index:    ",i,"     -     Simulation:    ",j,"\n"))
      if (!verbose) cat(j,"")

      msarrTmp <- msarr
      msarrTmp[2] <- sprintf('%i', nRep)

      # altering island specification
      I <- which(msarrTmp=="-I")
      islSize <- as.integer(msarrTmp[I+1+i])
      nOriIsl <- as.integer(msarrTmp[I+1])

      if (method=="partition") {
        nNonSweeping <- floor( islSize * P[[i]]$recRate[j] )
        nSweeping <- islSize - nNonSweeping
        msarrTmp[I+1+i] <- nSweeping
        msarrTmp[I+1+nIsl] <- paste(msarrTmp[I+1+nIsl], nNonSweeping, 0)
        msarrTmp[I+1] <- 1L + as.integer(msarrTmp[I+1])
      } else if (method=="nMigrants") {
        msarrTmp[I+1+i] <- islSize
        msarrTmp[I+1+nIsl] <- paste(msarrTmp[I+1+nIsl], 0, 0)
        msarrTmp[I+1] <- 1L + as.integer(msarrTmp[I+1])
      } else {
# 14012016
# NADA
      }

if (method=="partition" || method=="nMigrants") {
      # sweep age (absolute count of generations)
      Ts <- P[[i]]$sweepAge[j]

      # demographic events for the new twin subisland
      twinMS <- get_events(G$ms, i, Ts/(4*G$No))

      # near k
      islK <- get_size(G$ms, i, Ts/(4*G$No))

      # opts list
      addMigr <- ifelse(method=="nMigrants", P[[i]]$recRate[j]*G$No*4, 0)

      cmdList <- paste(paste(msarrTmp[-c(1:2)],collapse=" "),
                       twinMS,
                       # bottleneck to mimick multicoalescence
                       "-en", Ts/(4*G$No), i, Ib,
                       # migration
                       "-m", i, nOriIsl+1, addMigr,
                       # merge sub-islands after multicoalescence
                       ###### TO IMPROVE TO AVOID CONFLICT OF ej AGES
                       "-en", (Ts+islSize-1)/(4*G$No), i, islK,
                       "-ej", (Ts+islSize)/(4*G$No), nOriIsl+1, i)
      cmdList <- gsub("\\s+", " ", cmdList)
} else {
      # sweep age (absolute count of generations)
      Ts <- P[[i]]$sweepAge[j]

      # near k
      islK <- get_size(G$ms, i, Ts/(4*G$No))

      cmdList <- paste(paste(msarrTmp[-c(1:2)],collapse=" "),
                       # bottleneck to mimick multicoalescence
                       "-en", Ts/(4*G$No), i, Ib,
                       # merge sub-islands after multicoalescence
                       ###### TO IMPROVE TO AVOID CONFLICT OF ej AGES
                       "-en", (Ts+islSize-1)/(4*G$No), i, islK)
      cmdList <- gsub("\\s+", " ", cmdList)
}

      if (execute) {
        if (phyclust) {
          if (verbose) cat(paste(cmdList,"\n"))
          phyclust::ms(nsam = as.integer(msarrTmp[1]),
                       nreps = nRep,
                       opts = cmdList,
                       temp.file = paste(tempDir,"/segsites.txt",sep=""),
                       tbs.matrix = NULL)
        } else {
          # ms command
          stop("bad")
          cmd <- paste(msPath,paste(msarrTmp,collapse=" "),
                       twinMS,
                       # migration
                       "-m",i,nOriIsl+1,m,
                       # bottleneck to mimick multicoalescence
                       "-en",Ts/(4*G$No),i,Ib,
                       # merge sub-islands after multicoalescence
                       ###### TO IMPROVE TO AVOID CONFLICT OF ej AGES
                       "-en",(Ts+islSize-1)/(4*G$No),islK,
                       "-ej",(Ts+islSize)/(4*G$No),nOriIsl+1,i,
                       "> temp/segsites.txt")
          cmd <- gsub("\\s+", " ", cmd)
          if (verbose) cat(paste(cmd,"\n"))
          if (execute) system(cmd, intern=F)
        }
      }

      #cmdPy <- paste("python ",pyPath," -i ",tempDir,"/segsites.txt -o ",tempDir,"/sfs.txt -m ",nOriIsl+1," ",i, sep="")
      #if (execute) system(cmdPy)
      #if (execute) file.create(paste(tempDir,"/segsites.txt",sep=""), showWarnings=FALSE)
      #if (execute) SFS[[paste("i",i,sep="")]] <- rbind(SFS[[paste("i",i,sep="")]], apply(read.table(paste(tempDir,"/sfs.txt",sep=""), sep="\t", header=T), 2, mean))

    }

if (method=="partition" || method=="nMigrants") {
    cmdPy <- paste(pythonPath," ",pyPath," -i ",tempDir,"/segsites.txt -o ",tempDir,"/sfs.txt -m ",nOriIsl+1," ",i, sep="")
} else {
    cmdPy <- paste0(pythonPath," ",pyPath," -i ",tempDir,"/segsites.txt -o ",tempDir,"/sfs.txt")
}
    if (G$folded==TRUE) cmdPy <- paste(cmdPy,"--fold")
    if (execute) system(cmdPy)
    if (execute) SFS[[paste("i",i,sep="")]] <- as.matrix(read.table(paste(tempDir,"/sfs.txt",sep=""), header=T, sep="\t"))

	if (!verbose) cat("\n\n")
  }

  x$SFS <- SFS
  x$GENERAL$call.coalesce <- match.call()
  return(x)
}



