

#' @title Extract the deme effective size juste before an event time-scaled age
#' @description v1: only takes into account "-en" (perspectives: "-eN")
#' @param ms the ms-formatted demographic history
#' @param islandIndex the index of the focal island
#' @param eventAge the age of the event, scaled in 4No
#' @return the size coefficient in reference to No
#' #@export
#' @keywords internal
get_size <- function(ms, islandIndex, eventAge) {
print("size")
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


#' @title Extract all demographic events before an event time-scaled age and involving some specific island
#' @description does not handle -es -eM -eG -eN -N -M -ema
#' @param ms the ms-formatted demographic history
#' @param islandIndex the index of the focal island
#' @param eventAge the age of the event, scaled in 4No
#' @return a string of ms-formatted history events affecting the twin subisland
#' #@export
#' @keywords internal
get_events <- function(ms, islandIndex, eventAge) {

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


#' @title Performs coalescence and output multidimensional site frequency spectra (SFS)
#' @param x an initialized reference table produced by \code{\link{generate_priors}}
#' @param nRep (int) number of repetitions for each parameter combination (only for debug mode, otherwise prefer the default value)
#' @param execute (logical) if TRUE, execute the coalescence
#' @param verbose (logical) if TRUE, print all messages
#' @return A reference table, of class \code{referenceTable}, with all simulated multidimensional site frequency spectra.
#' @export
coalesce_migr <- function(x,
                     nRep = 1,
                     SAA = 1,
                     execute = TRUE,
                     verbose = FALSE) {

  if (class(x)!="referenceTable") stop("x is not an object of class referenceTable")

  G <- x$GENERAL
  P <- x$PRIORS
  SFS <- list()

  nIsl <- length(G$islandSizes)

  msarr <- unlist(strsplit(G$msDemography, " "))

  # altering theta by incorporating windowSize specification
  msarr[4] <- sprintf('%f', as.numeric(msarr[4]) * G$windowSize)

  if (execute) file.create(paste(tempDir,"/segsites.txt",sep=""), showWarnings=FALSE)

  ##############################
  ##############################
  # NEUTRAL KINGMAN
  cat("\nNeutral coalescence\n")
  if (execute) {
    if (verbose) cat(paste(c(msarr[-c(1:2)],"\n"), collapse=" "))
    phyclust::ms(nsam = msarr[1],
                 nreps = G$nSimul,
                 opts = paste(msarr[-c(1:2)], collapse=" "),
                 temp.file = paste(tempDir,"/segsites.txt",sep=""),
                 tbs.matrix = NULL)
  }
  cmd <- paste(pythonPath," ",pyPath," -i ",tempDir,"/segsites.txt -o ",tempDir,"/sfs.txt", sep="")
  if (execute) system(cmd, intern=F)
  if (execute) file.create(paste(tempDir,"/segsites.txt",sep=""), showWarnings=FALSE)
  if (execute) SFS[["i0"]] <- as.matrix(read.table(paste(tempDir,"/sfs.txt",sep=""), header=T, sep="\t"))

  ##############################
  ##############################
  # MULTICOALESCENT
  cat("\nMulticoalescence\n\n")

  # bottleneck intensity to mimick multicoalescence
  Ib <- 1/G$No

  for ( i in seq_along(G$islandSizes) ) {

    for ( j in 1:G$nSimul ) {
      if (verbose) cat(paste(">     islandIndex:    ",i,"     -     simulation:    ",j,"\n"))

      msarrTmp <- msarr
      msarrTmp[2] <- sprintf('%i', nRep)

      # altering island specification
      I <- which(msarrTmp=="-I")
      islSize <- as.integer(msarrTmp[I+1+i])
      nOriIsl <- as.integer(msarrTmp[I+1])
      msarrTmp[I+1+nIsl] <- paste(msarrTmp[I+1+nIsl],"0")
      msarrTmp[I+1] <- 1L + as.integer(msarrTmp[I+1])

      # sweep age (absolute count of generations)
      Ts <- P[[i]]$sweepAge[j]

      # migration rate

      #alpha <- P[[i]]$recRate[j] / SAA * log(2*G$No)
      #p <- 1 - exp(-alpha*G$windowSize)
      #m <- 4 * G$No * p

      #alpha <- P[[i]]$recRate[j]
      p <- 1 - exp(-1* P[[i]]$recRate[j] * G$windowSize * log(2*G$No))
      #p <- P[[i]]$recRate[j] * G$windowSize
      m <- 4 * G$No * P[[i]]$recRate[j]

      m <- 4 * G$No * P[[i]]$recRate[j] * G$windowSize / Ts

      m <- 4 * G$No * p/Ts

      #4Nm
      m <- 4*G$No*P[[i]]$recRate[j]

      #4Nexp(rd)
      p <- 1 - exp(-1* P[[i]]$recRate[j] * log(2*G$No))
      m <- 4*G$No*p

      # demographic events for the new twin subisland
      twinMS <- get_events(G$ms, i, Ts/(4*G$No))

      # near k
      islK <- get_size(G$ms, i, Ts/(4*G$No))

      # opts list
      cmdList <- paste(paste(msarrTmp[-c(1:2)],collapse=" "),
                       twinMS,
                       # migration
                       "-m",i,nOriIsl+1,m,
                       # bottleneck to mimick multicoalescence
                       "-en",Ts/(4*G$No),i,Ib,
                       # merge sub-islands after multicoalescence
                       ###### TO IMPROVE TO AVOID CONFLICT OF ej AGES
                       "-en",(Ts+islSize-1)/(4*G$No),i,islK,
                       "-ej",(Ts+islSize)/(4*G$No),nOriIsl+1,i)
      cmdList <- gsub("\\s+", " ", cmdList)

      if (execute) {
        if (verbose) cat(paste(cmdList,"\n"))
        phyclust::ms(nsam = as.integer(msarrTmp[1]),
                     nreps = nRep,
                     opts = cmdList,
                     temp.file = paste(tempDir,"/segsites.txt",sep=""),
                     tbs.matrix = NULL)
      }

      if (execute) system(paste(pythonPath," ",pyPath," -i ",tempDir,"/segsites.txt -o ",tempDir,"/sfs.txt -m ",nOriIsl+1," ",i, sep=""))
      if (execute) file.create(paste(tempDir,"/segsites.txt",sep=""), showWarnings=FALSE)
      if (execute) SFS[[paste("i",i,sep="")]] <- rbind(SFS[[paste("i",i,sep="")]],
                                                       apply(read.table(paste(tempDir,"/sfs.txt",sep=""), sep="\t", header=T), 2, mean))

    }

  }

  x$SFS <- SFS
  return(x)
}


#' @title Performs multicoalescence SIMPLE PARTITION WITHOUT MIGRATION, given a set of parameter priors
#' @description NA
#' @details Past resizings are always set in reference to  \emph{No}.
#' @param x an initialized reference table (class \emph{referenceTable}).
#' @param nRep (optional) number of repetitions to approximate the sweep (default 100)
#' @return A reference table with all simulated site frequency spectra.
#' @export
coalesce <- function(x, method = "partition", nRep = 1,
                          SAA = 1,
                          execute = TRUE,
                          phyclust = TRUE,
                          verbose = FALSE) {
  # method: partition || nMigrants
NNY <- 2
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

  cat("\nNeutral coalescence\n")
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
  cat("\nMulticoalescence\n\n")

  # bottleneck intensity to mimick multicoalescence
  Ib <- (1/2) /(G$No)

  for ( i in seq_along(G$islandSizes) ) {

    if (execute) file.create(paste(tempDir,"/segsites.txt",sep=""), showWarnings=FALSE)
    if (execute) file.create(paste(tempDir,"/sfs",sep=""), showWarnings=FALSE)

    for ( j in 1:G$nSimul ) {
      if (verbose) cat(paste(">     islandIndex:    ",i,"     -     simulation:    ",j,"\n"))

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
      } else {
        msarrTmp[I+1+i] <- islSize
        msarrTmp[I+1+nIsl] <- paste(msarrTmp[I+1+nIsl], 0, 0)
        msarrTmp[I+1] <- 1L + as.integer(msarrTmp[I+1])
      }

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

    cmdPy <- paste(pythonPath," ",pyPath," -i ",tempDir,"/segsites.txt -o ",tempDir,"/sfs.txt -m ",nOriIsl+1," ",i, sep="")
    if (G$folded==TRUE) cmdPy <- paste(cmdPy,"--fold")
    if (execute) system(cmdPy)
    if (execute) SFS[[paste("i",i,sep="")]] <- as.matrix(read.table(paste(tempDir,"/sfs.txt",sep=""), header=T, sep="\t"))

  }

  x$SFS <- SFS
  return(x)
}



