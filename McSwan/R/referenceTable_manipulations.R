
#' @title rbind two referenceTables
#' @export
rbind_reftb <- function(V, v) {
  # R and r two referenceTables
  if (class(V)!="referenceTable" || class(v)!="referenceTable") stop("one of your argument is not a referenceTable")
  
  if (V$GENERAL$msDemography != v$GENERAL$msDemography) stop("msDemography does not match")
  if (V$GENERAL$No != v$GENERAL$No) stop("No does not match")
  if (V$GENERAL$windowSize != v$GENERAL$windowSize) stop("windowSize does not match")
  if (!all.equal(V$GENERAL$islandSizes, v$GENERAL$islandSizes)) stop("islandSizes does not match")
  
  if (!all.equal(names(V$PRIORS), names(v$PRIORS))) stop("PRIORS elements do not match")
  if (!all.equal(names(V$SFS), names(v$SFS))) stop("SFS elements do not match")
  
  Vv <- list()
  class(Vv) <- "referenceTable"
  Vv$GENERAL <- V$GENERAL
  Vv$GENERAL$nSimul <- V$GENERAL$nSimul + v$GENERAL$nSimul
  
  SFS <- lapply(seq_along(V$SFS), function(x) rbind(V$SFS[[x]], v$SFS[[x]]))
  names(SFS) <- names(V$SFS)
  
  PRIORS <- lapply(seq_along(V$PRIORS), function(x) rbind(V$PRIORS[[x]], v$PRIORS[[x]]))
  names(PRIORS) <- names(V$PRIORS)
  
  Vv$PRIORS <- PRIORS
  Vv$SFS <- SFS
  
  rm(list=c("SFS","PRIORS","V","v")); invisible(gc())
  return(Vv)
}


#' @title rbind two validationTables
#' @export
rbind_valtb <- function(V, v) {
  # V and v two validationTables
  if (class(V)!="validationTable" || class(v)!="validationTable") stop("one of your argument is not a validationTable")
  
  if (V$GENERAL$msDemography != v$GENERAL$msDemography) stop("msDemography does not match")
  if (V$GENERAL$No != v$GENERAL$No) stop("No does not match")
  if (V$GENERAL$windowSize != v$GENERAL$windowSize) stop("windowSize does not match")
  if (V$GENERAL$islandSizes != v$GENERAL$islandSizes) stop("islandSizes does not match")
  
  if (!all.equal(names(V$PRIORS), names(v$PRIORS))) stop("PRIORS elements do not match")
  if (!all.equal(names(V$SFS), names(v$SFS))) stop("SFS elements do not match")
  
  Vv <- list()
  class(Vv) <- "validationTable"
  Vv$GENERAL <- V$GENERAL
  Vv$GENERAL$nSimul <- V$GENERAL$nSimul + v$GENERAL$nSimul
  
  SFS <- lapply(seq_along(V$SFS), function(x) rbind(V$SFS[[x]], v$SFS[[x]]))
  names(SFS) <- names(V$SFS)
  
  PRIORS <- lapply(seq_along(V$PRIORS), function(x) rbind(V$PRIORS[[x]], v$PRIORS[[x]]))
  names(PRIORS) <- names(V$PRIORS)
  
  Vv$PRIORS <- PRIORS
  Vv$SFS <- SFS
  
  rm(list=c("SFS","PRIORS","V","v")); invisible(gc())
  return(Vv)
}



