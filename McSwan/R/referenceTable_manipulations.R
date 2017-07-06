
#' @title Combine reference tables prior to dimension reduction


#' @title Combine McSwan objects
#' @export
#' @keywords internal
combine <- function(x, ...) UseMethod("combine")

#' @title Combine referenceTable objects
#' @description Combine two \code{referenceTable} objects prior to machine learning (i.e. prior to the \code{\link{dim_reduction}} call).
#' @param V a first \code{referenceTable} object with empty \code{DIMREDUC} slot
#' @param v a second \code{referenceTable} object with empty \code{DIMREDUC} slot
#' @return A \code{referenceTable} object combining the two original \code{referenceTable} objects.
#' @seealso \code{\link{generate_priors}}, \code{\link{coalesce}}, \code{\link{dim_reduction}}
#' @export
combine.referenceTable <- function(V, v) {
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






