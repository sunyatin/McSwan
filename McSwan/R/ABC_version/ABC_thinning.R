
#' @title Merge contiguous sweeps
#' @description This function merges contiguous regions detected by \code{\link{gscan}} as having experienced a selective sweep.
#' @param x a dataframe returned by \code{\link{gscan}}
#' @param reftb an initialized \code{referenceTable} object
#' @param bwadjust (numeric) bandwidth coefficient, the algorithm will extend contiguity as far as \code{bwadjust * windowSize} from any focal window
#' @param minWindows (integer) the minimum number of contiguous windows required to perform contiguity search
#' @param summary_stat (string) the central statistic to average the parameter estimates of merged contiguous windows ("mean" recommended)
#' @param plot_thinning (logical) plot a representation of the genomic fragment with the contiguous sweep regions
#' @return Returns a dataframe of the merged contiguous sweeps. One sweep per line, and in columns:
#' \itemize{ 
#' \item \code{\bold{sweep.center}} center of the contiguous sweep region
#' \item \code{\bold{sweep.lbound}} first position of the contiguous sweep region
#' \item \code{\bold{sweep.rbound}} last position of the contiguous sweep region
#' \item \code{\bold{BF}} aggregated Bayes Factor of the contiguous sweep region
#' \item \code{\bold{deme}} population detected to have experienced a selective sweep (given as "i" & index of the population in the \emph{MS} command)
#' \item \code{\bold{sweepAge}} posterior estimate for the sweep age
#' \item \code{\bold{sweepAge.IC.low}} lower boundary of the 95% credible interval for the sweep age estimation
#' \item \code{\bold{sweepAge.IC.up}} upper boundary of the 95% credible interval for the sweep age estimation
#' }
#' @seealso \code{\link{gscan}}
#' @export
#' @examples Please refer to the vignette.
ABC_thin <- function(x, reftb, bwadjust = 1, minWindows = 1, summary_stat = "mean", plot_thinning = FALSE) {
  
  # internally set options
  stat = "mean" # param estimates and BFs are averaged over the merged contiguous regions (otherwise we extract the central window estimates)
  epsilon = 1e-15
  method = "contig" # contiguity detection method, either "contig" (recommended) or "density"
  
  Y <- thin_internal(x, reftb, method = method, stat = stat, bwadjust = bwadjust, minWindows = minWindows, plot_thinning = plot_thinning, epsilon = epsilon, summary_stat = summary_stat)
  
  return(Y)
  
  template <- x$template
  PAC <- x$obsData
  
  H <- Y$estimation
  if (is.null(H) || nrow(H)==0) return(Y)
  for (r in 1:nrow(H)) {
    # get local observed SFS
    oSFS <- template
    localPAC <- subset(PAC, POS >= H[r,"sweep.lbound"] & POS < H[r,"sweep.rbound"]) # < only
    tmp <- table(subset(localPAC, nSNP == 1)$PAC)
    oSFS[names(tmp)] <- tmp
    tmp <- table(subset(localPAC, nSNP == .5)$PAC)
    oSFS[names(tmp)] <- oSFS[names(tmp)] + tmp / 2
    
    deme <- as.character(H[r,"deme"])
    
    lb <- t(apply(reftb$PRIORS[[deme]], 2, range))
    ss <- reftb$DIMREDUC$PLS[[deme]]$scores
    tproj <- project_target(reftb, oSFS, method="PLS", focalIsland = deme)
    colnames(ss) <- colnames(tproj) <- 1:ncol(ss)
    
    par <- reftb$PRIORS[[deme]]
    
    aa <- suppressWarnings(abc::abc(target = tproj,
                                param = par,
                                sumstat = ss,
                                tol = .005,
                                method = "loclinear",
                                hcorr = TRUE,
                                transf = "logit",
                                logit.bounds = lb))
    
    H[r,"sweepAge"] <- getmode(aa$adj.values[,1])
    H[r,"sweepAge.IC.low"] <- quantile(aa$adj.values[,1], .025, na.rm=T)
    H[r,"sweepAge.IC.up"] <- quantile(aa$adj.values[,1], .975, na.rm=T)
  }
  Y$estimation <- H
  
  return(Y)
}

#' @title Thinning - internal procedure
#' @keywords internal
ABC_thin_internal <- function(x, reftb, method = "contig", stat = "mean", bwadjust = 1, minWindows = 1, plot_thinning = F, epsilon = 1e-15, summary_stat = "mean") {
  ## output==list(density=NULL, estimation=NULL) means ' all is "i0" ' (by Ho-conservatism)

  AVALL <- x$adjVals
  x <- x$res
  
  RET <- list("density"=c(), "estimation"=c())
  demes <- names(reftb$SFS)
  demes <- demes[demes!="i0"]
  for (d in seq_along(demes)) {
    dn <- demes[d]
  
    #y <- subset(x, !is.na(deme.under.sel) & deme.under.sel == dn)
    INDICES <- with(x, !is.na(deme.under.sel) & deme.under.sel == dn)
    y <- x[INDICES,]
    AV <- AVALL[INDICES]
    
    gr <- c(min(x$start.pos), max(x$end.pos))
    ymid <- (y$start.pos+y$end.pos)/2
    L <- reftb$GENERAL$windowSize
    
    if (is.null(nrow(y))||nrow(y)<minWindows) next()
    
    if (method=="density") {
      
      if (is.null(bwadjust)) {
        D <- density( ymid,
                      from = gr[1],
                      to = gr[2])
      } else {
        D <- density( ymid,
                      bw = reftb$GENERAL$windowSize * bwadjust,
                      from = gr[1],
                      to = gr[2])
      }
      
      if (plot_thinning) {
        plot(D, main="Gaussian spline of sweeped regions")
        points(ymid, rep(0, length(ymid)), pch="|")
      }
      
      D$y[D$y<epsilon] <- 0
      
      # derivate
      dD <- list(x=D$x[-1], y=diff(D$y) / diff(D$x))
      
      # find near-zero intersection
      dDy <- dD$y
      zI <- sapply(2:(length(dDy)-1), function(j) dDy[j-1]>0&&dDy[j]>=0&&dDy[j+1]<0 ) # local maximum
      zI <- c(FALSE, zI, FALSE)
      
      xleft <- dD$x[zI]
      xright <- dD$x[c(FALSE, zI[-length(zI)])]
      xP <- ( xleft + xright ) / 2
      
      if (length(xP)==0) {
        adret <- list("density"=cbind("x"=NA, "y"=NA), "estimation"=rep(NA, 7))
      } else {
      
        # discardOrphanWindow? => still missing in the code
        idx <- sapply(xP, function(p) {
          which.min( abs(ymid - p) )
        })
        
        dy <- D$y
        lims <- sapply(which(c(F,zI)), function(p) {
          lm <- which(dy[1:p]<=epsilon)
          lm <- lm[length(lm)]
          
          rm <- which(dy[p:length(dy)]<=epsilon) + (p-1)
          rm <- rm[1]
          
          L <- which.min( abs(ymid - D$x[lm]) )
          L <- ifelse(length(L)==0, min(y$start.pos), y$start.pos[L])
          R <- which.min( abs(ymid - D$x[rm]) )
          R <- ifelse(length(R)==0, max(y$end.pos), y$end.pos[R])
          
          return(c(L, R))
        })
        lims <- t(lims)
        
        ret <- data.frame(sweep.center=xP, 
                          sweep.lbound = lims[,1],
                          sweep.rbound = lims[,2],
                          deme=dn,
                          sweepAge=y$sweepAge[idx], 
                          sweepAge.IC.low=y$sweepAge.IC.low[idx],
                          sweepAge.IC.up=y$sweepAge.IC.up[idx])
        #return(list("density"=cbind("x"=D$x, "y"=D$y), "estimation"=ret))
  
        adret <- list("density"=cbind("x"=D$x, "y"=D$y), "estimation"=ret)
        
      }
      
    } else if (method=="contig") {
      
      store <- TRUE
      startp <- y[1,"start.pos"]
      age <- age.low <- age.up <- bf <- XXX <- c()
      ret <- c()
      for (i in 1:nrow(y)) {
        
        if (i>1) store <- ifelse( y[i,"start.pos"] - y[i-1,"start.pos"] > bwadjust * L, FALSE, TRUE)
        
        if (!store) {
          if (stat=="mean") {
            pe.Age <- do.call(summary_stat, list(age))
            pe.lAge <- do.call(summary_stat, list(age.low))
            pe.rAge <- do.call(summary_stat, list(age.up))
            pe.BF <- mean(bf)
            pe.XXX <- do.call(summary_stat, list(XXX))
          } else {
            #ix <- which.min( abs(y$start.pos+y$end.pos)/2 - (startp+endp)/2 )
            ix <- order(abs(y$start.pos+y$end.pos)/2 - (startp+endp)/2)
            #ix <- ix[1:3]
            ix <- ix[1] #first closer
            pe.Age <- do.call(summary_stat, list(age[ix]))
            pe.lAge <- do.call(summary_stat, list(age.low[ix]))
            pe.rAge <- do.call(summary_stat, list(age.up[ix]))
            pe.BF <- do.call(summary_stat, list(bf[ix]))
            pe.XXX <- pe.Age
          }
          add <- data.frame(sweep.center = (startp+endp)/2, 
                            sweep.lbound = startp,
                            sweep.rbound = endp,
                            BF = pe.BF,
                            deme=dn, 
                            sweepAge = pe.Age,
                            sweepAgeSynth = pe.XXX,
                            sweepAge.IC.low = pe.lAge,
                            sweepAge.IC.up = pe.rAge)
          ret <- rbind(ret, add)
          startp <- y[i,"start.pos"]
          age <- age.low <- age.up <- bf <- XXX <- c()
        }
        
        endp <- y[i,"end.pos"]
        age <- c(age, y[i,"sweepAge"])
        age.low <- c(age.low, y[i,"sweepAge.IC.low"])
        age.up <- c(age.up, y[i,"sweepAge.IC.up"])
        bf <- c(bf, y[i,"BayesFactor"])
        XXX <- c(XXX, unlist(AV[[i]]))
      }
      
      if (store) {
        if (stat=="mean") {
          pe.Age <- do.call(summary_stat, list(age))
          pe.lAge <- do.call(summary_stat, list(age.low))
          pe.rAge <- do.call(summary_stat, list(age.up))
          pe.BF <- mean(bf)
          pe.XXX <- do.call(summary_stat, list(XXX))
        } else {
          #ix <- which.min( abs(y$start.pos+y$end.pos)/2 - (startp+endp)/2 )
          ix <- order(abs(y$start.pos+y$end.pos)/2 - (startp+endp)/2)
          #ix <- ix[1:3]
          ix <- ix[1] #first closer
          pe.Age <- do.call(summary_stat, list(age[ix]))
          pe.lAge <- do.call(summary_stat, list(age.low[ix]))
          pe.rAge <- do.call(summary_stat, list(age.up[ix]))
          pe.BF <- do.call(summary_stat, list(bf[ix]))
          pe.XXX <- pe.Age
        }
        add <- data.frame(sweep.center = (startp+endp)/2, 
                          sweep.lbound = startp,
                          sweep.rbound = endp,
                          BF = pe.BF,
                          deme=dn, 
                          sweepAge = pe.Age,
                          sweepAgeSynth = pe.XXX,
                          sweepAge.IC.low = pe.lAge,
                          sweepAge.IC.up = pe.rAge)
        ret <- rbind(ret, add)
      }
      
      rgg <- c(x$start.pos[1], x$end.pos[nrow(x)])
      ymid <- (y$start.pos+y$end.pos)/2
      if (plot_thinning) {
        plot(ymid, rep(0, length(ymid)), pch="|", type="p", ylim=c(0,1), xlim=rgg)
        segments(x0=ret$sweep.center, y0=0, y1=1, col="gray", lwd=1.5)
        rect(xleft=ret$sweep.lbound, ybottom=.1, xright=ret$sweep.rbound, ytop=.2, col="lightgray")
        abline(v=diff(rgg)*.5, col="red")
      }
      
      #return(list("density"=cbind("x"=ret$sweep.center, "y"=ret$BF), "estimation"=ret))
      adret <- list("density"=cbind("x"=ret$sweep.center, "y"=ret$BF), "estimation"=ret)
      
    } else {
      stop("the method you provided does not exist")
    }
   
    RET$density <- rbind(RET$density, adret$density)
    RET$estimation <- rbind(RET$estimation, adret$estimation)
  }

  if (method=="contig" && !is.null(RET$estimation)) {
    euRows <- with(RET$estimation, which(sweep.rbound-sweep.lbound+1 >= minWindows*L))
    RET$estimation <- RET$estimation[euRows,]
    RET$density <- RET$density[euRows,]
  }

RET$estimation$sweepAgeSynth <- NULL # 24 01 2017
  return(RET$estimation)
}

