
accuracy_f_dist_full <- function(reftb, nSimul, cutoff = 10, minSNP, tolABC = .01, 
                                 mean = NULL,
                                 sweepingIsl = NULL, sweepPos, recRate, ...) {
  
  if (is.null(sweepingIsl)||length(sweepingIsl)!=1) stop("only one single island index!")
  if (recRate[[1]]!="discr") stop("recRate must be _discr_")
  
  M <- simulate_obs(reftb = reftb,
                    nSimul = nSimul,
                    sweepingIsl = sweepingIsl,
                    sweepPos = sweepPos,
                    recRate = recRate, ...)
  M <<- M
  cat("\n\nScan\n\n")
  Ms <- scann(M, reftb, cutoff = cutoff, minSNP = minSNP, tolABC = tolABC)
  
  if (F) {
    # mean on BFs
    est <- Ms$estimates
    est[est$sweepingIsland!=paste("i",sweepingIsl,sep=""),"BayesFactor"] <- .1
    df <- data.frame(pos=Ms$trues$sweepPos * reftb$GENERAL$windowSize, BF=est$BayesFactor, recRate=Ms$trues$recRate)
    df <- aggregate(df, list(df$pos, df$recRate), mean); df$Group.1 <- df$Group.2 <- NULL
    df$recRate <- as.factor(df$recRate)
  } else {
    tmp <- data.frame(t=Ms$trues$sweepingIsland, e=Ms$estimates$sweepingIsland, bf=Ms$est$BayesFactor,
                      stringsAsFactors=F)
    tmp$y <- 0
    tmp[tmp$t==tmp$e & !is.na(tmp$bf) & tmp$bf>cutoff,"y"] <- 1
    df <- data.frame(pos=Ms$trues$sweepPos * reftb$GENERAL$windowSize, BF=tmp$y, recRate=Ms$trues$recRate)
    df <- aggregate(df, list(df$pos, df$recRate), mean); df$Group.1 <- df$Group.2 <- NULL
    df$recRate <- as.factor(df$recRate)
  }
  
  # ggplot
  g <- ggplot(data=df, aes(x=pos, y=BF, group=recRate, col=recRate)) +  
    geom_hline(yintercept = .1, size = .7, col = "black") +
    geom_point(size = 3, alpha=.7) + 
    geom_line(size = 1.2, alpha=.7) +
    scale_y_log10() +
    xlab("Distance (bp)") + ylab("") + ggtitle(paste("Focal island:",sweepingIsl))
  
  
  # ESTIMATION
  df2 <- data.frame(pos=Ms$trues$sweepPos * reftb$GENERAL$windowSize,
                    TM=Ms$trues$sweepingIsland,
                    EM=Ms$estimates$sweepingIsland,
                    t=Ms$trues$sweepAge, 
                    e=Ms$estimates$sweepAge, 
                    recRate=Ms$trues$recRate, stringsAsFactors=F)
  df2 <- na.omit(df2)
  df2 <- subset(df2, TM!="i0" & TM==EM)
  df2$recRate <- as.factor(df2$recRate)
  
  require(plyr)
  DF2 <- ddply(df2, .(pos,recRate), function(x) est_accuracy(x[,c("t","e")], mean=mean))
  
  q <- ggplot(data=DF2, aes(x=pos, y=V1, group=recRate, col=recRate)) +  
    geom_point(size = 3, alpha=.7) + 
    geom_line(size = 1.2, alpha=.7) +
    xlab("Distance (bp)") + ylab("NRMSError") + ggtitle(paste("Focal island:",sweepingIsl))
  
  return(list("raw"=Ms, "detection"=g, "estimation"=q, "det_DF"=df, "est_DF"=DF2))
}

