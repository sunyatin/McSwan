#plotting functions

plow <- function(pos.mcswan, score.mcswan, log10_mcswan=F,
                 pos.ref, score.ref, ref.show.above = NULL,
                 gff.start.pos.CDS = NULL,
                 gff.end.pos.CDS = NULL,
                 mcswan.signif.circles = NULL,
                 ref.signif.circles = NULL,                 
                 n_rays = 20,
                 do_radial = TRUE) {
  # CAREFUL! REFERENCE should be formatted so that higher value of scores
  # mean "most significant sweep"
  
  min.y <- -1.2
  max.y <- 1.5
  
  # automatically fits ref pos range in reference to McSwan pos range
  
  # prepare McSwan
  M <- data.frame(pos=pos.mcswan, score=score.mcswan)
  M <- na.omit(M)
  max.mcswan <- max(M$score[M$score!=Inf], na.rm=T)
  M$score[M$score==Inf] <- max.mcswan + 10
  max.mcswan <- max(M$score, na.rm=T)
  if (log10_mcswan==TRUE) {
    M$score[M$score<1] <- 1
    M$score <- log10(M$score)
    max.mcswan <- max(M$score, na.rm=T)
    if (!is.null(mcswan.signif.circles)) {
      mcswan.signif.circles <- log10(mcswan.signif.circles) / max.mcswan
    }
  }
  M$score <- M$score / max.mcswan
  
  # prepare ref
  R <- data.frame(pos=pos.ref, score=score.ref)
  R <- na.omit(R)
  R <- subset(R, pos>=min(M$pos) & pos<=max(M$pos))
  max.ref <- max(R$score, na.rm=T)
  R$score <- R$score / max.ref
  ref.signif.circles <- ref.signif.circles / max.ref
  
  print(summary(R$score))
  
  # min max pos
  pos.minmax <- range(c(M$pos, R$pos), na.rm=T)
  
  # vertical grids
  grid.x <- seq(pos.minmax[1], pos.minmax[2], length.out=n_rays)
  gridd <- data.frame(x=grid.x, xend=grid.x, y=min.y, yend=max.y)
  
  # position ticks
  pos.minor <- seq(pos.minmax[1], pos.minmax[2]+1, by=1e6)
  height.pos.minor <- .02
  pos.major <- seq(pos.minmax[1], pos.minmax[2]+1, by=1e7)
  height.pos.major <- .05
  # transform the tick dataframes
  dminor <- data.frame(position=pos.minor, height=height.pos.minor)
  dmajor <- data.frame(position=pos.major, height=height.pos.major, 
                       #label=paste(1:length(pos.major),"M",sep=""))
                       label=paste(round(pos.major/1e6, digits=0),"",sep="")) #!!!!!!
  
  # subset ref
  if (!is.null(ref.show.above)) R <- subset(R, score >= ref.show.above)
  
  print("ggplot!")
  g=ggplot(data=M) + 
    geom_segment(aes(x=pos, xend=pos, y=0, yend=score), col="#4D9FC8") +
    geom_segment(data=R, 
                 aes(x=pos, xend=pos, y=0, yend=-score), col="#B84C53") +
    geom_rect(xmin=pos.minmax[1], xmax=pos.minmax[2], 
              ymin=-.01, ymax=.01, col="white", fill="black", size=.2) +
    #geom_line(data=subset(R, val>.1), aes(x=pos, y=-val), col="orange") +
    theme_bw() +
    
    geom_segment(data=gridd, aes(x=x, xend=xend, y=y, yend=yend), 
                 col="gray", size=.1, linetype=2) +
    
    scale_x_continuous(limits=c(pos.minmax[1], pos.minmax[2]+1e7)) +
    
    scale_y_continuous(limits=c(min.y, max.y)) +
    
    theme(axis.text.x=element_blank()) +
    theme(axis.line.x=element_blank()) +
    theme(axis.title.x=element_blank()) +
    theme(panel.border=element_blank()) +
    theme(axis.ticks=element_blank()) +
    theme(legend.position="none") +
    theme(panel.grid.major.x=element_blank()) +
    theme(panel.grid.minor.x=element_blank()) +
    theme(axis.line.y=element_blank()) +
    theme(axis.text.y=element_blank()) +
    theme(axis.ticks=element_blank()) + xlab("") + ylab("")
  
  # position ruler: MINOR ticks
  y_ruler <- 1.31
  g <- g + geom_segment(data=dminor, 
                        aes(x=position, xend=position, 
                            y=y_ruler, yend=height+y_ruler))
  
  # position ruler: MAJOR ticks
  g <- g + geom_segment(data=dmajor, 
                        aes(x=position, xend=position, 
                            y=y_ruler, yend=height+y_ruler)) +
    geom_text(data=dmajor, size=2, alpha=.8, family="Times", 
              aes(x=position, y=height+.04+y_ruler, label=label))
  
  # significativity circles
  if (!is.null(mcswan.signif.circles)) {
    for ( h in seq_along(mcswan.signif.circles)) {
      g <- g + geom_segment(x=pos.minmax[1], xend=pos.minmax[2],
                            y=mcswan.signif.circles[h], yend=mcswan.signif.circles[h],
                            col="#698999", alpha=.6, size=.08, linetype=6)
    }
  }
  if (!is.null(ref.signif.circles)) {
    for ( h in seq_along(ref.signif.circles)) {
      g <- g + geom_segment(x=pos.minmax[1], xend=pos.minmax[2],
                            y=-ref.signif.circles[h], yend=-ref.signif.circles[h],
                            col="#780A0D", alpha=.6, size=.08, linetype=6)
    }
  }
  
  # GFF!
  if (!is.null(gff.start.pos.CDS)) {
    print("GFF")
    GFF <- data.frame(start=gff.start.pos.CDS, end=gff.start.pos.CDS)
    GFF <- subset(GFF, start>=pos.minmax[1] & end<=pos.minmax[2])
    g <- g +
      geom_rect(data=GFF, aes(xmin=start, xmax=end), 
                ymin=1.1, ymax=1.3, col="#767676", fill="#808080", size=.1)
  }
  
  if (do_radial) g <- g + coord_polar()
  
  return(g)
}




plot_barcodeRing <- function(mcswan.DF,
                 focal.deme, 
                 log10_mcswan=F,
                 
                 pybus.pos,
                 pybus.score,
                 pybus.log10 = F,
                 
                 sweed.pos,
                 sweed.score,
                 sweed.log10 = F,
                 
                 mcswan.width = 50,
                 yfactor = 100,
                 min.y = -200,
                 max.y = 100,
                 
                 plot.ticks = F,
                 
                 ref.show.above = NULL,
                 gff.start.pos.CDS = NULL,
                 gff.end.pos.CDS = NULL,
                 mcswan.signif.circles = NULL,
                 ref.signif.circles = NULL,                 
                 n_rays = 20,
                 do_radial = TRUE) {
  
  # min max pos
  pos.minmax <- range(c(mcswan.DF$sweep.lbound, mcswan.DF$sweep.rbound, pybus.pos, sweed.pos), na.rm=T)
  
  # vertical grids
  grid.x <- seq(pos.minmax[1], pos.minmax[2], length.out=n_rays)
  gridd <- data.frame(x=grid.x, xend=grid.x, y=min.y, yend=max.y)
  
  # position ticks
  pos.minor <- seq(pos.minmax[1], pos.minmax[2]+1, by=1e6)
  height.pos.minor <- yfactor * .01
  pos.major <- seq(pos.minmax[1], pos.minmax[2]+1, by=1e7)
  height.pos.major <- yfactor * .02
  # transform the tick dataframes
  dminor <- data.frame(position=pos.minor, height=height.pos.minor)
  dmajor <- data.frame(position=pos.major, height=height.pos.major, 
                       #label=paste(1:length(pos.major),"M",sep=""))
                       label=paste(round(pos.major/1e6, digits=0),"",sep="")) #!!!!!!
  
  # prepare dataframes
  
  # McSwan
  M <- subset(mcswan.DF, deme==focal.deme)
  # Pybus
  pybus <- data.frame(x=pybus.pos, y=pybus.score)
  if (pybus.log10) pybus$y <- log10(pybus$y)
  pybus$y <- (-mcswan.width/2) + -1 * pybus$y/max(pybus$y) * yfactor
  # SweeD
  sweed <- data.frame(x=sweed.pos, y=sweed.score)
  if (sweed.log10) sweed$y <- log10(sweed$y)
  sweed$y <- (mcswan.width/2) + sweed$y/max(sweed$y) * yfactor
  
  YMAX <- ifelse(max.y<max(sweed$y), max(sweed$y), max.y)

  # PLOOOOOT!
  print("ggplot!")
  g=ggplot(data=M) + 
    # pybus
# 705D90
    geom_line(data=pybus, aes(x=x, y=y), size=.5, col="#EE9A00") +
    # sweed
    geom_line(data=sweed, aes(x=x, y=y), col="#B84C53") +
    # McSwan
    geom_rect(aes(xmin=sweep.lbound, xmax=sweep.rbound, ymin=-mcswan.width/2+.013*yfactor, ymax=mcswan.width/2-.013*yfactor, fill=BF)) +

    theme_bw() +
    
    geom_segment(data=gridd, aes(x=x, xend=xend, y=y, yend=yend), 
                 col="gray", size=.1, linetype=2) +
    
    scale_x_continuous(limits=c(pos.minmax[1], pos.minmax[2]+1e7)) +
    
    scale_y_continuous(limits=c(min.y, YMAX+1)) +
    
    scale_fill_distiller("BF", palette="Blues", direction=1) +
    
    theme(axis.text.x=element_blank()) +
    theme(axis.line.x=element_blank()) +
    theme(axis.title.x=element_blank()) +
    theme(panel.border=element_blank()) +
    theme(axis.ticks=element_blank()) +
    #theme(legend.position="none") +
    theme(panel.grid.major.x=element_blank()) +
    theme(panel.grid.minor.x=element_blank()) +
    theme(axis.line.y=element_blank()) +
    theme(axis.text.y=element_blank()) +
    theme(axis.ticks=element_blank()) + xlab("") + ylab("")
  
  # position ruler: MINOR ticks
  if (plot.ticks) {
  y_ruler <- max(sweed$y)
  g <- g + geom_segment(data=dminor, 
                        aes(x=position, xend=position, 
                            y=y_ruler, yend=height+y_ruler))

  # position ruler: MAJOR ticks
  g <- g + geom_segment(data=dmajor, 
                        aes(x=position, xend=position, 
                            y=y_ruler, yend=height+y_ruler))
  }
  
  # position annotation
  g <- g + geom_text(data=dmajor, size=2, alpha=.8, 
              aes(x=position, y=height-.1*yfactor+y_ruler, label=label))
  
  # significativity circles
  if (!is.null(mcswan.signif.circles)) {
    for ( h in seq_along(mcswan.signif.circles)) {
      g <- g + geom_segment(x=pos.minmax[1], xend=pos.minmax[2],
                            y=mcswan.signif.circles[h], yend=mcswan.signif.circles[h],
                            col="#698999", alpha=.6, size=.08, linetype=6)
    }
  }
  if (!is.null(ref.signif.circles)) {
    for ( h in seq_along(ref.signif.circles)) {
      g <- g + geom_segment(x=pos.minmax[1], xend=pos.minmax[2],
                            y=-ref.signif.circles[h], yend=-ref.signif.circles[h],
                            col="#780A0D", alpha=.6, size=.08, linetype=6)
    }
  }
  
  # GFF!
  if (!is.null(gff.start.pos.CDS)) {
    print("GFF")
    GFF <- data.frame(start=gff.start.pos.CDS, end=gff.start.pos.CDS)
    GFF <- subset(GFF, start>=pos.minmax[1] & end<=pos.minmax[2])
    GFF.ymin <- max(sweed$y) + .1*yfactor
    GFF.ymax <- GFF.ymin + .1*yfactor
    g <- g +
      geom_rect(data=GFF, aes(xmin=start, xmax=end), 
                ymin=GFF.ymin, ymax=GFF.ymax, col="#767676", fill="#808080", size=.1)
  }
  
  if (do_radial) g <- g + coord_polar()
  
  return(g)
}



plot_1bubbleRing <- function(mcswan.DF,
                      focal.deme, 
                      log10_mcswan=F,
                      
                      pybus.pos,
                      pybus.score,
                      pybus.log10 = F,
                      pybus.signif = NULL,
                      
                      sweed.pos,
                      sweed.score,
                      sweed.log10 = F,
                      sweed.signif = NULL,
                      
                      mcswan.width = 50,
                      yfactor = 100,
                      min.y = -200,
                      max.y = 100,
                      
                      plot.ticks = F,
                      
                      ref.show.above = NULL,
                      gff.start.pos.CDS = NULL,
                      gff.end.pos.CDS = NULL,
                      mcswan.signif.circles = NULL,
                      ref.signif.circles = NULL,                 
                      n_rays = 20,
                      do_radial = TRUE) {
  
  # min max pos
  pos.minmax <- range(c(mcswan.DF$sweep.lbound, mcswan.DF$sweep.rbound, pybus.pos, sweed.pos), na.rm=T)
  
  # vertical grids
  grid.x <- seq(pos.minmax[1], pos.minmax[2], length.out=n_rays)
  gridd <- data.frame(x=grid.x, xend=grid.x, y=min.y, yend=max.y)
  
  # position ticks
  pos.minor <- seq(pos.minmax[1], pos.minmax[2]+1, by=1e6)
  height.pos.minor <- yfactor * .01
  pos.major <- seq(pos.minmax[1], pos.minmax[2]+1, by=1e7)
  height.pos.major <- yfactor * .02
  # transform the tick dataframes
  dminor <- data.frame(position=pos.minor, height=height.pos.minor)
  dmajor <- data.frame(position=pos.major, height=height.pos.major, 
                       #label=paste(1:length(pos.major),"M",sep=""))
                       label=paste(round(pos.major/1e6, digits=0),"",sep="")) #!!!!!!
  
  # prepare dataframes
  
  # McSwan
  M <- subset(mcswan.DF, deme==focal.deme)
  # Pybus
  pybus <- data.frame(x=pybus.pos, y=pybus.score)
  if (pybus.log10) pybus$y <- log10(pybus$y)
  maxPybus <- max(pybus$y)
  pybus$y <- (-mcswan.width/2) + -1 * pybus$y/maxPybus * yfactor
  # SweeD
  sweed <- data.frame(x=sweed.pos, y=sweed.score)
  maxSweeD <- max(sweed$y)
  if (sweed.log10) sweed$y <- log10(sweed$y)
  sweed$y <- (mcswan.width/2) + sweed$y/maxSweeD * yfactor
  
  YMAX <- ifelse(max.y<max(sweed$y), max(sweed$y), max.y)
  
  # PLOOOOOT!
  print("ggplot!")
  g=ggplot(data=M) + 
    
    # sweed
    geom_line(data=sweed, aes(x=x, y=y), col="#B84C53")
    
    # pybus
    # 705D90
  #if (!any(is.null(pybus.pos, pybus.score))) {
    g <- g + geom_line(data=pybus, aes(x=x, y=y), size=.5, col="#EE9A00")
  #} else {
    # nothing
  #}

    # McSwan bubbles
    #geom_point(aes(x=sweep.lbound/2+sweep.rbound/2, y=0, size=BF, col=BF), alpha=.6) +
    g <- g + geom_point(aes(x=sweep.lbound/2+sweep.rbound/2, y=0, size=BF), col="#005A9C", alpha=.6) +
    
    # McSwan sweep range
    geom_rect(aes(xmin=sweep.lbound, xmax=sweep.rbound, ymin=-mcswan.width/2+.013*yfactor, ymax=mcswan.width/2-.013*yfactor), fill="black", alpha=.8) +
    
    theme_bw() +
    
    geom_segment(data=gridd, aes(x=x, xend=xend, y=y, yend=yend), 
                 col="gray", size=.1, linetype=2) +
    
    scale_x_continuous(limits=c(pos.minmax[1], pos.minmax[2]+1e7)) +
    
    scale_y_continuous(limits=c(min.y, YMAX+1)) +
    
    scale_fill_distiller("BF", palette="Blues", direction=1) +
    
    theme(axis.text.x=element_blank()) +
    theme(axis.line.x=element_blank()) +
    theme(axis.title.x=element_blank()) +
    theme(panel.border=element_blank()) +
    theme(axis.ticks=element_blank()) +
    #theme(legend.position="none") +
    theme(panel.grid.major.x=element_blank()) +
    theme(panel.grid.minor.x=element_blank()) +
    theme(axis.line.y=element_blank()) +
    theme(axis.text.y=element_blank()) +
    theme(axis.ticks=element_blank()) + xlab("") + ylab("")
  
  # position ruler: MINOR ticks
  if (plot.ticks) {
    y_ruler <- max(sweed$y)
    g <- g + geom_segment(data=dminor, 
                          aes(x=position, xend=position, 
                              y=y_ruler, yend=height+y_ruler))
    
    # position ruler: MAJOR ticks
    g <- g + geom_segment(data=dmajor, 
                          aes(x=position, xend=position, 
                              y=y_ruler, yend=height+y_ruler))
  }
  
  # position annotation
  g <- g + geom_text(data=dmajor, size=2, alpha=.8, 
                     aes(x=position, y=height-.1*yfactor+y_ruler, label=label))
  
  # significativity circles
  if (!is.null(mcswan.signif.circles)) {
    for ( h in seq_along(mcswan.signif.circles)) {
      g <- g + geom_segment(x=pos.minmax[1], xend=pos.minmax[2],
                            y=mcswan.signif.circles[h], yend=mcswan.signif.circles[h],
                            col="#698999", alpha=.6, size=.08, linetype=6)
    }
  }
  if (!is.null(pybus.signif)) {
    for ( h in seq_along(pybus.signif)) {
      H <- (-mcswan.width/2) + -1 * pybus.signif[h]/maxPybus * yfactor
      g <- g + geom_segment(x=pos.minmax[1], xend=pos.minmax[2],
                            y=H, yend=H,
                            col="darkgray", alpha=.6, size=.08, linetype=6)
    }
  }
  if (!is.null(sweed.signif)) {
    for ( h in seq_along(sweed.signif)) {
      H <- (mcswan.width/2) + sweed.signif[h]/maxSweeD * yfactor
      g <- g + geom_segment(x=pos.minmax[1], xend=pos.minmax[2],
                            y=H, yend=H,
                            col="darkgray", alpha=.6, size=.08, linetype=6)
    }
  }
  
  # GFF!
  if (!is.null(gff.start.pos.CDS)) {
    print("GFF")
    GFF <- data.frame(start=gff.start.pos.CDS, end=gff.start.pos.CDS)
    GFF <- subset(GFF, start>=pos.minmax[1] & end<=pos.minmax[2])
    GFF.ymin <- max(sweed$y) + .1*yfactor
    GFF.ymax <- GFF.ymin + .1*yfactor
    g <- g +
      geom_rect(data=GFF, aes(xmin=start, xmax=end), 
                ymin=GFF.ymin, ymax=GFF.ymax, col="#767676", fill="#808080", size=.1)
  }
  
  if (do_radial) g <- g + coord_polar()
  
  return(g)
}



