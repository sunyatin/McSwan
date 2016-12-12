# demonstration kit
#require(caret); require(MASS); require(abc); require(pls); require(ggplot2); require(reshape2); require#(phyclust); require(RColorBrewer); require(McSwan); require(scales); require(sitools)


setwd("C:/Users/Windows/Desktop/McSwan_public/")

set_session(tempDir = "temp", javaPath = "C:/Program Files (x86)/Java/jre1.8.0_101/bin/javaw.exe")


mut <- 2.5e-8
No <- 4106
ms <- paste("80 1 -t",4*No*mut,
            "-I 2 40 40",
            "-en 0 2",58307/No,
            "-en",70/(4*No),2,23738/No,
            "-ej",2000/(4*No),"1 2",
            "-en",2683/(4*No),2,7083/No,
            "-en",2783/(4*No),2,7963/No)


reftb <- generate_priors(msDemography = ms, fold = F,
                         No = No,
                         windowSize = 2e4,
                         nSimul = 200) #recRate c'est la fraction des inds qui sweepent PAS


tt <- Sys.time()
    reftb <- coalesce(reftb, method="NULL", verbose=T)
print(Sys.time() - tt)


reftb <- dim_reduction(reftb, doPLSrecRate=F, nonvarTolerance = 1e-4, QR_multicollinearity=T,
                       relativeSFS = TRUE, minSDlda = 0, tolLDA = 1e-9)


# obs

#convert_VCF("../data/set the ancestral allele/CEU-20_LWK-20_POLARIZED.vcf", "../data/pops_CEU-20_LWK-20.txt", reftb, outPath="../temp/UNFOLD_POLARIZED.txt", chromosome=2, minQUAL=30)

obs <- get_SFS("../temp/UNFOLD_POLARIZED.txt", ms)

minPos <- NULL
maxPos <- NULL

print(Sys.time())
OA <- gscan(obs, reftb,
            minSNP = 20, # 20
            startPos = minPos,
            lastPos = maxPos,
            cutoff = 10, # 10 -- 2209_Homo2gfit0 = 5
            tolABC = c(.01, .005), #c(.01, .02)
            windowSlide = NULL,
            plot_simCloud = F,
            tolGFIT = 0, # 0.025 -- 2209_Homo2gfit0 = 0
            doESTIMATION = T)
print(Sys.time())

OR <- OA$res

save(OA, file="0810_UNFOLDED_obsHomo_VCFpolarized_relativeSFS_gfit0-tolLDA(1e-9)_allChr2")

OB <- thin(OA, reftb, method="contig", stat = "mean", bwadjust = 12, minWindows = 3, plot_thinning = F, summary_stat = "mean", obs = obs) # bw = 12 / minWindows = 3


### PLOT!

# SweeD
sweed <- read.table("../SweeD-30092016/SweeD_Report.CEU_5000_derived", skip=2, header=T)
if (!is.null(minPos)) sweed <- subset(sweed, Position > minPos & Position < maxPos)
sweedMax <- max(sweed$Likelihood, na.rm=T)
# Pybus
HB <- read.table("../data/whole_genome.CEU.Complete.boosting.scores", header=T)
HB <- subset(HB, chromosome==2)
HB$Position <- with(HB, (start+end)/2)
if (!is.null(minPos)) HB <- subset(HB, Position > minPos & Position < maxPos)
HB$score[HB$score<0] <- 0
HBMax <- max(HB$score, na.rm=T)
# GFF
gff <- read.table("../data/Homo_sapiens.GRCh37.75_chr2.gtf", comment.char="#", sep="\t", stringsAsFactors=F)
gff <- subset(gff, V3=="gene")
gene_name <- sapply(gff$V9, function(x) {
    z <- unlist(strsplit(x, ";"))
    z <- z[grepl("gene_name",z)]
    z <- gsub(" ", "", z)
    z <- gsub("gene_name", "", z)
    return(z)
})
gff$gene_name <- gene_name


# PLOT
g <- plot_1bubbleRing(mcswan.DF = OB$estimation,
                      focal.deme = "i1",
                      log10_mcswan = F,

                      pybus.pos = HB$Position,
                      pybus.score = HB$score,
                      pybus.signif = 0.401999708647892,

                      sweed.pos = sweed$Position,
                      sweed.score = sweed$Likelihood,

                      mcswan.width = 40,
                      yfactor = 100,
                      min.y = -200,
                      max.y = 140,
                      plot.ticks = T,

                      gff.start.pos.CDS = gff[,4],
                      gff.end.pos.CDS = gff[,5],
                      do_radial = T)

ggsave(g, filename="0810_Homo_UNFOLDED_VCFpolarizedOK_bubbles.pdf")


if (F) {

    fact <- 300
    g=ggplot(data=OB$estimation) +
        geom_rect(aes(xmin=136545410, xmax=136594750, ymin=0, ymax=300), col="dodgerblue2", alpha=.5, size=1) +
        geom_line(data=sweed, aes(x=Position, y=Likelihood/sweedMax*fact), col="#4794D0", size=.5) +
        geom_point(aes(x=sweep.lbound, y=BF, col=deme, group=deme), pch="+") +
        theme_minimal()  +
        ylim(c(-500,300)) +
        coord_polar() +
        #geom_segment(aes(x=sweep.lbound, xend=sweep.lbound, yend=0, y=BF, col=deme, group=deme)) +
        geom_line(data=HB, aes(x=Position, y=-1*score/HBMax*fact), size=.5, col="#1C3E5D") +
        geom_hline(yintercept=-.40/HBMax*fact, col="green", size=.3)


    fact <- 300
    g=ggplot(data=OB$estimation) +
        geom_rect(aes(xmin=136545410, xmax=136594750, ymin=0, ymax=300), col="dodgerblue2", alpha=.5, size=1) +
        geom_line(data=sweed, aes(x=Position, y=Likelihood/sweedMax*fact), col="#4794D0") +
        geom_point(aes(x=sweep.lbound, y=BF, col=deme, group=deme), pch="+") +
        theme_minimal()  +
        ylim(c(-500,300)) +
        coord_polar() +
        geom_rect(aes(xmin=sweep.lbound, xmax=sweep.rbound, ymin=0, ymax=BF, fill=deme, group=deme), size=.1) +
        geom_line(data=HB, aes(x=Position, y=-1*score/HBMax*fact), size=.5, col="#1C3E5D") +
        geom_hline(yintercept=-.40/HBMax*fact, col="green", size=.3)

    ggsave(g, filename="2909_Homo_UNFOLDED2.png")

}


OR <- OB$estimation

# time densities
plot(density(subset(ob, deme=="i1")$sweepAge), col="red", xlim=c(0,2000))
lines(density(subset(ob, deme=="i2")$sweepAge), col="blue")

tg <- 25
g=ggplot(data=OB$estimation) +
    geom_violin(aes(factor(deme), y=sweepAge*tg, group=deme, fill=deme), col="gray", trim=F, adjust=.5, alpha=.7) +
    #scale="count"
    geom_jitter(aes(factor(deme), sweepAge*tg, size=BF), col="black", alpha=.8, pch=21, width=.3) +
    scale_fill_brewer(palette = "Greens") +
    #, position="jitter"
    theme_minimal() +
    xlab("Population") + ylab("Sweep age (yrsBP)") + ylim(0, max(OB$estimation$sweepAge)*tg)
#+coord_flip()
ggsave(g, filename="0810_Homo_UNFOLDED_sweepAgeViolin.pdf")


######################################
######################################
######################################
# VALIDATION

pd <- generate_pseudoobs(reftb, nSimul=70, sweepingIsl=NULL, L=1e6, recRate = list("rlogunif", 1e-8, 2e-7), sweepAge = NULL, sweepPos = .7, verbose=T, doSFS=T, default_sweepAge_prior_func="runif")

pd <- generate_pseudoobs(reftb, nSimul=1, sweepingIsl=NULL, L=1e6, recRate = list("rlogunif", mut/10, mut*10), sweepAge = NULL, sweepPos = .5, verbose=T, doSFS=T, default_sweepAge_prior_func="runif")
# should increase heap space around 5 Go

pd <- generate_pseudoobs(reftb, nSimul=100, sweepingIsl=NULL, L=1e6, recRate = list("rlogunif", mut/10, mut*10), sweepAge = NULL, sweepPos = .5, verbose=T, doSFS=T, default_sweepAge_prior_func="runif", Smu=NULL)
#save(pd, file="1310_100POD_UNFOLDED_L(1e6)_p(05)_r(.1MUT-10MUT)_T(runif)_msCorrig_new_SmuNULL_NoTimes1")

Ts <- 1800
pd <- generate_pseudoobs(reftb, nSimul=100, sweepingIsl=NULL, L=1e6, recRate = list("rlogunif", 2e-7, 2e-7), sweepAge = list(list("runif", Ts, Ts), list("runif", Ts, Ts)), sweepPos = 0, verbose=T, doSFS=T, default_sweepAge_prior_func="runif", Smu=NULL, save_each_file=F)
#save(pd, file="0610_100POD_UNFOLDED_L(1e7)_p(05)_r(.4MUT)_T(1800)_msCorrig_SmuNULL_ANCIENT")

#save(pd, file="2609_POD_UNFOLDED_L(1e6)_p(05)_r(.1MUT_10MUT)_T(runif)_msCorrig")
#load("0909_POD_L(9e5)_p(07)_r(5e-9_2e-7)")
load("1009_POD_L(1e6)_p(07)_r(1e-8_2e-7)_T(runif)")
load("2009_POD_L(1e6)_p(05)_r(.1MUT_10MUT)_T(runif)")
load("2009_POD_L(2e6)_p(05)_r(.1MUT_10MUT)_T(runif)")


load("../WDIR/2609_POD_UNFOLDED_L(1e6)_p(05)_r(.1MUT_10MUT)_T(runif)_msCorrig")


load("../WDIR/2609_UNFOLDED_demo_10000sim_coalesceNoMigr_L(20000)_mscorrig_RELATIVEsfs_minSD(0)_tolLDA(1e-9).RData")

load("../WDIR/2609_POD_UNFOLDED_L(1e6)_p(05)_r(.1MUT_10MUT)_T(runif)_msCorrig")

#load("0610_50POD_UNFOLDED_L(1e6)_p(05)_r(.4MUT)_T(runif)_msCorrig")

print(reftb$GENERAL$windowSize)

reftb$GENERAL$windowSize <- 2e4

load("0610_50POD_UNFOLDED_L(5e5)_p(05)_r(.1MUT-10MUT)_T(runif)_msCorrig_new_SmuNULL_NoTimes3")

Ts <- 1800
pd <- generate_pseudoobs(reftb, nSimul=5, sweepingIsl=NULL, L=1e6, recRate = list("rlogunif", 2e-7, 2e-7), sweepAge = list(list("runif", Ts, Ts), list("runif", Ts, Ts)), sweepPos = .5, verbose=T, doSFS=T, default_sweepAge_prior_func="runif", Smu=NULL, save_each_file=F)


vj <- sliding_validation(pd, reftb, method="contig", stat = "mean", minSNP = 20, cutoff = 10, bwadjust = 12, minWindows = 3, tolABC = c(.01, .005), tolGFIT = 0, plot_thinning = F, summary_stat = "mean", windowSlide = NULL)

svj <- summary_cv(vj, file=paste(tempDir,"/1310TEST_new2_4040_nullSmu_VALIDATION_UNFOLD_relativeSFS_L2e4_modeNew_ABC005_bw15_gfit0.001_hiNO.pdf",sep=""))


with(svj$param.estimation, plot(log10(est.sweepEnd-est.sweepStart), abs(est.sweepAge-true.sweepAge)))
with(svj$param.estimation, plot( true.sweepAge, log10(est.sweepEnd-est.sweepStart)))

#save(svj, file="0110_VALIDATION_UNFOLD_relativeSFS")


###

load("0110_VALIDATION_UNFOLD_relativeSFS")
svj$performance.rates
svj$NRMSE
svj$number.of.sweeps
svj$sweep.inclusiveness







