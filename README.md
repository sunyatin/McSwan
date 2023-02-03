# McSwan

## Synopsis

McSwan is a R package to infer the extent and the timing of positive selection in the genomes of single or multiple populations.  

## Installation

```{r}
install.packages("devtools")
library(devtools)
install_github("sunyatin/McSwan/McSwan")
```

## Aims

The **M**ultiple-collision **C**oalescent **SW**eep **AN**alyser (McSwan) is a genome scan approach which aims at simultaneously detecting and dating hard selective sweeps from genome sequence data of multiple individuals from single or multiple populations.


## Model inputs and assumptions

McSwan comes as a cross-platform R package with Python modules (tested on Windows 8/10, OS X & Ubuntu). The typical inputs of a McSwan analysis are (i) a VCF file containing polarized or unpolarized individual genotypes; (ii) a neutral demographic model given as an *ms*-formatted string comprising any possible demographic switch natively available in *ms* (Hudson 2002). 

In its current implementation, McSwan discards SNPs having at least one missing genotype and applies to well-covered re-sequencing datasets. By default, McSwan will discard variants which are not biallelic SNPs.

## Initializing the working session

The working session must be initialized by specifying a temporary directory and the paths to external softwares.

```{r message=F, echo=F, eval=F}
setwd("C:/Users/Windows/Desktop/McSwan_public/")
```
```{r eval=F}
library(McSwan)

set_session(tempDir = "temp", # (MANDATORY) path to a temporary directory (will be created if it does not exist; if the directory already exists, its content will be overridden)
            pythonPath = "python", # (MANDATORY) path to your Python executable or, if Python is callable as an environment variable, simply set "python"
            javaPath = "C:/Program Files (x86)/Java/jre1.8.0_101/bin/javaw.exe" # (OPTIONAL) path to your Java executable (e.g. javaw.exe in Windows) or, if Java is installed as an environment variable, simply set "java"
            )
```

Note that `javaPath` argument is only mandatory if you set out to simulate pseudo-observed datasets for validation.

If you encounter bugs while setting the python path, please refer to the package help `?set_session`.

## Construction of the reference table

This step allows the construction of a *reference table* which contains a finite set of model realizations summarized into joint site frequency spectra. These simulations will allow learning the expected signatures of positive selection conditionally to your demographic scenario.

#### Specification of the neutral demographic model

For the tutorial, we will consider the joint demographic model of the Luhyan and North-Western European populations published by Excoffier *et al.* (2013). We will assume 20 diploid individuals (= 40 chromosomes) sampled in each population.

Coalescent simulations in McSwan are performed using Hudson's *ms* software. The demographic history must therefore be formatted as a string of switches according to *MS* conventions. If you don't know the *ms* syntax, please refer to the [manual](http://home.uchicago.edu/rhudson1/source/mksamples/msdir/msdoc.pdf). As for most coalescent simulators, demographic events are defined backward in time (time indices increase towards the past).

Note that an appropriate *ms* command can also be generated using the R package `coala` and visualized a posteriori using the Java software [PopPlanner](http://www.mabs.at/ewing/msms/popplanner.shtml).


```{r}
mut <- 2.5e-8 # mutation rate per base per generation
No <- 4106 # diploid effective size for your reference population
ms <- paste("80", # total number of sampled chromosomes
            "1", # number of simulation replicates -- MUST BE "1"
            "-t",4*No*mut, # population mutation rate $\theta$ = 4*No*mu
            "-I 2 40 40", # population structure
            "-en 0 2",58307/No, # population size rescaling
            "-en",70/(4*No),2,23738/No,
            "-ej",2000/(4*No),"1 2",
            "-en",2683/(4*No),2,7083/No,
            "-en",2783/(4*No),2,7963/No)
```

We can now check the *ms* command string:

```{r echo=F}
print(ms)
```

This command says that we have sampled 80 chromosomes at present (= 40 diploid individuals), 40 of which belong to a first population (= Europeans) and the remaining 40 to a second population (= Luhyans). Both populations exchange no gene flow. The reference effective size $N_o$ is the effective size of the European population at present. The fundamental population mutation rate is $\theta$ and equals $4N_o\mu$ with $\mu$ the mutation rate per base per generation. The effective size of the Luhyan population at present is 58,307. At time 70 generations before present, the Luhyan population is resized to 23,728 individuals. At time 2,000 generations, the European populations merge with the Luhyan population (i.e. forward in time, it corresponds to a split or divergence). At 2,683 generations BP, the ancestral population experiences a bottleneck ($N_e = 7,083$) until 2,783 generations BP where the $N_e$ increases to 7,963. Time is scaled in units of $4N_o$ and population size are rescaled in reference to the invariant $N_o$.

#### Parameter prior distributions

A **referenceTable** object is initialized using the `generate_priors()` function which takes the following (mandatory) arguments:

* `msDemography` The *ms*-formatted demographic history as defined previously.
* `No` The effective population size of reference ($N_o$).
* `fold` A switch specifying the nature of the joint spectra: either folded (`fold = TRUE`) if you want the spectra to be based on **minor** allele frequencies; or unfolded (`fold = FALSE`) if on **derived** allele frequencies. Please be very careful with this switch as it can significantly bias the results if misspecified. Use only `fold=TRUE` if your VCF is defined in terms of ANCESTRAL/DERIVED alleles instead of REF/ALT.
* `windowSize` The length of the sequences (in bp) to simulate, the larger the better (> 100k).
* `nSimul` The number of simulations to perform **for each** model, the larger the better (> 10k).

The default behaviour when you have $\geq 2$ populations is to set the sweep ages automatically relatively to the divergence times. If, however, you have a single population or if you want to manually force period boundaries, use the `sweepAge` argument. Please refer to the `?generate_priors` help for how to format the argument value.

```{r eval=F}
reftb <- generate_priors(msDemography = ms, 
                         fold = FALSE,
                         No = No,
                         windowSize = 20000,
                         nSimul = 2000)
```

Parameter prior values will be stored, for each model, in the `PRIORS` slot of the `reftb` object.

In our example, we will have three models: `i0`, the neutral model; `i1`, the selective model of European-specific sweeps; `i2`, the selective model of Luhyan specific sweeps. McSwan names the models as `i` & the index of the sweeping population. **Note that `i0` will always refer to the neutral model. Note also that *MS* starts population indexing at 1 (included).**

#### Generate joint site frequency spectra across the models

Based on the initialized prior parameter values, joint spectra can be simulated across all models (in the present case: `i0`, `i1` and `i2`) using the function `coalesce()`. Spectra will be stored into the `SFS` slot of the `reftb` object.

```{r eval=F}
reftb <- coalesce(reftb)
```

#### Machine learning

This step allows (i) the fitting of a linear discriminant model to classify neutral and selective demographic models, based on the expected SFS signatures; (ii) the fitting of partial least square regression models to infer the age of the sweep conditionally to the SFS simulated previously. These supervised learning algorithms will be used to infer the extent and dynamics of the positive selection along the genome.

```{r eval=F}
reftb <- dim_reduction(reftb, PLS_plotCV = TRUE)
```
  
> **IMPORTANT.**
Given the computation load required for the construction of the reference table (simulation of joint spectra and dimension reduction), it is highly recommended to save the reference table, e.g.:
`save(reftb, file = paste0(tempDir,"/reftb"))`.
> The reference table can be re-imported afterwards by executing `load(paste0(tempDir,"/reftb"))`.

> **IMPORTANT.**
You can have access to a summary of your referenceTable object at any time, using `print(X)` or simply by typing `X` with $X$ the name of your referenceTable object.

## Assessment of the statistical performance

To assess the accuracy of McSwan to (i) detect selective sweeps and (ii) estimate sweep ages, conditionally to your demographic history, we provide functions to perform a validation assessment using pseudo-observed datasets (PODs) simulated with the [MSMS](http://www.mabs.at/ewing/msms/) software.

#### Simulation of the pseudo-observed datasets (PODs)

We will simulate pseudo-observed chromosomes with varying sweep ages and recombination rates. The next function calls the *MSMS* simulator which is written in Java. To run *MSMS*, you must have set the path to the Java executable in `set_session(javaPath=...)`.

```{r eval=F}
POD <- generate_pseudoobs(reftb, 
                         nSimul = 5, # number of PODs to simulate for each model
                         L = 1e6, # length of the chromosomes to simulate
                         recRate = list("rlogunif", mut/10, mut*10), # distribution of recombination rates
                         sweepAge = NULL, # distribution of sweep ages -- automatic selection if set to NULL
                         sweepPos = .5, # relative position of the beneficial mutation in the simulated chromosome
                         verbose = T)
```

#### POD analysis

We can perform the genome scan over each POD using `gscan`:

```{r eval=F}
vj <- gscan(POD, # the pseudo-observed dataset generated with generate_pseudoobs()
             reftb, # the reference table
             firstPos = NULL, # starting position of the genome scan (NULL to automatically start at the first known position)
             lastPos = NULL, # last position of the genome scan (NULL to automatically end at the last known position)
             minSNP = 20, # minimum number of SNPs in a window to consider it a valid genomic region to work with
             windowSizes = seq(5e3, 2e5, length.out = 20), # a vector of window lengths over which genome scans will be performed iteratively (the more elements, the better), see ?gscan
             nSteps = 20, # number of overlapping shifts (the higher, the better)
             discard_extraRange = TRUE) # whether to discard sweep estimations going beyong the prior range (TRUE is recommended)
```

Infer the tiling paths of the sweeps. **Don't forget this step.** Note that you can build receiver operating characteristic (ROC) curves by iterating the `thin()` function over increasing `signif.threshold` values, extracting the TPR and TNR each time with a call to the `summary()` function.
```{r eval=F}
vjt <- thin(scanResult = vj, 
           reftb = reftb, 
           X = POD, # the pseudo-observed dataset used in the gscan() call
           max_L = 1e6, # the "horizon distance" (in base pairs), ie. the maximum distance from any position above which loci detected as under selection cannot be assumed to be related to the same sweep contig
           signif.threshold = .5) # (single or vector of numeric values) (between 0 and 1) a selection score cut-off above which we decide that a SNP is positive selected; if unique value, this cutoff will be used for all populations; alternatively, you can give a vector of cutoff values which will be specifically used for each population (e.g. signif.threshold=c(0.5, 0.2))
```

To calculate the performance statistics and output graphics to `file`:

```{r eval=F}
res <- summary(X = vjt, valtb = POD, file = paste0(tempDir,"/McSwan_perf_statistics.pdf"))
```

The function `summary()` outputs a list of statistics. You can have access to each of them by typing:

```{r eval=F}
res$performance.rates # the estimated performance rates (specificity = TNR, sensitivity = TPR, false discovery rates = FDR)
res$NRMSE # normalized/standardized root mean square error
res$number.of.sweeps # number of sweep contigs detected in each POD (ideally, should be 0 for the neutral model and 1 for the selective models)
res$sweep.inclusiveness # number of times the detected sweeps encompass the beneficial mutation
```

## Inference on observed datasets

#### Import the genotypes stored in a VCF file

* *OPTION 1*: for package testing only

For a basic test, you can load a ready-to-use observed dataset (named `obs`):
```{r eval=F}
load(system.file("data", "UNFOLDED_POLARIZED_2dSFS_CEU-LWK_20-20_chr2", package = "McSwan"))
```

This dataset corresponds to the derived allele counts for the SNPs along chromosome 2 genotyped for 20 & 20 diploid individuals within the Luhyan and North-Western European populations respectively.

* *OPTION 2*: import observed data

This procedure enables importating a VCF file for a specific chromosome/scaffold/contig.

First, convert the VCF to a McSwan-readable file containing per-SNP allele counts:

```{r eval=F}
convert_VCF(vcfPath = "../data/set the ancestral allele/CEU-20_LWK-20_POLARIZED.vcf", # path to the VCF input file
            pops = "../data/pops_CEU-20_LWK-20.txt", # a dataframe or alternatively a file giving the sample names and their population index, see below
            reftb, 
            outPath = "../temp/UNFOLD_POLARIZED.txt", # the file in which results will be exported
            chromosome = "2", # name of the chromosome (or contig or scaffold) to analyze
            minQUAL = 30) # minimum QUAL value of the SNPs
```

The `pops` dataframe (or file) must have two columns: (i) the names of the samples, matching the names in the VCF file; (ii) indices (integers) of the population to which each sample belongs, these indices must correspond to the indices of the populations specified in the *ms*-formatted demographic history. For instance, if the sample B101 in the VCF belongs to the third-indexed population of your demographic history, the first line will read: `B101    3`. Note that *ms* starts population indexing at 1 (included). If `pops` is an external file, it should be formatted similarly, with tab-separated columns, no header, no rownames and without quote marks around the sample names.

When the conversion is completed (can take some time), import the generated file:

```{r eval=F}
obs <- get_SFS("../temp/UNFOLD_POLARIZED.txt", reftb)
```

#### Get a summary of the SNP density

To have an idea about the distribution of the SNP density given a window size (useful to set the minimum SNP threshold or to choose an appropriate window size range, in bp), simply execute:

```{r eval=F}
summary(obs, windowSize = 10000)
```

#### Observed dataset analysis

We can now execute the genome scan with the `gscan()` function (in the example below, we focus on the 130-140 Mbp region of the chromosome 2 which contains a positively selected mutation near the LCT gene, allowing lactase persistence in the Europeans).

```{r eval=F}
OA <- gscan(obs, # the observed dataset generated with get_SFS()
            reftb, # the reference table
            firstPos = 130e6, # starting position of the genome scan (NULL to automatically start at the first known position)
             lastPos = 140e6, # last position of the genome scan (NULL to automatically end at the last known position)
             minSNP = 20, # minimum number of SNPs in a window to consider it a valid genomic region to work with
             windowSizes = seq(5e3, 2e5, length.out = 20), # a vector of window lengths over which genome scans will be performed iteratively (the more elements, the better), see ?gscan
             nSteps = 20, # number of overlapping shifts (the higher, the better)
            discard_extraRange = TRUE) # whether to discard sweep estimations going beyong the prior range (TRUE is recommended)
```

Next, we infer the tiling path of the selective sweeps with the `thin()` function. **This is a mandatory step!**

```{r eval=F}
OB <- thin(OA, 
           reftb, 
           X = obs, # the pseudo-observed dataset used in the gscan() call
           max_L = 1e6, # the "horizon distance" (in base pairs), ie. the maximum distance from any position above which loci detected as under selection cannot be assumed to be related to the same sweep contig
           signif.threshold = .5) # (single or vector of numeric values) (between 0 and 1) a selection score cut-off above which we decide that a SNP is positive selected; if unique value, this cutoff will be used for all populations; alternatively, you can give a vector of cutoff values which will be specifically used for each population (e.g. signif.threshold=c(0.5, 0.2))
```

The dataframe contains all the detected selective sweeps with one sweep on each line. In columns, you have:

* `sweep.center` the center of the sweep region
* `sweep.lbound` the first position of the sweep region
* `sweep.rbound` the last position of the sweep region
* `RPP` ratio of the LDA-derived posterior probabilities between the selective vs. neutral models, ie the highest ratio among the SNPs within the sweep region, should be >1
* `deme` the population detected to have experienced a selective sweep (given as `i` & index of the population, here `i1` = Europeans and `i2` = Luhyans)
* `sweepAge` the point sweep age estimate
* `sweepAge.IC.low` lower boundary of the 95% confidence interval for the sweep age estimation
* `sweepAge.IC.up` upper boundary of the 95% confidence interval for the sweep age estimation

We can plot the temporal density of population-specific sweeps:

```{r eval=F}
require(ggplot2)
ggplot(data=OB, aes(x=sweepAge, group=deme, fill=deme)) + stat_density(trim = FALSE)
```

or the distribution of sweeps along the genome:

```{r eval=F}
ggplot(data=OB, aes(x=sweep.center, y=BF, col=deme)) + geom_point() + geom_segment(aes(xend=sweep.center, yend=0)) + xlim(c(0, 250e6))
```

## Global functions

At any time, you can access a summary of the `referenceTable`, `observedDataset` and `validationTable` using `print(X)` with *X* the name of your objects.


If you want to parallelize the construction of the `referenceTable`, you can bind multiple `referenceTable` objects a posteriori (always two by two) using the `combine()` function.

## Troubleshoot

If you encounter any bug or problem, please contact remi (dot) tournebize (at) gmail (dot) com

## What is the optimal `signif.threshold` (score cutoff) to use for the `thin()` function?

To estimate the optimal `signif.threshold` for the CEU population (model `i1`) that maximizes TPR while minimizing TNR, you may use the following script:

```{r eval=F}
    test <- seq(.01, 1, by = .02)
    RES <- c()
    for (i in seq_along(test)) {
        ob <- thin(OA, reftb, X = obs, max_L = 1e6, signif.threshold = c(test[i], 0.5))
        res <- summary(X = ob, valtb = POD, file = paste0(tempDir,"/McSwan_perf_statistics.pdf"))
        RES <- rbind(RES, data.frame(signif.thresh = test[i], res$performance.rates$binary[1, c("TPR.sensitivitiy", "TNR.specificity")]))
    }
    
    ggplot(RES) + 
        geom_line(aes(x=signif.thresh, y=TPR.sensitivity), col="red") +
        geom_line(aes(x=signif.thresh, y=TNR.specificity), col="blue")
```

## Citation

Tournebize, R., Poncet, V., Jakobsson, M., Vigouroux, Y., & Manel, S. (2019). McSwan: A joint site frequency spectrum method to detect and date selective sweeps across multiple population genomes. Molecular ecology resources, 19(1), 283-295.

