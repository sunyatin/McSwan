

#' @title Filters and converts a VCF into a file of per-population per-SNP allele counts
#' @description This function calls a Python script which reads in a VCF file, extracts chromosome-specific variants, filters the variants to retain only high-quality biallelic SNPs, calculates the population-specific allele counts for each retained SNP (minor or derived allele counts, depending on your \code{referenceTable} specifications) and writes the results into a new external file.
#' @param vcfPath path to the VCF input file
#' @param pops a dataframe of sample names and population indices \emph{or}, alternatively, an external file with the same format (see \emph{Details})
#' @param reftb an initialized \code{referenceTable} object
#' @param outPath path to the file in which will be written the per-population per-SNP allele counts
#' @param chromosome name of the chromosome to analyze (should match a chromosome name in the VCF file)
#' @param minQUAL any SNP with \code{QUAL < minQUAL} will be filtered out (set \code{minQUAL = NULL} to disable the filtering)
#' @param haploidize set this option to TRUE to randomly draw a single allele out of the diploid genotypes (useful to mitigate the impact of low-coverage sequencing), note that if \code{haploidize=TRUE} the number of chromosomes in your MS command must be adjusted!
#' @details Non biallic SNPs in the VCF file will be automatically filtered out. If there are more samples in the VCF file than samples specified in \emph{pops}, those supernumary samples will be ignored.
#' 
#' The \code{pops} dataframe must have two columns: (i) the names of the samples, matching the names in the VCF file; (ii) the indices (integers) of the population to which each sample belongs, these indices must correspond to the indices of the populations specified in the \emph{MS}-formatted demographic history. For instance, if the sample B101 in the VCF belongs to the third-indexed population of your demographic history, the first line will read: first column:B101, second column: 3. Note that \emph{MS} starts population indexing at 1 (included). If \code{pops} is an external file, it should be formatted similarly, with tab-separated columns, no header, no rownames and without quote marks around the sample names.
#' @return No internal return. The function will write an external file. This file has a header summarizing the conversion parameters (starting with "#"), followed by the SNP-wise allele counts. For each line, there are three space-delimited fields:\itemize{
#'\item the SNP position on the chromosome; \item a SNP coefficient: always 1 if derived allele counts; can be 1 or 0.5 if minor allele counts -- 0.5 corresponding to a SNP which can have two possible alternative states; \item the per-population allele counts, formatted \code{ac.i1.i2.[...].in} with (\emph{i1}, \emph{i2}, ..., \emph{in}) the allele counts in populations (1, 2, ..., \emph{n}). }
#' @examples Please refer to the vignette.
#' @seealso \code{\link{get_SFS}}
#' @export
convert_VCF <- function(vcfPath, pops, reftb, outPath, chromosome, haploidize = FALSE, minQUAL = NULL) {

  ms <- reftb$GENERAL$msDemography
  nIslands <- ifelse(!grepl("-I", ms), 1, NA)
  if (is.na(nIslands)) {
    ms <- unlist(strsplit(ms, " "))
    nIslands <- ms[which(ms=="-I")+1]
  }
  
  # export populations
  if (is.data.frame(pops)) {
    write.table(pops, paste(tempDir,"/pops.txt",sep=""), row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
    # converter command
    cmd <- paste(pythonPath," ",pyVCF2PACPath," -vcf \"",vcfPath,"\" -o \"",outPath,"\" -chr ",chromosome," -I ",nIslands," -p ",tempDir,"/pops.txt",sep="")
  } else {
    if (!file.exists(pops)) stop("the pops file you specified has not been found")
    # converter command
    cmd <- paste(pythonPath," ",pyVCF2PACPath," -vcf \"",vcfPath,"\" -o \"",outPath,"\" -chr ",chromosome," -I ",nIslands," -p ",pops,sep="")
  }
  
  # add folding?
  doFold <- reftb$GENERAL$folded
  if (doFold) cmd <- paste(cmd,"-fold")
  
  # haploidize
  if (haploidize) cmd <- paste(cmd,"-haploidize")
  
  # min QUAL?
  if (!is.null(minQUAL)) cmd <- paste(cmd,"-minQ",minQUAL)
  
  cat("Converting...\n")
  system(cmd)
}

#' @title Import a file of per-population per-SNP allele counts
#' @param inputFile a file of per-population per-SNP allele counts, previously generated with \code{\link{convert_VCF}}
#' @param reftb an initialized \code{referenceTable} object
#' @details The function needs the \emph{MS}-formatted demographic model contained in the \code{referenceTable} object to generate the template of the joint allele frequency spectra.
#' @return Returns a list object of class \code{observedDataset}, containing three elements:\itemize{\item \code{template} the template of the joint allele frequency spectra \item \code{obsData} a dataframe of SNP positions and allele counts \item \code{folded} whether those are minor (\code{TRUE}) or  derived (\code{FALSE}) allele counts }
#' @examples Please refer to the vignette.
#' @seealso \code{\link{convert_VCF}}
#' @export
get_SFS <- function(inputFile, reftb) {
  
  ms <- reftb$GENERAL$msDemography
  
  # get template
  write.table(ms, paste(tempDir,"/template.txt",sep=""), row.names=F, col.names=F, quote=F)
  template_cmd <- paste(pythonPath," ",pyPath," -i ",tempDir,"/template.txt -o ",tempDir,"/template_sfs.txt",sep="")
  system(template_cmd)
  template <- unlist(read.table(paste(tempDir,"/template_sfs.txt",sep=""), sep="\t", header=T))
  
  # get data
  G <- read.table(inputFile, sep=" ", header=F, colClasses = c("integer","factor","factor"))
  
# 26.01.2017 remove within-sample monomorphic bins (for compatibility with the relativization as computed on the simulated refTable)
wMono <- names(template)[c(1,length(template))]
cat("Non-variant markers among the samples are discarded, i.e. those with allele counts: ",wMono,"\n")
G <- G[!G[,3]%in%wMono,]

  # is folded
  fold <- !grepl("unfolded", unlist(read.table(inputFile, header=F, skip = 3, nrows = 1, comment.char = "", sep = "~")))
  
  # combine
  out <- list(template=template, obsData=data.frame("POS"=G[,1], "nSNP"=G[,2], "PAC"=G[,3]), folded=fold)
  
  rm("G"); invisible(gc(F))
  
  class(out) <- "observedDataset"
  return(out)
}
  
  
  