% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/observed_manipulations.R
\name{convert_VCF}
\alias{convert_VCF}
\title{Filters and converts a VCF into a file of per-population per-SNP allele counts}
\usage{
convert_VCF(vcfPath, pops, reftb, outPath, chromosome, haploidize = FALSE,
  minQUAL = NULL)
}
\arguments{
\item{vcfPath}{path to the VCF input file}

\item{pops}{a dataframe of sample names and population indices \emph{or}, alternatively, an external file with the same format (see \emph{Details})}

\item{reftb}{an initialized \code{referenceTable} object}

\item{outPath}{path to the file in which will be written the per-population per-SNP allele counts}

\item{chromosome}{name of the chromosome to analyze (should match a chromosome name in the VCF file)}

\item{haploidize}{set this option to TRUE to randomly draw a single allele out of the diploid genotypes (useful to mitigate the impact of low-coverage sequencing), note that if \code{haploidize=TRUE} the number of chromosomes in your MS command must be adjusted!}

\item{minQUAL}{any SNP with \code{QUAL < minQUAL} will be filtered out (set \code{minQUAL = NULL} to disable the filtering)}
}
\value{
No internal return. The function will write an external file. This file has a header summarizing the conversion parameters (starting with "#"), followed by the SNP-wise allele counts. For each line, there are three space-delimited fields:\itemize{
\item the SNP position on the chromosome; \item a SNP coefficient: always 1 if derived allele counts; can be 1 or 0.5 if minor allele counts -- 0.5 corresponding to a SNP which can have two possible alternative states; \item the per-population allele counts, formatted \code{ac.i1.i2.[...].in} with (\emph{i1}, \emph{i2}, ..., \emph{in}) the allele counts in populations (1, 2, ..., \emph{n}). }
}
\description{
This function calls a Python script which reads in a VCF file, extracts chromosome-specific variants, filters the variants to retain only high-quality biallelic SNPs, calculates the population-specific allele counts for each retained SNP (minor or derived allele counts, depending on your \code{referenceTable} specifications) and writes the results into a new external file.
}
\details{
Non biallic SNPs in the VCF file will be automatically filtered out. If there are more samples in the VCF file than samples specified in \emph{pops}, those supernumary samples will be ignored.

The \code{pops} dataframe must have two columns: (i) the names of the samples, matching the names in the VCF file; (ii) the indices (integers) of the population to which each sample belongs, these indices must correspond to the indices of the populations specified in the \emph{MS}-formatted demographic history. For instance, if the sample B101 in the VCF belongs to the third-indexed population of your demographic history, the first line will read: first column:B101, second column: 3. Note that \emph{MS} starts population indexing at 1 (included). If \code{pops} is an external file, it should be formatted similarly, with tab-separated columns, no header, no rownames and without quote marks around the sample names.
}
\examples{
Please refer to the vignette.
}
\seealso{
\code{\link{get_SFS}}
}

