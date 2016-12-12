

#' @title convert a VCF to a proper format
#' @param pops a dataframe containing two columns: \emph{sample.ID} & \emph{population.ID}. Note that \emph{population.ID} needs to be an integer and that \emph{population.ID} must correspond to the index of the island specified in your \emph{ms} command
#' @export
convert_VCF <- function(vcfPath, pops, reftb, outPath, chromosome, minQUAL = NULL) {
  
  # export populations
  if (is.data.frame(pops)) {
    write.table(pops, paste(tempDir,"/pops.txt",sep=""), row.names=F, col.names=F, sep="\t", quote=F)
  } else {
    if (!file.exists(pops)) stop("the pops file you specified has not been found")
  }
  
  # converter command
  cmd <- paste("python ",tempDir,"/vcf2pac.py -vcf \"",vcfPath,"\" -o \"",outPath,"\" -chr ",chromosome," -p ",tempDir,"/pops.txt",sep="")
  
  # add folding?
  doFold <- reftb$GENERAL$folded
  if (doFold) cmd <- paste(cmd,"-fold")
  
  # min QUAL?
  if (!is.null(minQUAL)) cmd <- paste(cmd,"-minQ",minQUAL)
  
  system(cmd)
}

#' @title import the observed dataset
#' @export
get_SFS <- function(inputFile, ms) {
  
  # get template
  write.table(ms, paste(tempDir,"/template.txt",sep=""), row.names=F, col.names=F, quote=F)
  template_cmd <- paste("python ",pyPath," -i ",tempDir,"/template.txt -o ",tempDir,"/template_sfs.txt",sep="")
  system(template_cmd)
  template <- unlist(read.table(paste(tempDir,"/template_sfs.txt",sep=""), sep="\t", header=T))
  
  # get data
  G <- read.table(inputFile, sep=" ", header=F, colClasses = c("integer","factor","factor"))
  
  # is folded
  fold <- !grepl("unfolded", unlist(read.table(inputFile, header=F, skip = 3, nrows = 1, comment.char = "", sep = "~")))
  
  # combine
  out <- list(template=template, obsData=data.frame("POS"=G[,1], "nSNP"=G[,2], "PAC"=G[,3]), folded=fold)
  
  rm("G"); invisible(gc(F))
  
  class(out) <- "observedDataset"
  return(out)
}
  
  
  