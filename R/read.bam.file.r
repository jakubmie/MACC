#' A function to count GC content in the selected bins.
#'
#' @param fname a character showing a path to files to be read. 
#' @param mc.cores a number of cores for parallel computing.
#' @param tl an integer for tag shifting.
#' @keywords bam
#' @export
#' @examples
#' read.bam.file()

read.bam.file <- function(fname,mc.cores=2,tl=50) {  #require(Rsamtools) 
  bam <- Rsamtools::scanBam(fname,param=Rsamtools::ScanBamParam(what=c("flag","rname","pos","strand"),flag=Rsamtools::scanBamFlag(isUnmappedQuery=FALSE)))[[1]];
  strm <- as.integer(bam$strand=="+")
  tmp <- tolower(unique(as.vector(bam$rname)))
  if(all(sapply(tmp,function(a) (substr(a,1,3)!="chr")))) bam$rname <- paste("chr",as.vector(bam$rname),sep="")
  tags <- tapply(1:length(bam$pos),bam$rname,function(ii) 
				as.numeric(na.omit(strm[ii]*bam$pos[ii]-(1-strm[ii])*(bam$pos[ii]+tl-1)))
				) 
  return(tags)
}
