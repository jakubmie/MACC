#' A function to read genomic data.
#'
#' @param genome a string corresponding to genome of interest.
#' @param chrn a vector of chromosome names or NULL for all chromosomes.
#' @param bin a size of a bin.
#' @param mc.cores a number of cores for parallel computing.
#' @keywords genome
#' @export
#' @examples
#' read.genomic.data()

read.genomic.data <-  function(genome="",chrn=NULL,bin=300,mc.cores=2){
			avg <- c("dm3","mm9","hg19")
			if(!genome%in%avg){stop("Only genomic data for ",paste(avg[-length(avg)],collapse=", ")," and ",avg[length(avg)]," are available\n")}
			if (genome=="dm3"){
			      data(dm3.genomic)
			      if(is.null(chrn)){ 
					  chrn <- dm3.genomic.data$dchrn
			      }else{
				  if(any(!chrn%in%dm3.genomic.data$dchrn)){
					  stop("Not all of the selected chromosomes are present in the selected genome\n")
				  }else{
					  names(chrn) <- chrn
				  }
			      }
                             CpG <- dm3.genomic.data$dCpG
                             gc.cont.all <- dm3.genomic.data$gc.cont.all[chrn]
			     chr.stat <- dm3.genomic.data$chr.stat
			}
			if (genome=="mm9"){
			      data(mm9.genomic)
			      if(is.null(chrn)){ 
					  chrn <- mm9.genomic.data$dchrn
			      }else{
				  if(any(!chrn%in%mm9.genomic.data$dchrn)){
					  stop("Not all of the selected chromosomes are present in the selected genome\n")
				  }else{
					  names(chrn) <- chrn
				  }
			      }			
			      CpG <- mm9.genomic.data$dCpG			 
			      gc.cont.all <- mm9.genomic.data$gc.cont.all[chrn]
                              chr.stat <- mm9.genomic.data$chr.stat 
			}			
			if (genome=="hg19"){
 			      data(hg19.genomic)
			      if(is.null(chrn)){ 
					  chrn <- hg19.genomic.data$dchrn
			      }else{
				  if(any(!chrn%in%hg19.genomic.data$dchrn)){
					  stop("Not all of the selected chromosomes are present in the selected genome\n")
				  }else{
					  names(chrn) <- chrn
				  }
			      }			
			      CpG <- hg19.genomic.data$dCpG
			      gc.cont.all <- hg19.genomic.data$gc.cont.all[chrn]
                              chr.stat <- hg19.genomic.data$chr.stat 
			}
		       gc.cont     <- count.gc.cont(gcc=gc.cont.all,bin=bin,mc.cores=mc.cores,chrn=chrn)
		       if(!is.null(CpG)){
			CpG   <- CpG[which(CpG[,1]%in%chrn),] 
			CpG <- GenomicRanges::GRanges(seqnames=as.vector(CpG[,1]),ranges=IRanges::IRanges(start=as.vector(CpG[,2]),end=as.vector(CpG[,3])),strand="*")			
		       }
		       if(any(!chrn%in%chr.stat[,1])){
			     tmp <- paste(chrn[!chrn%in%chr.stat[,1]],collapse=" ")
			     stop("Genomic data for ",tmp," are not available\n")
		       }else{ chr.stat <- chr.stat[chr.stat[,1]%in%chrn,]}  
		       rm(gc.cont.all)
		       cat(paste("\ndone for ",bin,"bp bins\n",sep=""))
		      return(list(chr.stat=chr.stat,chrn=chrn,gc.cont=gc.cont,CpG=CpG))
		      }

