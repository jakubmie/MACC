#' A function to read positions of aligned tags.
#'
#' @param setlist a list of vectors with names of files contains positions of aligned reads. Each vector should correspond to one sample, and order of files should correspond to ascending order of MNase concentrations used in the experiment.
#' @param path2files a path to files to be read.
#' @param chrn a vector of chromosome names or NULL for all chromosomes.
#' @param mc.cores a number of cores for parallel computing.
#' @param filter.anomalies logical. Should tags associated with anomalously high read counts be removed?
#' @keywords tags
#' @export
#' @examples
#' read.tags()

read.tags <- function(setlist=NULL,path2files=NULL,chrn=NULL,mc.cores=2,filter.anomalies=TRUE){
	      if(is.null(setlist)){stop("File names were not provided")}
	      if(is.null(path2files)){stop("A path to files was not provided")}
	      if(is.null(chrn)){stop("Chromosome names were not provided")}
              if(any(names(chrn)!=chrn)) {names(chrn) <- chrn}
              tt <- list()
    for (nsets in names(setlist)) {#nsets <- names(setlist)[1]
        sets <- setlist[[nsets]]
        names(sets) <- sets
        for (set in sets) {#sets[2]->set
            fname <- paste(path2files, set, sep = "")
            cat("Reading ", fname, "\t")
            tags <- read.bam.file(fname = fname)
            if (any(!chrn %in% names(tags))) {
                tmp <- paste(chrn[!chrn %in% names(tags)], collapse = " ")
                stop(paste(tmp, " not in ", set, " chromosome names", 
                  sep = ""))
            }
            tags <- tags[chrn]
	    if(filter.anomalies){
		    tags <- parallel::mclapply(tags,function(c){
			      inp <- MACC:::rta(tv=c[which(c>0)])
			      inn <- MACC:::rta(tv=c[which(c<0)])
			      pp  <- which(!inp)
			      np  <- which(!inn)
			      p   <- unique(c(c((pp*2-1),pp*2),c((np*2-1),np*2)))
			      return(c[which(!((1:length(c)) %in% p))])
		    },mc.cores=mc.cores)
	    }
	    tt[[nsets]][[set]] <- tags
            cat("done\n")
        }
    }
    return(tt)
}

