#' A function to compute number of tags in the selected bins.
#'
#' @param profs a list of signle positions representing used reads. 
#' @param chr.stat a data.frame with sizes of analyzed chromosmes.
#' @param chrn a vector of chromosome names or NULL for all chromosomes.
#' @param mc.cores a number of cores for parallel computing.
#' @param bin a size of a bin.
#' @param scale.genome.to.100Mb logical. Should the values be normalized to genome size?
#' @param normalize.perBinSize logical. Should the values be normalized to been size?
#' @keywords counts bins
#' @export
#' @examples
#' count.tags.in.bins()



count.tags.in.bins <- function(profs=NULL,chr.stat=NULL,chrn=NULL,mc.cores=2,bin=NULL, scale.genome.to.100Mb=TRUE,normalize.perBinSize=TRUE){
  if(is.null(profs)){stop("Profiles were not provided")}
  if(is.null(chr.stat)){stop("Lengths of chromosoms were not provided")}	      
  if(is.null(chrn)){stop("Chromosome names were not provided")}
  if(any(names(chrn)!=chrn)) {names(chrn) <- chrn}
  ll <- lapply(profs,function(fcen){#fcen <- profs[[1]] 
	    f <- list(); 
	    sets <- names(fcen); names(sets) <- sets
	    f$counts <- list(); 
	    for (s in sets){
	      f$counts[[s]] <- list();
	    }
	    for (chr in chrn){# chr="chr2L"
	      poss <- lapply(sets,function(s){ceiling(fcen[[s]][[chr]]/bin)}) 
	      rpos <- c(1,ceiling(chr.stat$length[which(chr.stat$chr==chr)]/bin))
	      f$pos[[chr]] <- (rpos[1]:rpos[2])*bin-bin/2
	      for (s in sets){
		f$counts[[s]][[chr]] <- rep(0,times=(diff(rpos)+1));
	      }
	      frnw  <- parallel::mclapply(sets,function(s){descr::freq(poss[[s]],plot=F) },mc.cores=mc.cores)
	      tposs <- lapply(sets,function(s){frnw[[s]][,1][!(names(frnw[[s]][,1])%in%c("Total","NA's"))] }) 
	      rm(poss,frnw);gc()
	      for (s in sets){#sets[1]->s
		pos <- as.numeric(names(tposs[[s]]))
		ind <- which(!is.na(match((rpos[1]:rpos[2]), pos)))
		f$counts[[s]][[chr]][ind] <- tposs[[s]]
	      }
	      cat(chr," ")
	    };
	    f$values <- parallel::mclapply(f$counts, function(v) {#f$counts[[1]]-> v
						lib.frac <- 1e+06/sum(unlist(v))
						genome.frac <- sum(chr.stat[which(chr.stat$chr %in%   chrn), 2]) * (1/ifelse(scale.genome.to.100Mb,1e+08,1))
						frac <- lib.frac * genome.frac
						lapply(v,function(vv) vv * frac * ifelse(normalize.perBinSize,200/bin,1))
			}, mc.cores = mc.cores)

	    cat("\n")
	    return(f)
 })
 names(ll) <- names(setlist)
 return(ll)
}
