#' A function to count GC content in the selected bins.
#'
#' @param gcc a vector containing information about GC content. 
#' @param chrn a vector of chromosome names or NULL for all chromosomes.
#' @param bin a size of a bin.
#' @param mc.cores a number of cores for parallel computing.
#' @keywords GC content
#' @export
#' @examples
#' count.gc.cont()


count.gc.cont <-  function(gcc=NULL,chrn=NULL,bin=300,mc.cores=2){
			if(is.null(gcc)){stop("GC content was not provided")
			}else{
			  gc.count <- parallel::mclapply(1:length(gcc),function(ii){
						size <- bin/100
						ll   <- gcc[[ii]]
						ll   <- c(ll,rep(0,bin-length(ll)%%(bin)))
						n    <- length(ll)
						sl   <- seq(50,n,by=100)
						ll2  <- ll[sl]
						gc   <- sapply(1:ceiling(length(ll)/(size*100)),function(i) {
								     mean(ll2[((i-1)*size+1):(size*i)])
							})
						gc		
				      },mc.cores=mc.cores)
			}
			names(gc.count) <- names(gcc)
			gc.count <- gc.count[chrn]
			return(gc.count)
		}
