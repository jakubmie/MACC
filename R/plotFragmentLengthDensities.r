#' A function to plot distributions of fragment lengths.
#'
#' @param lengths a list of tag lengths. This may be retuned by read.tags function.
#' @param col a vector of colors.
#' @param xlim a numeric vector of length 2 giving range for x-axis.
#' @keywords lengths
#' @export
#' @examples
#' plotFragmentLengthDensities()

plotFragmentLengthDensities <- function(lengths=ll,col= c("#CDCD00", "#A0B43C", "#739B78", "#4682B4"),xlim=NULL){
   ld <- lapply(lengths,function(l) density(unlist(l),bw=15))
   ry <- range(sapply(ld,function(x) range(x$y)))
   if(is.null(xlim)) rx <- range(sapply(ld,function(x) range(x$x)))
   plot(rx,ry,type="n",ylab="density",xlab="tag lengths",xlim=rx,ylim=ry)
   for(ii in 1:length(ld)){
     lines(ld[[ii]],col=col[ii],lwd=2)
   }
   legend("topleft",names(lengths),lty=1,bty="n",col=col,lwd=2)
   }

