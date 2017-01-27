#' A function to plot scatter plot showing association between MACC scores and GC content 
#'
#' @param data a table which typically is one of the element of a list return by macc function.
#' @param sname a character string giving the header of the plot.
#' @param ylim a vector containing ranges for the y-axis.
#' @param gc.cont a results of ount.gc.count function containing GC content in the selected bins.
#' @keywords GC correction
#' @export
#' @examples
#' plotGCvsMACC()


plotGCvsMACC <- function(data,sname,ylim=NULL,gc.cont=gc.cont){
	for.slope <- which(is.finite(data$MACC))
	gc.c      <- unlist(gc.cont)[for.slope]
	sl        <- data$Slope[for.slope]
	lF        <- sl-data$MACC[for.slope]
	rr        <- quantile(abs(sl),prob=.999)
	if(is.null(ylim)){ylim <- c(-rr,rr)}
	smoothScatter(gc.c,sl,ylim=ylim,ylab="",xlab="")
	ind <- which(is.finite(gc.c) & is.finite(sl-lF))
	lines(lowess(gc.c[ind],(sl-lF)[ind],f=1),lwd=2,col="red4")
	lines(lowess((sl)~gc.c[ind]),lwd=2,col="blue4")
	mtext(side=1,"GC content",line=2.2,cex=1.15);mtext(side=2,"MACC",line=2.2,cex=1.15)
	mtext(side=3,sname,cex=1.4,line=.2)
	legend("topright",c("non-GC corrected","GC corrected"),bty="n",lwd=2,col=c("blue4","red4"))
}

