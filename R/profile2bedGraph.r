#' A wrapper function to compute MNase accessibility (MACC) scores with one command.
#'
#' @param tab a matrix with four columns corresponding to chromosome names, start of the bins, end of the bins and interesting values respectively.
#' @param path2files a character string giving the name of directory where the new file will be saved.
#' @param fname a character string giving the name of the file to be created.
#' @param nn integer indicating the number of decimal places to be used.
#' @keywords bedGraph
#' @export
#' @examples
#' profile2bedGraph()


profile2bedGraph <- function(tab,path="./",fname="profile",nn=4){
		  tab <- tab[which(is.finite(tab[,4])),]
		  tab[,4] <- round(tab[,4],nn)
		  nname <- fname
		  fname <- paste(path,fname,".bedGraph",sep="") 
		  write(paste('track type=bedGraph name=',nname,' visibility=dense color=255,51,51 altColor=0,100,200 priority=20',sep=""),file=fname)
		  write.table(tab,file=fname,row.names=FALSE,col.names=FALSE,quote=FALSE,sep=" ",append=TRUE) 
		  }

