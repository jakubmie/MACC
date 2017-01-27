#' A function to compute positions of sequenced tags.
#'
#' @param cp a vector with tag coordinates.
#' @param max.length  an integer corresponding to maximal tag length.
#' @param shift  an integer corresponding to shift of location applied to tags with length longer than max.length.
#' @keywords tags
#' @export
#' @examples
#' tag.centers()

tag.centers <-  function(cp=NULL,max.length=160,shift=73){#c <- tags[[1]]
			if(is.null(cp)){stop("Alignment positions were not provided")
			}else{
				tp <- which(cp>0); 
				tn <- which(cp<0); 
				cen <- ((abs(cp[tn])+abs(cp[tp]))/2 )
				w   <- which( (abs(cp[tn])-abs(cp[tp])  )  >max.length)
				if(length(w)>0){
				  b   <- sort(sample(1:length(w),floor(length(w)/2)))
				  cen[w][b]  <- abs(cp[tn][w[b]]+shift)
				  cen[w][-b] <- abs(cp[tp][w[-b]]+shift)
				  }
				return( cen )
			}
      }

