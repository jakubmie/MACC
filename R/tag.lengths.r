#' A function to compute lengths of sequenced tags.
#'
#' @param cp a vector with tag coordinates.
#' @keywords tags
#' @export
#' @examples
#' tag.lengths()

tag.lengths <-  function(cp=NULL){#c <- tags[[1]]
			if(is.null(cp)){stop("Alignment positions were not provided")
			}else{
				tp <- which(cp>0); 
				tn <- which(cp<0); 
				len<- abs(cp[tn]+cp[tp] )
				return( len )
			}
      }

