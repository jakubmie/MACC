#' A function to remove tag anomalies.
#'
#' @param tv a vector with tag coordinates from the same strand.
#' @keywords tags 
#' @examples
#' rta()


rta <- function(tv) {                                                                                    
  tt <- table(tv);                                                                                                                                                                                                                                            
  stt <- sort(as.numeric(tt));                                                                                                                     
  stt <- stt[1:(length(stt))];                                                                                                   
  mstt <- mean(stt);                                                                                                                               
  if(mstt > 5){mtc <- mstt}else{mtc <- 5}; tcd <- sqrt(var(stt));                                                                                                                                                                                                                                                 
  thr <- ceiling(mtc+7*tcd);                                                                                                                       
  thr.o <- ceiling(mtc+10*tcd);                                                                                                                                                                                                                                                      
  tt <- tt[tt > thr]                                                                                                                                                                                                                                                         
  tp <- as.numeric(names(tt));                                                                                                                     
  pti <- tp>0;                                                                                                                                     
  it <- intersect(tp[pti],(-1)*tp[!pti]);                                                                                                                                                                                                                                 
  it <- unique(c(it,tp[tt > thr.o]));                                                                                                              
  sit <- c(it,(-1)*it);                                                                                                                                                                                                                                                                              
  return(!tv %in% sit);                                                                                                                                                                                                                                                                                                                                                                                         
} 

