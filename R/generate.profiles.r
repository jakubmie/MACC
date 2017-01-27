#' A function to read positions of aligned tags.
#'
#' @param tags a list of signle positions representing used reads. 
#' @param chrn a vector of chromosome names or NULL for all chromosomes.
#' @param mc.cores a number of cores for parallel computing.
#' @keywords tags
#' @export
#' @examples
#' generate.profiles()

generate.profiles <- function(tags,chrn=chrn,mc.cores=2){
    fcen <- flen  <- list()
    for (nsets in names(tags)) {#nsets <- names(setlist)[1]
        sets <- setlist[[nsets]]
        names(sets) <- sets
        for (set in sets) {#sets[1]->set
            ttags <- tags[[nsets]][[set]]
            fcen[[nsets]][[set]] <- parallel::mclapply(chrn, 
                function(c) {
                  tag.centers(cp = ttags[[c]], max.length = 160, 
                    shift = 73)
                }, mc.cores = mc.cores)
            flen[[nsets]][[set]] <- parallel::mclapply(chrn, 
                function(c) {
                  tag.lengths(cp = ttags[[c]])
                }, mc.cores = mc.cores)
            names(fcen[[nsets]][[set]]) <- names(flen[[nsets]][[set]])  <-chrn
            rm(ttags)
            gc()
            cat("done\n")
        }
    }
    return(list(frag.pos = fcen, lengths = flen))
}

