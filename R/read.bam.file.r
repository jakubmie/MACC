#' A function to count GC content in the selected bins.
#'
#' @param fname a character showing a path to files to be read. 
#' @param mc.cores a number of cores for parallel computing.
#' @param tl an integer for tag shifting.
#' @keywords bam
#' @export
#' @examples
#' read.bam.file()

read.bam.file <- function(fname,mc.cores=2,tl=50) {  #require(Rsamtools) 
    bam1 <- Rsamtools::scanBam(fname, param = Rsamtools::ScanBamParam(what = c("flag", 
        "rname", "pos", "mpos", "strand", "qwidth", "qname"), 
        flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)))[[1]]
    indp <- which((bam1$flag == 99 | bam1$flag == 163))
    if (length(indp) > 0) {
        bam1$rname <- as.vector(bam1$rname)
        revi <- which((bam1$flag == 83 | bam1$flag == 147))
        tl <- mean(bam1$qwidth[revi], na.rm = T)
        bam <- list()
        bam$pos <- rep(0, length(indp) * 2)
        bam$pos[seq(1, length(bam$pos) - 1, by = 2)] <- bam1$pos[indp]
        bam$pos[seq(2, length(bam$pos), by = 2)] <- (-1 * (bam1$mpos[indp] + 
            tl))
        bam$rname <- rep("", length(indp) * 2)
        bam$rname[seq(1, length(bam$pos) - 1, by = 2)] <- as.character(bam1$rname[indp])
        bam$rname[seq(2, length(bam$pos), by = 2)] <- as.character(bam1$rname[indp])
        chrnn <- names(table(bam$rname))
        names(chrnn) <- chrnn
        tags <- parallel::mclapply(chrnn, function(n) {
            ii <- which(bam$rname == n)
            return(bam$pos[ii])
        }, mc.cores = mc.cores)
        return(tags)
    }
    else {
        stop("Input BAM file does not contain valid alignment")
    }
}
