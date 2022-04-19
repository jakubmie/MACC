#' A function to compute MNase accessibility (MACC) scores.
#'
#' @param tags.in.bins an list returned by count.tags.in.bins function and contatining localization of bins and number of tags in each of the selected bins.
#' @param tit.points a vector of MNase concentrations used in the experiment. The concentrations should be provided in ascending order.
#' @param gc.cont a list of vectors containing information about GC content. Each vector corresponds to different chromosome.
#' @param chrn a vector of chromosome names or NULL for all chromosomes.
#' @param chr.stat a data.frame with sizes of analyzed chromosmes.
#' @param mc.cores a number of cores for parallel computing.
#' @param bin a size of a bin.
#' @param CpG a GenomeRanges object with localization of CpG island.
#' @param scale.genome.to.100Mb logical. Should the values be normalized to genome size?
#' @param normalize.perBinSize logical. Should the values be normalized to been size?#' @keywords MACC
#' @export
#' @examples
#' macc()


macc <- function (tags.in.bins = NULL, tit.points = NULL, gc.cont = NULL, 
    chrn = NULL, chr.stat = NULL, mc.cores = 2, bin = 300, CpG = NULL) 
{
    if (is.null(tags.in.bins)) {
        stop("File names were not provided")
    }
    if (is.null(tit.points)) {
        stop("MNase titration concentrations were not provided")
    }
    if (is.null(gc.cont)) {
        stop("GC contents were not provided")
    }
    if (is.null(chr.stat)) {
        stop("Lengths of chromosoms were not provided")
    }
    if (is.null(chrn)) {
        stop("Chromosome names were not provided")
    }
    if (!all(chrn %in% chr.stat$chr)) {
        stop("Some chromosom lengths are missing")
    }
    if (any(names(chrn) != chrn)) {
        names(chrn) <- chrn
    }
    if (any(sapply(tags.in.bins, function(x) length(x$values)) != 
        length(tit.points))) {
        stop(paste("Number of MNase profiles in are different than number of provided MNase concentrations", 
            sep = ""))
    }
    if (any(sapply(tags.in.bins, function(x) sapply(x$values, 
        function(y) length(unlist(y)))) != length(unlist(gc.cont)))) {
        stop(paste("Number of bins in are different than in provided gc.cont object", 
            sep = ""))
    }
    Res <- lapply(1:length(tags.in.bins), function(ijk) {#ijk<-1
        frg.bins <- tags.in.bins[[ijk]]
        Pos <- lapply(names(frg.bins$pos), function(nn) {
            vec <- frg.bins$pos[[nn]]
            data.frame(chr = rep(nn, length(vec)), start = as.integer(vec - 
                (bin/2) + 1), stop = as.integer(vec + bin/2))
        })
        Pos <- do.call("rbind", Pos)
        val <- parallel::mclapply(1:length(names(frg.bins$values)), 
            function(ii) {
                unlist(lapply(chrn, function(c) {
                  frg.bins$values[[ii]][[c]]
                }))
            }, mc.cores = mc.cores)
        names(val) <- paste("mn", tit.points, sep = "")
        cntsms <- do.call(cbind, val)
        coefs <- rep(0, nrow(cntsms))
        no.data <- which(rowSums(cntsms) == 0)
        for.slope <- setdiff(1:nrow(cntsms), no.data)
        coefs[for.slope] <- unlist(parallel::mclapply(for.slope,function(jj) {
                coef(lm(cntsms[jj, ] ~ (log2(tit.points))))[2] * (-1)
            }, mc.cores = mc.cores))
        if (length(no.data) > 0) 
            coefs[no.data] <- NA
        Mat <- cbind(Pos, cntsms, rowMeans(cntsms),coefs)
        colnames(Mat)[c(1:3, ncol(Mat)-1,ncol(Mat))] <- c("chr", "start", 
            "end","Pooled" ,"Slope")
        cat("slope for ", paste(names(tags.in.bins)[ijk], " was computed", 
            sep = ""), "\n")
        if (is.null(CpG)) {
            lFit <- limma::loessFit(coefs[for.slope], unlist(gc.cont)[for.slope], 
                span = 0.15)$fitted
            m <- coefs
            m[for.slope] <- coefs[for.slope] - lFit
        }
        else {
            cat("---->\t\tBeen here\t\t<---------")
            Vgr <- GenomicRanges::GRanges(seqnames = as.vector(Pos[, 
                1]), ranges = IRanges::IRanges(start = as.vector(Pos[, 
                2]), end = as.vector(Pos[, 3])), strand = "*")
            cpg <- as.vector(as.numeric(countOverlaps(Vgr, CpG, 
                minoverlap = bin) != 0))
            sum(cpg)
            no.in.cpg <- sort(intersect(for.slope, which(cpg != 
                1)))
            in.cpg <- sort(intersect(for.slope, which(cpg == 
                1)))
            lFit1 <- limma::loessFit(coefs[no.in.cpg], unlist(gc.cont)[no.in.cpg], 
                span = 0.25)$fitted
            lFit2 <- limma::loessFit(coefs[in.cpg], unlist(gc.cont)[in.cpg], 
                span = 0.25)$fitted
            m <- coefs
            m[no.in.cpg] <- coefs[no.in.cpg] - lFit1
            m[in.cpg] <- coefs[in.cpg] - lFit2
        }
        Mat <- cbind(Mat, m); colnames(Mat)[ncol(Mat)] <- c("MACC")
        cat("MACC for ", paste(names(tags.in.bins)[ijk], " was computed", 
            sep = ""), "\n")
        return(Mat)
    })
    return(Res)
} 
