#' A wrapper function to compute MNase accessibility (MACC) scores with one command.
#'
#' @param setlist a list of vectors with names of files contains positions of aligned reads. Each vector should correspond to one sample, and order of files should correspond to ascending order of MNase concentrations used in the experiment.
#' @param path2files a path to files to be read.
#' @param bin a size of a bin.
#' @param tit.points a vector of MNase concentrations used in the experiment. The concentrations should be provided in ascending order.
#' @param genome a string corresponding to genome of interest.
#' @param chrn a vector of chromosome names or NULL for all chromosomes.
#' @param mc.cores a number of cores for parallel computing.
#' @param scale.genome.to.100Mb logical. Should the values be normalized to genome size?
#' @param normalize.perBinSize logical. Should the values be normalized to been size?
#' @keywords MACC
#' @export
#' @examples
#' getMACCscores()


getMACCscores <- function(setlist=NULL,path2files=path2files,bin=300,tit.points=NULL,genome=NULL,chrn=NULL,mc.cores=2,filter.anomalies=TRUE,scale.genome.to.100Mb=TRUE,normalize.perBinSize=TRUE){
		gen.dat      <- read.genomic.data(genome=genome,chrn=chrn,mc.cores=mc.cores)
		tags         <- read.tags(setlist=setlist,path2files=path2files,chrn=gen.dat$chrn,mc.cores=mc.cores,filter.anomalies=filter.anomalies)
		profiles     <- generate.profiles(tags,chrn=gen.dat$chrn,mc.cores=mc.cores)
		tags.in.bins <- count.tags.in.bins(profs=profiles$frag.pos,chr.stat=gen.dat$chr.stat,chrn=gen.dat$chrn,mc.cores=mc.cores,bin=bin, 				  scale.genome.to.100Mb=scale.genome.to.100Mb,normalize.perBinSize=normalize.perBinSize)
		Res         <- macc(tags.in.bins=tags.in.bins,tit.points=tit.points,gc.cont=gen.dat$gc.cont,chr.stat=gen.dat$chr.stat, chrn=gen.dat$chrn,bin=bin, mc.cores=mc.cores,CpG=gen.dat$CpG)
		Res
		}


		