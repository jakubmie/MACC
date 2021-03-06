% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getMACCscores.r
\name{getMACCscores}
\alias{getMACCscores}
\title{A wrapper function to compute MNase accessibility (MACC) scores with one command.}
\usage{
getMACCscores(setlist = NULL, path2files = path2files, bin = 300,
  tit.points = NULL, genome = NULL, chrn = NULL, mc.cores = 2,
  filter.anomalies = TRUE, scale.genome.to.100Mb = TRUE,
  normalize.perBinSize = TRUE)
}
\arguments{
\item{setlist}{a list of vectors with names of files contains positions of aligned reads. Each vector should correspond to one sample, and order of files should correspond to ascending order of MNase concentrations used in the experiment.}

\item{path2files}{a path to files to be read.}

\item{bin}{a size of a bin.}

\item{tit.points}{a vector of MNase concentrations used in the experiment. The concentrations should be provided in ascending order.}

\item{genome}{a string corresponding to genome of interest.}

\item{chrn}{a vector of chromosome names or NULL for all chromosomes.}

\item{mc.cores}{a number of cores for parallel computing.}

\item{scale.genome.to.100Mb}{logical. Should the values be normalized to genome size?}

\item{normalize.perBinSize}{logical. Should the values be normalized to been size?}
}
\description{
A wrapper function to compute MNase accessibility (MACC) scores with one command.
}
\examples{
getMACCscores()
}
\keyword{MACC}

