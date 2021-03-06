% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.genomic.data.r
\name{read.genomic.data}
\alias{read.genomic.data}
\title{A function to read genomic data.}
\usage{
read.genomic.data(genome = "", chrn = NULL, bin = 300, mc.cores = 2)
}
\arguments{
\item{genome}{a string corresponding to genome of interest.}

\item{chrn}{a vector of chromosome names or NULL for all chromosomes.}

\item{bin}{a size of a bin.}

\item{mc.cores}{a number of cores for parallel computing.}
}
\description{
A function to read genomic data.
}
\examples{
read.genomic.data()
}
\keyword{genome}

