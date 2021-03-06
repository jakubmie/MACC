% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.tags.r
\name{read.tags}
\alias{read.tags}
\title{A function to read positions of aligned tags.}
\usage{
read.tags(setlist = NULL, path2files = NULL, chrn = NULL, mc.cores = 2,
  filter.anomalies = TRUE)
}
\arguments{
\item{setlist}{a list of vectors with names of files contains positions of aligned reads. Each vector should correspond to one sample, and order of files should correspond to ascending order of MNase concentrations used in the experiment.}

\item{path2files}{a path to files to be read.}

\item{chrn}{a vector of chromosome names or NULL for all chromosomes.}

\item{mc.cores}{a number of cores for parallel computing.}

\item{filter.anomalies}{logical. Should tags associated with anomalously high read counts be removed?}
}
\description{
A function to read positions of aligned tags.
}
\examples{
read.tags()
}
\keyword{tags}

