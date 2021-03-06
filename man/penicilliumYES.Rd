% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{penicilliumYES}
\alias{penicilliumYES}
\title{Data set of three species of Penicillium fungi, from \emph{sparseLDA}}
\format{This data set contains the following matrices:
\describe{
  \item{X}{A matrix with 36 columns and 3754 rows. The training and test data. The first 12
   rows are \emph{P. Melanoconidium} species, rows 13-24 are \emph{P. Polonicum} species,
    and the last 12 rows are \emph{P. Venetum} species. The samples are ordered so that
    each pair of three is from the same isolate.}
  \item{Y}{A matrix of dummy variables for the training data.}
  \item{Z}{Z matrix of probabilities for the subcalsses of the training data.}
}}
\source{
\url{http://www.imm.dtu.dk/~lhc}.
}
\usage{
data(penicilliumYES)
}
\description{
The data set \code{penicilliumYES} has \emph{36} rows and \emph{3754} columns. The variables
are 1st order statistics from multi-spectral images of three species of Penicillium fungi:
 \emph{Melanoconidium}, \emph{Polonicum}, and \emph{Venetum}. These are the data used in
 the Clemmemsen et al "Sparse Discriminant Analysis" paper.
}
\details{
The X matrix is not normalized.
}
\references{
{
  \insertRef{clemmensen2007method}{ddsPLS}
}
}
\keyword{datasets}
