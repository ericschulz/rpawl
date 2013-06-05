\docType{class}
\name{tuningparameters}
\alias{show,tuningparameters-method}
\alias{tuningparameters}
\alias{tuningparameters,ANY-method}
\alias{tuningparameters-class}
\title{MCMC Tuning Parameters}
\description{
  This class holds tuning parameters for the
  Metropolis-Hastings and Wang-Landau algorithms.
}
\section{Objects from the Class}{
  Objects can be created by calls of the function
  \code{"tuningparameters"}.
}
\examples{
showClass("tuningparameters")
mhparameters <- tuningparameters(nchains = 10, niterations = 1000, adaptiveproposal = TRUE)
}
\author{
  Luke Bornn <bornn@stat.harvard.edu>, Pierre E. Jacob
  <pierre.jacob.work@gmail.com>
}
\seealso{
  \code{\link{adaptiveMH}} \code{\link{preexplorationAMH}}
  \code{\link{pawl}}
}
\keyword{classes}

