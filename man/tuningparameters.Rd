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

\section{Slots}{
  \describe{ \item{\code{nchains}:}{Object of class
  \code{"numeric"}: it should be an integer representing
  the desired number of parallel chains.}
  \item{\code{niterations}:}{Object of class
  \code{"numeric"}: it should be an integer representing
  the desired number of iterations.}
  \item{\code{computemean}:}{Object of class
  \code{"logical"}: specifies whether the mean of all
  chains should be computed at each iteration (useful if
  the chains are not to be stored).}
  \item{\code{computemeanburnin}:}{Object of class
  \code{"numeric"}: if \code{computemean} is set to TRUE,
  specifies after which iteration the mean of the chain has
  to be computed. Default is 0 (no burnin).}
  \item{\code{saveeverynth}:}{Object of class
  \code{"numeric"}: specifies when the chains are to be
  stored: for instance at every iteration (=1), every 10th
  iteration (=10), etc. Default is -1, meaning the chains
  are not stored.} }
}

\section{Methods}{
  \describe{ \item{show}{\code{signature(object =
  "tuningparameters")}: provides a little summary of a
  binning object when called (or when \code{print} is
  called).} }
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

