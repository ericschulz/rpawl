\name{pawl}
\alias{pawl}
\title{Parallel Adaptive Wang-Landau}
\usage{
  pawl(target, binning, AP, proposal, verbose = TRUE)
}
\arguments{
  \item{target}{Object of class \code{\link{target}}:
  specifies the target distribution.  See the help of
  \code{\link{target}}. If the target is discrete, target
  must contain the slots \code{dproposal}, \code{rproposal}
  and \code{proposalparam} that specify the proposal kernel
  in the Metropolis-Hastings step. Otherwise the default is
  an adaptive gaussian random walk.}

  \item{binning}{Object of class \code{\link{binning}},
  defining the initial bins used by the Wang-Landau
  algorithm.  The binning object also contains some
  parameters specifying if the automatic binning mechanism
  is active or not, for instance.}

  \item{AP}{Object of class \code{\link{tuningparameters}}:
  specifies the number of chains, the number of iterations,
  and what should be stored during along the run. See the
  help of \code{\link{tuningparameters}}.}

  \item{proposal}{Object of class \code{\link{proposal}}:
  specifies the proposal distribution to be used to propose
  new values and to compute the acceptance rate. See the
  help of \code{\link{proposal}}. If this is not specified
  and the target is continuous, then the default is an
  adaptive gaussian random walk.}

  \item{verbose}{Object of class \code{"logical"}: if TRUE
  (default) then prints some indication of progress in the
  console.}
}
\value{
  The function returns a list holding various information:

  and other quantities, that you can browse by calling
  \code{"names(results)"} where \code{"results"} is the
  output of the function. \itemize{ \item finalchains The
  last point of each chain. \item acceptrates The vector of
  acceptance rates at each step. \item sigma The vector of
  the standard deviations used by the MH kernel along the
  iterations. If the proposal was adaptive, this allows to
  check how the adaptation behaved. \item allchains If
  asked in the tuning parameters, the chain history. \item
  alllogtarget If asked in the tuning parameters, the
  associated log density evaluations. \item meanchains If
  asked in the tuning parameters, the mean (component-wise)
  of each chain. \item logthetahistory If asked in the
  tuning parameters, all the log theta penalties. }
}
\description{
  Implements the Parallel Adaptive Wang-Landau algorithm.
}
\author{
  Luke Bornn <bornn@stat.harvard.edu>, Pierre E. Jacob
  <pierre.jacob.work@gmail.com>
}
\seealso{
  \code{\link{adaptiveMH}, \link{binning}}
}
\keyword{~kwd1}
\keyword{~kwd2}

