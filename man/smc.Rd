\name{smc}
\alias{smc}
\title{Sequential Monte Carlo}
\usage{
  smc(target, AP, verbose = TRUE)
}
\arguments{
  \item{target}{Object of class \code{\link{target}}:
  specifies the target distribution.  See the help of
  \code{\link{target}}. The target must be defined on a
  continuous state space.}

  \item{AP}{Object of class \code{"smcparameters"}:
  specifies the number of particles, the ESS threshold, the
  sequence of distributions, etc. See the help of
  \code{\link{smcparameters}}.}

  \item{verbose}{Object of class \code{"logical"}: if TRUE
  (default) then prints some indication of progress in the
  console.}
}
\value{
  The function returns a list holding various information:
  \itemize{ \item particles a matrix with rows representing
  particles and columns components of each particle. \item
  weights a vector of weights associated to each particle.
  See also the convenience function
  \code{\link{normalizeweight}}. \item ESSarray a vector of
  the ESS computed at each iteration. \item resamplingtimes
  a vector indicating the iteration at which resampling was
  performed. }
}
\description{
  Sequential Monte Carlo samplers, using a sequence of
  tempered distributions.
}
\author{
  Luke Bornn <bornn@stat.harvard.edu>, Pierre E. Jacob
  <pierre.jacob.work@gmail.com>
}
\seealso{
  \code{\link{smcparameters}}
}
\keyword{~kwd1}
\keyword{~kwd2}

