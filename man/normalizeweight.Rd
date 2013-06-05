\name{normalizeweight}
\alias{normalizeweight}
\title{Normalize weights}
\usage{
  normalizeweight(log_weights)
}
\arguments{
  \item{log_weights}{Object of class \code{"numeric"}: a
  real-valued vector}
}
\value{
  The function returns an object of class
  \code{"data.frame"}, with columns for the chain indices,
  the chain values, the iteration indices, and the
  associated log density values.
}
\description{
  Simple function that normalize vectors (ie takes log
  weights and returns normalized weights, in the SMC
  context).
}
\details{
  Simple function that takes log weights (ie any
  real-valued vector), computes the exponential of it,
  divides it by its sum and returns it.
}
\author{
  Luke Bornn <bornn@stat.harvard.edu>, Pierre E. Jacob
  <pierre.jacob.work@gmail.com>
}
\seealso{
  \code{\link{smc}}
}
\keyword{~kwd1}
\keyword{~kwd2}

