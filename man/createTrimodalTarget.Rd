\name{createTrimodalTarget}
\alias{createTrimodalTarget}
\title{Trimodal target distribution}
\usage{
  createTrimodalTarget()
}
\value{
  The function returns an object of class
  \code{\link{target}}, with a name, a dimension, a
  function giving the log density, a function to generate
  sample from the distribution, parameters of the
  distribution, and a function to draw init points for the
  MCMC algorithms.
}
\description{
  Create the trimodal distribution as in Liang, Liu and
  Caroll, 2007: Stochastic approximation in Monte Carlo
  computation.
}
\details{
  This distribution is a mixture of three bivariate
  Gaussian distributions. Their covariance matrices are
  such that an adaptive MCMC algorithm which proposal
  variance adapts to one of the component, will likely fail
  to explore the others.
}
\author{
  Luke Bornn <bornn@stat.harvard.edu>, Pierre E. Jacob
  <pierre.jacob.work@gmail.com>
}
\references{
  Liang, Liu and Caroll: Stochastic approximation in Monte
  Carlo computation. Published in JASA, 2007.
}
\seealso{
  \code{\link{target}}, \code{\link{createMixtureTarget}}
}

