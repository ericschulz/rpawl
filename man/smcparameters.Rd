\docType{class}
\name{smcparameters}
\alias{show,smcparameters-method}
\alias{smcparameters}
\alias{smcparameters,ANY-method}
\alias{smcparameters-class}
\title{SMC Tuning Parameters}
\description{
  This class holds parameters for the Sequential Monte
  Carlo sampler.
}
\section{Objects from the Class}{
  Objects can be created by calls of the function
  \code{"smcparameters"}.
}

\section{Slots}{
  \describe{ \item{\code{nparticles}:}{Object of class
  \code{"numeric"}: an integer representing the desired
  number of particles.} \item{\code{temperatures}:}{Object
  of class \code{"numeric"}: a vector of temperatures,
  default being \code{"seq(from = 0.01, to = 1, length.out
  = 100)"}.} \item{\code{nmoves}:}{Object of class
  \code{"numeric"}: number of move steps to be performed
  after each resampling step, default being 1.}
  \item{\code{ESSthreshold}:}{Object of class
  \code{"numeric"}: resampling occurs when the Effective
  Sample Size goes below \code{"ESSthreshold"} multiplied
  by the number of particles \code{"nparticles"}.}
  \item{\code{movetype}:}{Object of class
  \code{"character"}: type of Metropolis-Hastings move step
  to be performed; can be either set to
  \code{"independent"} or \code{"randomwalk"}, default
  being \code{"independent"}.}
  \item{\code{movescale}:}{Object of class
  \code{"numeric"}: if \code{movetype} is set to
  \code{"randomwalk"}, this parameter specifies the amount
  by which the estimate of the standard deviation of the
  target distribution is multiplied; the product being used
  to propose new points in the random-walk MH step. Default
  is 10\%, ie a new point is proposed from a Normal
  distribution, centered on the latest point, with standard
  deviation equal to 10\% of the standard deviation of the
  already-generated chain.}
  \item{\code{resamplingscheme}:}{Object of class
  \code{"character"}: type of resampling to be used; either
  "multinomial", "residual" or "systematic", the default
  being "systematic".} }
}

\section{Methods}{
  \describe{ \item{show}{\code{signature(object =
  "smcparameters")}: provides a little summary of a binning
  object when called (or when \code{print} is called).} }
}
\examples{
showClass("smcparameters")
smcparam<- smcparameters(nparticles=5000,
                       temperatures = seq(from = 0.0001, to = 1, length.out= 100),
                       nmoves = 5, ESSthreshold = 0.5, movetype = "randomwalk",
                       movescale = 0.1)
}
\author{
  Luke Bornn <bornn@stat.harvard.edu>, Pierre E. Jacob
  <pierre.jacob.work@gmail.com>
}
\seealso{
  \code{\link{smc}}
}
\keyword{classes}

