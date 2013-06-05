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

