\docType{class}
\name{binning}
\alias{binning}
\alias{binning,ANY-method}
\alias{binning-class}
\alias{binning-method}
\alias{binning-methods}
\alias{show,binning-method}
\title{Class \code{"binning"}}
\description{
  This class holds all the parameters of the Parallel
  Adaptive Wang-Landau algorithm that are related to the
  bins: it includes the functions that take points and
  return the point locations with respect to the bins,
  parameters related to the number of bins, the split
  mechanism, the adaptation rate of the stochastic
  approximation schedule, etc.
}
\section{Objects from the Class}{
  Objects should created by calls of the function
  \code{binning}. Examples are provided that should help
  understanding this class. Essentially it is a list of
  parameters, most of which have a reasonable default value
  so you do not need to think about it too much.
}
\examples{
showClass("binning")
 getPos <- function(points, logdensity) points
 positionbinning <- binning(position = getPos,
                       name = "position",
                       binrange = c(-4, 0),
                       ncuts = 4,
                       autobinning = TRUE,
                       useLearningRate = TRUE)
}
\author{
  Luke Bornn <bornn@stat.harvard.edu>, Pierre E. Jacob
  <pierre.jacob.work@gmail.com>
}
\keyword{classes}

