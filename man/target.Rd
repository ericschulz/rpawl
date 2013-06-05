\docType{class}
\name{target}
\alias{show,target-method}
\alias{target}
\alias{target,ANY-method}
\alias{target-class}
\title{Class: target distribution}
\description{
  This class represents target distributions, that is,
  probability distributions from which we want to sample
  using MCMC or Wang-Landau.
}
\section{Objects from the Class}{
  Objects should created by calls of the function
  \code{target}. Examples are provided that should help
  implementing any continuous probability distributions. %%
  ~~ describe objects here ~~
}
\examples{
showClass("target")
 # starting points for MCMC algorithms
 rinit <- function(size) rnorm(size)
 # target log density function: a gaussian distribution N(mean = 2, sd = 3)
 parameters <- list(mean = 2, sd = 3)
 logdensity <- function(x, parameters) dnorm(x, parameters$mean, parameters$sd, log = TRUE)
 # creating the target object
 gaussiantarget <- target(name = "gaussian", dimension = 1,
                   rinit = rinit, logdensity = logdensity,
                   parameters = parameters)
 print(gaussiantarget)
}
\author{
  Luke Bornn <bornn@stat.harvard.edu>, Pierre E. Jacob
  <pierre.jacob.work@gmail.com>
}
\keyword{classes}

