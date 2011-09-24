\name{target}
\Rdversion{1.1}
\docType{class}
\alias{target-class}
\alias{target}
\alias{target,ANY-method}
\alias{show,target-method}

\title{Class: target distribution}
\description{
  This class represents target distributions, that is, probability distributions from which we want to sample using MCMC or Wang-Landau.
}
\section{Objects from the Class}{
Objects should created by calls of the function \code{target}.
Examples are provided that should help implementing any
continuous probability distributions.
%%  ~~ describe objects here ~~ 
}
\section{Important slots}{
  \describe{
    \item{\code{dimension}:}{Object of class \code{"numeric"}: should be an integer specifying the dimension of the state space on which the target distribution is defined.}
    \item{\code{logdensity}:}{Object of class \code{"function"} : should be a function taking n points in the state space and parameters, and returning a vector of n real values. See the example below. This function is in most cases the most time-consuming part in a MCMC algorithm, so make sure it runs reasonably fast!}
    \item{\code{rinit}:}{Object of class \code{"function"} : this function should take an integer as argument, say n. Then the function should return a matrix of dimension n times d (where d is the dimension of the state space), representing n points in the state space. These n points will be used as starting points of a parallel MCMC algorithm.}
  }
}
\section{Optional slots}{
  \describe{
    \item{\code{parameters}:}{Object of class \code{"list"} : you can put anything in that list (and nothing, which is the defaults), the important thing is that calls to \code{logdensity(x, parameters)}} return sensible values. For example, for a gaussian target distribution, you can put the mean and the variance in the \code{parameters} list (see example below). If need be, you can put a whole data set in there.
    \item{\code{type}:}{Object of class \code{"character"} : could be "continuous" or "discrete"; default is "continuous". }
    \item{\code{name}:}{Object of class \code{"character"} : ... if you want to name your distribution (default is "unspecified").}
    \item{\code{generate}:}{Object of class \code{"function"} : does not have to be specified, but if it is specified it should be a function to generate from the distribution (like rnorm is to the standard normal distribution).}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "target")}: provides a little summary of a target object when called (or when \code{print} is called).}
	 }
}
%\references{
%%  ~~put references to the literature/web site here~~
%}
\author{
    Luke Bornn <l.bornn@stat.ubc.ca>,
    Pierre Jacob <pierre.jacob@ensae.fr>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
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
\keyword{classes}
