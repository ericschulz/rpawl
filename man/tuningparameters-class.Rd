\name{tuningparameters-class}
\Rdversion{1.1}
\docType{class}
\alias{tuningparameters-class}
\alias{tuningparameters}
\alias{show,tuningparameters-method}

\title{Class \code{"tuningparameters"}}
\description{
    This class holds tuning parameters for the Metropolis-Hastings and Wang-Landau algorithms.
}
\section{Objects from the Class}{
Objects can be created by calls of the function \code{"tuningparameters"}.
}
\section{Slots}{
  \describe{
    \item{\code{nchains}:}{Object of class \code{"numeric"}: it should be an integer
    representing the desired number of parallel chains.}
    \item{\code{niterations}:}{Object of class \code{"numeric"}: it should be an integer
    representing the desired number of iterations.}
    \item{\code{adaptiveproposal}:}{Object of class \code{"logical"}: specifies whether an adaptive
    proposal should be used, if available.}
    \item{\code{adaptationrate}:}{Object of class \code{"function"}: specifies the rate at
    which the adaptation of the proposal is performed. It should be a function defined on [0, + infty[ such that it is not integrable but its square is integrable, e.g. t -> 1/t for instance. The default is t -> t^{-0.6}.
    }
    \item{\code{sigma_init}:}{Object of class \code{"numeric"}: it should be a positive real number specifying the standard deviation of the proposal distribution at the first iteration. If the proposal is adaptive, it acts as a starting point for the adaptation. If it is not adaptive, then this value is used throughout all the iterations.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "tuningparameters")}: provides a little summary of a binning object when called (or when \code{print} is called).}
	 }
}
%\references{
%%  ~~put references to the literature/web site here~~
%}
\author{
    Pierre Jacob <pierre.jacob@ensae.fr>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{adaptiveMH}}
\code{\link{preexplorationAMH}}
}
\examples{
showClass("tuningparameters")
mhparameters <- tuningparameters(nchains = 10, niterations = 1000, adaptiveproposal = TRUE) 
}
\keyword{classes}
