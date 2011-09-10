\name{tuningparameters}
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
    \item{\code{computemean}:}{Object of class \code{"logical"}: specifies whether the mean of all chains should
        be computed at each iteration (useful if the chains are not to be stored).}
    \item{\code{saveeverynth}:}{Object of class \code{"numeric"}: specifies when the chains are to be stored:
    for instance at every iteration (=1), every 10th iteration (=10), etc. Default is -1, meaning the chains are not stored.}
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
  Luke Bornn <l.bornn@stat.ubc.ca>,
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
