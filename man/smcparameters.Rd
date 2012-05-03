\name{smcparameters}
\Rdversion{1.1}
\docType{class}
\alias{smcparameters-class}
\alias{smcparameters}
\alias{smcparameters,ANY-method}
\alias{show,smcparameters-method}

\title{Class \code{"smcparameters"}}
\description{
    This class holds parameters for the Sequential Monte Carlo sampler.
}
\section{Objects from the Class}{
Objects can be created by calls of the function \code{"smcparameters"}.
}
\section{Slots}{
  \describe{
    \item{\code{nparticles}:}{Object of class \code{"numeric"}: an integer
    representing the desired number of particles.}
    \item{\code{temperatures}:}{Object of class \code{"numeric"}: }
    \item{\code{nmoves}:}{Object of class \code{"numeric"}: }
    \item{\code{ESSthreshold}:}{Object of class \code{"numeric"}: }
    \item{\code{movetype}:}{Object of class \code{"numeric"}: }
    \item{\code{movescale}:}{Object of class \code{"numeric"}: }
    \item{\code{resamplingscheme}:}{Object of class \code{"numeric"}: }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "smcparameters")}: provides a little summary of a binning object when called (or when \code{print} is called).}
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
showClass("smcparameters")
smcparam<- smcparameters(nparticles=5000, 
                        temperatures = seq(from = 0.0001, to = 1, length.out= 100),
                        nmoves = 5, ESSthreshold = 0.5, movetype = "randomwalk",
                        movescale = 0.1)
}
\keyword{classes}

