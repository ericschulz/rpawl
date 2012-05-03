\name{preexplorationAMH}
\alias{preexplorationAMH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Pre exploration Adapative Metropolis-Hastings
}
\description{
  This function takes a target distribution, an integer representing the number of parallel chains,
  and an integer representing a number of iterations, and runs adaptive Metropolis-Hastings algorithm
  using them. The chains are then used to create a range called SuggestedRange, to be used to bin
  the state space according to the energy levels. The energy is here defined as minus the log density of the 
  target distribution.
}
\usage{
preexplorationAMH(target, nchains, niterations, proposal, verbose)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{target}{Object of class \code{"target"}: this argument describes the target distribution.
        See \code{\link{target}} for details.
    }
    \item{nchains}{Object of class \code{"numeric"}: specifies the number of parallel chains.
    }
    \item{niterations}{Object of class \code{"numeric"}: specifies the number of iterations.
    }
    \item{proposal}{Object of class \code{"proposal"}: specifies the proposal distribution to be used to 
        propose new values and to compute the acceptance rate. See the help of \code{\link{proposal}}. If this
        is not specified and the target is continuous, then the default is an adaptive gaussian random walk.
    }
    \item{verbose}{
        Object of class \code{"logical"}: if TRUE (default) then prints some indication of progress
        in the console.
    }
}
\details{
    The adaptive Metropolis-Hastings algorithm used in the function is described in more details 
    in the help page of \code{\link{adaptiveMH}}
}
\value{
    The function returns a list holding the following entries:
    \item{LogEnergyRange}{This holds the minimum and maximum energy values seen by the chains during the exploration.}
    \item{LogEnergyQtile}{Returns the first 10\% quantile of the energy values seen by the chains during the exploration.}
    \item{SuggestedRange}{This holds the suggested range, that is, the first 10\% quantile and the maximum value
        of the energy values seen during the exploration. This can be passed as the \code{binrange} argument of the
        \code{binning} class, see the \code{trimodal} example.}
    \item{finalchains}{The last point of each chain.}
}

%\references{
%% ~put references to the literature/web site here ~
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
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
