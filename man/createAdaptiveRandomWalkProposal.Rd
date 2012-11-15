\name{createAdaptiveRandomWalkProposal}
\alias{createAdaptiveRandomWalkProposal}
\title{
    Adaptive Random Walk proposal distribution for MCMC algorithms
}
\description{
Create the adaptive gaussian random walk proposal that is used as a default in \code{\link{adaptiveMH}} and
\code{\link{pawl}}, whenever the target distribution is continuous.
}
\usage{
    createAdaptiveRandomWalkProposal(nchains, targetdimension, adaptiveproposal, adaptationrate, sigma_init)
}
\arguments{
    \item{nchains}{Object of class \code{"numeric"}: it should be an integer
        representing the desired number of parallel chains.}
    \item{targetdimension}{Object of class \code{"numeric"}: it should be an integer
        representing the dimension of the target distribution.} 
    \item{adaptiveproposal}{Object of class \code{"logical"}: specifies whether an adaptive
        proposal (Robbins-Monroe type of adaptation) should be used. Default is FALSE.}
    \item{adaptationrate}{Object of class \code{"function"}: specifies the rate at
        which the adaptation of the proposal is performed. It should be a function defined on [0, + infty[ such that it is not integrable but its square is integrable, e.g. t -> 1/t for instance. The default is t -> t^{-0.6}.}
    \item{sigma_init}{Object of class \code{"numeric"}: it should be a positive real number specifying the standard deviation of the proposal distribution at the first iteration. If the proposal is adaptive, it acts as a starting point for the adaptation. If it is not adaptive, then this value is used throughout all the iterations. Default is 1.}
}
\value{
    The function returns an object of class \code{\link{proposal-class}}, to be used in calls to \code{\link{adaptiveMH}} and
\code{\link{pawl}}.
}
\author{
    Luke Bornn <bornn@stat.harvard.edu>,
    Pierre E. Jacob <pierre.jacob.work@gmail.com>
}
%\note{
%%%  ~~further notes~~
%}
%
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{proposal-class}}, \code{\link{adaptiveMH}},
\code{\link{pawl}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.



