\name{adaptiveMH}
\alias{adaptiveMH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Adaptive Metropolis-Hastings
}
\description{
    Adaptive Metropolis-Hastings algorithm, with parallel chains. The adaptation is such that
    it targets an acceptance rate.
}
\usage{
adaptiveMH(target, AP, proposal, verbose)
}
\arguments{
    \item{target}{Object of class \code{\link{target}}: specifies the target distribution.
        See the help of \code{\link{target}}. If the target is discrete, target must contain the slots \code{dproposal}, \code{rproposal} and \code{proposalparam} that specify the proposal kernel in the Metropolis-Hastings step. Otherwise the default is an adaptive gaussian random walk.
    }
    \item{AP}{Object of class \code{"tuningparameters"}: specifies the number of chains, the number
        of iterations, and what should be stored during along the run. See the help of \code{\link{tuningparameters}}.
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
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{
    The function returns a list holding various information:
    \item{finalchains}{The last point of each chain.}
    \item{acceptrates}{The vector of acceptance rates at each step.}
    \item{sigma}{The vector of the standard deviations used by the MH kernel along the iterations. If the proposal was adaptive, this allows to check how the adaptation behaved.}
    \item{allchains}{If asked in the tuning parameters, the chain history.}
    \item{alllogtarget}{If asked in the tuning parameters, the associated log density evaluations.}
    \item{meanchains}{If asked in the tuning parameters, the mean (component-wise) of each chain.}
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
    Luke Bornn <l.bornn@stat.ubc.ca>,
    Pierre Jacob <pierre.jacob@ensae.fr>
}
%\note{
%%%  ~~further notes~~
%}
%
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{preexplorationAMH}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
