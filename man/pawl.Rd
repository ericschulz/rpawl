\name{pawl}
\alias{pawl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Parallel Adaptive Wang-Landau
}
\description{
    Implements the Parallel Adaptive Wang-Landau algorithm.
}
\usage{
    pawl(target, binning, AP, proposal)
}
\arguments{
    \item{target}{Object of class \code{\link{target}}: specifies the target distribution.
        See the help of \code{\link{target}}. If the target is discrete, target must contain the slots \code{dproposal}, \code{rproposal} and \code{proposalparam} that specify the proposal kernel in the Metropolis-Hastings step. Otherwise the default is an adaptive gaussian random walk.
    }
    \item{binning}{Object of class \code{\link{binning}}, defining the initial bins used by the Wang-Landau algorithm.
        The binning object also contains some parameters specifying if the automatic binning mechanism is active or not, for instance.
    }
    \item{AP}{Object of class \code{\link{tuningparameters}}: specifies the number of chains, the number
        of iterations, and what should be stored during along the run. See the help of \code{\link{tuningparameters}}.
    }
    \item{proposal}{Object of class \code{\link{proposal}}: specifies the proposal distribution to be used to 
        propose new values and to compute the acceptance rate. See the help of \code{\link{proposal}}. If this
        is not specified and the target is continuous, then the default is an adaptive gaussian random walk.
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
    \item{logthetahistory}{If asked in the tuning parameters, all the log theta penalties.}
    and other quantities, that you can browse by calling \code{"names(results)"} where \code{"results"} is the output
    of the function.
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
    \code{\link{adaptiveMH}, \link{binning}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

