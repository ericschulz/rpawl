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
adaptiveMH(Target, AP)
}
\arguments{
  \item{Target}{Object of class \code{"target"}: specifies the target distribution.
      See the help of \code{target}. If the target is discrete, target must contain the slots \code{dproposal}, \code{rproposal} and \code{proposalparam} that specify the proposal kernel in the Metropolis-Hastings step. Otherwise the default is an adaptive gaussian random walk.
}
  \item{AP}{Object of class \code{"tuningparameters"}: specifies the number of chains, the number
      of iterations, if the proposal should be adaptive or not. See the help of \code{tuningparameters}.
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{
    The function returns a list holding various information:
    \item{allchains}{The full chain history.}
    \item{alllogtarget}{The associated log density evaluations.}
    \item{acceptrates}{The vector of acceptance rates at each step.}
    \item{sigma}{The vector of the standard deviations used by the MH kernel along the iterations. If the proposal was adaptive, this allows to check how the adaptation behaved.}
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
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
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
