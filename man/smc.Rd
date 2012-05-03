\name{smc}
\alias{smc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Sequential Monte Carlo
}
\description{
   Sequential Monte Carlo samplers, using a sequence of tempered distributions. 
}
\usage{
smc(target, AP, verbose)
}
\arguments{
    \item{target}{Object of class \code{\link{target}}: specifies the target distribution.
        See the help of \code{\link{target}}. The target must be defined on a continuous state space.
    }
    \item{AP}{Object of class \code{"smcparameters"}: specifies the number of particles, the ESS 
        threshold, the sequence of distributions, etc. See the help of \code{\link{smcparameters}}.
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
    \item{particles}{}
    \item{weights}{}
    \item{ESSarray}{}
    \item{resamplingtimes}{}
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
    \code{\link{smcparameters}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

