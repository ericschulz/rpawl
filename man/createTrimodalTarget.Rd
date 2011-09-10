\name{createTrimodalTarget}
\alias{createTrimodalTarget}
\title{
    Trimodal target distribution
}
\description{
    Create the trimodal distribution as in Liang, Liu and Caroll, 2007: Stochastic approximation in Monte Carlo computation.
}
\usage{
    createTrimodalTarget()
}
%\arguments{
%}
\details{
    This distribution is a mixture of three bivariate Gaussian distributions. Their covariance matrices are such that
    an adaptive MCMC algorithm which proposal variance adapts to one of the component, will likely fail to explore the others.
}
\value{
    The function returns an object of class \code{\link{target}}, with a name, a dimension, a function giving the log density,
    a function to generate sample from the distribution, parameters of the distribution, and a function to draw init points for
    the MCMC algorithms.
}
\references{
    Liang, Liu and Caroll: Stochastic approximation in Monte Carlo computation. Published in JASA, 2007.
}
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
    \code{\link{target}}, \code{\link{createMixtureTarget}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

