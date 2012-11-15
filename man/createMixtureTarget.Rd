\name{createMixtureTarget}
\alias{createMixtureTarget}
\title{
    Mixture target distribution
}
\description{
    Create the posterior distribution of the parameters of a mixture of univariate gaussian distributions, with a fixed (known)
    number of components.
}
\usage{
    createMixtureTarget(mixturesample, mixturesize, ncomponents, mixtureparameters)
}
\arguments{
    All the arguments are optional, since if none is given, a mixture distribution with 4 components will be
    created, as in Jasra, Holmes, Stephens, "MCMC and label switching problem in Bayesian mixture models", published
    in Statistical Science (2005).
    \item{mixturesample}{Object of class \code{"vector"}: data set to be used. If not provided, a synthetic data set is generated.} 
    \item{mixturesize}{Object of class \code{"numeric"}: represents the data set size if a data set is to be generated.}
    \item{ncomponents}{Object of class \code{"numeric"}: represents the fixed number of components to be used.}
    \item{mixtureparameters}{Object of class \code{"list"}: provides the parameters to be used if a data set has to be generated.
    The parameters include the number of components, the component weights, means and variances.}
}
\value{
    The function returns an object of class \code{\link{target-class}}, with a name, a dimension, a function giving the log density,
    a function to generate sample from the distribution, parameters of the distribution, and a function to draw init points for
    the MCMC algorithms. The log density involves a likelihood and a prior, and the prior is as in Richardson and Green, "On Bayesian 
    analysis of mixtures with an unknown number of components", published in JRSS B, 1997.
}
\references{
Jasra, Holmes, Stephens, "MCMC and label switching problem in Bayesian mixture models", published
    in Statistical Science (2005).
Richardson and Green, "On Bayesian 
    analysis of mixtures with an unknown number of components", published in JRSS B, 1997.
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
    \code{\link{target-class}}, \code{\link{createTrimodalTarget}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.


