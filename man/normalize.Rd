\name{normalizeweight}
\alias{normalizeweight}
\title{
    Normalize weights 
}
\description{
    Simple function that normalize vectors (ie takes log weights and returns normalized weights,
    in the SMC context).}
\usage{
    normalizeweight(log_weights)
}
\arguments{
    \item{log_weights}{
        Object of class \code{"numeric"}: a real-valued vector}.
    }
\details{
    Simple function that takes log weights (ie any real-valued vector), computes
    the exponential of it, divides it by its sum and returns it.}
\value{
    The function returns an object of class \code{"data.frame"}, with columns for the chain indices, the chain values,
    the iteration indices, and the associated log density values.
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
\code{\link{smc}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line



