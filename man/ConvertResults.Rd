\name{ConvertResults}
\alias{ConvertResults}
\title{
    Convert Results
}
\description{
    Convert results from \code{\link{pawl}} and \code{\link{adaptiveMH}}. The result is a
    data set that is more convenient to use with \code{"ggplot2"} functions.}
\usage{
    ConvertResults(results, verbose)
}
\arguments{
    \item{results}{
        Object of class \code{"list"}: either the output of \code{\link{pawl}} or of \code{\link{adaptiveMH}}.
    }
    \item{verbose}{
        Object of class \code{"logical"}: if TRUE (default) then prints some indication of progress
        in the console.
    }
}
\details{
Essentially
    it concatenates the parallel chains in a single column, and adds a column
    with the associated log density values.
    If more than 1000 parallel chains are used, the function can take some time to return
    its output.
}
\value{
    The function returns an object of class \code{"data.frame"}, with columns for the chain indices, the chain values,
    the iteration indices, and the associated log density values.
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
\code{\link{adaptiveMH}}, \code{\link{pawl}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


