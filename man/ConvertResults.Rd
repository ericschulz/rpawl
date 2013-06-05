\name{ConvertResults}
\alias{ConvertResults}
\title{Convert Results}
\usage{
  ConvertResults(results, verbose = TRUE)
}
\arguments{
  \item{results}{Object of class \code{"list"}: either the
  output of \code{\link{pawl}} or of
  \code{\link{adaptiveMH}}.}

  \item{verbose}{Object of class \code{"logical"}: if TRUE
  (default) then prints some indication of progress in the
  console.}
}
\value{
  The function returns an object of class
  \code{"data.frame"}, with columns for the chain indices,
  the chain values, the iteration indices, and the
  associated log density values.
}
\description{
  Convert results from \code{\link{pawl}} and
  \code{\link{adaptiveMH}}. The result is a data set that
  is more convenient to use with \code{"ggplot2"}
  functions.
}
\details{
  Essentially it concatenates the parallel chains in a
  single column, and adds a column with the associated log
  density values.  If more than 1000 parallel chains are
  used, the function can take some time to return its
  output.
}
\author{
  Luke Bornn <bornn@stat.harvard.edu>, Pierre E. Jacob
  <pierre.jacob.work@gmail.com>
}
\seealso{
  \code{\link{adaptiveMH}}, \code{\link{pawl}}
}
\keyword{~kwd1}
\keyword{~kwd2}

