\name{PlotNbins}
\alias{PlotNbins}
\title{Plot of the increase of the number of bins along the iterations}
\usage{
  PlotNbins(results)
}
\arguments{
  \item{results}{Object of class \code{"list"}: either the
  output of \code{\link{pawl}} or of
  \code{\link{adaptiveMH}}.}
}
\value{
  The function returns a ggplot2 object.
}
\description{
  This function takes the result of \code{\link{pawl}}, and
  draws a plot of the increase of the number of bins along
  the iterations.
}
\author{
  Luke Bornn <bornn@stat.harvard.edu>, Pierre E. Jacob
  <pierre.jacob.work@gmail.com>
}
\seealso{
  \code{\link{ggplot}}
}

