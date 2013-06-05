\name{PlotLogTheta}
\alias{PlotLogTheta}
\title{Plot of the log theta penalties}
\usage{
  PlotLogTheta(results)
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
  draws a trace plot of the log theta penalties along the
  iterations.
}
\author{
  Luke Bornn <bornn@stat.harvard.edu>, Pierre E. Jacob
  <pierre.jacob.work@gmail.com>
}
\seealso{
  \code{\link{ggplot}}
}

