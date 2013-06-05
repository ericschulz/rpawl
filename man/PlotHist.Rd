\name{PlotHist}
\alias{PlotHist}
\title{Plot a histogram of one component of the chains}
\usage{
  PlotHist(results, component)
}
\arguments{
  \item{results}{Object of class \code{"list"}: either the
  output of \code{\link{pawl}} or of
  \code{\link{adaptiveMH}}.}

  \item{component}{Object of class \code{"numeric"}:
  specifies the index of the component to plot on the
  x-axis.}
}
\value{
  The function returns a ggplot2 object.
}
\description{
  This function takes the result of
  \code{\link{adaptiveMH}} or of \code{\link{pawl}}, and a
  component index, and draws a histogram of it.
}
\author{
  Luke Bornn <bornn@stat.harvard.edu>, Pierre E. Jacob
  <pierre.jacob.work@gmail.com>
}
\seealso{
  \code{\link{ggplot}}
}

