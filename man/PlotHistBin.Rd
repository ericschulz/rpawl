\name{PlotHistBin}
\alias{PlotHistBin}
\title{Plot a histogram of the binning coordinate}
\usage{
  PlotHistBin(results, binning)
}
\arguments{
  \item{results}{Object of class \code{"list"}: either the
  output of \code{\link{pawl}} or of
  \code{\link{adaptiveMH}}.}

  \item{binning}{Object of class \code{\link{binning}},
  defining the initial bins used by the Wang-Landau
  algorithm.}
}
\value{
  The function returns a ggplot2 object.
}
\description{
  This function takes the result of
  \code{\link{adaptiveMH}} or of \code{\link{pawl}}, and a
  \code{\link{binning}} object, and draws a histogram of
  the chains according to the binning coordinate.
}
\author{
  Luke Bornn <bornn@stat.harvard.edu>, Pierre E. Jacob
  <pierre.jacob.work@gmail.com>
}
\seealso{
  \code{\link{ggplot}}
}

