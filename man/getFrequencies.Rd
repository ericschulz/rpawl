\name{getFrequencies}
\alias{getFrequencies}
\title{Observed Frequencies in each bin.}
\usage{
  getFrequencies(results, binning)
}
\arguments{
  \item{results}{Object of class \code{"list"}: either the
  output of \code{\link{pawl}} or of
  \code{\link{adaptiveMH}}.}

  \item{binning}{Object of class \code{\link{binning}}: the
  binning on which the frequencies have to be computed.}
}
\value{
  The function is supposed to be used for the lines that it
  prints in the console.  However it also returns a vector
  of sampling frequencies corresponding to the initial
  bins.
}
\description{
  This function provides a convenient way to check whether
  the target frequencies have been reached.  Since new bins
  can be created during the algorithm, this function
  aggregates them in the right way so that the user can
  know if the initial bins (on which the desired
  frequencies were specified) have been visited enough.
}
\author{
  Luke Bornn <bornn@stat.harvard.edu>, Pierre E. Jacob
  <pierre.jacob.work@gmail.com>
}
\seealso{
  \code{\link{pawl}}
}

