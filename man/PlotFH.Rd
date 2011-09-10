\name{PlotFH}
\alias{PlotFH}
\title{
  Plot of the Flat Histogram occurrences
}
\description{
    This function takes the result of \code{\link{pawl}}, and draws 
    a plot of the occurrences of the Flat Histogram criteria along the iterations.
}
\usage{
PlotFH(results)
}
\arguments{
    \item{results}{
        Object of class \code{"list"}: either the output of \code{\link{pawl}} or of \code{\link{adaptiveMH}}.
    }
}
\value{
    The function returns a ggplot2 object.
}

\author{
  Luke Bornn <l.bornn@stat.ubc.ca>,
  Pierre Jacob <pierre.jacob@ensae.fr>
}

\seealso{
\code{\link{ggplot}}
}


