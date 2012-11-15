\name{PlotAllVar}
\alias{PlotAllVar}
\title{
  Trace plot of all the variables
}
\description{
    This function takes the result of \code{\link{adaptiveMH}} or of \code{\link{pawl}}, and draws 
    a trace plot for each component of the chains
}
\usage{
PlotAllVar(results)
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
  Luke Bornn <bornn@stat.harvard.edu>,
  Pierre E. Jacob <pierre.jacob.work@gmail.com>
}

\seealso{
\code{\link{ggplot}}
}


