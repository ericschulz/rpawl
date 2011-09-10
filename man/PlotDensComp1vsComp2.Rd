\name{PlotDensComp1vsComp2}
\alias{PlotDensComp1vsComp2}

\title{
  Plot one component versus another in a density plot
}
\description{
    This function takes the result of \code{\link{adaptiveMH}} or of \code{\link{pawl}}, and
    component indices, and draws a 2D density plot with the first component on the x-axis and
    the second on the y-axis.
}
\usage{
PlotDensComp1vsComp2(results, comp1, comp2)
}
\arguments{
    \item{results}{
        Object of class \code{"list"}: either the output of \code{\link{pawl}} or of \code{\link{adaptiveMH}}.
    }
    \item{comp1}{Object of class \code{"numeric"}: specifies the index of the component to plot on the x-axis.
    }
    \item{comp2}{Object of class \code{"numeric"}: specifies the index of the component to plot on the y-axis.
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


