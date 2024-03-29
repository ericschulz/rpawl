\docType{package}
\name{PAWL-package}
\alias{PAWL}
\alias{PAWL-package}
\title{PARALLEL ADAPTIVE WANG-LANDAU}
\description{
  The package implements the Parallel Adaptive Wang-Landau
  algorithm on various examples.  The provided demos allow
  to reproduce the figures of the article.
}
\details{
  \tabular{ll}{ Package: \tab PAWL\cr Type: \tab Package\cr
  Version: \tab 1.0\cr Date: \tab 2011-08-11\cr License:
  \tab GPL (>= 2)\cr LazyLoad: \tab yes\cr Depends: \tab
  mvtnorm\cr Suggests: \tab ggplot2\cr } The main function
  is \code{pawl}. It takes algorithmic parameters in
  arguments (see the help of the \code{pawl} function), as
  well a target distribution. Look at the demos to learn
  how to specify a target distribution.
}
\examples{
demo(discreteexample)
 demo(gaussianexample)
 demo(mixture2kexample)
}
\author{
  Luke Bornn <bornn@stat.harvard.edu>, Pierre E. Jacob
  <pierre.jacob.work@gmail.com>
}
\keyword{package}

