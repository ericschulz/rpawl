\docType{class}
\name{proposal}
\alias{proposal}
\alias{proposal,ANY-method}
\alias{proposal-class}
\alias{proposal-method}
\alias{proposal-methods}
\alias{show,proposal-method}
\title{Class \code{"proposal"}}
\description{
  This class holds a proposal distribution to be used in a
  Metropolis-Hastings kernel.
}
\section{Objects from the class}{
  Objects should created by calls of the function
  \code{proposal}.
}

\section{Important slots}{
  \itemize{ \item \code{rproposal} Object of class
  \code{"function"} \item \code{dproposal} Object of class
  \code{"function"} }
}

\section{Optional slots}{
  \itemize{ \item \code{proposalparam} Object of class
  \code{"list"} \item \code{adaptiveproposal} Object of
  class \code{"logical"} \item \code{adaptationrate} Object
  of class \code{"function"} \item \code{sigma_init} Object
  of class \code{"numeric"} }
}

\section{Methods}{
  \describe{ \item{show}{\code{signature(object =
  "proposal")}: provides a little summary of a proposal
  object when called (or when \code{print} is called).} }
}
\author{
  Luke Bornn <bornn@stat.harvard.edu>, Pierre E. Jacob
  <pierre.jacob.work@gmail.com>
}
\keyword{classes}

