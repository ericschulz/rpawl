#'Class: target distribution
#'
#'This class represents target distributions, that is, probability
#'distributions from which we want to sample using MCMC or Wang-Landau.
#'
#'
#'@name target
#'@aliases target-class target target,ANY-method show,target-method
#'@docType class
#'@section Objects from the Class: Objects should created by calls of the
#'function \code{target}. Examples are provided that should help implementing
#'any continuous probability distributions.
#'@section Important slots:
#'  \describe{
#'    \item{\code{dimension}:}{Object of class \code{"numeric"}: should be an integer specifying the dimension of the state space on which the target distribution is defined.}
#'    \item{\code{logdensity}:}{Object of class \code{"function"} : should be a function taking n points in the state space and parameters, and returning a vector of n real values. See the example below. This function is in most cases the most time-consuming part in a MCMC algorithm, so make sure it runs reasonably fast!}
#'    \item{\code{rinit}:}{Object of class \code{"function"} : this function should take an integer as argument, say n. Then the function should return a matrix of dimension n times d (where d is the dimension of the state space), representing n points in the state space. These n points will be used as starting points of a parallel MCMC algorithm.}
#'  }
#'
#'@section Optional slots:
#'  \describe{
#'    \item{\code{parameters}:}{Object of class \code{"list"} : you can put anything in that list (and nothing, which is the defaults), the important thing is that calls to \code{logdensity(x, parameters)} return sensible values. For example, for a gaussian target distribution, you can put the mean and the variance in the \code{parameters} list (see example below). If need be, you can put a whole data set in there.}
#'    \item{\code{type}:}{Object of class \code{"character"} : could be "continuous" or "discrete"; default is "continuous". }
#'    \item{\code{name}:}{Object of class \code{"character"} : ... if you want to name your distribution (default is "unspecified").}
#'    \item{\code{generate}:}{Object of class \code{"function"} : does not have to be specified, but if it is specified it should be a function to generate from the distribution (like rnorm is to the standard normal distribution).}
#'  }
#'
#'@section Methods:
#'  \describe{
#'    \item{show}{\code{signature(object = "target")}: provides a little summary of a target object when called (or when \code{print} is called).}
#'   }
#'
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@keywords classes
#'@examples
#'
#'  showClass("target")
#'  # starting points for MCMC algorithms
#'  rinit <- function(size) rnorm(size)
#'  # target log density function: a gaussian distribution N(mean = 2, sd = 3)
#'  parameters <- list(mean = 2, sd = 3)
#'  logdensity <- function(x, parameters) dnorm(x, parameters$mean, parameters$sd, log = TRUE)
#'  # creating the target object
#'  gaussiantarget <- target(name = "gaussian", dimension = 1,
#'                    rinit = rinit, logdensity = logdensity,
#'                    parameters = parameters)
#'  print(gaussiantarget)
#'
NULL
setClass("target",
         representation(dimension="numeric", parameters="list",
                        type = "character", name = "character",
                        rinit = "function", dinit = "function",
                        generate = "function",
                        logdensity = "function", logdensityupdate =
                        "function", updateavailable = "logical"))
setGeneric("target", function(...)standardGeneric("target"))
target.constructor <- function(..., dimension, parameters, name,
                               type, rinit, dinit, logdensity, generate,
                               logdensityupdate){
  if (missing(dimension))
    stop(sQuote("dimension"), "has to be specified.")
  if (missing(parameters))
    parameters <- list()
  if (missing(name))
    name <- "unspecified"
  if (missing(type))
    type <- "continuous"
  if (missing(rinit))
    stop(sQuote("rinit"), "has to be specified to draw initial points of the Markov chains")
  if (missing(dinit)){
    cat("dinit has to be specified if you want to use SMC")
    dinit <- function(x) 0
  }
  if (missing(logdensity))
    stop(sQuote("logdensity"), "has to be specified (up to an (additive) normalizing constant)")    
  if (missing(generate))
      generate <- function(size, parameters, ...) stop(sQuote("generate"), "is not specified")
  if (missing(logdensityupdate)){
      logdensityupdate <- function(...) NULL
      updateavailable <- FALSE
  } else { updateavailable <- TRUE }

  new("target", dimension = dimension, parameters = parameters,
      name = name, type = type, rinit = rinit, dinit = dinit, logdensity = logdensity,
      generate = generate, logdensityupdate = logdensityupdate,
      updateavailable = updateavailable)
}
setMethod("target",
          definition = function(..., dimension, parameters, name, type,
                                rinit, dinit, logdensity, generate, logdensityupdate){
            target.constructor(dimension = dimension,
                               parameters = parameters, name = name,
                               type = type, rinit = rinit, dinit = dinit,
                               logdensity = logdensity,
                               generate = generate, logdensityupdate =
                               logdensityupdate)
          })

setMethod(f = "show", signature = "target",
          def = function(object){
            cat("Object of class ", class(object), ".\n", sep = "")
            cat("*name:", object@name, "\n")
            cat("*type:", object@type, "\n")
            cat("*dimension of the state space:", object@dimension, "\n")
            cat("*target parameters:", names(object@parameters), "\n")
            cat("*update available:", object@updateavailable, "\n")
          })



