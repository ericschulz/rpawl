setClass("target",
         representation(dimension="numeric", parameters="list",
                        type = "character", name = "character",
                        rinit = "function", generate = "function",
                        logdensity = "function", 
                        rproposal = "function",
                        dproposal = "function",
                        proposalparam = "list"))


setGeneric("target", function(...)standardGeneric("target"))
target.constructor <- function(..., dimension, parameters, name, 
                               type, rinit, logdensity, generate,
                               rproposal, dproposal, proposalparam){
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
  if (missing(logdensity))
    stop(sQuote("logdensity"), "has to be specified (up to an (additive) normalizing constant)")    
  if (missing(generate))
      generate <- function(size, parameters, ...) stop(sQuote("generate"), "is not specified")
  isproposalspecified <- !(missing(rproposal))
  if (!isproposalspecified){
      rproposal <- function(...) NULL
      dproposal <- function(...) NULL
      proposalparam <- list()
      cat("missing rproposal and dproposal: default is a gaussian random walk\n")
      if (type == "discrete")
          stop("error: missing rproposal and dproposal, while the target is discrete")
  }
  new("target", dimension = dimension, parameters = parameters, 
      name = name, type = type, rinit = rinit, logdensity = logdensity,
      generate = generate, rproposal = rproposal, dproposal = dproposal, 
      proposalparam = proposalparam)
}
setMethod("target",
          definition = function(..., dimension, parameters, name, type, 
                                rinit, logdensity, generate, rproposal, dproposal, proposalparam){
            target.constructor(dimension = dimension, 
                               parameters = parameters, name = name,
                               type = type, rinit = rinit,
                               logdensity = logdensity,
                               generate = generate, rproposal = rproposal,
                               dproposal = dproposal, proposalparam = proposalparam)
          })

setMethod(f = "show", signature = "target", 
          def = function(object){
            cat("Object of class ", class(object), ".\n", sep = "")
            cat("name:", object@name, "\n")
            cat("type:", object@type, "\n")
            cat("dimension of the state space:", object@dimension, "\n")
            cat("target parameters:", names(object@parameters), "\n")
          })