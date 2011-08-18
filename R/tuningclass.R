setClass("tuningparameters",
         representation(nchains = "numeric", niterations = "numeric", 
                        adaptiveproposal = "logical", adaptationrate = "function", 
                        sigma_init = "numeric"))

setGeneric("tuningparameters", function(...)standardGeneric("tuningparameters"))
tuningparameters.constructor <- function(..., nchains, niterations, 
                                         adaptiveproposal, adaptationrate,
                                         sigma_init){
    if (missing(nchains))
        nchains <- 1
    if (missing(niterations))
        stop(sQuote("niterations"), "has to be specified")
    if (missing(adaptiveproposal)){
        adaptiveproposal <- TRUE 
        cat("adaptiveproposal unspecified: set to TRUE\n")
    }
    if (missing(adaptationrate)){
        adaptationrate <- function(x) x^(-0.6)
        cat("adaptationrate unspecified: set to x -> x^(-0.6)\n")
    }
    if (missing(sigma_init)){
        sigma_init <- 1
        cat("sigma_init unspecified: set 1\n")
    }
    new("tuningparameters", 
        nchains = nchains, niterations = niterations, 
        adaptiveproposal = adaptiveproposal,
        adaptationrate = adaptationrate, 
        sigma_init = sigma_init)
}
setMethod("tuningparameters",
          definition = function(..., nchains, niterations, adaptiveproposal, adaptationrate, 
                                sigma_init){
              tuningparameters.constructor(
                               nchains = nchains, niterations = niterations, 
                               adaptiveproposal = adaptiveproposal,
                               adaptationrate = adaptationrate, 
                               sigma_init = sigma_init)
          })

setMethod(f = "show", signature = "tuningparameters", 
          def = function(object){
              cat("Object of class ", class(object), ".\n", sep = "")
              cat("*number of parallel chains:", object@nchains, "\n")
              cat("*number of iterations:", object@niterations, "\n")
              cat("*adaptive proposal:", object@adaptiveproposal, "\n")
          })


