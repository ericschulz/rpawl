setClass("tuningparameters",
         representation(nchains = "numeric", niterations = "numeric", 
                        computemean = "logical",
                        saveeverynth = "numeric"))

setGeneric("tuningparameters", function(...)standardGeneric("tuningparameters"))
tuningparameters.constructor <- function(..., nchains, niterations, storeall, 
                                         computemean, saveeverynth){
    if (missing(nchains))
        nchains <- 1
    if (missing(niterations))
        stop(sQuote("niterations"), "has to be specified")
    if (missing(saveeverynth)){
        if (missing(storeall)){
            cat("storeall unspecified: set to FALSE\n")
            storeall <- FALSE
            saveeverynth <- -1
        } else {
            if (storeall){
                saveeverynth <- 1
            } else {
                saveeverynth <- -1
            }
        }
    }
    if (missing(computemean)){
      computemean <- FALSE
      cat("computemean unspecified: set to FALSE\n")
    }
    new("tuningparameters", 
        nchains = nchains, niterations = niterations, 
        computemean = computemean, saveeverynth = saveeverynth)
}
setMethod("tuningparameters",
          definition = function(..., nchains, niterations, storeall, computemean, saveeverynth){
              tuningparameters.constructor(
                               nchains = nchains, niterations = niterations, 
                               storeall = storeall, computemean = computemean,
                               saveeverynth = saveeverynth)
          })

setMethod(f = "show", signature = "tuningparameters", 
          def = function(object){
              cat("Object of class ", class(object), ".\n", sep = "")
              cat("*number of parallel chains:", object@nchains, "\n")
              cat("*number of iterations:", object@niterations, "\n")
              cat("*compute mean:", object@computemean, "\n")
              cat("*save every nth iteration:", object@saveeverynth, "\n")
          })


