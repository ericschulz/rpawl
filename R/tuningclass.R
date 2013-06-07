###################################################
#    This file is part of RPAWL.
#
#    RPAWL is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    RPAWL is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with RPAWL.  If not, see <http://www.gnu.org/licenses/>.
###################################################
#'MCMC Tuning Parameters
#'
#'This class holds tuning parameters for the Metropolis-Hastings and
#'Wang-Landau algorithms.
#'
#'
#'@name tuningparameters
#'@aliases tuningparameters-class tuningparameters tuningparameters,ANY-method
#'show,tuningparameters-method
#'@docType class
#'@section Objects from the Class: Objects can be created by calls of the
#'function \code{"tuningparameters"}.
#'@section Slots:
#'  \describe{
#'    \item{\code{nchains}:}{Object of class \code{"numeric"}: it should be an integer
#'    representing the desired number of parallel chains.}
#'    \item{\code{niterations}:}{Object of class \code{"numeric"}: it should be an integer
#'    representing the desired number of iterations.}
#'    \item{\code{computemean}:}{Object of class \code{"logical"}: specifies whether the mean of all chains should
#'        be computed at each iteration (useful if the chains are not to be stored).}
#'    \item{\code{computemeanburnin}:}{Object of class \code{"numeric"}: 
#'        if \code{computemean} is set to TRUE, specifies after which iteration the mean of the chain has
#'        to be computed. Default is 0 (no burnin).}
#'    \item{\code{saveeverynth}:}{Object of class \code{"numeric"}: specifies when the chains are to be stored:
#'    for instance at every iteration (=1), every 10th iteration (=10), etc. Default is -1, meaning the chains are not stored.}
#'  }
#'
#'@section Methods:
#'  \describe{
#'    \item{show}{\code{signature(object = "tuningparameters")}: provides a little summary of a binning object when called (or when \code{print} is called).}
#'   }
#'

#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@seealso \code{\link{adaptiveMH}} \code{\link{preexplorationAMH}}
#'\code{\link{pawl}}
#'@keywords classes
#'@examples
#'
#'showClass("tuningparameters")
#'mhparameters <- tuningparameters(nchains = 10, niterations = 1000, adaptiveproposal = TRUE) 
#'
NULL
setClass("tuningparameters",
         representation(nchains = "numeric", niterations = "numeric", 
                        computemean = "logical", computemeanburnin = "numeric",
                        saveeverynth = "numeric"))

setGeneric("tuningparameters", function(...)standardGeneric("tuningparameters"))
tuningparameters.constructor <- function(..., nchains, niterations, storeall, 
                                         computemean, computemeanburnin, saveeverynth){
    if (missing(nchains))
        nchains <- 1
    if (missing(niterations))
        stop(sQuote("niterations"), "has to be specified")
    if (missing(saveeverynth)){
        if (missing(storeall)){
            #cat("storeall unspecified: set to FALSE\n")
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
      #cat("computemean unspecified: set to FALSE\n")
    }
    if (missing(computemeanburnin)){
        computemeanburnin <- 0
    }
    new("tuningparameters", 
        nchains = nchains, niterations = niterations, 
        computemean = computemean, computemeanburnin = computemeanburnin, 
        saveeverynth = saveeverynth)
}
setMethod("tuningparameters",
          definition = function(..., nchains, niterations, storeall, computemean, 
                                computemeanburnin, saveeverynth){
              tuningparameters.constructor(
                               nchains = nchains, niterations = niterations, 
                               storeall = storeall, computemean = computemean,
                               computemeanburnin = computemeanburnin,
                               saveeverynth = saveeverynth)
          })

setMethod(f = "show", signature = "tuningparameters", 
          def = function(object){
              cat("Object of class ", class(object), ".\n", sep = "")
              cat("*number of parallel chains:", object@nchains, "\n")
              cat("*number of iterations:", object@niterations, "\n")
              cat("*compute mean:", object@computemean, "\n")
              cat("*compute mean (burnin):", object@computemeanburnin, "\n")
              cat("*save every nth iteration:", object@saveeverynth, "\n")
          })


