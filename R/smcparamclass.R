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

#'SMC Tuning Parameters
#'
#'This class holds parameters for the Sequential Monte Carlo sampler.
#'
#'
#'@name smcparameters
#'@aliases smcparameters-class smcparameters smcparameters,ANY-method
#'show,smcparameters-method
#'@docType class
#'@section Objects from the Class: Objects can be created by calls of the
#'function \code{"smcparameters"}.
#'@section Slots:
#' \describe{
#' \item{\code{nparticles}:}{Object of class \code{"numeric"}: an integer
#'                           representing the desired number of particles.}
#' \item{\code{temperatures}:}{Object of class \code{"numeric"}: a vector of temperatures, default being
#'                             \code{"seq(from = 0.01, to = 1, length.out = 100)"}.}
#' \item{\code{nmoves}:}{Object of class \code{"numeric"}: number of move steps to be performed after each
#'                       resampling step, default being 1.}
#' \item{\code{ESSthreshold}:}{Object of class \code{"numeric"}: resampling occurs when the Effective Sample Size
#'                             goes below \code{"ESSthreshold"} multiplied by the number of particles \code{"nparticles"}.}
#' \item{\code{movetype}:}{Object of class \code{"character"}: type of Metropolis-Hastings move step to be performed;
#'                         can be either set to \code{"independent"} or \code{"randomwalk"}, default being \code{"independent"}.}
#' \item{\code{movescale}:}{Object of class \code{"numeric"}: if \code{movetype} is set to \code{"randomwalk"}, this parameter
#'                          specifies the amount by which the estimate of the standard deviation of the target distribution is multiplied; the product
#'                          being used to propose new points in the random-walk MH step. Default is 10\%, ie a new point is proposed from a Normal
#'                          distribution, centered on the latest point, with standard deviation equal to 10\% of the standard deviation of the already-generated
#'                          chain.}
#' \item{\code{resamplingscheme}:}{Object of class \code{"character"}: type of resampling to be used; either "multinomial", "residual"
#'                                 or "systematic", the default being "systematic".}
#' }
#'@section Methods:
#'\describe{
#' \item{show}{\code{signature(object = "smcparameters")}: provides a little summary of a binning object when called (or when \code{print} is called).}
#'}

#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@seealso \code{\link{smc}}
#'@keywords classes
#'@examples
#'
#'showClass("smcparameters")
#'smcparam<- smcparameters(nparticles=5000, 
#'                        temperatures = seq(from = 0.0001, to = 1, length.out= 100),
#'                        nmoves = 5, ESSthreshold = 0.5, movetype = "randomwalk",
#'                        movescale = 0.1)
#'
NULL

setClass("smcparameters",
         representation(nparticles = "numeric",  
                        temperatures = "numeric",
                        nmoves = "numeric",
                        ESSthreshold = "numeric",
                        movetype = "character",
                        movescale = "numeric",
                        resamplingscheme = "character"))

setGeneric("smcparameters", function(...)standardGeneric("smcparameters"))
smcparameters.constructor <- function(..., nparticles, temperatures, nmoves, ESSthreshold,
                                      movetype, movescale, resamplingscheme){
    if (missing(nparticles))
        nparticles <- 1
    if (missing(temperatures))
        temperatures <- seq(from = 0.01, to = 1, length.out = 100)
    if (missing(nmoves))
        nmoves <- 1
    if (missing(ESSthreshold))
        ESSthreshold <- 0.5
    if (missing(movetype)){
        movetype <- "independent"
    }
    if (movetype != "independent" && movetype != "randomwalk"){
        stop(sQuote("movetype"), " should be either 'independent' or 'randomwalk'!")
    }
    if (missing(movescale)){
        movescale <- 0.1
    }
    if (missing(resamplingscheme))
        resamplingscheme <- "systematic"
    if (resamplingscheme != "multinomial" && 
        resamplingscheme != "residual" && resamplingscheme != "systematic"){
        stop(sQuote("resamplingscheme"), 
             " should be either 'multinomial', 'residual' or 'systematic'!")
    }
    new("smcparameters", 
        nparticles = nparticles, temperatures = temperatures, nmoves = nmoves,
        ESSthreshold = ESSthreshold, movetype = movetype, movescale = movescale, 
        resamplingscheme = resamplingscheme)
}
setMethod("smcparameters",
          definition = function(..., nparticles, temperatures, nmoves, ESSthreshold, movetype, movescale,
                                resamplingscheme){
              smcparameters.constructor(
                                        nparticles = nparticles, temperatures = temperatures,
                                        nmoves = nmoves, ESSthreshold = ESSthreshold, movetype = movetype,
                                        movescale = movescale, resamplingscheme = resamplingscheme)
          })

setMethod(f = "show", signature = "smcparameters", 
          def = function(object){
            cat("Object of class ", class(object), ".\n", sep = "")
            cat("*number of particles:", object@nparticles, "\n")
            cat("*number of distributions:", length(object@temperatures), "\n")
            cat("*temperatures: ...", tail(object@temperatures), "\n")
            cat("*number of moves:", object@nmoves, "\n")
            cat("*ESS threshold:", object@ESSthreshold, "\n")
            cat("*move type:", object@movetype, "\n")
            if (object@movetype == "randomwalk") cat("*move scale:", object@movescale, "\n")
            cat("*resampling scheme:", object@resamplingscheme, "\n")
          })


