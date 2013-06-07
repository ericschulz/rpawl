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
#'Class \code{"proposal"}
#'
#'This class holds a proposal distribution to be used in a Metropolis-Hastings kernel.
#'@rdname proposal
#'@name proposal
#'@aliases proposal-class proposal proposal-methods proposal-method
#'proposal,ANY-method show,proposal-method
#'@docType class
#'@section Objects from the class:
#'Objects should created by calls of the function \code{proposal}.
#'
#' @section Important slots:
#' \itemize{
#'     \item \code{rproposal} Object of class \code{"function"}
#'     \item \code{dproposal} Object of class \code{"function"}
#'}
#' @section Optional slots:
#' \itemize{
#'     \item \code{proposalparam} Object of class \code{"list"}
#'     \item \code{adaptiveproposal} Object of class \code{"logical"}
#'     \item \code{adaptationrate} Object of class \code{"function"}
#'     \item \code{sigma_init} Object of class \code{"numeric"}
#'}
#'@section Methods:
#'  \describe{
#'    \item{show}{\code{signature(object = "proposal")}: provides a little summary of a proposal object when called (or when \code{print} is called).}
#'  }
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@keywords classes
#'
NULL



#### the standard dmvnorm is slow because it checks a lot of thing (like is Sigma symmetric)
### Fast rmvnorm
fastrmvnorm <- function(n, mu, sigma = diag(length(mu))){
    ev <- eigen(sigma, symmetric = TRUE)
    retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
    retval <- sweep(retval, 2, mu, "+")
    return(retval)
}

setClass("proposal",
         representation( 
                        rproposal = "function",
                        dproposal = "function",
                        proposalparam = "list",
                        adaptiveproposal = "logical", 
                        adaptationrate = "function", 
                        sigma_init = "numeric"))

setGeneric("proposal", function(...)standardGeneric("proposal"))
proposal.constructor <- function(..., 
                               rproposal, dproposal, proposalparam,
                               adaptiveproposal, adaptationrate, sigma_init){
    if (missing(rproposal))
        stop(sQuote("rproposal"), "has to be specified.")
    if (missing(dproposal))
        stop(sQuote("dproposal"), "has to be specified.")
    if (missing(proposalparam))
        proposalparam <- list()
    if (missing(adaptiveproposal)){
        adaptiveproposal <- FALSE 
        cat("adaptiveproposal unspecified: set to FALSE\n")
    }
    if (missing(adaptationrate)){
        if (adaptiveproposal){
            adaptationrate <- function(x) x^(-0.6)
            cat("adaptationrate unspecified: set to x -> x^(-0.6)\n")
        } else {
            adaptationrate <- function(x) NULL
        }
    }
    if (missing(sigma_init)){
        sigma_init <- 1
        cat("sigma_init unspecified: set 1\n")
    }
    new("proposal", rproposal = rproposal, dproposal = dproposal, 
        proposalparam = proposalparam, adaptiveproposal = adaptiveproposal,
        adaptationrate = adaptationrate, sigma_init = sigma_init)
}
setMethod("proposal",
          definition = function(..., rproposal, dproposal, proposalparam,
                                adaptiveproposal, adaptationrate, sigma_init){
              proposal.constructor(rproposal = rproposal, dproposal = dproposal, 
                                   proposalparam = proposalparam, 
                                   adaptiveproposal = adaptiveproposal,
                                   adaptationrate = adaptationrate, sigma_init = sigma_init)
          })

setMethod(f = "show", signature = "proposal", 
          def = function(object){
            cat("Object of class ", class(object), ".\n", sep = "")
            cat("adaptive proposal:", object@adaptiveproposal, "\n")
          })



#'Adaptive Random Walk proposal distribution for MCMC algorithms
#'
#'Create the adaptive gaussian random walk proposal that is used as a default
#'in \code{\link{adaptiveMH}} and \code{\link{pawl}}, whenever the target
#'distribution is continuous.
#'
#'
#'@param nchains Object of class \code{"numeric"}: it should be an integer
#'representing the desired number of parallel chains.
#'@param targetdimension Object of class \code{"numeric"}: it should be an
#'integer representing the dimension of the target distribution.
#'@param adaptiveproposal Object of class \code{"logical"}: specifies whether
#'an adaptive proposal (Robbins-Monroe type of adaptation) should be used.
#'Default is FALSE.
#'@param adaptationrate Object of class \code{"function"}: specifies the rate
#'at which the adaptation of the proposal is performed. It should be a function
#'defined on [0, + infty[ such that it is not integrable but its square is
#'integrable, e.g. t -> 1/t for instance. The default is t -> t^-0.6.
#'@param sigma_init Object of class \code{"numeric"}: it should be a positive
#'real number specifying the standard deviation of the proposal distribution at
#'the first iteration. If the proposal is adaptive, it acts as a starting point
#'for the adaptation. If it is not adaptive, then this value is used throughout
#'all the iterations. Default is 1.
#'@return The function returns an object of class \code{\link{proposal-class}},
#'to be used in calls to \code{\link{adaptiveMH}} and \code{\link{pawl}}.
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@seealso \code{\link{proposal-class}}, \code{\link{adaptiveMH}},
#'\code{\link{pawl}}
createAdaptiveRandomWalkProposal <- function(nchains, targetdimension, adaptiveproposal,
                                       adaptationrate, sigma_init){
    if (missing(adaptiveproposal)){
        adaptiveproposal <- FALSE 
        cat("adaptiveproposal unspecified: set to FALSE\n")
    }
    if (missing(adaptationrate)){
        adaptationrate <- function(x) x^(-0.6)
        cat("adaptationrate unspecified: set to x -> x^(-0.6)\n")
    }
    if (missing(sigma_init)){
        sigma_init <- 1
        cat("sigma_init unspecified: set 1\n")
    }
    proposalcovmatrix <- 1 / targetdimension * diag((targetdimension))
    proposalparam <- list(sigma = sigma_init, nchains = nchains, targetdimension = targetdimension,
                          proposalcovmatrix = proposalcovmatrix)
    rproposal <- function(currentstates, proposalparam){
        list(states = currentstates + proposalparam$sigma * fastrmvnorm(proposalparam$nchains, 
              mu = rep(0, proposalparam$targetdimension), sigma = proposalparam$proposalcovmatrix))
    }
    dproposal <- function(currentstates, proposedstates, proposalparam) 0
    proposalinstance <- proposal(rproposal = rproposal, 
                                 dproposal = dproposal, 
                                 proposalparam = proposalparam,
                                 adaptiveproposal = adaptiveproposal, 
                                 adaptationrate = adaptationrate,
                                 sigma_init = sigma_init)
    return(proposalinstance)
}


