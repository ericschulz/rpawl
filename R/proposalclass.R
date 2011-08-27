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


