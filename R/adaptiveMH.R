#### the standard dmvnorm is slow because it checks a lot of thing (like is Sigma symmetric)
#### so I replace it with a faster version without checks (it's dangerous though!!!)
fastdmvnorm <- function(x, mu, Sigma){
    distval <- mahalanobis(x, center = mu, cov = Sigma)
    logdet <- sum(log(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values))
    logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
    return(exp(logretval))
}
## Fast rmvnorm
fastrmvnorm <- function(n, mu, sigma = diag(length(mu))){
    ev <- eigen(sigma, symmetric = TRUE)
    retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
    retval <- sweep(retval, 2, mu, "+")
    return(retval)
}

# Metropolis-Hastings transition kernel
MHkernel <- function(currentChains, currentLogTarget, nchains, target, 
                     rproposal, dproposal, proposalparam){
    proposals <- rproposal(currentChains, proposalparam) 
    proposalLogTarget <- target@logdensity(proposals, target@parameters)
    loguniforms <- log(runif(nchains))
    accepts <- (loguniforms < (proposalLogTarget + dproposal(proposals, currentChains, proposalparam)
                               - currentLogTarget) - dproposal(currentChains, proposals, proposalparam))
    currentChains[accepts,] <- proposals[accepts,]
    currentLogTarget[accepts] <- proposalLogTarget[accepts]
    return(list(newchains = currentChains, newlogtarget = currentLogTarget, accepts = accepts))
}

## Adaptive Metropolis-Hastings 
## AP stands for Algorithmic Parameters
adaptiveMH <- function(target, AP, proposal){
    print("Launching Adaptive Metropolis-Hastings algorithm with parameters:")
    print(AP)
    acceptrates <- c()
    chains <- as.matrix(target@rinit(AP@nchains), ncol = target@dimension)
    alllogtarget <- matrix(NA, nrow = (AP@niterations + 1), ncol = AP@nchains)
    if (AP@storeall){
      allchains <- array(NA, dim = c(AP@niterations + 1, AP@nchains, target@dimension))
      allchains[1,,] <- chains
    }
    currentlogtarget <- target@logdensity(chains, target@parameters)
    alllogtarget[1,] <- currentlogtarget
    # Setting the proposal distribution for the MH kernel
    if (missing(proposal) & target@type == "continuous"){
        proposal <- createAdaptiveRandomWalkProposal(nchains = AP@nchains, 
                                               targetdimension = target@dimension,
                                               adaptiveproposal = TRUE)
    }
    dproposal <- proposal@dproposal
    rproposal <- proposal@rproposal
    proposalparam <- proposal@proposalparam
    if (proposal@adaptiveproposal){
        sigma <- rep(0, AP@niterations + 1)
        sigma[1] <- proposal@proposalparam$sigma
    }
    if (target@type == "discrete" & proposal@adaptiveproposal){
        proposal@adaptiveproposal <- FALSE
        cat("switching off adaptive proposal, because the target is discrete\n")
    }
    for (iteration in 1:AP@niterations){
        ## at each time...
        #cat("iteration:", iteration, "\n")
        if (!(iteration %% 100)){
            cat("Iteration", iteration, "\n")
        }
        ## sample new chains from the MH kernel ...
        mhresults <- MHkernel(chains, currentlogtarget, AP@nchains, target, 
                              rproposal, dproposal, proposalparam)
        chains <- mhresults$newchains
        currentlogtarget <- mhresults$newlogtarget
        acceptrates <- c(acceptrates, mean(mhresults$accepts))
        if (proposal@adaptiveproposal){
            sigma[iteration + 1] <- max(10^(-10 - target@dimension), sigma[iteration] + 
                proposal@adaptationrate(iteration) * (2 * (mean(mhresults$accepts) > 0.234) - 1))
            proposalparam$sigma <- sigma[iteration + 1]
        } 
        if (AP@storeall){
          allchains[iteration + 1,,] <- chains
        }
        alllogtarget[iteration + 1,] <- currentlogtarget
    }
    results <- list(acceptrates = acceptrates, alllogtarget = alllogtarget,
                    finalchains = chains)
    if (proposal@adaptiveproposal)
        results$sigma <- sigma
    if (AP@storeall)
        results$allchains <- allchains
    return(results)
}




