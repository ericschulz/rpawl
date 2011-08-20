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
#MHkernel <- function(currentChains, currentLogTarget, nchains, target, 
#                     rproposal, dproposal, proposalparam){
#    rproposalresults <- rproposal(currentChains, proposalparam)
#    proposals <- rproposalresults$states
#    if (target@updateavailable){
#        proposalLogTarget <- currentLogTarget + target@logdensityupdate(currentChains, 
#                                 target@parameters, rproposalresults$others)
#    } else {
#        proposalLogTarget <- target@logdensity(proposals, target@parameters)
#    }
#    loguniforms <- log(runif(nchains))
#    accepts <- (loguniforms < (proposalLogTarget + dproposal(proposals, currentChains, proposalparam)
#                               - currentLogTarget) - dproposal(currentChains, proposals, proposalparam))
#    currentChains[accepts,] <- proposals[accepts,]
#    currentLogTarget[accepts] <- proposalLogTarget[accepts]
#    return(list(newchains = currentChains, newlogtarget = currentLogTarget, accepts = accepts))
#}
#
## Adaptive Metropolis-Hastings 
## AP stands for Algorithmic Parameters
adaptiveMH <- function(target, AP, proposal){
    print("Launching Adaptive Metropolis-Hastings algorithm with parameters:")
    print(AP)
    acceptrates <- rep(0, AP@niterations)
    chains <- as.matrix(target@rinit(AP@nchains), ncol = target@dimension)
    currentlogtarget <- target@logdensity(chains, target@parameters)
    if (AP@saveeverynth > 0){
      nallchains <- 1 + floor((AP@niterations) / AP@saveeverynth)
      allchains <- array(NA, dim = c(nallchains, AP@nchains, target@dimension))
      alllogtarget <- matrix(NA, nrow = nallchains, ncol = AP@nchains)
      nstoredchains <- 1
      allchains[nstoredchains,,] <- chains
      alllogtarget[nstoredchains,] <- currentlogtarget
    }
    if (AP@computemean){
        sumchains <- chains
    }
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
        rproposalresults <- rproposal(chains, proposalparam)
        proposals <- rproposalresults$states
        if (target@updateavailable){
            proposalLogTarget <- currentlogtarget + target@logdensityupdate(chains, 
                                                                            target@parameters, rproposalresults$others)
        } else {
            proposalLogTarget <- target@logdensity(proposals, target@parameters)
        }
        loguniforms <- log(runif(AP@nchains))
        accepts <- (loguniforms < (proposalLogTarget + dproposal(proposals, chains, proposalparam)
                                   - currentlogtarget) - dproposal(chains, proposals, proposalparam))
        chains[accepts,] <- proposals[accepts,]
        currentlogtarget[accepts] <- proposalLogTarget[accepts]
        acceptrates[iteration] <- mean(accepts)
        if (proposal@adaptiveproposal){
            sigma[iteration + 1] <- max(10^(-10 - target@dimension), sigma[iteration] + 
                proposal@adaptationrate(iteration) * (2 * (mean(accepts) > 0.234) - 1))
            proposalparam$sigma <- sigma[iteration + 1]
        } 
        if (AP@saveeverynth > 0 & iteration %% AP@saveeverynth == 0){
            nstoredchains <- nstoredchains + 1
            allchains[nstoredchains,,] <- chains
            alllogtarget[nstoredchains,] <- currentlogtarget
        }
        if (AP@computemean){
            sumchains <- sumchains + chains
        }
    }
    results <- list(acceptrates = acceptrates,
                    finalchains = chains)
    if (proposal@adaptiveproposal)
        results$sigma <- sigma
    if (AP@saveeverynth > 0){
        results$allchains <- allchains
        results$alllogtarget <- alllogtarget
    }
    if (AP@computemean)
        results$meanchains <- sumchains / (AP@niterations + 1)
    return(results)
}




