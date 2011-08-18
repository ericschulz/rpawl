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

# Metropolis-Hastings transition kernel, with random walk
MHkernel <- function(currentChains, currentLogTarget, AP, target, currentsigma, proposalcovmatrix){
    #print("mhkernel start")
    proposals <- currentChains + currentsigma * fastrmvnorm(AP@nchains, 
                            mu = rep(0, target@dimension), sigma = proposalcovmatrix)
    #print("beforelogdens")
#    print("proposals:")
#    print(proposals)
#    print("weirdproposals:")
#    print(sum(is.infinite(proposals)))
#    print(sum(is.na(proposals)))
    proposalLogTarget <- target@logdensity(proposals, target@parameters)
#    print(proposalLogTarget)
#    print("afterlogdens")
#    print("a")
    loguniforms <- log(runif(AP@nchains))
#    print(loguniforms)
    accepts <- (loguniforms < (proposalLogTarget - currentLogTarget))
#    print(accepts)
#    print("aaa")
    currentChains[accepts,] <- proposals[accepts,]
    currentLogTarget[accepts] <- proposalLogTarget[accepts]
#    print("aaaaa")
#    print(currentChains)
    res <- list(newchains = currentChains, newlogtarget = currentLogTarget, accepts = accepts)
#    print("aaaaaaaa")
    #print(res)
    return(res)
}

## Adaptive Metropolis-Hastings 
## AP stands for Algorithmic Parameters
adaptiveMH <- function(target, AP){
    print("Launching Adaptive Metropolis-Hastings algorithm with parameters:")
    print(AP)
    acceptrates <- c()
    chains <- as.matrix(target@rinit(AP@nchains), ncol = target@dimension)
    alllogtarget <- matrix(NA, nrow = (AP@niterations + 1), ncol = AP@nchains)
    allchains <- array(NA, dim = c(AP@niterations + 1, AP@nchains, target@dimension))
    allchains[1,,] <- chains
    chainsize <- AP@nchains
    currentlogtarget <- target@logdensity(chains, target@parameters)
    alllogtarget[1,] <- currentlogtarget
    sigma <- rep(0,AP@niterations + 1)
    sigma[1] <- AP@sigma_init
    proposalcovmatrix <- 1 / target@dimension * diag((target@dimension))
    for (iteration in 1:AP@niterations){
        ## at each time...
        #cat("iteration:", iteration, "\n")
        if (!(iteration %% 100)){
            cat("Iteration", iteration, "\n")
        }
        ## sample new chains from the MH kernel ...
        #print("beforemh")
        mhresults <- MHkernel(chains, currentlogtarget, AP, target, 
                              sigma[iteration], proposalcovmatrix)
        #print("aftermh")
        chains <- mhresults$newchains
        currentlogtarget <- mhresults$newlogtarget
        acceptrates <- c(acceptrates, mean(mhresults$accepts))
        if (AP@adaptiveproposal){
            sigma[iteration + 1] <- max(10^(-10 - target@dimension), sigma[iteration] + 
                AP@adaptationrate(iteration) * (2 * (mean(mhresults$accepts) > 0.234) - 1))
        } else {
            sigma[iteration + 1] <- sigma[iteration]
        }
        allchains[iteration + 1,,] <- chains
        alllogtarget[iteration + 1,] <- currentlogtarget
    }
    return(list(acceptrates = acceptrates, allchains = allchains, alllogtarget = alllogtarget, 
                sigma = sigma))
}

