###############
## Metropolis-Hastings transition kernel targeting the biased distribution
MHkernelPawl <- function(currentChains, currentLogTarget, 
                         currentLocations, currentReaction, logTheta, 
                         nchains, binning, target, rproposal, dproposal, proposalparam){
    proposals <- rproposal(currentChains, proposalparam)
    #proposals <- currentChains + currentsigma * fastrmvnorm(nchains, 
    #                                                        mu = rep(0, target@dimension),
    #                                                        sigma = proposalcovmatrix)
    proposalLogTarget <- target@logdensity(proposals, target@parameters)
    proposalReaction <- binning@position(proposals, proposalLogTarget)
    proposalLocations <- binning@getLocations(binning@bins, proposalReaction)
    loguniforms <- log(runif(nchains))
    accepts <- (loguniforms < ((proposalLogTarget + dproposal(proposals, currentChains, proposalparam)
                                - logTheta[proposalLocations]) 
                               - (currentLogTarget + dproposal(currentChains, proposals, proposalparam)
                                  - logTheta[currentLocations ])))
    currentChains[accepts,] <- proposals[accepts,]
    currentLogTarget[accepts] <- proposalLogTarget[accepts]
    currentLocations[accepts] <- proposalLocations[accepts]
    currentReaction[accepts] <- proposalReaction[accepts]
    return(list(newchains = currentChains, newlogtarget = currentLogTarget, 
                newlocations = currentLocations, newreaction = currentReaction, accepts = accepts))
}

## Function to check if the flat histogram criterion is reached
## Here c is taken to be 1 / d.
checkFlatHistogram <- function(bincount, binning){
    difference <- abs((bincount / sum(bincount)) - binning@desiredfreq)
    #threshold <- binning@fhthreshold * (min(binning@desiredfreq))
    FHreached <- all(difference < binning@fhthreshold * binning@desiredfreq) &
                                                         all(bincount > 10)
#    if (FHreached){
#        print("/*inside FH")
#        print(bincount / sum(bincount))
#        print(binning@bins)
#        print(binning@desiredfreq)
#        print((binning@fhthreshold / length(binning@bins)))
#        print("*/")
#    }
    return(FHreached)
}


## Particle Wang-Landau function with flat histogram criterion
## AP stands for Algorithmic Parameters
pawl <- function(target, binning, AP, verbose = TRUE){
    print("Launching Particle Wang-Landau algorithm ...") 
    # Init some algorithmic parameters ...
    nbins <- length(binning@bins)
    # The bias is saved in a list (and not a matrix) because
    # the dimension might change, if bins are split
    logtheta <- list()
    logtheta[[1]] <- rep(0, nbins)
    # keep track of the count of visits to each bin
    bincount <- rep(0, nbins)
    # as well as the count of visits to each bin since last FH
    tempbincount <- rep(0, nbins)
    # We compute the bins' middles in case the automatic bin split
    # mechanism is enabled. We also keep the count of visits in the bins,
    # as well as the count of visits in the half left part of the bins ...
    binning@binmids <- getBinMiddles(binning)
    innerbinleftcount <- rep(0, nbins - 2)
    # k is the temperature of the stochastic approximation tempering,
    # that makes the update of the bias less and less important. k is
    # increased when the flat histogram criterion is met.
    k <- 1
    khistory <- c(k)
    FHtimes <- c()
    # Due to the bin split mechanism, we keep track of the bins ...
    splitTimes <- c()
    binshistory <- list()
    binshistory[[1]] <- binning@bins
    nbinsvector <- c(nbins)
    # We keep track of the chains, and of the history in "allchains"...
    chains <- as.matrix(target@rinit(AP@nchains))
    allchains <- array(NA, dim = c(AP@niterations + 1, AP@nchains, target@dimension))
    allchains[1,,] <- chains
    acceptrates <- c()
    # Cov matrix and standard dev of the adaptive proposal ...
    proposalcovmatrix <- 1 / target@dimension * diag((target@dimension))
    sigma <- rep(0, AP@niterations + 1)
    sigma[1] <- AP@sigma_init
    # Setting the proposal distribution for the MH kernel
    proposalNotSpecified <- (is.null(target@rproposal(chains, target@proposalparam)))
    if (proposalNotSpecified){
        if (target@type == "continuous"){ 
            proposalparam <- list(sigma = sigma[1])
            rproposal <- function(currentstates, proposalparam){
                currentstates + proposalparam$sigma * fastrmvnorm(AP@nchains, 
                                       mu = rep(0, target@dimension), sigma = proposalcovmatrix)
            }
            dproposal <- function(currentstates, proposedstates, proposalparam) 0
        }
        else {
            stop("missing proposal for the MH kernel in the discrete case\n")
        }
    } else {
        dproposal <- target@dproposal
        rproposal <- target@rproposal
        proposalparam <- target@proposalparam
    }
    if (target@type == "discrete" & AP@adaptiveproposal){
        AP@adaptiveproposal <- FALSE
        cat("switching off adaptive proposal, because the target is discrete\n")
    }
    # We keep track of the log densities computed along the iterations ...
    alllogtarget <- matrix(NA, nrow = (AP@niterations + 1), ncol = AP@nchains)
    currentlogtarget <- target@logdensity(chains, target@parameters)
    alllogtarget[1,] <- currentlogtarget
    # We compute the locations of the chains, that is, in which 
    # bins they are.
    currentreaction <- binning@position(chains, currentlogtarget)
    currentlocations <- binning@getLocations(binning@bins, currentreaction)
    #print(currentlocations)
    thisiterationcount <- tabulate(currentlocations, nbins = nbins)
    bincount <- bincount + thisiterationcount
    tempbincount <- tempbincount + thisiterationcount
    for (iteration in 1:AP@niterations){
        #print(iteration)
        if (iteration %% 100 == 0){
            cat("Iteration", iteration, "\n")
        }
        ## Sample new values from the MH kernel targeting the biased distribution ...
        mhresults <- MHkernelPawl(chains, currentlogtarget, currentlocations, currentreaction,
                                  logtheta[[iteration]], AP@nchains, binning, target, 
                                  rproposal, dproposal, proposalparam)
        #mhresults <- MHkernelPawl(chains, sigma[iteration], currentlogtarget, 
        #                          currentlocations, currentreaction, logtheta[[iteration]], 
        #                          AP, binning, target, proposalcovmatrix)
        chains <- mhresults$newchains
        allchains[iteration + 1,,] <- chains
        acceptrates <- c(acceptrates, mean(mhresults$accepts))
        currentlogtarget <- mhresults$newlogtarget
        alllogtarget[iteration + 1,] <- currentlogtarget
        currentreaction <- mhresults$newreaction
        #print("currentlocations:")
        #print(currentlocations)
        currentlocations <- mhresults$newlocations
        ## update the proportions of visit in each bin
        thisiterationcount <- tabulate(currentlocations, nbins = nbins)
        bincount <- bincount + thisiterationcount 
        tempbincount <- tempbincount + thisiterationcount
        ## update the bias using all the chains ...
        #currentproportions <- getProportions(length(binning@bins), currentlocations)
        currentproportions <- thisiterationcount / AP@nchains
        if (binning@useLearningRate){
            if (binning@useFH){
                logtheta[[iteration + 1]] <- logtheta[[iteration]] + 
                binning@learningrate(k) * (currentproportions - binning@desiredfreq)
            } else {
                logtheta[[iteration + 1]] <- logtheta[[iteration]] + 
                binning@learningrate(iteration) * (currentproportions - binning@desiredfreq)
            }
        } else {
            logtheta[[iteration + 1]] <- logtheta[[iteration]] + 
                    (currentproportions - binning@desiredfreq)
        }
        ## update the adaptive proposal standard deviation ...
        if (AP@adaptiveproposal){
            sigma[iteration + 1] <- max(10^(-10 - target@dimension), sigma[iteration] + 
                            AP@adaptationrate(iteration) * (2 * (mean(mhresults$accepts) > 0.234) - 1))
        } else {
            sigma[iteration + 1] <- sigma[iteration]
        }
        proposalparam$sigma <- sigma[iteration + 1]
        ## update the inner distribution of the chains in each bin
        ## (proportions in the left hand side of each bin)
        if (binning@autobinning){
#            print("aaah")
#            print(currentreaction)
#            print(currentlocations)
#            print(binning@bins)
#            print(slotNames(binning))
#            print("?")
            innerbinleftcount <- innerbinleftcount + getInnerLeftCounts(binning,
                                            currentreaction, currentlocations)
            #print("tchoum")
        } 
        ## check if flat histogram is reached ...
        if (binning@useFH){
            FHreached <- checkFlatHistogram(tempbincount, binning)
            lastFHtime <- ifelse(!is.null(FHtimes),
                                 FHtimes[length(FHtimes)],
                                 0)
            enoughElapsedTime <- (iteration >= lastFHtime + binning@minSimEffort)
            if (FHreached && enoughElapsedTime){
                if (verbose) print("Flat histogram criterion met!")
                ## if flat histogram is reached, update the temperature k ...
                k <- k + 1
                khistory <- c(khistory, k)
                FHtimes <- c(FHtimes, iteration)
                # If the automatic bin split mechanism is enabled ...
                if (binning@autobinning & iteration < AP@niterations){
                    # Find which bins would benefit from a split
                    foundbins <- findBinsToSplit(binning, innerbinleftcount, tempbincount)
                    #binsToSplit <- findBinsResults$binsToSplit
                    #newcuts <- findBinsResults$newcuts
                    if (!is.null(foundbins$binsToSplit)){
                        splitresults <- binsplitter(binning, foundbins, 
                                              logtheta[[iteration + 1]], binning@desiredfreq)
                        newbins <- splitresults$newbins
                        nbins <- length(newbins)
                        binning@bins <- newbins
                        binning@binmids <- getBinMiddles(binning)
                        binning@desiredfreq <- splitresults$newdesiredfreq
                        bincount <- rep(0, nbins)
                        # Reset temperature
                        k <- 1
                        # "For the record"...
                        khistory[length(khistory)] <- 1
                        nbinsvector <- c(nbinsvector, nbins)
                        binshistory[[length(binshistory) + 1]] <- newbins
                        splitTimes <- c(splitTimes, iteration + 1)
                        logtheta[[iteration + 1]] <- splitresults$newthetas
                    }
                }
                tempbincount <- rep(0, nbins)
                innerbinleftcount <- rep(0, nbins - 2)
            }
        }
    }
    return(list(chains = chains, acceptrates = acceptrates, logtheta = logtheta,
                finallocations = currentlocations, FHtimes = FHtimes, 
                allchains = allchains, sigma = sigma, bins = binning@bins, 
                alllogtarget = alllogtarget, splitTimes = splitTimes, nbins = nbinsvector,
                binshistory = binshistory, khistory = khistory, bincount = bincount))
}

getFrequencies <- function(results, binning){
    finalfrequencies <- results$bincount / sum(results$bincount)
    innerfinalbins <- results$bins
    innerfinalbins <- innerfinalbins[2:(length(innerfinalbins))]
    innerinitbins <- results$binshistory[[1]]
    innerinitbins <- innerinitbins[2:(length(innerinitbins))]
    samplefrequencies <- c(finalfrequencies[1])
    for (index in 1:(length(innerinitbins)-1)){
        indexstart <- which(innerinitbins[index] == innerfinalbins)
        indexstop <- which(innerinitbins[index+1] == innerfinalbins)
        samplefrequencies <- c(samplefrequencies, sum(finalfrequencies[(1+indexstart):indexstop]))
    }
    samplefrequencies <- c(samplefrequencies, finalfrequencies[length(finalfrequencies)])
    cat("Do the obtained frequencies match the desired frequencies?\n")
    cat("final bins:", results$bins, "\n")
    cat("corresponding frequencies:", finalfrequencies, "\n")
    cat("initial bins:", results$binshistory[[1]], "\n")
    cat("desired frequencies: ", binning@desiredfreq, "\n")
    cat("obtained frequencies:", samplefrequencies, "\n")
    return(samplefrequencies)
}



