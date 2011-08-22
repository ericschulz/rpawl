## Function to check if the flat histogram criterion is reached
## Here c is taken to be 1 / d.
checkFlatHistogram <- function(bincount, binning){
    difference <- abs((bincount / sum(bincount)) - binning@desiredfreq)
    FHreached <- all(difference < binning@fhthreshold * binning@desiredfreq) &
                                                         all(bincount > 10)
    return(FHreached)
}


## Particle Wang-Landau function with flat histogram criterion
## AP stands for Algorithmic Parameters
pawl <- function(target, binning, AP, proposal, verbose = TRUE){
    print("Launching Particle Wang-Landau algorithm ...") 
    # Init some algorithmic parameters ...
    nbins <- length(binning@bins)
    # The bias is saved in a list (and not a matrix) because
    # the dimension might change, if bins are split
    #logtheta <- list()
    #logtheta[[1]] <- rep(0, nbins)
    logtheta <- matrix(ncol = nbins, nrow = AP@niterations + 1)
    logtheta[1,] <- 0
    logthetahistory <- list()
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
    # We keep track of the log densities computed along the iterations ...
    currentlogtarget <- target@logdensity(chains, target@parameters)
    if (AP@saveeverynth > 0){
      nallchains <- 1 + floor((AP@niterations) / AP@saveeverynth)
      cat("saving every", AP@saveeverynth, "iterations\n")
      totalsize <- nallchains * AP@nchains * target@dimension
      cat("hence saving a vector of size", nallchains, "x", AP@nchains, "x", target@dimension, "=", totalsize, "\n")
      if (totalsize > 10^8){
        cat("which bigger than 10^8: you better have a lot of memory available!!\n")
        suggestedmaxnallchains <- floor((10^8) / (target@dimension * AP@nchains)) + 1
        cat("you can maybe set saveeverynth to something bigger than ", 
            floor(AP@niterations / suggestedmaxnallchains) + 1, "\n")
        cat("type anything to continue, or Ctrl-C to abort\n")
        y<-scan(n=1)
      }
      allchains <- array(NA, dim = c(nallchains, AP@nchains, target@dimension))
      nstoredchains <- 1
      allchains[nstoredchains,,] <- chains
      alllogtarget <- matrix(NA, nrow = nallchains, ncol = AP@nchains)
      alllogtarget[nstoredchains,] <- currentlogtarget
    }
    if (AP@computemean){
        sumchains <- chains
    }
    acceptrates <- rep(0, AP@niterations)
    # Setting the proposal distribution for the MH kernel
    if (missing(proposal) & target@type == "continuous"){
        proposal <- createAdaptiveRandomWalkProposal(nchains = AP@nchains, 
                                               targetdimension = target@dimension,
                                               adaptiveproposal = TRUE)
    }
    if (proposal@adaptiveproposal){
        sigma <- rep(0, AP@niterations + 1)
        sigma[1] <- proposal@proposalparam$sigma
    }
    dproposal <- proposal@dproposal
    rproposal <- proposal@rproposal
    proposalparam <- proposal@proposalparam
    # standard dev of the adaptive proposal ...
    if (proposal@adaptiveproposal){
        sigma <- rep(0, AP@niterations + 1)
        sigma[1] <- proposal@proposalparam$sigma
    }
    if (target@type == "discrete" & proposal@adaptiveproposal){
        proposal@adaptiveproposal <- FALSE
        cat("switching off adaptive proposal, because the target is discrete\n")
    }
    # We compute the locations of the chains, that is, in which 
    # bins they are.
    currentreaction <- binning@position(chains, currentlogtarget)
    currentlocations <- binning@getLocations(binning@bins, currentreaction)
    thisiterationcount <- tabulate(currentlocations, nbins = nbins)
    bincount <- bincount + thisiterationcount
    tempbincount <- tempbincount + thisiterationcount
    lastFHtime <- 0
    lastSplittime <- 0
    iterstep <- max(100, AP@niterations / 50)
    for (iteration in 1:AP@niterations){
        if (!(iteration %% iterstep)){
          cat("Iteration", iteration, "/", AP@niterations, "\n")
        }
        ## Sample new values from the MH kernel targeting the biased distribution ...
        rproposalresults <- rproposal(chains, proposalparam)
        proposals <- rproposalresults$states
        if (target@updateavailable){
            proposalLogTarget <- currentlogtarget + target@logdensityupdate(chains, 
                                              target@parameters, rproposalresults$others)
        } else {
            proposalLogTarget <- target@logdensity(proposals, target@parameters)
        }
        proposalReaction <- binning@position(proposals, proposalLogTarget)
        proposalLocations <- binning@getLocations(binning@bins, proposalReaction)
        loguniforms <- log(runif(AP@nchains))

        accepts <- (loguniforms < ((proposalLogTarget + dproposal(proposals, chains, proposalparam)
                                    - logtheta[iteration - lastSplittime,][proposalLocations]) 
        - (currentlogtarget + dproposal(chains, proposals, proposalparam)
           - logtheta[iteration - lastSplittime,][currentlocations ])))
        chains[accepts,] <- proposals[accepts,]
        currentlogtarget[accepts] <- proposalLogTarget[accepts]
        currentlocations[accepts] <- proposalLocations[accepts]
        currentreaction[accepts] <- proposalReaction[accepts]
        if (AP@saveeverynth > 0 & iteration %% AP@saveeverynth == 0){
            nstoredchains <- nstoredchains + 1
            allchains[nstoredchains,,] <- chains
            alllogtarget[nstoredchains,] <- currentlogtarget
        }
        if (AP@computemean){
            sumchains <- sumchains + chains
        }
        acceptrates[iteration] <- mean(accepts)
        ## update the proportions of visit in each bin
        thisiterationcount <- tabulate(currentlocations, nbins = nbins)
        bincount <- bincount + thisiterationcount 
        tempbincount <- tempbincount + thisiterationcount
        ## update the bias using all the chains ...
        #currentproportions <- getProportions(length(binning@bins), currentlocations)
        currentproportions <- thisiterationcount / AP@nchains
        if (binning@useLearningRate){
            if (binning@useFH){
                logtheta[iteration + 1 - lastSplittime,] <- logtheta[iteration - lastSplittime,] + 
                binning@learningrate(k) * (currentproportions - binning@desiredfreq)
            } else {
                logtheta[iteration + 1 - lastSplittime,] <- logtheta[iteration - lastSplittime,] + 
                binning@learningrate(iteration) * (currentproportions - binning@desiredfreq)
            }
        } else {
            logtheta[iteration + 1 - lastSplittime,] <- logtheta[iteration - lastSplittime,] + 
                    (currentproportions - binning@desiredfreq)
        }
        ## update the adaptive proposal standard deviation ...
        if (proposal@adaptiveproposal){
            sigma[iteration + 1] <- max(10^(-10 - target@dimension), sigma[iteration] + 
                proposal@adaptationrate(iteration) * (2 * (mean(accepts) > 0.234) - 1))
            proposalparam$sigma <- sigma[iteration + 1]
        } 
        ## update the inner distribution of the chains in each bin
        ## (proportions in the left hand side of each bin)
        if (binning@autobinning){
            innerbinleftcount <- innerbinleftcount + getInnerLeftCounts(binning,
                                            currentreaction, currentlocations)
        } 
        ## check if flat histogram is reached ...
        if (binning@useFH){
            FHreached <- checkFlatHistogram(tempbincount, binning)
            enoughElapsedTime <- (iteration >= lastFHtime + binning@minSimEffort)
            if (FHreached && enoughElapsedTime){
                if (verbose) print("Flat histogram criterion met!")
                ## if flat histogram is reached, update the temperature k ...
                k <- k + 1
                khistory <- c(khistory, k)
                FHtimes <- c(FHtimes, iteration)
                lastFHtime <- iteration
                # If the automatic bin split mechanism is enabled ...
                if (binning@autobinning & iteration < AP@niterations){
                    # Find which bins would benefit from a split
                    foundbins <- findBinsToSplit(binning, innerbinleftcount, tempbincount)
                    #binsToSplit <- findBinsResults$binsToSplit
                    #newcuts <- findBinsResults$newcuts
                    if (!is.null(foundbins$binsToSplit)){
                        cat("split at time", iteration, "\n")
                        splitresults <- binsplitter(binning, foundbins, 
                                              logtheta[iteration + 1 - lastSplittime,], binning@desiredfreq)
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
                        splitTimes <- c(splitTimes, iteration)
                        lastSplittime <- iteration
                        if (length(splitTimes) == 1){
                          logthetahistory[[1]] <- logtheta[1:(iteration + 1),]
                        } else {
                          diffSplitTimes <- splitTimes[length(splitTimes)] - splitTimes[length(splitTimes) - 1]
                          logthetahistory[[length(splitTimes)]] <- logtheta[1:(diffSplitTimes + 1),]
                        }
                        logtheta <- matrix(ncol = nbins, nrow = AP@niterations + 1 - iteration)
                        logtheta[1,] <- splitresults$newthetas
                    }
                }
                tempbincount <- rep(0, nbins)
                innerbinleftcount <- rep(0, nbins - 2)
            }
        }
    }
    results <- list(chains = chains, acceptrates = acceptrates, logtheta = logtheta,
                finallocations = currentlocations, FHtimes = FHtimes, 
                finalbins = binning@bins, finaldesiredfreq = binning@desiredfreq,
                splitTimes = splitTimes, nbins = nbinsvector,
                binshistory = binshistory, khistory = khistory, bincount = bincount)
    if (proposal@adaptiveproposal)
        results$sigma <- sigma
    if (AP@saveeverynth > 0){
        results$allchains <- allchains
        results$alllogtarget <- alllogtarget
    }
    if (AP@computemean){
        results$meanchains <- sumchains / (AP@niterations + 1)
    }
    logthetahistory[[length(logthetahistory) + 1]] <- logtheta
    results$logthetahistory <- logthetahistory
    return(results)
}

getFrequencies <- function(results, binning){
    finalfrequencies <- results$bincount / sum(results$bincount)
    innerfinalbins <- results$finalbins
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
    cat("final bins:", results$finalbins, "\n")
    cat("corresponding frequencies:", finalfrequencies, "\n")
    cat("initial bins:", results$binshistory[[1]], "\n")
    cat("desired frequencies: ", binning@desiredfreq, "\n")
    cat("obtained frequencies:", samplefrequencies, "\n")
    return(samplefrequencies)
}



