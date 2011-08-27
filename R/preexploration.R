preexplorationAMH <- function(target, nchains, niterations, proposalinstance){
    cat("Pre exploration to get log energy range\n")
    preexpparameters <- tuningparameters(nchains = nchains, niterations = niterations,
                                         adaptiveproposal = TRUE, storeall = FALSE)
    if (missing(proposalinstance))
      preexp <- adaptiveMH(target, preexpparameters)
    else 
      preexp <- adaptiveMH(target, preexpparameters, proposalinstance)
    burnin <- min(1000, niterations / 10)
    logtarget <- as.vector(preexp$alllogtarget[burnin:(niterations + 1),])
    LogEnergyRange <- range(-logtarget)
    LogEnergyQtile <- as.numeric(quantile(-logtarget, probs = 0.1))
    return(list(LogEnergyRange = LogEnergyRange, 
                LogEnergyQtile = LogEnergyQtile,
                SuggestedRange = c(LogEnergyQtile, LogEnergyRange[2]),
                finalchains = preexp$finalchains))
}
