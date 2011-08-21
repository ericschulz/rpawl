preexplorationAMH <- function(target, nchains, niterations){
    print("Pre exploration to get log energy range")
    preexpparameters <- tuningparameters(nchains = nchains, niterations = niterations,
                                         adaptiveproposal = TRUE, storeall = TRUE)
    preexp <- adaptiveMH(target, preexpparameters)
    burnin <- min(1000, niterations / 10)
    logtarget <- as.vector(preexp$alllogtarget[burnin:(niterations + 1),])
    LogEnergyRange <- range(-logtarget)
    LogEnergyQtile <- as.numeric(quantile(-logtarget, probs = 0.1))
    return(list(LogEnergyRange = LogEnergyRange, 
                LogEnergyQtile = LogEnergyQtile,
                SuggestedRange = c(LogEnergyQtile, LogEnergyRange[2])))
}
