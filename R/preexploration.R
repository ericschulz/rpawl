#'Pre exploration Adapative Metropolis-Hastings
#'
#'This function takes a target distribution, an integer representing the number
#'of parallel chains, and an integer representing a number of iterations, and
#'runs adaptive Metropolis-Hastings algorithm using them. The chains are then
#'used to create a range called SuggestedRange, to be used to bin the state
#'space according to the energy levels. The energy is here defined as minus the
#'log density of the target distribution.
#'
#'The adaptive Metropolis-Hastings algorithm used in the function is described
#'in more details in the help page of \code{\link{adaptiveMH}}
#'
#'@param target Object of class \code{"target"}: this argument describes the
#'target distribution.  See \code{\link{target}} for details.
#'@param nchains Object of class \code{"numeric"}: specifies the number of
#'parallel chains.
#'@param niterations Object of class \code{"numeric"}: specifies the number of
#'iterations.
#'@param proposal Object of class \code{"proposal"}: specifies the proposal
#'distribution to be used to propose new values and to compute the acceptance
#'rate. See the help of \code{\link{proposal}}. If this is not specified and
#'the target is continuous, then the default is an adaptive gaussian random
#'walk.
#'@param verbose Object of class \code{"logical"}: if TRUE (default) then
#'prints some indication of progress in the console.
#'@return The function returns a list holding the following entries:
#'\itemize{
#'\item LogEnergyRange This holds the minimum and maximum energy values
#'seen by the chains during the exploration.
#'\item LogEnergyQtile Returns the first 10\% quantile of the energy
#'values seen by the chains during the exploration.
#'\item SuggestedRange This holds the suggested range, that is, the first
#'10\% quantile and the maximum value of the energy values seen during the
#'exploration. This can be passed as the \code{binrange} argument of the
#'\code{binning} class, see the \code{trimodal} example.
#'\item finalchains The last point of each chain.
#'}
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@seealso \code{\link{adaptiveMH}}
#'@keywords ~kwd1 ~kwd2
preexplorationAMH <- function(target, nchains, niterations, proposal, verbose = TRUE){
    if (verbose) cat("Pre exploration to get log energy range\n")
    preexpparameters <- tuningparameters(nchains = nchains, niterations = niterations,
                                         adaptiveproposal = TRUE, storeall = FALSE)
    if (missing(proposal))
      preexp <- adaptiveMH(target, preexpparameters, verbose = verbose)
    else 
      preexp <- adaptiveMH(target, preexpparameters, proposal, verbose = verbose)
    burnin <- min(1000, niterations / 10)
    logtarget <- as.vector(preexp$alllogtarget[burnin:(niterations + 1),])
    LogEnergyRange <- range(-logtarget)
    LogEnergyQtile <- as.numeric(quantile(-logtarget, probs = 0.1))
    LogEnergyQtiles <- as.numeric(quantile(-logtarget, probs = c(0.1, 0.9)))
    return(list(LogEnergyRange = LogEnergyRange, 
                LogEnergyQtile = LogEnergyQtile,
                #SuggestedRange = c(LogEnergyQtile, LogEnergyRange[2]),
                SuggestedRange = c(LogEnergyQtiles[1], 
                                   LogEnergyQtiles[1] + 2 * (LogEnergyQtiles[2] - LogEnergyQtiles[1])),
                finalchains = preexp$finalchains))
}
