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
## Adaptive Metropolis-Hastings 
## AP stands for Algorithmic Parameters


#'Adaptive Metropolis-Hastings
#'
#'Adaptive Metropolis-Hastings algorithm, with parallel chains. The adaptation
#'is such that it targets an acceptance rate.
#'
#'
#'@param target Object of class \code{\link{target}}: specifies the target
#'distribution.  See the help of \code{\link{target}}. If the target is
#'discrete, target must contain the slots \code{dproposal}, \code{rproposal}
#'and \code{proposalparam} that specify the proposal kernel in the
#'Metropolis-Hastings step. Otherwise the default is an adaptive gaussian
#'random walk.
#'@param AP Object of class \code{"tuningparameters"}: specifies the number of
#'chains, the number of iterations, and what should be stored during along the
#'run. See the help of \code{\link{tuningparameters}}.
#'@param proposal Object of class \code{"proposal"}: specifies the proposal
#'distribution to be used to propose new values and to compute the acceptance
#'rate. See the help of \code{\link{proposal}}. If this is not specified and
#'the target is continuous, then the default is an adaptive gaussian random
#'walk.
#'@param verbose Object of class \code{"logical"}: if TRUE (default) then
#'prints some indication of progress in the console.
#'@return The function returns a list holding various information:
#'\itemize{
#'\item finalchains The last point of each chain.
#'\item acceptrates The vector of acceptance rates at each step.
#'\item sigma The vector of the standard deviations used by the MH kernel
#'along the iterations. If the proposal was adaptive, this allows to check how
#'the adaptation behaved.
#'\item allchains If asked in the tuning parameters, the chain history.
#'\item alllogtarget If asked in the tuning parameters, the associated
#'log density evaluations.
#'\item meanchains If asked in the tuning parameters, the mean
#'(component-wise) of each chain.
#'}
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@seealso \code{\link{preexplorationAMH}}
#'@keywords ~kwd1 ~kwd2
adaptiveMH <- function(target, AP, proposal, verbose = TRUE){
    if (verbose) print("Launching Adaptive Metropolis-Hastings algorithm with parameters:")
    if (verbose) print(AP)
    acceptrates <- rep(0, AP@niterations)
    chains <- as.matrix(target@rinit(AP@nchains), ncol = target@dimension)
    currentlogtarget <- target@logdensity(chains, target@parameters)
    if (AP@saveeverynth > 0){
      nallchains <- 1 + floor((AP@niterations) / AP@saveeverynth)
      if (verbose) cat("saving every", AP@saveeverynth, "iterations\n")
      totalsize <- nallchains * AP@nchains * target@dimension
      if (verbose) cat("hence saving a vector of size", nallchains, "x", AP@nchains, "x", target@dimension, "=", totalsize, "\n")
      if (totalsize > 10^8){
        if (verbose) cat("which bigger than 10^8: you better have a lot of memory available!!\n")
        suggestedmaxnallchains <- floor((10^8) / (target@dimension * AP@nchains)) + 1
        if (verbose) cat("you can maybe set saveeverynth to something bigger than ", 
            floor(AP@niterations / suggestedmaxnallchains) + 1, "\n")
        if (verbose) cat("type anything to continue, or Ctrl-C to abort\n")
        y<-scan(n=1)
      }
      allchains <- array(NA, dim = c(nallchains, AP@nchains, target@dimension))
      nstoredchains <- 1
      allchains[nstoredchains,,] <- chains
    }
    alllogtarget <- matrix(NA, nrow = AP@niterations + 1, ncol = AP@nchains)
    alllogtarget[1,] <- currentlogtarget
    if (AP@computemean){
        sumchains <- matrix(0, nrow = AP@nchains, ncol = target@dimension)
    }

    # Setting the proposal distribution for the MH kernel
    if (missing(proposal) & target@type == "continuous"){


#'Class \code{"proposal"}
#'
#'This class holds a proposal distribution to be used in a Metropolis-Hastings
#'kernel.
#'
#'
#'@name proposal
#'@aliases proposal-class proposal proposal,ANY-method show,proposal-method
#'@docType class
#'@section Objects from the Class: Objects should created by calls of the
#'function \code{proposal}.
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@keywords classes
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
        if (verbose) cat("switching off adaptive proposal, because the target is discrete\n")
    }
    iterstep <- max(100, AP@niterations / 50)
    for (iteration in 1:AP@niterations){
        ## at each time...
        #cat("iteration:", iteration, "\n")
        if (!(iteration %% iterstep)){
            if (verbose) cat("Iteration", iteration, "/", AP@niterations, "\n")
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
        }
        alllogtarget[iteration + 1,] <- currentlogtarget
        if (AP@computemean && iteration > AP@computemeanburnin){
            sumchains <- sumchains + chains
        }
    }
    results <- list(acceptrates = acceptrates,
                    finalchains = chains)
    if (proposal@adaptiveproposal)
        results$sigma <- sigma
    if (AP@saveeverynth > 0){
        results$allchains <- allchains
    }
    results$alllogtarget <- alllogtarget
    if (AP@computemean)
        results$meanchains <- sumchains / (AP@niterations - AP@computemeanburnin)
    return(results)
}




