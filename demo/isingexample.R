# remove all objects
rm(list = ls())
# setting a seed for the RNG
set.seed(17)
# try to detach the package if it was already loaded
try(detach(package:PAWL, unload = TRUE), silent = TRUE)
# load the package
library(PAWL)
library(fields)
data(icefloe)

imgsize <- dim(IceFloe)[1]
targetdimension <- imgsize * imgsize
targetparameters <- list(a = 1, b = 0.8, imgmatrix = IceFloe, targetdimension = targetdimension,
                         imagesize = imgsize)

##### initialization 
rInitDistribution <- function(size){
  matrix(sample(0:1,  targetdimension * size, replace = TRUE), nrow = size)
#  initpoints <- matrix(nrow = size, ncol = targetdimension)
#  for(i in 1:size){
#    initpoints[i,] <- c(IceFloe)
#  }
#  return(initpoints)
}

logdensity <- function(X, parameters){
    out <- .Call("IsingCounter", chains = X, dataimg = parameters$imgmatrix,
                 imagesize = parameters$imagesize)
    countData <- out$countData
    countSimilar <- out$countSimilarities
    return(parameters$a * countData + parameters$b * (countSimilar / 2))
}
logdensityupdate <- function(chains, parameters, updateparam){
    res <- .Call("IsingUpdate", chains = chains, dataimg = parameters$imgmatrix,
                imagesize = parameters$imagesize, flippedindex = updateparam)
    updateLikelihood <- 2*(1/2 - res$equalToData) * parameters$a
    updatePrior <- res$SimilarityChange * parameters$b / 2
    return(updateLikelihood + updatePrior)
}
isingtarget <- target(name = "ising", type = "discrete", dimension = targetdimension,
                         rinit = rInitDistribution, logdensity = logdensity,
                         parameters = targetparameters, logdensityupdate = logdensityupdate)
print(isingtarget)

proposalparam <- list(imagesize = imgsize, targetdimension = targetdimension)

rproposal <- function(states, proposalparam){
    nchains <- dim(states)[1]
    index_to_flip <- sample.int(n = proposalparam$targetdimension, 
                                size = nchains, replace = TRUE)
    states[nchains * (index_to_flip - 1) + 1:nchains] <- !(states[nchains *
                                                           (index_to_flip - 1) + 1:nchains])
  return(list(states = states, others = index_to_flip))
}
# function to compute the density of the proposal kernel
# (necessary to compute the acceptance rate)
dproposal <- function(states, ys, proposalparam){
  return(rep(0, dim(states)[1]))
}

proposalinstance <- proposal(rproposal = rproposal, 
                             dproposal = dproposal,
                             proposalparam = proposalparam)

######
# Adaptive Metropolis-Hastings
######
#mhparameters <- tuningparameters(nchains = 10, niterations = 2500, 
#                                 storeall = FALSE)
mhparameters <- tuningparameters(nchains = 100, niterations = 25000, 
                                 saveeverynth = 1000, computemean = TRUE) 
print(mhparameters)

Rprof(tmp <- tempfile())
amhresults <- adaptiveMH(isingtarget, mhparameters, proposalinstance)
Rprof()
# display profiling results
print(summaryRprof(tmp))
unlink(tmp)
#
#meanchains <- matrix(apply(amhresults$finalchains, 2, mean), ncol = imgsize)
#image.plot(1:40, 1:40, 1 - meanchains, zlim=c(0,1), 
#           col=gray((64:1)^2 / (64)^2), xlab=expression(X[1]), ylab=expression(X[2]))

meanchains <- matrix(apply(amhresults$meanchains, 2, mean), ncol = imgsize)
X11()
image.plot(1:40, 1:40, 1 - meanchains, zlim=c(0,1), 
           col=gray((64:1)^2 / (64)^2), xlab=expression(X[1]), ylab=expression(X[2]))

#meanchains <- apply(amhresults$allchains, c(2, 3), mean)
#meanchains <- matrix(apply(meanchains, 2, mean), ncol = imgsize)
#X11()
#image.plot(1:40, 1:40, 1 - meanchains, zlim=c(0,1), 
#           col=gray((64:1)^2 / (64)^2), xlab=expression(X[1]), ylab=expression(X[2]))

getLogEnergy <- function(points, logdensity) -logdensity
densitybinning <- binning(position = getLogEnergy,
                          name = "minus log target density",
                          binrange = c(-5900, -5600),
                          ncuts = 9,
                          autobinning = FALSE)
#
#
Rprof(tmp <- tempfile())
pawlresults <- pawl(isingtarget, densitybinning, mhparameters, proposalinstance)
Rprof()
# display profiling results
print(summaryRprof(tmp))
unlink(tmp)
#chains <- ConvertResults(pawlresults)
#hist(subset(chains, iterations > 50)$logdens)
#PlotHistBin(pawlresults, densitybinning)
getFrequencies(pawlresults, densitybinning)
X11();print(PlotLogTheta(pawlresults))
meanchains <- apply(pawlresults$allchains, c(2, 3), mean)
meanchains <- pawlresults$meanchains
meanchains <- matrix(apply(meanchains, 2, mean), ncol = imgsize)
X11()
image.plot(1:40, 1:40, 1 - meanchains, zlim=c(0,1), 
           col=gray((64:1)^2 / (64)^2), xlab=expression(X[1]), ylab=expression(X[2]))






