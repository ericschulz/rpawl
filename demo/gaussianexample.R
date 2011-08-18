# remove all objects
rm(list = ls())
# try to detach the package if it was already loaded
try(detach(package:PAWL, unload = TRUE), silent = TRUE)
# load the package
library(PAWL)
# starting points for MCMC algorithms
rinit <- function(size) rnorm(size)
# target log density function: a gaussian distribution N(mean = 2, sd = 3)
parameters <- list(mean = 2, sd = 3)
logdensity <- function(x, parameters) dnorm(x, parameters$mean, parameters$sd, log = TRUE)
# creating the target object
gaussiantarget <- target(name = "gaussian", dimension = 1,
                         rinit = rinit, logdensity = logdensity,
                         parameters = parameters)
# setting a seed for the RNG
set.seed(17)

######
# Adaptive Metropolis-Hastings
######
mhparameters <- tuningparameters(nchains = 10, niterations = 1000, adaptiveproposal = TRUE) 
amhresults <- adaptiveMH(gaussiantarget, mhparameters)
# check that it's working
PlotHist(results = amhresults, component = 1)
curve(dnorm(x, mean = gaussiantarget@parameters$mean,
            sd = gaussiantarget@parameters$sd), add = TRUE, lwd = 2)
######
# Parallel Adaptive Wang-Landau
######
N <- 10
T <- 1000
# first create a "binning" object
# here we bin according to the (only) dimension of the 
# state space
getPos <- function(points, logdensity) points
# we further specify some parameters, like the bins,
# the desired frequency in each bin
positionbinning <- binning(position = getPos,
                            name = "position",
                            binrange = c(-8, -3),
                            ncuts = 2,
                            autobinning = TRUE,
                            desiredfreq = c(0.2, 0.2, 0.6),
                            splitThreshold = 0.2)
# get a summary of the binning
print(positionbinning)
pawlparameters <- tuningparameters(nchains = N, niterations = T)
# get a summary of the tuning parameters
print(pawlparameters)
# launching the algorithm...
pawlresults <- pawl(gaussiantarget, binning = positionbinning, AP = pawlparameters)
# now some plotting:
# plot the log thetas
PlotLogTheta(pawlresults)
# trace plot of all the variables
# (here there is only one variable)
PlotAllVar(pawlresults)
# and finally a histogram of the binned coordinate
# (here the state space)
PlotHistBin(pawlresults, positionbinning)
getFrequencies(pawlresults, positionbinning)

