# remove all objects
graphics.off()
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

#######
## Adaptive Metropolis-Hastings
#######
#mhparameters <- tuningparameters(nchains = 10, niterations = 1000, storeall = TRUE)
#amhresults <- adaptiveMH(gaussiantarget, mhparameters)
#range(amhresults$alllogtarget)
## check that it's working
#PlotHist(results = amhresults, component = 1)
#curve(dnorm(x, mean = gaussiantarget@parameters$mean,
#            sd = gaussiantarget@parameters$sd), add = TRUE, lwd = 2)


######
# Parallel Adaptive Wang-Landau
######
N <- 20
T <- 5000
preexpresults <- preexplorationAMH(gaussiantarget, N, 1000)
binrange <- preexpresults$SuggestedRange
rinit <- function(size)
  preexpresults$finalchains
gaussiantarget@rinit <- rinit
# first create a "binning" object
# here we bin according to the (only) dimension of the 
# state space
#getPos <- function(points, logdensity) points
## we further specify some parameters, like the bins,
## the desired frequency in each bin
#ncuts <- 2
#positionbinning <- binning(position = getPos,
#                            name = "position",
#                            binrange = c(-20, -3),
#                            ncuts = ncuts,
#                            desiredfreq = c(0.1, rep(0.8/(ncuts-1), ncuts-1), 0.1),
#                            useLearningRate = FALSE,
#                            #desiredfreq = c(0.1, rep(0.2/18, 18), 0.7),
#                            autobinning = FALSE,
#                            splitThreshold = 0.17,
#                            diagnose = TRUE)


getPos <- function(points, logdensity) - logdensity 
# we further specify some parameters, like the bins,
# the desired frequency in each bin
ncuts <- 20 
energybinning <- binning(position = getPos,
                            name = "energy",
                            binrange = binrange,
                            ncuts = ncuts,
                            useLearningRate = TRUE,
                            #desiredfreq = c(0.1, rep(0.2/18, 18), 0.7),
                            autobinning = FALSE,
                            useFH = TRUE,
                            #splitThreshold = 0.17,
                            diagnose = TRUE)



proposal <- createAdaptiveRandomWalkProposal(N, gaussiantarget@dimension,
                                             adaptiveproposal = FALSE)

# get a summary of the binning
print(energybinning)
pawlparameters <- tuningparameters(nchains = N, niterations = T, storeall = TRUE)
# get a summary of the tuning parameters
print(pawlparameters)
# launching the algorithm...
#Rprof(tmp <- tempfile())
pawlresults <- pawl(gaussiantarget, binning = energybinning, AP = pawlparameters,
                    proposal = proposal)
#Rprof()
getFrequencies(pawlresults, energybinning)
print(pawlresults$bincount)
# display profiling results
#print(summaryRprof(tmp))
#unlink(tmp)
# now some plotting:
# plot the log thetas
X11();print(PlotLogTheta(pawlresults))
# trace plot of all the variables
# (here there is only one variable)
#PlotAllVar(pawlresults)
# and finally a histogram of the binned coordinate
# (here the state space)
X11()
#PlotHistBin(pawlresults, energybinning)
#getFrequencies(pawlresults, positionbinning)
chains <- ConvertResults(pawlresults)
Xnames <- grep("X", names(chains), value = TRUE)
positions <- energybinning@position(chains[,Xnames], chains$logdens)
hist(positions, nclass = 1000, 
     main = "Histogram of the binned coordinate", 
     xlab = paste("binned coordinate", sep = ""), prob = TRUE,
     col = "orange")
#
allreac <- (pawlresults$allreaction)
bins = pawlresults$binshistory[[1]]
allproportions <- tabulate(energybinning@getLocations(bins, 
                                c(allreac)), nbins = length(bins))

print(allproportions / sum(allproportions))
