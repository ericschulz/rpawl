rm(list = ls())
try(detach(package:PAWL, unload = TRUE), silent = TRUE)
library(PAWL)

### toy example: the state space is made of three states
# target distribution:
targetpdf <- c(0.5, 0.4, 0.1)
# target log probability:
parameters <- list(logpi = log(targetpdf))
logdensity <- function(x, parameters){
    parameters$logpi[x]
}
# we have to specify a proposal mechanism
# for the Metropolis-Hastings kernel
# since the (default) gaussian random walk
# is not applicable here
transitionmatrix <- t(matrix(c(0.5, 0.5, 0.0,
                               0.3, 0.5, 0.2,
                               0.4, 0.2, 0.4), ncol = 3))

proposalparam <- list(transitionmatrix = transitionmatrix, card = 3)
# function that generates proposal:
rproposal <- function(states, proposalparam){
    for (index in 1:length(states)){
        states[index] <- sample(x = 1:proposalparam$card, 
                                size = 1, prob = proposalparam$transitionmatrix[states[index],])
    }
  return(states)
}
# function to compute the density of the proposal kernel
# (necessary to compute the acceptance rate)
dproposal <- function(states, ys, proposalparam){
    for (index in 1:(length(states))){
        states[index] <- log(transitionmatrix[states[index], ys[index]])
    }
  return(states)
}
# function to draw starting points for the MCMC algorithms:
rinit <- function(size) return(rep(1, size))
# define the target
discretetarget <- target(name = "discrete toy example", dimension = 1, type = "discrete",
                         rinit = rinit, logdensity = logdensity, parameters = parameters,
                         rproposal = rproposal, dproposal = dproposal, proposalparam = proposalparam)

# specify Metropolis-Hastings tuning parameters:
mhparameters <- tuningparameters(nchains = 10, niterations = 100)
Rprof(tmp <- tempfile())
amhresults <- adaptiveMH(discretetarget, mhparameters)
Rprof()
print(summaryRprof(tmp))
unlink(tmp)
chains <- ConvertResults(amhresults)

cat("target probabilities:", targetpdf, "\n")
amhcount <- tabulate(chains$X1, nbins = 3)
cat("obtained frequencies:", amhcount / sum(amhcount), "\n")

# we bin such that states 1 and 2 are in bin 1, and state 3 is in bin 2
getPos <- function(points, logdensity) 2 - (points <= 2)
# we further specify some parameters, like the bins,
# the desired frequency in each bin...
positionbinning <- binning(position = getPos,
                            name = "position",
                            bins = c(1, 2),
                            desiredfreq = c(0.8, 0.2),
                            useLearningRate = FALSE)
pawlresults <- pawl(discretetarget, binning = positionbinning, AP = mhparameters)
pawlchains <- ConvertResults(pawlresults)
cat("desired frequencies:", positionbinning@desiredfreq, "\n")
pawlcount <- tabulate(getPos(pawlchains$X1, pawlchains$logdens), nbins = 2)
cat("obtained frequencies:", pawlcount / sum(pawlcount), "\n")
# show the trace plot of log theta:
PlotLogTheta(pawlresults)

