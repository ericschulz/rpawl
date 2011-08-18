rm(list = ls())
try(detach(package:PAWL, unload = TRUE), silent = TRUE)
library(PAWL)

set.seed(17)
trimodal <- createTrimodalTarget()
N <- 5
Tprelim <- 1000
preexp = preexplorationAMH(target = trimodal, nchains = N, niterations = Tprelim)
print("Suggesting this energy range:")
print(preexp$SuggestedRange)

T <- 1000
getLogEnergy <- function(points, logdensity) -logdensity
densitybinning <- binning(position = getLogEnergy,
                            name = "minus log target density",
                            binrange = preexp$SuggestedRange,
                            autobinning = TRUE)

print(densitybinning)
pawlparameters <- tuningparameters(nchains = N, niterations = T)
print(pawlparameters)

### Launching the algorithm...
pawlresults <- pawl(trimodal, binning = densitybinning, AP = pawlparameters)
getFrequencies(pawlresults, densitybinning)

PlotAllVar(pawlresults)

PlotDensComp1vsComp2(pawlresults, "X1", "X2")
X11()
print(PlotComp1vsComp2(pawlresults, "X1", "X2"))

PlotLogTheta(pawlresults)

X11()
par(mfrow = c(3, 1))
PlotHist(pawlresults, 1)
PlotHist(pawlresults, 2)
PlotHistBin(pawlresults, densitybinning)
par(mfrow = c(1, 1))

PlotFH(pawlresults)


### And now the adaptive MCMC with the same number of target density
### evaluations. Since we don't use a preliminary exploration here,
### the chains are run for N + Nprelim iterations.
#
mhparameters <- tuningparameters(nchains = N, niterations = T + Tprelim, adaptiveproposal = TRUE) 
#
#### launching the algorithm...
amhresults <- adaptiveMH(trimodal, mhparameters)

PlotAllVar(amhresults)
PlotDensComp1vsComp2(amhresults, "X1", "X2")
X11()
print(PlotComp1vsComp2(amhresults, "X1", "X2"))

