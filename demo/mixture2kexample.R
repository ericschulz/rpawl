rm(list = ls())
try(detach(package:PAWL, unload = TRUE), silent = TRUE)
library(PAWL)

set.seed(17)
mixture <- createMixtureTarget()
N <- 10
T <- 2000
betaindex <- mixture@dimension
getBeta <- function(points, logdensity) exp(points[,betaindex])
betabinning <- binning(position = getBeta,
                       name = "beta",
                       binrange = c(2, 16),
                       autobinning = FALSE,
                       fhthreshold = 0.5,
                       useLearningRate = FALSE)

print(betabinning)
pawlparameters <- tuningparameters(nchains = N, niterations = T)
print(pawlparameters)

### Launching the algorithm...
pawlresults <- pawl(mixture, binning = betabinning, AP = pawlparameters)
chains <- ConvertResults(pawlresults)
getFrequencies(pawlresults, betabinning)
print(PlotHistBin(pawlresults, betabinning))

##X11()
print(PlotLogTheta(pawlresults))
plot(pawlresults$sigma, type = "l")
##X11()
#PlotAllVar(chains)
#X11()
print(PlotComp1vsComp2(chains, "X3", "X4"))

#X11()
print(PlotHistBin(pawlresults, betabinning))
#print(PlotHist(pawlresults, component = 7))
#X11(); print(PlotHist(pawlresults, component = 3))

#PlotFH(pawlresults)
#binincrease <- data.frame(cbind(c(1, pawlresults$splitTimes, T), 
#                                c(pawlresults$nbins, max(pawlresults$nbins))))
#g <- ggplot(binincrease, aes(x = X1, y = X2)) + geom_step()
#g <- g + ylim(0, 1.5 * max(pawlresults$nbins)) + ylab("Number of bins") + xlab("iterations")
#g <- g + opts(title = "Number of bins along the iterations")
#print(g)


### Adaptive MH for comparison
mhparameters <- tuningparameters(nchains = N, niterations = T, adaptiveproposal = TRUE) 
#### launching the algorithm...
amhresults <- adaptiveMH(mixture, mhparameters)

#PlotAllVar(amhresults)
#X11()
print(PlotComp1vsComp2(amhresults, "X3", "X4"))
#X11(); print(PlotHist(amhresults, component = 3))
#
