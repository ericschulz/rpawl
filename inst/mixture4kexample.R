rm(list = ls())
try(detach(package:PAWL, unload = TRUE), silent = TRUE)
library(PAWL)

set.seed(17)
MP <- list(
           ncomponents = 4,
           componentweights = rep(0.25, 4), 
           componentmeans = c(-3, 0, 3, 6),
           componentvariances = rep(0.55^2, 4))
mixture <- createMixtureTarget(mixturesize = 100, ncomponents = 4, mixtureparameters = MP)
print(mixture)

N <- 10
T <- 1000
betaindex <- mixture@dimension
getBeta <- function(points, logdensity) exp(points[,betaindex])
betabinning <- binning(position = getBeta,
                       name = "beta",
                       binrange = c(1, 16),
                       autobinning = TRUE,
                       fhthreshold = 0.5,
                       useLearningRate = TRUE)

print(betabinning)
pawlparameters <- tuningparameters(nchains = N, niterations = T)
print(pawlparameters)

### Launching the algorithm...
Rprof(tmp <- tempfile())
pawlresults <- pawl(mixture, binning = betabinning, AP = pawlparameters)
Rprof()
print(summaryRprof(tmp))
unlink(tmp)

getFrequencies(pawlresults, betabinning)

chains <- ConvertResults(pawlresults)
pawlresults$allchains <- NULL
pawlresults$alllogtarget <- NULL
gc()
X11(); print(PlotHistBin(chains, betabinning))
X11(); print(PlotLogTheta(pawlresults))


allchains <- subset(chains, select = c("X5", "X6", "X13", "indexchain", "iterations", "logdens"))
rm(chains); gc(); 

names(allchains) <- c("Mu1", "Mu2", "beta", "indexchain", "iterations", "logdens")
alllocations <- betabinning@getLocations(pawlresults$bins, exp(allchains$beta))
finaltheta <- exp(pawlresults$logtheta[[T]])
finaltheta <- finaltheta / sum(finaltheta)
allchains$importanceweights <- finaltheta[alllocations]
burnin <- min(1000, T / 10)
subchains <- subset(allchains, iterations > burnin)
maxnumberpoints <- max(20000, T / 50)
iterstep <- max(floor(T / maxnumberpoints), 1)
subchains <- subset(subchains, iterations %% iterstep == 0)


X11(); plot(pawlresults$sigma, type = "l")

library(ggplot2)
schematictargets <- data.frame(cbind(c(rep(-3, 3), rep(0, 3), rep(3, 3), rep(6, 3)), 
                                     c(0, 3, 6, -3, 3, 6, -3, 0, 6, -3, 0, 3)))
g <- ggplot(schematictargets, aes(x = X1, y = X2))
g <- g + geom_point(size = 20, colour = "red")
g <- g + geom_point(size = 15, colour = "white")
g <- g + geom_point(size = 10, colour = "red")
g <- g + geom_point(size = 5, colour = "white")
g <- g + xlim(-5, 8) + ylim(-5, 8)
g <- g + xlab(expression(mu[i])) + ylab(expression(mu[j]))
X11(); print(g)


g <- ggplot(subchains, aes(x = Mu1, y = Mu2))
g <- g + geom_point(aes(alpha = importanceweights))
g <- g + xlab(expression(mu[1])) + ylab(expression(mu[2]))
g <- g + opts(legend.position = "none")
g <- g + scale_alpha(to=c(0.005, 0.1))
g <- g + xlim(-5, 8) + ylim(-5, 8)
X11(); print(g)



