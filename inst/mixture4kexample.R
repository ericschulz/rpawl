library(ggplot2)
theme_update(
axis.title.x = theme_text(size=25),
axis.title.y = theme_text(size=25, angle = 90),
axis.text.x = theme_text(size=25),
axis.text.y = theme_text(size=25),
strip.text.x = theme_text(size=25),
strip.text.y = theme_text(size=25),
plot.title = theme_text(size=25),
legend.text = theme_text(size=25),
legend.title = theme_text(size=25),
strip.background = theme_rect(fill = "whitesmoke"))

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
T <- 10^5
betaindex <- mixture@dimension
getBeta <- function(points, logdensity) exp(points[,betaindex])
betabinning <- binning(position = getBeta,
                       name = "beta",
                       binrange = c(1, 16),
                       ncuts = 20,
                       autobinning = TRUE,
                       fhthreshold = 0.5,
                       useLearningRate = TRUE)

print(betabinning)
pawlparameters <- tuningparameters(nchains = N, niterations = T, storeall = TRUE)
print(pawlparameters)

### Launching the algorithm...
#Rprof(tmp <- tempfile())
pawlresults <- pawl(mixture, binning = betabinning, AP = pawlparameters)
#Rprof()
#print(summaryRprof(tmp))
#unlink(tmp)

getFrequencies(pawlresults, betabinning)
print(PlotLogTheta(pawlresults))
X11();print(PlotHistBin(pawlresults, betabinning))

chains <- ConvertResults(pawlresults)
allchains <- subset(chains, select = c("X5", "X6", "X13", "indexchain", "iterations", "logdens"))
rm(chains); gc(); 

names(allchains) <- c("Mu1", "Mu2", "beta", "indexchain", "iterations", "logdens")
alllocations <- betabinning@getLocations(pawlresults$finalbins, exp(allchains$beta))
logtheta <- pawlresults$logthetahistory[[length(pawlresults$logthetahistory)]]
finaltheta <- exp(logtheta[dim(logtheta)[1],])
finaltheta <- finaltheta / sum(finaltheta)
allchains$importanceweights <- finaltheta[alllocations]
burnin <- min(1000, T / 10)

subchains <- subset(allchains, iterations > burnin)
totalnpoints <- dim(subchains)[1]
subchains$index <- 1:totalnpoints
maxnumberpoints <- min(totalnpoints, 500000)
subchains <- subset(subchains, index %in% sample(1:totalnpoints, maxnumberpoints, replace = FALSE))
schematictargets <- data.frame(cbind(c(rep(-3, 3), rep(0, 3), rep(3, 3), rep(6, 3)), 
                                     c(0, 3, 6, -3, 3, 6, -3, 0, 6, -3, 0, 3)))
g <- ggplot(schematictargets, aes(x = X1, y = X2))
g <- g + geom_point(size = 20, colour = "red")
g <- g + geom_point(size = 15, colour = "white")
g <- g + geom_point(size = 10, colour = "red")
g <- g + geom_point(size = 5, colour = "white")
g <- g + xlim(-5, 8) + ylim(-5, 8)
g <- g + xlab(expression(mu[1])) + ylab(expression(mu[2]))
ggsave(g, file = "Mixture4kTargets.png")


g <- ggplot(subchains, aes(x = Mu1, y = Mu2))
g <- g + geom_point(aes(alpha = importanceweights))
g <- g + xlab(expression(mu[1])) + ylab(expression(mu[2]))
g <- g + opts(legend.position = "none")
g <- g + scale_alpha(to=c(0.005, 0.1))
g <- g + xlim(-5, 8) + ylim(-5, 8)
ggsave(g, file = "Mixture4kWeightedPoints.png")

g <- ggplot(subchains, aes(x = Mu1, y = Mu2))
g <- g + geom_point(alpha = 1/20)
g <- g + xlab(expression(mu[1])) + ylab(expression(mu[2]))
g <- g + opts(legend.position = "none")
g <- g + scale_alpha(to=c(0.005, 0.1))
g <- g + xlim(-5, 8) + ylim(-5, 8)
ggsave(g, file = "Mixture4kPAWLPoints.png")

mhparameters <- tuningparameters(nchains = N, niterations = T, 
                                 storeall = TRUE) 
#
#### launching the algorithm...
amhresults <- adaptiveMH(mixture, mhparameters)
chains <- ConvertResults(amhresults)
allchains <- subset(chains, select = c("X5", "X6", "X13", "indexchain", "iterations", "logdens"))
rm(chains); gc(); 

names(allchains) <- c("Mu1", "Mu2", "beta", "indexchain", "iterations", "logdens")
burnin <- min(1000, T / 10)
subchains <- subset(allchains, iterations > burnin)
totalnpoints <- dim(subchains)[1]
subchains$index <- 1:totalnpoints
maxnumberpoints <- min(totalnpoints, 500000)
subchains <- subset(subchains, index %in% sample(1:totalnpoints, maxnumberpoints, replace = FALSE))

g <- ggplot(subchains, aes(x = Mu1, y = Mu2))
g <- g + geom_point(alpha = 1/20)
g <- g + xlab(expression(mu[1])) + ylab(expression(mu[2]))
g <- g + opts(legend.position = "none")
g <- g + scale_alpha(to=c(0.005, 0.1))
g <- g + xlim(-5, 8) + ylim(-5, 8)
ggsave(g, file = "Mixture4kAMHPoints.png")


