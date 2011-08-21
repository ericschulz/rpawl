rm(list = ls())
try(detach(package:PAWL, unload = TRUE), silent = TRUE)
library(PAWL)

set.seed(17)
trimodal <- createTrimodalTarget()
N <- 10
Tprelim <- 1000
preexp = preexplorationAMH(target = trimodal, nchains = N, niterations = Tprelim)
print("Suggesting this energy range:")
print(preexp$SuggestedRange)

T <- 2000
getLogEnergy <- function(points, logdensity) -logdensity
densitybinning <- binning(position = getLogEnergy,
                            name = "minus log target density",
                            binrange = preexp$SuggestedRange,
                            ncuts = 5,
                            useFH = TRUE,
                            autobinning = TRUE,
                            splitThreshold = 0.3)

print(densitybinning)
pawlparameters <- tuningparameters(nchains = 1000, niterations = T, storeall = TRUE)
print(pawlparameters)

### Launching the algorithm...
pawlresults <- pawl(trimodal, binning = densitybinning, AP = pawlparameters)
getFrequencies(pawlresults, densitybinning)

chains <- ConvertResults(pawlresults)

# 2D density plot of the components
PlotDensComp1vsComp2(chains, "X1", "X2")
# cloud of points with colour representing the density values
PlotComp1vsComp2(chains, "X1", "X2")

logtheta <- pawlresults$logtheta
names(logtheta) <- 1:length(logtheta)
st <- pawlresults$splitTimes
st <- c(0, st, length(logtheta))
library(foreach)
library(reshape)
allmdata <- list()
logtheta[match(paste(1:2), names(logtheta))]
matrix(unlist(logtheta[match(paste(1:2), names(logtheta))]), ncol = 6, byrow = TRUE)
df <- foreach (i= 1:(length(st)-1), .combine = rbind) %do% {
    substart <- st[i] + 1
    substop  <- st[i+1] - 1
    sublogtheta <- logtheta[match(paste(substart:substop), names(logtheta))]
    sublogtheta <- matrix(unlist(sublogtheta), ncol=length(sublogtheta[[1]]), byrow=TRUE)
    theta <- exp(sublogtheta) / apply(exp(sublogtheta), 1, sum)
    thetaDF <- data.frame(theta)
    names(thetaDF) <- paste("theta", seq(1, pawlresults$nbins[i]))
    thetaDF$iterations <- substart:substop
    mdata <- melt(thetaDF, id = c("iterations"))
    names(mdata) <- c("iterations", "estimator", "value")
    mdata
}

# trace plot of log theta around the first split
g <- ggplot(df, aes(x = iterations, y = value, colour = estimator))
g <- g + geom_line() + scale_y_log()
g <- g + geom_vline(xintercept = pawlresults$splitTimes, linetype = 1)
g <- g + xlim(0, 500)
print(g)

### We can get precise estimates
## of the true thetas, which are equal (in bin i) to:
##  psi_i / phi_i
## (renormalized)
## where psi_i is the integral of the target over bin i
## and phi_i is the desired frequency of the bin
proposedvalues <- trimodal@generate(10^6, trimodal@parameters)
logtde <- trimodal@logdensity(proposedvalues, trimodal@parameters)
proposedvalues <- cbind(proposedvalues, logtde)
BINS <- pawlresults$finalbins
locations <- densitybinning@getLocations(BINS, -proposedvalues[,"logtde"])
truethetas <- tabulate(locations)
truethetas <- truethetas / pawlresults$finaldesiredfreq
truethetas <- truethetas / sum(truethetas)
print(truethetas)

# trace plot of log theta between the last
# bin split and the final iteration
g <- ggplot(df, aes(x = iterations, y = value, colour = estimator))
g <- g + geom_line() + scale_y_log()
g <- g + geom_hline(yintercept = truethetas)
g <- g + xlim(st[2], T)
print(g)

X11()
par(mfrow = c(3, 1))
PlotHist(chains, 1)
PlotHist(pawlresults, 2)
PlotHistBin(chains, densitybinning)
par(mfrow = c(1, 1))

PlotFH(pawlresults)

Xnames <- grep("X", names(chains), value = TRUE)
positions <- data.frame(densitybinning@position(chains[,Xnames], chains$logdens))
names(positions) <- c("energy")
positions$index <- 1:(dim(positions)[1])
g <- ggplot(data = subset(positions, index < 50000), aes(x = energy))
g <- g + geom_histogram(binwidth = 0.025, aes(y = ..density..))
g <- g + geom_vline(xintercept = densitybinning@bins, linesize = 2)
g <- g + geom_vline(xintercept = pawlresults$finalbins, linetype = 2)
print(g)



hist(positions, nclass = 100, 
     main = "Histogram of the binned coordinate", 
     xlab = paste("binned coordinate", sep = ""), prob = TRUE,
     col = "orange")
abline(v = binning@bins, lwd = 2)
abline(v = results$finalbins, lwd = 2, lty = 3)

#### And now the adaptive MCMC with the same number of target density
#### evaluations. Since we don't use a preliminary exploration here,
#### the chains are run for N + Nprelim iterations.
##
#mhparameters <- tuningparameters(nchains = N, niterations = T + Tprelim, adaptiveproposal = TRUE,
#                                 storeall = TRUE) 
##
##### launching the algorithm...
#amhresults <- adaptiveMH(trimodal, mhparameters)
#
#PlotAllVar(amhresults)
#PlotDensComp1vsComp2(amhresults, "X1", "X2")
#X11()
#print(PlotComp1vsComp2(amhresults, "X1", "X2"))
#
