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
T <- max(chains$iterations)
burnin <- min(1000, T / 10)
subchains <- subset(chains, iterations > burnin)
totalnpoints <- dim(subchains)[1]
subchains$index <- 1:totalnpoints
maxnumberpoints <- 50000
subchains <- subset(subchains, index > totalnpoints - maxnumberpoints)
g <- ggplot(subchains, aes(x = X1, y = X2))
g <- g + stat_bin2d() + geom_density2d()
g <- g + opts(legend.position = "none")
g <- g + xlab(expression(X[1])) + ylab(expression(X[2]))
pdf(file = "Trimodal2Ddensity.pdf")
print(g)
dev.off()
# cloud of points with colour representing the density values
g <- ggplot(data = subchains, aes(x = X1, y = X2))
g <- g + geom_point(aes(alpha = logdens, size = logdens, colour = logdens))  
g <- g + xlab(expression(X[1])) + ylab(expression(X[2]))
g <- g + opts(legend.position = "none")
ggsave(g, file = "TrimodalCloud.png")

# trace plot of log theta
pdf(file = "TrimodalLogThetasSplit.pdf")
print(PlotLogTheta(pawlresults))
dev.off()

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
g <- g + geom_hline(yintercept = truethetas, linetype = 3)
g <- g + opts(legend.position = "none")
g <- g + xlim(st[2], T)
pdf(file = "TrimodalLogThetasStable.pdf")
print(g)
dev.off()

# histogram of the energy values
Xnames <- grep("X", names(chains), value = TRUE)
positions <- data.frame(densitybinning@position(chains[,Xnames], chains$logdens))
names(positions) <- c("energy")
npoints <- dim(positions)[1]
positions$index <- 1:npoints
maxnumberpoints <- 500000
g <- ggplot(data = subset(positions, index > npoints - maxnumberpoints), aes(x = energy))
g <- g + geom_histogram(binwidth = 0.025, aes(y = ..density..))
g <- g + geom_vline(xintercept = densitybinning@bins, size = 2)
g <- g + geom_vline(xintercept = pawlresults$finalbins, linetype = 2, size = 2)
g <- g + xlim(0, 15)
pdf(file = "TrimodalHistogramBins.pdf", width = 21, height = 7)
print(g)
dev.off()


#### And now the adaptive MCMC with the same number of target density
#### evaluations. Since we don't use a preliminary exploration here,
#### the chains are run for N + Nprelim iterations.
##
#mhparameters <- tuningparameters(nchains = N, niterations = T + Tprelim,
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
