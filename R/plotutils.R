###################################################
#    This file is part of RPAWL.
#
#    RPAWL is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    RPAWL is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with RPAWL.  If not, see <http://www.gnu.org/licenses/>.
###################################################

## convert results to a more proper format for plotting
ConvertResults <- function(results){
    if (!("allchains" %in% names(results))){
        stop("error while converting results: allchains not available\n")
    }
    allchains <- results$allchains
    alllogtarget <- results$alllogtarget
    nbiterations <- dim(allchains)[1]
    nbchains <- dim(allchains)[2]
    targetdim <- dim(allchains)[3]
    library(foreach)
    print("converting chains before plotting...")
    df <- foreach (indexchain = 1:(nbchains), .combine = rbind) %do% {
        subchains <- as.data.frame(allchains[,indexchain,])
        names(subchains) <- paste("X", 1:targetdim, sep = "")
        subchains$indexchain <- indexchain
        subchains$iterations <- 1:nbiterations
        subchains$logdens <- alllogtarget[,indexchain]
        subchains
    }
    df$indexchain <- factor(df$indexchain)
    print("...done")
    return(df)
}

## Plot of one component against another
PlotComp1vsComp2 <- function(results, comp1, comp2 ){
    if ("indexchain" %in% names(results)){
        chains <- results
    } else {
        chains <- ConvertResults(results)
    }
    T <- max(chains$iterations)
    burnin <- min(1000, T / 10)
    subchains <- subset(chains, iterations > burnin)
    maxnumberpoints <- max(5000, T / 50)
    iterstep <- max(floor(T / maxnumberpoints), 1)
    library(ggplot2)
    subchains <- subset(subchains, iterations %% iterstep == 0)
    g <- ggplot(data = subchains, aes_string(x = comp1, y = comp2))
    g <- g + geom_point(aes(alpha = logdens, size = logdens, colour = logdens))  
    g <- g + xlab(comp1) + ylab(comp2)
    g <- g + opts(title = paste(comp1, "versus", comp2))
    return(g)
}

## 2D density plot of one component against another
PlotDensComp1vsComp2 <- function(results, comp1, comp2){
    if ("indexchain" %in% names(results)){
        chains <- results
    } else {
        chains <- ConvertResults(results)
    }
    T <- max(chains$iterations)
    burnin <- min(1000, T / 10)
    subchains <- subset(chains, iterations > burnin)
    maxnumberpoints <- max(5000, T / 50)
    iterstep <- max(floor(T / maxnumberpoints), 1)
    subchains <- subset(subchains, iterations %% iterstep == 0)
    library(ggplot2)
    g <- ggplot(subchains, aes_string(x = comp1, y = comp2))
    g <- g + stat_bin2d() + geom_density2d()
    g <- g + opts(title = "2D density")
    return(g)
}

## Plot of the log theta since the last bin split
#### (THIS SHOULD BE MODIFIED TO SHOW THE BIN SPLITS LIKE
#### IN THE ARTICLE)
PlotLogTheta <- function(results){
    T <- length(results$acceptrates)
    maxnumberpoints <- max(5000, T / 50)
    iterstep <- max(floor(T / maxnumberpoints), 1)
    if (is.null(results$splitTimes)){
        lastsplit <- 1
    } else {
        lastsplit <- results$splitTimes[length(results$splitTimes)]
    }
    logtheta <- results$logtheta[lastsplit:(length(results$logtheta))]
    logtheta <- matrix(unlist(logtheta), ncol = length(logtheta[[1]]), byrow = TRUE)
    logthetaDF <- data.frame(logtheta)
    names(logthetaDF) <- paste("logtheta", seq(1, results$nbins[length(results$nbins)]))
    logthetaDF$iterations <- lastsplit:(lastsplit + dim(logthetaDF)[1] - 1)
    logthetaDF <- subset(logthetaDF, iterations %% iterstep == 0)
    library(ggplot2)
    mthetaDF <- melt(logthetaDF, id = c("iterations"))
    g <- ggplot(mthetaDF, aes(x = iterations, y = value, colour = variable))
    g <- g + geom_line()
    g <- g + opts(title = "Log theta estimators")
    return(g)
}     

## Plot of the temperature k
PlotFH <- function(results){
    T <- length(results$acceptrates)
    times <- c(1, results$FHtimes, T)
    ks <- results$khistory
    ks <- c(ks, ks[length(ks)])
    kincrease <- data.frame(cbind(times, ks))
    library(ggplot2)
    g <- ggplot(kincrease, aes(x = times, y = ks)) + geom_step()
    g <- g + xlab("iterations") + ylab("k")
    g <- g + opts(title = "Number of flat histogram criteria met along the iterations")
    return(g)
}

## Plot of the number of bins
#### WHAT THE HELL GOES WRONG HERE??
PlotNbins <- function(results){
    if (!(is.null(results$splitTimes))){
        binincrease <- data.frame(cbind(c(1, results$splitTimes, T), 
                                        c(results$nbins, results$nbins[length(results$nbins)])))
        g <- ggplot(binincrease, aes(x = X1, y = X2)) + geom_step()
        g <- g + ylim(0, 1.5 * results$nbins[length(results$nbins)])
        g <- g + ylab("Number of bins") + xlab("iterations")
        g <- g + opts(title = "Number of bins along the iterations")
        return(g)
    } else {
        return(paste("number of bins always was equal to", length(results$bins)))
    }
}

## Trace plot of all the variables
PlotAllVar <- function(results){
    if ("indexchain" %in% names(results)){
        chains <- results
    } else {
        chains <- ConvertResults(results)
    }
    T <- max(chains$iterations)
    burnin <- min(1000, T / 10)
    chains <- subset(chains, iterations > burnin)
    chains[,!(names(chains) %in% c("logdens"))]
    maxnumberpoints <- max(5000, T / 50)
    iterstep <- floor(T / maxnumberpoints) + 1
    library(ggplot2)
    meltedchains <- melt(chains[,!(names(chains) %in% c("logdens"))], id = c("iterations", "indexchain"))
    g <- ggplot(subset(meltedchains, iterations %% iterstep == 0), aes(x = iterations, y = value, colour = indexchain))
    g <- g + geom_line() + facet_wrap( ~ variable) + opts(axis.text.x = theme_text(angle = -45, size = 20))
    g <- g + opts(legend.position = "none")
    g <- g + opts(title = "Trace plot of the chains")
    return(g)
}

## Simple histogram of one component
PlotHist <- function(results, component){
    if ("indexchain" %in% names(results)){
        chains <- results
    } else {
        chains <- ConvertResults(results)
    }
    hist(chains[, component], nclass = 100, 
         main = "Histogram of the chains", 
         xlab = paste("x[", component, "]", sep = ""), prob = TRUE,
         col = "orange")
}

## Simple histogram of the binned component 
PlotHistBin <- function(results, binning){
    if ("indexchain" %in% names(results)){
        chains <- results
    } else {
        chains <- ConvertResults(results)
    }
    Xnames <- grep("X", names(chains), value = TRUE)
    positions <- binning@position(chains[,Xnames], chains$logdens)
    hist(positions, nclass = 100, 
         main = "Histogram of the binned coordinate", 
         xlab = paste("binned coordinate", sep = ""), prob = TRUE,
         col = "orange")
    abline(v = binning@bins, lwd = 2)
    abline(v = results$bins, lwd = 2, lty = 3)
}






