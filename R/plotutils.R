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


#'Convert Results
#'
#'Convert results from \code{\link{pawl}} and \code{\link{adaptiveMH}}. The
#'result is a data set that is more convenient to use with \code{"ggplot2"}
#'functions.
#'
#'Essentially it concatenates the parallel chains in a single column, and adds
#'a column with the associated log density values.  If more than 1000 parallel
#'chains are used, the function can take some time to return its output.
#'
#'@param results Object of class \code{"list"}: either the output of
#'\code{\link{pawl}} or of \code{\link{adaptiveMH}}.
#'@param verbose Object of class \code{"logical"}: if TRUE (default) then
#'prints some indication of progress in the console.
#'@return The function returns an object of class \code{"data.frame"}, with
#'columns for the chain indices, the chain values, the iteration indices, and
#'the associated log density values.
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@seealso \code{\link{adaptiveMH}}, \code{\link{pawl}}
#'@keywords ~kwd1 ~kwd2
ConvertResults <- function(results, verbose = TRUE){
    if (!("allchains" %in% names(results))){
        stop("error while converting results: allchains not available\n")
    }
    allchains <- results$allchains
    alllogtarget <- results$alllogtarget
    nbiterations <- dim(allchains)[1]
    nbchains <- dim(allchains)[2]
    targetdim <- dim(allchains)[3]
    indexchain <- 0
    library(foreach)
    if (verbose) print("converting chains before plotting...")
    df <- foreach (indexchain = 1:(nbchains), .combine = rbind) %do% {
        subchains <- as.data.frame(allchains[,indexchain,])
        names(subchains) <- paste("X", 1:targetdim, sep = "")
        subchains$indexchain <- indexchain
        subchains$iterations <- 1:nbiterations
        subchains$logdens <- alllogtarget[,indexchain]
        subchains
    }
    df$indexchain <- factor(df$indexchain)
    if (verbose) print("...done")
    return(df)
}

## Plot of one component against another


#'Plot one component versus another in a scatter plot
#'
#'This function takes the result of \code{\link{adaptiveMH}} or of
#'\code{\link{pawl}}, and component indices, and draws a cloud of points with
#'the first component on the x-axis and the second on the y-axis.
#'
#'
#'@param results Object of class \code{"list"}: either the output of
#'\code{\link{pawl}} or of \code{\link{adaptiveMH}}.
#'@param comp1 Object of class \code{"numeric"}: specifies the index of the
#'component to plot on the x-axis.
#'@param comp2 Object of class \code{"numeric"}: specifies the index of the
#'component to plot on the y-axis.
#'@return The function returns a ggplot2 object.
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@seealso \code{\link{ggplot}}
PlotComp1vsComp2 <- function(results, comp1, comp2 ){
    if ("indexchain" %in% names(results)){
        chains <- results
    } else {
        chains <- ConvertResults(results)
    }
    iterations <- 1; logdens <- 1; index <- 1
    T <- max(chains$iterations)
    burnin <- floor(min(1000, T / 10))
    subchains <- subset(chains, iterations > burnin)
    totalnpoints <- dim(subchains)[1]
    subchains$index <- 1:totalnpoints
    maxnumberpoints <- 50000
    subchains <- subset(subchains, index > totalnpoints - maxnumberpoints)
    library(ggplot2)
    g <- ggplot(data = subchains, aes_string(x = comp1, y = comp2))
    g <- g + geom_point(aes(alpha = logdens, size = logdens, colour = logdens))  
    g <- g + xlab(comp1) + ylab(comp2)
    g <- g + labs(title = paste(comp1, "versus", comp2))
    return(g)
}

## 2D density plot of one component against another


#'Plot one component versus another in a density plot
#'
#'This function takes the result of \code{\link{adaptiveMH}} or of
#'\code{\link{pawl}}, and component indices, and draws a 2D density plot with
#'the first component on the x-axis and the second on the y-axis.
#'
#'
#'@param results Object of class \code{"list"}: either the output of
#'\code{\link{pawl}} or of \code{\link{adaptiveMH}}.
#'@param comp1 Object of class \code{"numeric"}: specifies the index of the
#'component to plot on the x-axis.
#'@param comp2 Object of class \code{"numeric"}: specifies the index of the
#'component to plot on the y-axis.
#'@return The function returns a ggplot2 object.
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@seealso \code{\link{ggplot}}
PlotDensComp1vsComp2 <- function(results, comp1, comp2){
    if ("indexchain" %in% names(results)){
        chains <- results
    } else {
        chains <- ConvertResults(results)
    }
    iterations <- 1; index <- 1
    T <- max(chains$iterations)
    burnin <- min(1000, T / 10)
    subchains <- subset(chains, iterations > burnin)
    totalnpoints <- dim(subchains)[1]
    subchains$index <- 1:totalnpoints
    maxnumberpoints <- 50000
    subchains <- subset(subchains, index > totalnpoints - maxnumberpoints)
    library(ggplot2)
    g <- ggplot(subchains, aes_string(x = comp1, y = comp2))
    g <- g + stat_bin2d() + geom_density2d()
    g <- g + theme(title = "2D density")
    return(g)
}

## Plot of the log theta since the last bin split


#'Plot of the log theta penalties
#'
#'This function takes the result of \code{\link{pawl}}, and draws a trace plot
#'of the log theta penalties along the iterations.
#'
#'
#'@param results Object of class \code{"list"}: either the output of
#'\code{\link{pawl}} or of \code{\link{adaptiveMH}}.
#'@return The function returns a ggplot2 object.
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@seealso \code{\link{ggplot}}
PlotLogTheta <- function(results){
  st <- results$splitTimes
  T <- length(results$acceptrates)
  st <- c(0, st, T)
  library(foreach)
  library(reshape)
  library(ggplot2)
  i <- 1; iterations <- 1; value <- 1; estimator <- 1
  df <- foreach (i= 1:(length(st) - 1), .combine = rbind) %do% {
    substart <- st[i] + 1
    substop  <- st[i+1] + 1
    sublogtheta <- data.frame(results$logthetahistory[[i]])
    #thetaDF <- exp(sublogtheta) / apply(exp(sublogtheta), 1, sum)
    thetaDF <- sublogtheta
    names(thetaDF) <- paste("theta", seq(1, results$nbins[i]))
    thetaDF$iterations <- substart:substop
    mdata <- melt(thetaDF, id = c("iterations"))
    names(mdata) <- c("iterations", "estimator", "value")
    mdata
  }
  # trace plot of log theta around the first split
  df$i <- 1:(dim(df)[1])
  maxnumberpoints <- 1000
  iterstep <- floor(dim(df)[1] / maxnumberpoints) + 1
  g <- ggplot(subset(df, i %% iterstep == 0), aes(x = iterations, y = value, colour = estimator))
  g <- g + geom_line()
  g <- g + geom_vline(xintercept = results$splitTimes, linetype = 1)
  g <- g + theme(legend.position = "none")
  return(g)
}     

## Plot of the temperature k


#'Plot of the Flat Histogram occurrences
#'
#'This function takes the result of \code{\link{pawl}}, and draws a plot of the
#'occurrences of the Flat Histogram criteria along the iterations.
#'
#'
#'@param results Object of class \code{"list"}: either the output of
#'\code{\link{pawl}} or of \code{\link{adaptiveMH}}.
#'@return The function returns a ggplot2 object.
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@seealso \code{\link{ggplot}}
PlotFH <- function(results){
    T <- length(results$acceptrates)
    times <- c(1, results$FHtimes, T)
    ks <- results$khistory
    ks <- c(ks, ks[length(ks)])
    kincrease <- data.frame(cbind(times, ks))
    library(ggplot2)
    g <- ggplot(kincrease, aes(x = times, y = ks)) + geom_step()
    g <- g + xlab("iterations") + ylab("k")
    g <- g + theme(title = "Number of flat histogram criteria met along the iterations")
    return(g)
}

## Plot of the number of bins
#### WHAT THE HELL GOES WRONG HERE??


#'Plot of the increase of the number of bins along the iterations
#'
#'This function takes the result of \code{\link{pawl}}, and draws a plot of the
#'increase of the number of bins along the iterations.
#'
#'
#'@param results Object of class \code{"list"}: either the output of
#'\code{\link{pawl}} or of \code{\link{adaptiveMH}}.
#'@return The function returns a ggplot2 object.
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@seealso \code{\link{ggplot}}
PlotNbins <- function(results){
    if (!(is.null(results$splitTimes))){
        binincrease <- data.frame(cbind(c(1, results$splitTimes, T), 
                                        c(results$nbins, results$nbins[length(results$nbins)])))
        g <- ggplot(binincrease, aes_string(x = "X1", y = "X2")) + geom_step()
        g <- g + ylim(0, 1.5 * results$nbins[length(results$nbins)])
        g <- g + ylab("Number of bins") + xlab("iterations")
        g <- g + theme(title = "Number of bins along the iterations")
        return(g)
    } else {
        return(paste("number of bins always was equal to", length(results$finalbins)))
    }
}

## Trace plot of all the variables


#'Trace plot of all the variables
#'
#'This function takes the result of \code{\link{adaptiveMH}} or of
#'\code{\link{pawl}}, and draws a trace plot for each component of the chains
#'
#'
#'@param results Object of class \code{"list"}: either the output of
#'\code{\link{pawl}} or of \code{\link{adaptiveMH}}.
#'@return The function returns a ggplot2 object.
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@seealso \code{\link{ggplot}}
PlotAllVar <- function(results){
    if ("indexchain" %in% names(results)){
        chains <- results
    } else {
        chains <- ConvertResults(results)
    }
    iterations <- 1; value <- 1; indexchain <- 1
    T <- max(chains$iterations)
    burnin <- min(1000, T / 10)
    chains <- subset(chains, iterations > burnin)
    chains[,!(names(chains) %in% c("logdens"))]
    maxnumberpoints <- max(5000, T / 50)
    iterstep <- floor(T / maxnumberpoints) + 1
    library(ggplot2)
    meltedchains <- melt(chains[,!(names(chains) %in% c("logdens"))], id = c("iterations", "indexchain"))
    g <- ggplot(subset(meltedchains, iterations %% iterstep == 0), aes(x = iterations, y = value, colour = indexchain))
    g <- g + geom_line() + facet_wrap( ~ variable) + theme(axis.text.x = element_text(angle = -45, size = 20))
    g <- g + theme(legend.position = "none")
    g <- g + labs(title = "Trace plot of the chains")
    return(g)
}

## Simple histogram of one component


#'Plot a histogram of one component of the chains
#'
#'This function takes the result of \code{\link{adaptiveMH}} or of
#'\code{\link{pawl}}, and a component index, and draws a histogram of it.
#'
#'
#'@param results Object of class \code{"list"}: either the output of
#'\code{\link{pawl}} or of \code{\link{adaptiveMH}}.
#'@param component Object of class \code{"numeric"}: specifies the index of the
#'component to plot on the x-axis.
#'@return The function returns a ggplot2 object.
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@seealso \code{\link{ggplot}}
PlotHist <- function(results, component){
    if ("indexchain" %in% names(results)){
        chains <- results
    } else {
        chains <- ConvertResults(results)
    }
    hist(chains[, component], nclass = 100, 
         main = "Histogram of the chains", 
         xlab = paste("x[", component, "]", sep = ""), prob = TRUE,
         col = "black")
}

## Simple histogram of the binned component 


#'Plot a histogram of the binning coordinate
#'
#'This function takes the result of \code{\link{adaptiveMH}} or of
#'\code{\link{pawl}}, and a \code{\link{binning}} object, and draws a histogram
#'of the chains according to the binning coordinate.
#'
#'
#'@param results Object of class \code{"list"}: either the output of
#'\code{\link{pawl}} or of \code{\link{adaptiveMH}}.
#'@param binning Object of class \code{\link{binning}}, defining the initial
#'bins used by the Wang-Landau algorithm.
#'@return The function returns a ggplot2 object.
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@seealso \code{\link{ggplot}}
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
         col = "black")
    abline(v = binning@bins, lwd = 2, col = "red")
    abline(v = results$finalbins, lwd = 2, lty = 3, col = "red")
}






