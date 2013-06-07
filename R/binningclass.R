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
#' Class \code{"binning"}
#'
#'This class holds all the parameters of the Parallel Adaptive Wang-Landau
#'algorithm that are related to the bins: it includes the functions that take
#'points and return the point locations with respect to the bins, parameters
#'related to the number of bins, the split mechanism, the adaptation rate of
#'the stochastic approximation schedule, etc.
#'
#'
#'@name binning
#'@rdname binning
#'@aliases binning-class binning binning-methods binning-method
#'binning,ANY-method show,binning-method
#'@docType class
#'@section Objects from the Class: Objects should created by calls of the
#'function \code{binning}. Examples are provided that should help understanding
#'this class. Essentially it is a list of parameters, most of which have a
#'reasonable default value so you do not need to think about it too much.
#'@section Important slots:
#'  \describe{
#'    \item{\code{position}:}{Object of class \code{"function"}: should be a function taking points and associated log density values, and returning a "reaction coordinate", that is, a value that will be associated with bins. Typically, it can be the log density itself, or one component of a d-dimensional point. See the example below.}
#'    \item{\code{binrange}:}{ Object of class \code{"numeric"}: it should be a vector of size 2, holding the minimum and the maximum on the reaction coordinate scale. The bins are going to be between those two (inner bins), while a bin will go from - infinity to the minimum, and a bin will go from maximum to + infinity (outer bins).
#'    }
#'    \item{\code{ncuts}:}{ Object of class \code{"numeric"}: how many cuts will be made in the bin range specified by the previous argument. This induce the number of initial bins. Bins are automatically created by the following line: 
#'
#'      \code{
#'      bins <- c(-Inf, seq(from = binrange[1], to = binrange[2], length.out = ncuts))}
#'
#'    There are then (ncuts +1) bins. The default for ncuts is 9, resulting in 10 bins.
#'    }
#'  }
#'
#'@section Optional slots:
#'  \describe{
#'    \item{\code{bins}:}{ Object of class \code{"numeric"}: you can specify the bins directly, in which case you do not need to specify \code{binrange}.
#'    }
#'    \item{\code{name}:}{ Object of class \code{"character"}: ... if you want to name the instance (default is "unspecified").
#'    }
#'    \item{\code{autobinning}:}{ Object of class \code{"logical"}: activate or not the splitting mechanism, to create new inner bins automatically. This does not create new bins outside the specified bin range, it just add new bins inside to help reaching the Flat Histogram criteria more quickly.
#'    }
#'    \item{\code{desiredfreq}:}{ Object of class \code{"numeric"}: you can specify the desired frequency of each bin. The default is 1 / nbins in each bin, where nbins is the number of bins. Note that if autobinning is enable, when a bin is split into two bins, the desired frequencies of the new bins are equal to half of the desired frequency of the former bin.
#'    }
#'    \item{\code{useLearningRate}:}{ Object of class \code{"logical"}: active or not the stochastic approximation schedule. That is, if it is not activated, then no schedule are used in the update of theta (the penalty associated to the bins). Default is TRUE.
#'    }
#'    \item{\code{useFH}:}{ Object of class \code{"logical"}: active or not the Flat Histogram checks. If it is not activated, then the stochastic approximation decreases at each step. Default is TRUE, unless useLearningRate is FALSE, in which case there is no point checking for Flat Histograms.
#'    }
#'    \item{\code{fhthreshold}:}{ Object of class \code{"numeric"}: specifies the threshold to accept Flat Histogram. The default is 0.5. Smaller values make the Flat Histogram criterion harder to reach. 
#'    }
#'    \item{\code{minSimEffort}:}{ Object of class \code{"numeric"}: specifies the minimum number of iterations after a Flat Histogram, for a new Flat Histogram criterion to be accepted. It prevents the criterion to be accepted at every iteration when using a large number of parallel chains. Default is 200.
#'    }
#'    \item{\code{learningrate}:}{ Object of class \code{"function"}: specifies the learning rate, that is, the rate at which the stochastic schedule decreases. It should be a function defined on [0, + infty[ such that it is not integrable but its square is integrable, e.g. t -> 1/t for instance. The default is t -> t^{-0.6}.
#'    }
#'    \item{\code{splitThreshold}:}{ Object of class \code{"numeric"}: specifies the threshold to split a bin into two new bins. The default is 0.1 (read 10\%), which means that a bin is split if at least 90\% of the points in that bin are on the half right (or left) side of the bin. Larger values (e.g. 25\%) result in more splits, and hence more final bins.
#'    }
#'  }
#'
#'@section Methods:
#'  \describe{
#'    \item{show}{\code{signature(object = "binning")}: provides a little summary of a binning object when called (or when \code{print} is called).}
#'   }
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@keywords classes
#'@examples
#'
#'  showClass("binning")
#'  getPos <- function(points, logdensity) points
#'  positionbinning <- binning(position = getPos,
#'                        name = "position",
#'                        binrange = c(-4, 0),
#'                        ncuts = 4,
#'                        autobinning = TRUE,
#'                        useLearningRate = TRUE)
#'
NULL

setClass("binning",
         representation(position="function", getLocations="function",
                        name = "character", 
                        bins = "numeric",
                        desiredfreq = "numeric", fhthreshold = "numeric",
                        splitThreshold = "numeric",
                        minSimEffort = "numeric", 
                        learningrate = "function", useLearningRate = "logical",
                        useFH = "logical", binmids = "numeric",
                        autobinning = "logical",
                        smoothbinning = "logical",
                        alongenergy = "logical"))

setGeneric("binning", function(...)standardGeneric("binning"))
binning.constructor <- function(position, name, 
                                autobinning, ncuts, binrange, bins,
                                desiredfreq, fhthreshold, splitThreshold,
                                minSimEffort, learningrate, useLearningRate, 
                                useFH, alongenergy, smoothbinning){
  if (missing(name))
    name <- "unspecified"
  if (missing(position))
    stop(sQuote("position"), "has to be specified")
  if (missing(bins)){
      if (missing(binrange)){
          stop("error: ", sQuote("bins"), " is not specified, so you have to specify at least ",
sQuote("binrange"), " (and maybe ", sQuote("ncuts"), " too)\n")
      }
      if (missing(ncuts)){
          ncuts <- 9
          cat("number of cuts (ncuts) unspecified: set to 9\n")
      }
      bins <- c(-Inf, seq(from = binrange[1], to = binrange[2], length.out = ncuts))
  }
  if (missing(desiredfreq)){
      desiredfreq <- rep(1 / length(bins), length(bins))
      cat("desiredfreq unspecified: set to 1 / nbins for each bin\n")
  } else {
      if (length(desiredfreq) != length(bins))
          stop(sQuote("desiredfreq"), " provided but not the right length")
      if (sum(desiredfreq) != 1)
          stop(sQuote("desiredfreq"), " should sum to 1")
  }
  if (missing(fhthreshold)){
      fhthreshold <- 0.5 
      cat("fhthreshold unspecified: set to 0.5 \n")
  }
  if (missing(splitThreshold)){
      splitThreshold <- 0.1
      cat("splitThreshold unspecified: set to 0.1\n")
  }
  if (missing(minSimEffort)){
      minSimEffort <- 200
      cat("minSimEffort unspecified: set 200\n")
  }
  if (missing(learningrate)){
      learningrate <- function(x) x^(-0.6)
      cat("learningrate unspecified: set to x -> x^(-1)\n")
  }
  if (missing(useLearningRate)){
      useLearningRate <- TRUE
      cat("useLearningRate unspecified: set to TRUE\n")
  }
  if (missing(useFH)){
      useFH <- useLearningRate
      cat("useFH unspecified: set to", useLearningRate, "\n")
  }
  if (missing(autobinning)){
      autobinning <- FALSE
      cat("autobinning unspecified: set to FALSE\n")
  }
  if (missing(alongenergy)){
      alongenergy <- FALSE
  }
  if (missing(smoothbinning)){
      smoothbinning <- FALSE
  }
  if (!useLearningRate & useFH){
    stop("error: you specified useFH to be TRUE and useLearningRate to be FALSE; check the help files")      
  }
  getLocations <- function(bins, somevector) findInterval(somevector, bins)
  new("binning", position = position, getLocations = getLocations, 
      name = name, autobinning = autobinning, 
      bins = bins, desiredfreq = desiredfreq, fhthreshold = fhthreshold,
      splitThreshold = splitThreshold, minSimEffort = minSimEffort,
      learningrate = learningrate, useLearningRate = useLearningRate,
      useFH = useFH, alongenergy = alongenergy, smoothbinning = smoothbinning)
}
setMethod("binning",
          definition = function(position, name, 
                                autobinning, ncuts, binrange, bins, desiredfreq, fhthreshold,
                                splitThreshold, minSimEffort,
                                learningrate, useLearningRate, useFH, alongenergy, smoothbinning){
            binning.constructor(position = position,
                                name = name, autobinning = autobinning, 
                                ncuts = ncuts, binrange = binrange, bins = bins,
                                desiredfreq = desiredfreq, fhthreshold = fhthreshold,
                                splitThreshold = splitThreshold, minSimEffort = minSimEffort,
                                learningrate = learningrate, useLearningRate = useLearningRate,
                                useFH = useFH, alongenergy = alongenergy, smoothbinning = smoothbinning)
          })

setMethod(f = "show", signature = "binning", 
          def = function(object){
            cat("Object of class ", class(object), ".\n", sep = "")
            cat("*name:", object@name, "\n")
            cat("*autobinning:", object@autobinning, "\n")
            cat("*bins:", object@bins, "\n")
            cat("*Flat Histogram threshold:", object@fhthreshold, "\n")
            cat("*split threshold:", object@splitThreshold, "\n")
            cat("*use Flat Histogram criterion:", object@useFH, "\n")
            cat("*use learning rate:", object@useLearningRate, "\n")
            cat("*min number of iterations between FH:", object@minSimEffort, "\n")
          })

## Function taking bins and giving middles of inner bins ...
getBinMiddles <- function(object){
    nbins <- length(object@bins)
    binmids <- vector(length = nbins - 2)
    for (indexbin in 1:(nbins - 2)){
        binmids[indexbin] <- object@bins[indexbin + 1] + 
        (object@bins[indexbin + 2] - object@bins[indexbin + 1]) / 2
    }
    return(binmids)
}
getInnerLeftCounts <- function(object, currentreaction, currentlocations){
    nbins <- length(object@bins)
    innerbinleftcount <- rep(0, nbins - 2)
    for (indexinnerbin in 1:(nbins - 2)){
        localreaction <- currentreaction[currentlocations == indexinnerbin + 1]
        innerbinleftcount[indexinnerbin] <- sum(localreaction < object@binmids[indexinnerbin])
    }
    return(innerbinleftcount)
}
# Find which bins would benefit from a split
findSkewedBins <- function(object, innerbinleftcount, bincount){
    nbins <- length(object@bins)
    innerbincount <- bincount[2:(nbins-1)]
    skewedbins <- c(); newcuts <- c()
    for (indexbin in 1:(nbins - 2)){
        proportionleft <- innerbinleftcount[indexbin] / innerbincount[indexbin]
        if (is.na(proportionleft)){
            # in this case there is no point in that bin
            next
        }
        if (proportionleft <= object@splitThreshold | 
            (1 - proportionleft) <= object@splitThreshold){
            skewedbins <- c(skewedbins, indexbin)
            newcuts <- c(newcuts, object@binmids[indexbin])
        }
    }
    return(list(newcuts = newcuts, skewedbins = skewedbins))
}

binsplitter <- function(object, foundbins, oldthetas, olddesiredfreq, alongenergy){
    binsToSplit <- foundbins$binsToSplit
    newcuts <- foundbins$newcuts
    oldbins <- object@bins
    oldnbins <- length(oldbins)
    newbins <- sort(c(oldbins, newcuts))
    # update the bias: if bin i is split, the former bias theta(i) is cut into
    # two biases theta_1(i) and theta_2(i), each equal to theta(i) / 2
    newthetas <- c(oldthetas[1])
    newdesiredfreq <- c(olddesiredfreq[1])
    for (i in 2:(oldnbins - 1)){
        if (sum(i == binsToSplit)){
            if (alongenergy){
                newthetas <- c(newthetas, c(oldthetas[i] + log(2), 
                                            oldthetas[i] - log(2)))
            } else {
                newthetas <- c(newthetas, rep(oldthetas[i] - log(2), 2))
            }
            newdesiredfreq <- c(newdesiredfreq, rep(olddesiredfreq[i] / 2, 2))
        } else {
            newthetas <- c(newthetas, oldthetas[i])
            newdesiredfreq <- c(newdesiredfreq, olddesiredfreq[i])
        }
    }
    newthetas <- c(newthetas, oldthetas[oldnbins])
    newdesiredfreq <- c(newdesiredfreq, olddesiredfreq[oldnbins])
    return(list(newbins = newbins, newthetas = newthetas, newdesiredfreq = newdesiredfreq))
}

#getPos <- function(points, logdensity) points
#positionbinning <- binning(position = getPos,
#                            name = "position",
#                            binrange = c(-4, 0),
#                            ncuts = 4,
#                            autobinning = TRUE,
#                            useLearningRate = TRUE)
#
#getBinMiddles(positionbinning)
#positionbinning@binmids <- getBinMiddles(positionbinning)
#print(positionbinning)
#nbins <- length(positionbinning@bins)
#tempbincount <- c(5, 20, 4, 5, 19)
#innerbincount <- tempbincount[2:(nbins-1)]
#innerbinleftcount <- innerbincount- c(2, 2, 0.1)
#print(innerbinleftcount / innerbincount)
#foundbins <- findBinsToSplit(positionbinning, innerbinleftcount, tempbincount)
#splitttt <- binsplitter(positionbinning, foundbins, c(1, 1, 1, 1, 1), c(0.2, 0.2, 0.2, 0.2, 0.2))
#print(splitttt)
#
#
#
#
