setClass("binning",
         representation(position="function", getLocations="function",
                        name = "character", autobinning = "logical",
                        bins = "numeric",
                        desiredfreq = "numeric", fhthreshold = "numeric",
                        splitThreshold = "numeric",
                        minSimEffort = "numeric", 
                        learningrate = "function", useLearningRate = "logical",
                        useFH = "logical", binmids = "numeric"))

setGeneric("binning", function(...)standardGeneric("binning"))
binning.constructor <- function(position, name, 
                                autobinning, ncuts, binrange, bins,
                                desiredfreq, fhthreshold, splitThreshold,
                                minSimEffort, learningrate, useLearningRate, 
                                useFH){
  if (missing(name))
    name <- "unspecified"
  if (missing(position))
    stop(sQuote("position"), "has to be specified")
  #if (missing(getLocations))
  #  stop(sQuote("getLocations"), "has to be specified")
  if (missing(autobinning)){
      autobinning <- FALSE
      cat("autobinning unspecified: set to FALSE\n")
  }
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
  if (!useLearningRate & useFH){
    stop("error: you specified useFH to be TRUE and useLearningRate to be FALSE; check the help files")      
  }
  getLocations <- function(bins, somevector) findInterval(somevector, bins)
  new("binning", position = position, getLocations = getLocations, 
      name = name, autobinning = autobinning, 
      bins = bins, desiredfreq = desiredfreq, fhthreshold = fhthreshold,
      splitThreshold = splitThreshold, minSimEffort = minSimEffort,
      learningrate = learningrate, useLearningRate = useLearningRate,
      useFH = useFH)
}
setMethod("binning",
          definition = function(position, name, 
                                autobinning, ncuts, binrange, bins, desiredfreq, fhthreshold,
                                splitThreshold, minSimEffort,
                                learningrate, useLearningRate, useFH){
            binning.constructor(position = position,
                                name = name, autobinning = autobinning, 
                                ncuts = ncuts, binrange = binrange, bins = bins,
                                desiredfreq = desiredfreq, fhthreshold = fhthreshold,
                                splitThreshold = splitThreshold, minSimEffort = minSimEffort,
                                learningrate = learningrate, useLearningRate = useLearningRate,
                                useFH = useFH)
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
findBinsToSplit <- function(object, innerbinleftcount, tempbincount){
    nbins <- length(object@bins)
    innerbincount <- tempbincount[2:(nbins-1)]
    binsToSplit <- c(); newcuts <- c()
    for (indexbin in 1:(nbins - 2)){
        #print(innerbinleftcount)
        #print(innerbincount)
        proportionleft <- innerbinleftcount[indexbin] / innerbincount[indexbin]
        if (proportionleft <= object@splitThreshold | 
            (1 - proportionleft) <= object@splitThreshold){
            #cat("inner bin", indexbin, ", proportion on the left side", proportionleft, "\n")
            #cat("->hence we split inner bin", indexbin, "\n")
            binsToSplit <- c(binsToSplit, indexbin)
            newcuts <- c(newcuts, object@binmids[indexbin])
        }
    }
    if (is.null(newcuts))
        print("no bin to split")
    else
        cat("need to split inner bins:", binsToSplit, "\n")
    return(list(newcuts = newcuts, binsToSplit = binsToSplit))
}

binsplitter <- function(object, foundbins, oldthetas, olddesiredfreq){
    binsToSplit <- foundbins$binsToSplit
    newcuts <- foundbins$newcuts
    oldbins <- object@bins
    oldnbins <- length(oldbins)
    newbins <- sort(c(oldbins, newcuts))
    # update the bias: if bin i is split, the former bias theta(i) is cut into
    # two biases theta_1(i) and theta_2(i), each equal to theta(i) / 2
    newthetas <- c(oldthetas[1])
    newdesiredfreq <- c(olddesiredfreq[1])
    for (i in 1:(oldnbins - 2)){
        if (sum(i == binsToSplit)){
            newthetas <- c(newthetas, rep(oldthetas[i + 1] - log(2), 2))
            newdesiredfreq <- c(newdesiredfreq, rep(olddesiredfreq[i + 1] / 2, 2))
        } else {
            newthetas <- c(newthetas, oldthetas[i + 1])
            newdesiredfreq <- c(newdesiredfreq, olddesiredfreq[i + 1])
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
