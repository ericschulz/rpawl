

#'Class \code{"binning"}
#'
#'This class holds all the parameters of the Parallel Adaptive Wang-Landau
#'algorithm that are related to the bins: it includes the functions that take
#'points and return the point locations with respect to the bins, parameters
#'related to the number of bins, the split mechanism, the adaptation rate of
#'the stochastic approximation schedule, etc.
#'
#'
#'@name binning
#'@aliases binning-class binning binning-methods binning-method
#'binning,ANY-method show,binning-method
#'@docType class
#'@section Objects from the Class: Objects should created by calls of the
#'function \code{binning}. Examples are provided that should help understanding
#'this class. Essentially it is a list of parameters, most of which have a
#'reasonable default value so you do not need to think about it too much.
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





#'Image of ice floes
#'
#'This data represents a binary matrix, representing an image of ice floes.
#'
#'
#'@name IceFloe
#'@aliases IceFloe icefloe
#'@docType data
#'@format A matrix containing 40 rows and 40 columns
#'@source Banfield, J. and Raftery, A. (1992). Ice floe identification in
#'satellite images using mathematical morphology and clustering about principal
#'curves. Journal of the American Statistical Association, 87(417):7-16.
#'@keywords datasets
NULL





#'PARALLEL ADAPTIVE WANG-LANDAU
#'
#'The package implements the Parallel Adaptive Wang-Landau algorithm on various
#'examples.  The provided demos allow to reproduce the figures of the article.
#'
#'\tabular{ll}{ Package: \tab PAWL\cr Type: \tab Package\cr Version: \tab
#'1.0\cr Date: \tab 2011-08-11\cr License: \tab GPL (>= 2)\cr LazyLoad: \tab
#'yes\cr Depends: \tab mvtnorm\cr Suggests: \tab ggplot2\cr } The main function
#'is \code{pawl}. It takes algorithmic parameters in arguments (see the help of
#'the \code{pawl} function), as well a target distribution. Look at the demos
#'to learn how to specify a target distribution.
#'
#'@name PAWL-package
#'@aliases PAWL-package PAWL
#'@docType package
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@keywords package
#'@examples
#'
#'  demo(discreteexample)
#'  demo(gaussianexample)
#'  demo(mixture2kexample)
#'
NULL





#'Pollution Data
#'
#'This data contains 1 response (mortality, normalized to have mean zero) along
#'with 15 pollution-related explanatory variables.
#'
#'
#'@name Pollution
#'@aliases Pollution pollution
#'@docType data
#'@format A matrix containing 60 rows and 16 columns
#'@source McDonald, G.C. and Schwing, R.C. (1973) 'Instabilities of regression
#'estimates relating air pollution to mortality', Technometrics, vol.15,
#'463-482.
#'@keywords datasets
NULL





#'SMC Tuning Parameters
#'
#'This class holds parameters for the Sequential Monte Carlo sampler.
#'
#'
#'@name smcparameters
#'@aliases smcparameters-class smcparameters smcparameters,ANY-method
#'show,smcparameters-method
#'@docType class
#'@section Objects from the Class: Objects can be created by calls of the
#'function \code{"smcparameters"}.
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@seealso \code{\link{smc}}
#'@keywords classes
#'@examples
#'
#'showClass("smcparameters")
#'smcparam<- smcparameters(nparticles=5000, 
#'                        temperatures = seq(from = 0.0001, to = 1, length.out= 100),
#'                        nmoves = 5, ESSthreshold = 0.5, movetype = "randomwalk",
#'                        movescale = 0.1)
#'
NULL





#'Class: target distribution
#'
#'This class represents target distributions, that is, probability
#'distributions from which we want to sample using MCMC or Wang-Landau.
#'
#'
#'@name target
#'@aliases target-class target target,ANY-method show,target-method
#'@docType class
#'@section Objects from the Class: Objects should created by calls of the
#'function \code{target}. Examples are provided that should help implementing
#'any continuous probability distributions. %% ~~ describe objects here ~~
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@keywords classes
#'@examples
#'
#'  showClass("target")
#'  # starting points for MCMC algorithms
#'  rinit <- function(size) rnorm(size)
#'  # target log density function: a gaussian distribution N(mean = 2, sd = 3)
#'  parameters <- list(mean = 2, sd = 3)
#'  logdensity <- function(x, parameters) dnorm(x, parameters$mean, parameters$sd, log = TRUE)
#'  # creating the target object
#'  gaussiantarget <- target(name = "gaussian", dimension = 1,
#'                    rinit = rinit, logdensity = logdensity,
#'                    parameters = parameters)
#'  print(gaussiantarget)
#'
NULL





#'MCMC Tuning Parameters
#'
#'This class holds tuning parameters for the Metropolis-Hastings and
#'Wang-Landau algorithms.
#'
#'
#'@name tuningparameters
#'@aliases tuningparameters-class tuningparameters tuningparameters,ANY-method
#'show,tuningparameters-method
#'@docType class
#'@section Objects from the Class: Objects can be created by calls of the
#'function \code{"tuningparameters"}.
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@seealso \code{\link{adaptiveMH}} \code{\link{preexplorationAMH}}
#'\code{\link{pawl}}
#'@keywords classes
#'@examples
#'
#'showClass("tuningparameters")
#'mhparameters <- tuningparameters(nchains = 10, niterations = 1000, adaptiveproposal = TRUE) 
#'
NULL



