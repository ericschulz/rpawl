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
#rm(list = ls())
#library(PAWL)


#'Mixture target distribution
#'
#'Create the posterior distribution of the parameters of a mixture of
#'univariate gaussian distributions, with a fixed (known) number of components.
#'
#'
#'@param mixturesample Object of class \code{"vector"}: data set to be used. If
#'not provided, a synthetic data set is generated.
#'@param mixturesize Object of class \code{"numeric"}: represents the data set
#'size if a data set is to be generated.
#'@param ncomponents Object of class \code{"numeric"}: represents the fixed
#'number of components to be used.
#'@param mixtureparameters Object of class \code{"list"}: provides the
#'parameters to be used if a data set has to be generated.  The parameters
#'include the number of components, the component weights, means and variances.
#'@return The function returns an object of class \code{\link{target-class}},
#'with a name, a dimension, a function giving the log density, a function to
#'generate sample from the distribution, parameters of the distribution, and a
#'function to draw init points for the MCMC algorithms. The log density
#'involves a likelihood and a prior, and the prior is as in Richardson and
#'Green, "On Bayesian analysis of mixtures with an unknown number of
#'components", published in JRSS B, 1997.
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@seealso \code{\link{target-class}}, \code{\link{createTrimodalTarget}}
#'@references Jasra, Holmes, Stephens, "MCMC and label switching problem in
#'Bayesian mixture models", published in Statistical Science (2005). Richardson
#'and Green, "On Bayesian analysis of mixtures with an unknown number of
#'components", published in JRSS B, 1997.
createMixtureTarget <- function(mixturesample, mixturesize, ncomponents,
                                mixtureparameters){
  if (missing(ncomponents)){
    cat("ncomponents unspecified: setting it to 2\n")
    ncomponents <- 2
  }
  if (!missing(mixturesample)){
    mixturesize <- length(mixturesample)
    cat("sample provided, of length", mixturesize, "\n")
    rmixture <- function(n, MP) stop(sQuote("generate"), "is not specified")
  } else {
    if (missing(mixtureparameters)){
      cat("mixtureparameters unspecified: putting some default parameters\n")
      mixtureparameters <- list(
                                ncomponents = ncomponents,
                                componentweights = rep(0.5, ncomponents), 
                                componentmeans = seq(from = 0, to = ncomponents * 1.25, 
                                                     length.out = ncomponents), 
                                componentvariances = rep(0.3, ncomponents))
      print(mixtureparameters)
    }
    # Function to generate n-sample from a mixture model with parameters MP
    rmixture <- function(n, MP){
      tempmixturesample <- c()
      repartition <- rmultinom(1, n, MP$componentweights)
      for (j in 1:MP$ncomponents){
        tempmixturesample <- c(tempmixturesample, 
                               rnorm(repartition[j], MP$componentmeans[j], 
                                     sqrt(MP$componentvariances[j])))
      }
      return(sample(tempmixturesample))
    }
    if (missing(mixturesize)){
      cat("no sample provided and no sample size: setting size to 100\n")
      mixturesize <- 100
    }
    cat("generating sample of size", mixturesize, "\n")
    mixturesample <- rmixture(mixturesize, mixtureparameters)
  }
  # Function to generate hyper parameters using the sample
  makehyperparameters <- function(mixturesample){
    Range = max(mixturesample) - min(mixturesample)
    hyperparameters <- list(
                            delta = 3,
                            gammadelta = gamma(3), 
                            xi = mean(mixturesample),
                            tau = Range**2 / 4,
                            alpha = 2,
                            gammaalpha = gamma(2),
                            g = 0.2,
                            gammag = gamma(0.2),
                            h = 100 * 0.2 / (2 * (Range**2)))
    return(hyperparameters)
  }
  # Generate hyperparameters for the loaded sample
  parameters <- list(ncomponents = ncomponents, 
                     hyperparameters = makehyperparameters(mixturesample),
                     mixturesample = mixturesample)

  ###########################
  ### Specify prior distributions
  #
  # Function to generate from the prior distribution
  # e.g. to get starting values for algorithms
  rprior <- function(size, ncomponents, hyperparameters){
    weightsprior <- matrix(rgamma(ncomponents * size, 
                                  shape = hyperparameters$delta, rate = 1),
                           ncol = ncomponents, nrow = size)
    musprior <- matrix(rnorm(ncomponents * size, 
                             hyperparameters$xi, sqrt(hyperparameters$tau)),
                       ncol = ncomponents, nrow = size)
    betasprior <- matrix(rgamma(size, shape = hyperparameters$g, rate = hyperparameters$h),
                         ncol = 1, nrow = size)
    precisionsprior <- matrix(NA, ncol = ncomponents, nrow = size)
    for (indexcomponent in 1:ncomponents){
      precisionsprior[, indexcomponent] <- rgamma(size, shape = hyperparameters$alpha, rate = 1 / betasprior)
    }
    return(cbind(log(weightsprior), musprior, log(precisionsprior), log(betasprior)))
  }
  # Define the "target" dictionary to be used by the MCMC algorithms
  dimension <- parameters$ncomponents * 3 + 1 
  # Function to generate initial values
  # Here we use the prior distribution
  rinit <- function(size){
    initvalues <- rprior(size, parameters$ncomponents, parameters$hyperparameters)
    return(initvalues)
  }
  dinit <- function(x, TP){
    resultsprior <- .Call("mixturelogpriordensity", x, as.numeric(TP$hyperparameters))
    return(resultsprior)
  }
  logdensity <- function(x, TP){
    resultsprior <- .Call("mixturelogpriordensity", x, as.numeric(TP$hyperparameters))
    resultslikelihood <- .Call("mixtureloglikelihood", x, TP$mixturesample)
    results <- resultsprior + resultslikelihood
    results[is.infinite(results)] = -(10**150)
    return(results)
  }
  targetinstance <- target(name = paste("mixture with", ncomponents, "components"), dimension = dimension,
                           rinit = rinit, dinit = dinit, logdensity = logdensity,
                           parameters = parameters, generate = rmixture)
  return(targetinstance)
}

#mixture <- createMixtureTarget()
#slotNames(mixture)
#show(mixture)
#
#a <- mixture@rinit(4)
#mixture@logdensity(a, mixture@parameters)


