#### the standard dmvnorm is slow because it checks a lot of thing (like is Sigma symmetric)
#### so I replace it with a faster version without checks (it's dangerous though!!!)
fastdmvnorm <- function(x, mu, Sigma){
    distval <- mahalanobis(x, center = mu, cov = Sigma)
    logdet <- sum(log(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values))
    logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
    return(exp(logretval))
}



#'Trimodal target distribution
#'
#'Create the trimodal distribution as in Liang, Liu and Caroll, 2007:
#'Stochastic approximation in Monte Carlo computation.
#'
#'This distribution is a mixture of three bivariate Gaussian distributions.
#'Their covariance matrices are such that an adaptive MCMC algorithm which
#'proposal variance adapts to one of the component, will likely fail to explore
#'the others.
#'
#'@return The function returns an object of class \code{\link{target}}, with a
#'name, a dimension, a function giving the log density, a function to generate
#'sample from the distribution, parameters of the distribution, and a function
#'to draw init points for the MCMC algorithms.
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@seealso \code{\link{target}}, \code{\link{createMixtureTarget}}
#'@references Liang, Liu and Caroll: Stochastic approximation in Monte Carlo
#'computation. Published in JASA, 2007.
createTrimodalTarget <- function(){
  rinit <- function(size){
  fastrmvnorm(n = size, mu = c(0, 0), sigma =  diag(c(0.1, 0.1)))
  }
  logdensity <- function(X, TP){
    component1 <- TP$componentw[1] * fastdmvnorm(X, TP$mean1, TP$cov1)
    component2 <- TP$componentw[2] * fastdmvnorm(X, TP$mean2, TP$cov2)
    component3 <- TP$componentw[3] * fastdmvnorm(X, TP$mean3, TP$cov3)
    return(log(component1 + component2 + component3))
  }   
  generate <- function(n, TP){
      n1 <- floor(n * TP$componentw[1])
      n2 <- floor(n * TP$componentw[2])
      n3 <- floor(n * TP$componentw[3])
      targetsample <- matrix(0, nrow = n1 + n2 + n3, ncol = 2)
      targetsample[1:n1,] <- rmvnorm(n1, TP$mean1, TP$cov1)
      targetsample[(n1 + 1):(n1 + n2),]<- rmvnorm(n2, TP$mean2, TP$cov2)
      targetsample[(n1 + n2 + 1):(n1 + n2 + n3),] <- rmvnorm(n3, TP$mean3, TP$cov3)
      return(targetsample)
  }
  parameters <- parameters <- list(mean1 = c(-8, -8), mean2 = c(6, 6), mean3 = c(0, 0),
        cov1 = matrix(c(1, .9, .9, 1), ncol = 2), cov2 = matrix(c(1, -.9, -.9, 1), ncol = 2), 
        cov3 = matrix(c(1, 0, 0, 1), ncol = 2), componentw = c(1/3, 1/3, 1/3))
  targetinstance <- target(name = "trimodal", dimension = 2,
                      rinit = rinit, logdensity = logdensity,
                      generate = generate, parameters = parameters)
  return(targetinstance)
}
#trimodal <- createTrimodalTarget()
#slotNames(trimodal)
#show(trimodal)
