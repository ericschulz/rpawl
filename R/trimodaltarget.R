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
