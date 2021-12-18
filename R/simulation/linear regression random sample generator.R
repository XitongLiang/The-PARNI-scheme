# linear regression random sample generator
# library(MASS)

lrrsg <- function(n, p, a = 0,
                 beta_basis = c(2,-3,2,2,-3,3,-2,3,-2,3), 
                 rho = 0,
                 SNR = 1, sigma2 = 1, seed = NULL,
                 hyper_par = FALSE){
  
  beta <- rep(0, p)
  b_n <- length(beta_basis)
  
  set.seed(seed)
  
  position <- 1:b_n
  
  beta[position] <- beta_basis
  
  beta <- SNR * sqrt(sigma2 * log(p)/n) * beta
  
  b <- sqrt(1-rho^2)
  
  X <- matrix(rnorm(n*p), nrow = n, ncol = p)
  
  for (j in 2:p) {
    X[,j] <- rho * X[,j-1] + b * X[,j]
  }
  
  
  # responses
  y <- a + X %*% beta + rnorm(n)
  
  if(hyper_par){
    
    y <- y - mean(y)
    X <- t(t(X)-colMeans(X))
    
    return(list(y = y,
                X = X,
                p = p,
                n = n,
                yty = sum(y^2),
                XtX = t(X) %*% X,
                ytX = t(y) %*% X,
                h = b_n/p,
                g = 9))
  }
  else{
  return(list(a = a,
              beta = beta, 
              y = y, 
              X = X, 
              position = position))
  }
}




