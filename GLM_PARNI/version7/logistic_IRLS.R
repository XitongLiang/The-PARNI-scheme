logistic_IRLS <- function(y, X, g, inv_variances, n_trials, p_z, p_gam, max_iter = 100){
  
  # length_y <- length(y)
  # X <- as.matrix(X)
  
  # calculate the starting value of beta
  pi <- sapply(y/n_trials, bounded_x)
  eta <- logit(pi)
  
  W <- n_trials*pi*(1-pi)          
  z <- eta  
  XW <- t(W*X)  
  XWX <- XW %*% X       
  XWz <- XW %*% z                        
  
  if (p_gam > 0) {
    diag(XWX) <- diag(XWX) + c(inv_variances, rep(1/g, p_gam))
  }
  else {
    diag(XWX) <- diag(XWX) + inv_variances
  }
  
  inv_RXWX <- solve(XWX)
  
  # calucation of startval  betahat
  betahat <- inv_RXWX%*%XWz
  
  # initialisation
  diff_beta <- 10  
  iter <- 0
  
  while((diff_beta > 1e-6) & (iter <= max_iter)) {
    
    
    # calculation of linear predictors and means
    eta <- as.vector(X%*%betahat) 
    mu <- inv_logit(eta)            
    
    # calculation of diagonal elements of W 
    W <- n_trials*mu*(1-mu)        
    
    # calculation of adjusted dependent variate z
    z <- eta +  (y-n_trials*mu)/(n_trials*mu*(1-mu)) 
    
    if (any(is.nan(z))) {
      print(y[is.nan(z)])
      print(mu[is.nan(z)])
      print(which(is.nan(z)))
      print(W[which(is.nan(z))])
      nan_idx <- which(is.nan(z))
      z[nan_idx] <- 0 # eta[nan_idx] + (-1)^(2*y[nan_idx])
      
      
    }
    
    # calculation of X'W
    XW <- t(W*X)
    
    # calculation of [X'WX]^-1, X'Wz and U
    XWX <- XW %*% X
    XWz <- XW %*% z 
    
    # U <- XW %*% (z - eta)
    
    if (p_gam > 0) {
      diag(XWX) <- diag(XWX) + c(inv_variances, rep(1/g, p_gam))
    }
    else {
      diag(XWX) <- diag(XWX) + inv_variances
    }
    
    inv_RXWX <- solve(XWX)
    
    # update betahat and go back
    betahat_new <- inv_RXWX%*%XWz 
    
    diff_beta <- max(abs(betahat_new - betahat))
    
    betahat <- betahat_new
    
    # update the counter
    iter <- iter + 1
    
    # print(y[1:10])
    # print(betahat)
    # print(z[1:10])
    
    
  }
  
  # inv_Fisher <- inv_RXWX
  
  
  return(list(betahat = betahat,
              inv_hessian = inv_RXWX,
              iter = iter, 
              eta = eta))
  
}