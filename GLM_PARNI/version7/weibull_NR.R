weibull_NR <- function(t, X, d, k, g, inv_variances, p_z, p_gam, max_iter = 100){
  
  # length_y <- length(y)
  # X <- as.matrix(X)
  
  # calculate the starting value of beta
  lambda <- 1
  eta <- log(lambda)
  
  tilda_y <- k*(t^k - d)
  W <- k^2*t^k       
  
  XW <- t(W*X)  
  XWX <- XW %*% X       
  
  if (p_gam > 0) {
    diag(XWX) <- diag(XWX) + c(inv_variances, rep(1/g, p_gam))
  }
  else {
    diag(XWX) <- diag(XWX) + inv_variances
  }
  
  inv_RXWX <- solve(XWX)
  
  # calucation of startval  betahat
  betahat <- inv_RXWX %*% (t(X) %*% tilda_y)
  
  # initialisation
  diff_beta <- 10  
  iter <- 0
  
  while((diff_beta > 1e-6) & (iter <= max_iter)) {
    
    # calculation of linear predictors and means
    eta <- as.vector(X%*%betahat) 
    lambda <- exp(eta)            
    
    # calculation of diagonal elements of W 
    W <- k^2 * (t*exp(eta))^k   
    tilda_y <- k*((t*exp(eta))^k - d)
    
    
    # calculation of X'W
    XW <- t(W*X)
    
    # calculation of [X'WX]^-1, X'Wz and U
    XWX <- XW %*% X
      
    # U <- XW %*% (z - eta)
    
    if (p_gam > 0) {
      diag(XWX) <- diag(XWX) + c(inv_variances, rep(1/g, p_gam))
    }
    else {
      diag(XWX) <- diag(XWX) + inv_variances
    }
    
    inv_RXWX <- solve(XWX)
    
    # update betahat and go back
    betahat_new <- betahat - inv_RXWX %*% (t(X) %*% tilda_y + c(inv_variances, rep(1/g, p_gam))*betahat) 
    
    diff_beta <- max(abs(betahat_new - betahat))
    
    betahat <- betahat_new
    
    # update the counter
    iter <- iter + 1
    
  }
  
  # inv_Fisher <- inv_RXWX
  
  
  return(list(betahat = as.vector(betahat),
              inv_hessian = inv_RXWX,
              iter = iter, 
              eta = as.vector(eta)))
  
}



