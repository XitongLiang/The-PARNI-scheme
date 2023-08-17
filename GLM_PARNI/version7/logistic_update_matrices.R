logistic_update_matrices <- function(ALL, hyper_par, ...) {
  
  
  kappa <- hyper_par$kappa
  Z <- hyper_par$Z
  p_z <- hyper_par$p_z
  p <- hyper_par$p
  inv_variances <- hyper_par$inv_variances
  
  curr <- ALL$curr
  p_gam <- ALL$p_gam
  omega <- ALL$omega
  g <- ALL$g
  
  
  if (p_gam >= 1) {
    
    X_gamma <- as.matrix(hyper_par$X[, which(curr==1)], ncol = p_gam)
    J_gamma <- cbind(hyper_par$Z, X_gamma)
    # print(J_gamma)
    A <- t(J_gamma) %*% kappa
    
    tilda_J <- J_gamma * sqrt(omega)
    
    B <- t(tilda_J) %*% tilda_J
    diag(B) <- diag(B) + c(inv_variances, rep(1/g, p_gam))
    
    inv_B <- solve(B)
    L_invB <- chol(inv_B)
    
    ALL$A <- A
    ALL$L_invB <- L_invB
    
    
    # L_B <- chol(B)
    # inv_L_B <- backsolve(L_B, diag(p_z + p_gam))
    # ALL$inv_L_B <- inv_L_B
    
    # log_det_L_B <- sum(log(diag(inv_L_B)))
    # ALL$llh <- - p_gam*log(g)/2 + log_det_L_B + sum((t(inv_L_B) %*% A)^2)/2
    # cat(log_det_L_B, ",", -determinant(B)$modulus/2, "\n")
    
    log_det_L_B <- sum(log(diag(L_invB)))
    ALL$llh <- - p_gam*log(g)/2 + log_det_L_B + sum((L_invB %*% A)^2)/2
    
    
  }
  else if (p_gam==0){
    
    J_gamma <- hyper_par$Z
    A <- t(J_gamma) %*% kappa
    tilda_J <- J_gamma * sqrt(omega)
    B <- t(tilda_J) %*% tilda_J
    diag(B) <- diag(B) + inv_variances
    # L_B <- chol(B)
    # inv_L_B <- backsolve(L_B, diag(p_z))
    
    # ALL$A <- A
    # ALL$inv_L_B <- inv_L_B
    
    
    
    inv_B <- solve(B)
    L_invB <- chol(inv_B)
    
    ALL$A <- A
    ALL$L_invB <- L_invB
    
    
    
    if (p_z == 1) {
      log_det_L_B <- log(L_invB)
    }
    else {
      log_det_L_B <- sum(log(diag(L_invB)))
    }
    ALL$llh <- log_det_L_B + sum((L_invB %*% A)^2)/2
    
  }
  
  ALL$approx_llh <- ALL$llh
  
  ALL$J_gamma <- J_gamma
  ALL$lmp <- hyper_par$log_m_prior(p_gam, hyper_par$h, p)
  
  
  return(ALL)
  
}



