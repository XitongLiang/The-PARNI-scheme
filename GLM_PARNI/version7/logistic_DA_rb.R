logistic_DA_rb <- function(ALL, hyper_par) {
  
  
  p <- hyper_par$p
  kappa <- hyper_par$kappa
  X <- hyper_par$X
  p_z <- hyper_par$p_z
  h <- hyper_par$h_exp
  
  curr <- ALL$curr
  omega <- ALL$omega
  A <- ALL$A # t(J_gamma) %*% hyper_par$kappa
  L_invB <- ALL$L_invB
  g <- ALL$g
  
  includes <- which(curr == 1)
  
  tilda_X <- X * sqrt(omega)
  kappa_X <- t(kappa) %*% X
  # kappa_tildaX <- t(kappa) %*% tilda_X
  tilda_diag_XtX <- colSums(tilda_X^2)
  
  # X_gamma <- as.matrix(X[, which(curr==1)], ncol = p_gam)
  # J_gamma <- cbind(hyper_par$Z, X_gamma)
  
  J_gamma <- ALL$J_gamma
  tilda_J <- J_gamma * sqrt(omega)
  inv_B <- t(L_invB) %*% L_invB
  
  kappa_Jg_invB <- t(kappa) %*% J_gamma %*% inv_B
  kappa_Jg_invL <- as.numeric(t(kappa) %*% J_gamma %*% t(L_invB))
  tXT_tJg_invL <- t(tilda_X) %*% tilda_J %*% t(L_invB)
  
  d_vec <- tilda_diag_XtX + 1/g - rowSums(tXT_tJg_invL^2)
  
  d_vec[includes] <- 1
  Bayes_factor <- as.numeric(d_vec^{-1/2}*g^{-1/2}*exp((kappa_X - colSums(kappa_Jg_invL * t(tXT_tJg_invL)))^2/(2*d_vec)))
  q_j <- 0
  
  for (j in includes) {
    
    q_j <- q_j + 1
    d_j <- 1/inv_B[(p_z+q_j),(p_z+q_j)]
    Bayes_factor[j] <- d_j^{-1/2}*g^{-1/2}*exp(d_j*kappa_Jg_invB[(p_z+q_j)]^2/2)
    
  }
  
  
  PIPs <- h*Bayes_factor/(1-h+h*Bayes_factor)
  
  return(PIPs)
  
  
  
}

