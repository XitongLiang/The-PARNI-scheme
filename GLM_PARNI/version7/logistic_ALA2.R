# fixed n_tails
logistic_ALA2 <- function(ALL, hyper_par, eta) {
  
  
  curr <- ALL$curr
  p_gam <- ALL$p_gam
  g <- ALL$g
  
  p_z <- hyper_par$p_z
  # kappa <- hyper_par$kappa
  y <- hyper_par$y
  inv_variances <- hyper_par$inv_variances
  
  if (p_gam >= 1) {
    X_gamma <- as.matrix(hyper_par$X[, which(curr==1)], ncol = p_gam)
    J_gamma <- cbind(hyper_par$Z, X_gamma)
  }
  else {
    J_gamma <- hyper_par$Z
  }
  
  
  mu <- inv_logit(eta)
  
  W <- mu*(1-mu)
  
  z <- eta + (y-mu)/W
  
  # B <- t(J_gamma) %*% J_gamma
  # print(B)
  
  JgT_W <- t(W * J_gamma)
  B <- JgT_W %*% J_gamma
  
  # tilda_J <- J_gamma * sqrt(omega)
  # B <- t(tilda_J) %*% tilda_J
  diag(B) <- diag(B) + c(inv_variances, rep(1/g, p_gam))
  
  inv_B <- solve(B)
  # L_invB <- chol(inv_B)
  
  JgT_W_z <- JgT_W %*% z
  # eta_hat <- inv_B %*% JgT_ty
  
  tilde_theta <- inv_B %*% JgT_W_z
  ALL$tilde_betahat <- tilde_theta
  # print(tilde_theta)
  
  eta <- as.vector(J_gamma %*% tilde_theta) 
  mu <- inv_logit(eta)            
  new_W <- mu*(1-mu)  
  
  JgT_W <- t(new_W * J_gamma)
  B <- JgT_W %*% J_gamma
  
  diag(B) <- diag(B) + c(inv_variances, rep(1/g, p_gam))
  inv_B <- solve(B)
  L_invB <- chol(inv_B)
  
  d <- p_z + p_gam
  
  if (d == 1) {
    log_det <- log(L_invB)
  }
  else {
    log_det <- sum(log(diag(L_invB)))
  }
  
  # print(sum((L_invB %*% t(J_gamma) %*% (y-mu))^2)/2)
  
  grad_theta <- t(J_gamma) %*% (y-mu) - c(inv_variances, rep(1/g, p_gam))*tilde_theta
  
  ALL$approx_llh <- log_det + d*log(2*pi)/2 - 
    logistic_log_post_beta(tilde_theta, y, J_gamma, g, inv_variances, 1, p_z, p_gam) + 
    # sum((L_invB %*% JgT_W_z)^2)/2
    sum((L_invB %*% grad_theta)^2)/2
  # + (p_gam/2)*log(4/(n_trials*g)) + n_trials / 8 *sum((L_invB %*% JgT_ty)^2)

  ALL$lmp <- hyper_par$log_m_prior(p_gam, hyper_par$h, hyper_par$p)
  
  return(ALL)
  
}


