logistic_ALA <- function(ALL, hyper_par, ...) {
  
  curr <- ALL$curr
  p_gam <-ALL$p_gam
  g <- ALL$g
  
  p_z <- hyper_par$p_z
  # n_trials <- hyper_par$n_trials
  tilde_y <- 4*(hyper_par$y - 1/2)
  inv_variances <- hyper_par$inv_variances
  
  if (p_gam >= 1) {
    X_gamma <- as.matrix(hyper_par$X[, which(curr==1)], ncol = p_gam)
    J_gamma <- cbind(hyper_par$Z, X_gamma)
  }
  else {
    J_gamma <- hyper_par$Z
  }
  
  B <- t(J_gamma) %*% J_gamma
  # print(B)
  
  diag(B) <- diag(B) + 4 * c(inv_variances, rep(1/g, p_gam))
  
  inv_B <- solve(B)
  L_invB <- chol(inv_B)
  
  JgT_ty <- (t(J_gamma) %*% tilde_y)
  # eta_hat <- inv_B %*% JgT_ty
  
  d <- p_z + p_gam
  
  if (d == 1) {
    log_det <- log(L_invB)
  }
  else {
    log_det <- sum(log(diag(L_invB)))
  }
  
  
  ALL$approx_llh <- log_det + (p_gam/2)*log(4/g) + 1 / 8 *sum((L_invB %*% JgT_ty)^2)
  
  ALL$lmp <- hyper_par$log_m_prior(p_gam, hyper_par$h, hyper_par$p)
  
  # print((n_trials/8) * sum((L_invJTJ %*% JgT_ty)^2))
  # print((n_trials/8) * t(eta_hat) %*% JgTJg %*% eta_hat)
  
  return(ALL)
  
  
}





