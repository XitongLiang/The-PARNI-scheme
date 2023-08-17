# laplace approximation for logistic regression
weibull_laplace <- function(ALL, hyper_par, ...) {
  
  curr <- ALL$curr
  p_gam <-ALL$p_gam
  k <- ALL$k
  g <- ALL$g
  
  p_z <- hyper_par$p_z
  t <- hyper_par$t
  d <- hyper_par$d
  inv_variances <- hyper_par$inv_variances
  
  
  # print(which(curr==1))
  # print(k)
  
  if (p_gam >= 1) {
    
    X_gamma <- as.matrix(hyper_par$X[, which(curr==1)], ncol = p_gam)
    J_gamma <- cbind(hyper_par$Z, X_gamma)
    
  }
  else {
    J_gamma <- hyper_par$Z
  }
  
  
  optimals <- weibull_NR(t, J_gamma, d, k, g, inv_variances, p_z, p_gam)
  
  # print(optimals)
  # 
  # nlm_results <- nlm(log_post_beta_ind, rep(0, p_z+p_gam), hessian = TRUE,
  #                    y, J_gamma, g, n_trials, p_z, p_gam)
  # 
  # print(nlm_results)
  
  optimal_betahat <- optimals$betahat
  optimal_inv_hessian <- optimals$inv_hessian 
  
  d_beta <- length(optimal_betahat)
  
  if (d_beta == 1) {
    
    det_inv_hessian <- log(optimal_inv_hessian)
    
  }
  else {
    
    L_hessian <- chol(optimal_inv_hessian)
    det_inv_hessian <- sum(log(diag(L_hessian)))
    # cat(- det_inv_hessian, -log(det(optimal_inv_hessian))/2, "\n")
  }
  
  
  ALL$llh <- - weibull_log_post_beta(optimal_betahat, t, J_gamma, d, k, g, inv_variances, p_z, p_gam) + det_inv_hessian + d_beta*log(2*pi)/2
  ALL$approx_llh <- ALL$llh
  
  ALL$lmp <- hyper_par$log_m_prior(p_gam, hyper_par$h, hyper_par$p)
  
  ALL$eta <- optimals$eta
  
  return(ALL)
  
}


