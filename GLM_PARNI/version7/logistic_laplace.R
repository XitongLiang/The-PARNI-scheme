# laplace approximation for logistic regression
logistic_laplace <- function(ALL, hyper_par, ...) {
  
  curr <- ALL$curr
  p_gam <-ALL$p_gam
  g <- ALL$g
  
  p_z <- hyper_par$p_z
  n_trials <- 1 # hyper_par$n_trials
  y <- hyper_par$y
  inv_variances <- hyper_par$inv_variances
  
  if (p_gam >= 1) {
    X_gamma <- as.matrix(hyper_par$X[, which(curr==1)], ncol = p_gam)
    J_gamma <- cbind(hyper_par$Z, X_gamma)
  }
  else {
    J_gamma <- hyper_par$Z
  }
  
  # print(which(curr==1))
  
  optimals <- logistic_IRLS(y, J_gamma, g, inv_variances, n_trials, p_z, p_gam)
  
  # print(optimals)
  # 
  # nlm_results <- nlm(log_post_beta_ind, rep(0, p_z+p_gam), hessian = TRUE,
  #                    y, J_gamma, g, n_trials, p_z, p_gam)
  # 
  # print(nlm_results)
  
  optimal_betahat <- optimals$betahat
  optimal_inv_hessian <- optimals$inv_hessian 
  
  d <- length(optimal_betahat)
  
  if (d == 1) {
    
    det_inv_hessian <- log(optimal_inv_hessian)
    
  }
  else {
    
    L_hessian <- chol(optimal_inv_hessian)
    det_inv_hessian <- sum(log(diag(L_hessian)))
    # cat(- det_inv_hessian, -log(det(optimal_inv_hessian))/2, "\n")
  }
  
  ALL$llh <- - logistic_log_post_beta(optimal_betahat, y, J_gamma, g, inv_variances, n_trials, p_z, p_gam) + det_inv_hessian + d*log(2*pi)/2
  ALL$approx_llh <- ALL$llh
  
  ALL$lmp <- hyper_par$log_m_prior(p_gam, hyper_par$h, hyper_par$p)
  
  ALL$eta <- optimals$eta
  
  return(ALL)
  
}
