# laplace approximation for logistic regression
cox_laplace <- function(ALL, hyper_par, ...) {
  
  
  curr <- ALL$curr
  p_gam <-ALL$p_gam
  g <- ALL$g
  
  
  t <- hyper_par$t
  d <- hyper_par$d
  
  Z <- hyper_par$Z
  p_z <- hyper_par$p_z
  inv_variances <- hyper_par$inv_variances
  
  # print(which(curr==1))
  # print(k)
  
  p_theta <- p_z + p_gam
  
  
  if (p_theta >= 1) {
    
    X_gamma <- as.matrix(hyper_par$X[, which(curr==1)], ncol = p_gam)
    
    optimals <- cox_NR_coxph(t, Z, X_gamma, d, p_gam, p_z, g, inv_variances)
    
    optimal_betahat <- optimals$betahat
    optimal_inv_hessian <- optimals$inv_hessian 
    
    if (p_theta == 1) {
      
      det_inv_hessian <- log(optimal_inv_hessian)
      
    }
    else {
      
      L_hessian <- chol(optimal_inv_hessian)
      det_inv_hessian <- sum(log(diag(L_hessian)))
      # cat(- det_inv_hessian, -log(det(optimal_inv_hessian))/2, "\n")
    }
    
    ALL$llh <- - cox_log_post_beta(optimal_betahat, t, optimals$J, d, 
                                   g, inv_variances, p_z, p_gam) + det_inv_hessian + p_theta*log(2*pi)/2
    
    ALL$eta <- optimals$eta
    
  }
  else {
    
    
    ALL$llh <- - cox_log_post_beta(NULL, t, NULL, d)
    ALL$eta <- rep(0, hyper_par$n)
    
  }
  
  ALL$approx_llh <- ALL$llh
  ALL$lmp <- hyper_par$log_m_prior(p_gam, hyper_par$h, hyper_par$p)
  
  
  return(ALL)
  
}


