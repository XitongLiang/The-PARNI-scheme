# negative log-posterior distribution of beta
cox_log_post_beta <- function(beta, t, X, d, g = Inf, inv_variances = 0, p_z = 0, p_gam = 0, t2 = NULL) {
  
  if ((p_z + p_gam) > 0) {
    
    log_prior <- 0
    
    if (p_z > 0) {
      
      log_prior <- log_prior + sum(beta[1:p_z]^2*inv_variances)/2 + p_z*log(2*pi)/2 - sum(log(inv_variances))/2
      
      if (p_gam > 0){
        
        log_prior <- log_prior + sum(beta[-(1:p_z)]^2)/(2*g) + p_gam*log(2*pi)/2 + p_gam*log(g)/2
        
      }
      
    }
    else  {
      
      # if (p_gam > 0)
      log_prior <- log_prior + sum(beta^2)/(2*g) + p_gam*log(2*pi)/2 + p_gam*log(g)/2
      
    }
    
    eta <- as.vector(X %*% beta)
    
    if (is.null(t2)) {
      exp_eta <- exp(eta)
      t2 <- cox_compute_t2(exp_eta, t) # sapply(t, cox_sum_R_j, t, exp_eta)
    }
    
    log_llh <- - sum(d*(eta - log(t2)))
    
  }
  else {
    
    exp_eta <- rep(1, length(t)) # exp(eta)
    t2 <- cox_compute_t2(exp_eta, t)
    
    log_prior <- 0
    log_llh <- sum(d*log(t2))
  }
  
  
  log_pi <- log_prior + log_llh
  
  # print(log_pi)
  
  return(log_pi)
  
}



