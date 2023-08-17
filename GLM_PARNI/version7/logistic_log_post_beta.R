logistic_log_post_beta <- function(beta, y, X, g, inv_variances, n_trials, p_z, p_gam) {
  
  if (p_gam != 0) {
    log_prior <- sum(beta[-(1:p_z)]^2)/(2*g) + p_gam*log(2*pi)/2 + p_gam*log(g)/2
  }
  else {
    log_prior <- 0
  }
  
  if (inv_variances[1] != 0) {
    log_prior <- log_prior + beta[1]^2*inv_variances[1]/2 + log(2*pi)/2 - log(inv_variances[1])/2
  }
  
  if (p_z > 1) {
    if (sum(inv_variances[-1]) != 0) {
      log_prior <- log_prior + sum(beta[2:p_z]^2*inv_variances[-1])/2 + (p_z-1)*log(2*pi)/2 - sum(log(inv_variances[-1]))/2
    }
  }
  
  eta <- X %*% beta
  
  # mu <- inv_logit(eta)
  # log_llh <- - sum(y * eta + n_trials * log(1-mu))
  
  log_llh <- - sum(y * eta - n_trials * log(1+exp(eta)))
  log_pi <- log_prior + log_llh
  # print(log_pi)
  return(log_pi)
  
}










