# negative log-posterior distribution of beta
weibull_log_post_beta <- function(beta, t, X, d, k, g, inv_variances, p_z, p_gam) {
  
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
  
  log_llh <- - sum( d*(log(k) + k*eta + (k-1)*log(t)) - (t*exp(eta))^k )
  # plot(eta, t)
  # print(sum(d*k*eta))
  # print(sum(d))
  # print(sum(d*(k-1)*log(t)))
  # print(sum((t*exp(eta))^k))
  
  log_pi <- log_prior + log_llh
  # print(log_pi)
  
  return(log_pi)
  
}







