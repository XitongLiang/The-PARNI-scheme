compute_LA <- function(gamma, hyper_par, t = 1){
  
  p_gam <- sum(gamma)
  
  log_llh <- hyper_par$log_llh(gamma, hyper_par)
  
  log_m_prior <- hyper_par$log_m_prior(p_gam, hyper_par$h, hyper_par$p)
  
  return(list(# curr = gamma,
              llh = log_llh,
              lmp = log_m_prior,
              log_post = t*log_llh + log_m_prior,
              p_gam = p_gam))
}