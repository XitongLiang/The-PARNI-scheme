compute_LA_DAG <- function(gamma, hyper_par){
  
  p_gam <- sum(gamma)
  
  LA <- list(curr = gamma,
             p_gam = p_gam)
  
  LA <- hyper_par$log_llh(LA, hyper_par)
  
  # cat(LA$llh, log_llh_DAG(LA, hyper_par)$llh, "\n")
  
  log_m_prior <- hyper_par$log_m_prior(p_gam, hyper_par$h, hyper_par$p)
  
  LA$lmp <- log_m_prior
  LA$log_post <- LA$llh + log_m_prior
  
  return(LA)
}