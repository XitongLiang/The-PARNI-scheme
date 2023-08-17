hc_update_g <- function(ALL, hyper_par, g_rw_var, iter, n_chain) {
  
  g <- ALL$g
  epsilon <- rnorm(1, 0, g_rw_var)
  
  g_new <- g * exp(epsilon)
  
  ALL_new <- ALL
  ALL_new$g <- g_new
  
  
  ALL_new <- try(hyper_par$log_post(ALL_new, hyper_par), silent = TRUE)
  
  
  if (is.list(ALL_new)) {
      
    acc_rate <- ALL_new$llh + log_half_cauchy_wj(g_new) -
      ALL$llh - log_half_cauchy_wj(g)
    
  }
  else {
    acc_rate <- -Inf
  }
  
  # cat(g_new, ALL_new$llh, acc_rate, "\n")
  
  if (log(runif(1)) < acc_rate) {
    ALL <- ALL_new
    g_acc <- TRUE
  }
  else {
    g_acc <- FALSE
  }
  
  g_rw_var <- g_rw_var * exp(iter^(-0.7)/n_chain*(min(1, exp(acc_rate)) - 0.25))
  
  return(list(ALL = ALL,
              g_acc = g_acc,
              g_rw_var = g_rw_var))
  
}
