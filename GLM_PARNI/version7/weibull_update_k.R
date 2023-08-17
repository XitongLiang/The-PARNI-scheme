weibull_update_k <- function(ALL, hyper_par, k_rw_var, iter, n_chain) {
  
  k <- ALL$k
  epsilon <- rnorm(1, 0, k_rw_var)
  
  k_new <- k*exp(epsilon)
  
  
  # cat(k, k_new, epsilon, "\n")
  
  ALL_new <- ALL
  ALL_new$k <- k_new
  
  # if (k_new >=2) {
  #   ALL_new$llh <- -Inf
  # }
  # else {
  #   ALL_new <- log_post(ALL_new, hyper_par)
  # }
  
  ALL_new <- try(hyper_par$log_post(ALL_new, hyper_par), silent = TRUE)
  
  
  # print(is.list(ALL_new))
  # print(k_new)
  
  if (is.list(ALL_new)) {
    acc_rate <- ALL_new$llh + dnorm(log(k_new), 0, hyper_par$log_k_var, log = TRUE) -
      ALL$llh - dnorm(log(k), 0, hyper_par$log_k_var, log = TRUE)
  }
  else {
    acc_rate <- -Inf
  }
  
  # print(acc_rate)
  
  if (log(runif(1)) < acc_rate) {
    ALL <- ALL_new
    k_acc <- TRUE
  }
  else {
    k_acc <- FALSE
  }
  
  k_rw_var <- exp(log(k_rw_var) + iter^(-0.7)/n_chain*(min(1, exp(acc_rate)) - 0.25))
  
  return(list(ALL = ALL,
              k_acc = k_acc,
              k_rw_var = k_rw_var))
  
}
