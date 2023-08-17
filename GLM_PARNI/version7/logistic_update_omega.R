logistic_update_omega_r <- function(ALL, hyper_par, initn = FALSE) {
  
  n <- hyper_par$n
  n_trials <- 1
  # y <- hyper_par$y
  omega <- rep(NA, n)
  
  if (initn) {
    
    for (i in 1:n) {
      # omega[i] <- rpolyagamma(n_trials[i], logit(y[i]/n_trials[i]))
      omega[i] <- rpolyagamma(n_trials, 0)
    }
  }
  else {
    
    J_gamma <- ALL$J_gamma
    
    A <- ALL$A
    L_invB <- ALL$L_invB
    
    # cat((inv_L_B %*% (t(inv_L_B) %*% A))[1], (t(L_invB) %*% (L_invB %*% A))[1], "\n")
    
    theta_gamma <- t(L_invB) %*% (L_invB %*% A) + t(L_invB) %*% rnorm(n = hyper_par$p_z + ALL$p_gam)
    eta <- J_gamma %*% theta_gamma
    for (i in 1:n) {
      omega[i] <- rpolyagamma(n_trials, eta[i])

    }
  }
  
  ALL$omega <- omega
  
  return(ALL)
  
}




library(pgdraw)
logistic_update_omega_cpp <- function(ALL, hyper_par, initn = FALSE) {

  n <- hyper_par$n
  n_trials <- 1
  # y <- hyper_par$y

  if (initn) {
    # omega <- pgdraw(n_trials, rep(0,n))
    omega <- pgdraw(n_trials, rep(0,n))
  }
  else {

    J_gamma <- ALL$J_gamma

    A <- ALL$A
    L_invB <- ALL$L_invB

    theta_gamma <- t(L_invB) %*% (L_invB %*% A) + t(L_invB) %*% rnorm(n = hyper_par$p_z + ALL$p_gam)
    
    omega <- pgdraw(n_trials, J_gamma %*% theta_gamma)
  }

  ALL$omega <- omega
  
  return(ALL)

}

