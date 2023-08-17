# correlated pseudo-marginal for logistic regression
weibull_cpm <- function(ALL, hyper_par, ...) {
  
  curr <- ALL$curr
  p_gam <-ALL$p_gam
  v_indx <- ALL$v_indx
  v <- ALL$v
  u <- ALL$u
  k <- ALL$k
  g <- ALL$g
  
  p_z <- hyper_par$p_z
  t <- hyper_par$t
  d <- hyper_par$d
  inv_variances <- hyper_par$inv_variances
  var_intercept <- hyper_par$var_intercept
  var_fixed_covariate <- hyper_par$var_fixed_covariate
  
  n_particles <- hyper_par$n_particles
  rho <- hyper_par$rho
  
  new_v_indx <- which(curr == 1)
  
  initn <- is.null(u)
  
  if (initn) {
    if (p_gam == 0) {
      v <- NULL
    }
    else {
      v <- matrix(rnorm(p_gam*n_particles), ncol = n_particles)
    }
    u <- matrix(rnorm(p_z*n_particles), ncol = n_particles)
  }
  else {
    
    if (length(v) != 0) {
      if (p_gam > 0) {
        # both models aren't empty
        common_old <- v_indx %in% new_v_indx
        common_new <- new_v_indx %in% v_indx
        
        p_com <- sum(common_old)
        new_v <- matrix(NA, nrow = p_gam, ncol = n_particles)
        
        # print("g")
        # print(v_indx[common_old])
        # print(new_v_indx[common_new])
        
        if (p_com == 0) {
          new_v <- matrix(rnorm(p_gam*n_particles), ncol = n_particles)
        }
        else {
          new_v[common_new,] <- rho * v[common_old,] + 
            sqrt(1-rho^2) * matrix(rnorm(p_com*n_particles), ncol = n_particles)
          new_v[!common_new,] <- matrix(rnorm((p_gam - p_com)*n_particles), ncol = n_particles)
          
          # for (i in 1:p_com) {
          #   print(cor(v[which(common_old)[i],], new_v[which(common_new)[i],]))
          # }
          
          
        }
        
        v <- new_v
      }
      else {
        # new mode is empty
        v <- NULL
      }
    }
    else {
      if (p_gam > 0) {
        # old model is empty {
        v <- matrix(rnorm(p_gam*n_particles), ncol = n_particles)
      }
      else {
        # both models are empty
        v <- NULL
      }
    }
    u <- rho * u + sqrt(1-rho^2) * matrix(rnorm(p_z*n_particles), ncol = n_particles)
  }
  
  
  # update particles
  ALL$v <- v
  ALL$u <- u
  ALL$v_indx <- new_v_indx
  
  if (p_gam >= 1) {
    X_gamma <- as.matrix(hyper_par$X[, new_v_indx], ncol = p_gam)
    J_gamma <- cbind(hyper_par$Z, X_gamma)
  }
  else {
    J_gamma <- hyper_par$Z
  }
  
  optimals <- weibull_NR(t, J_gamma, d, k, g, inv_variances, p_z, p_gam)
  optimal_betahat <- optimals$betahat
  optimal_inv_hessian <- optimals$inv_hessian 
  
  working_vs <- mvnorm_trans(rbind(u,v), optimal_betahat, optimal_inv_hessian)
  
  working_v <- working_vs$v
  density_v <- working_vs$density
  # print(density_v)
  # numerators of improtance sampling
  # beta prior
  
  
  if (p_gam > 1) {
    if (n_particles == 1) {
      log_prior <- - sum(working_v[-(1:p_z),]^2)/(2*g) - p_gam*log(2*pi)/2 - p_gam*log(g)/2
    }
    else {
      log_prior <- - colSums(working_v[-(1:p_z),]^2)/(2*g) - p_gam*log(2*pi)/2 - p_gam*log(g)/2
    }
    
    # log_prior <- - colSums(working_v[-(1:p_z),]^2)/(2*g) - (p_gam+1)*log(2*pi)/2 - (p_gam+1)*log(g)/2
  }
  else if (p_gam == 0) {
    log_prior <- 0
  }
  else {
    log_prior <- - working_v[-(1:p_z),]^2/(2*g) - log(2*pi)/2 - log(g)/2
  }
  
  
  # if (inv_variances[1] != 0) {
  #   log_prior_fixed <- - working_v[1,]^2*inv_variances[1]/2 - log(2*pi)/2 + log(inv_variances[1])/2
  #   log_prior <- log_prior + log_prior_fixed
  # }
  
  if (!is.null(var_intercept)) {
    log_prior_fixed <- - working_v[1,]^2/(2*var_intercept) - log(2*pi)/2 - log(var_intercept)/2
    log_prior <- log_prior + log_prior_fixed
  }
  
  if (!is.null(var_fixed_covariate)) {
    if (n_particles == 1) {
      log_prior_fixed <- - sum(working_v[2:p_z,]^2)/(2*var_fixed_covariate) - (p_z-1)*log(2*pi)/2 - (p_z-1)*log(var_fixed_covariate)/2 
    }
    else if (p_z == 2){
      log_prior_fixed <- - working_v[2,]^2/(2*var_fixed_covariate) - log(2*pi)/2 - (p_z-1)*log(var_fixed_covariate)/2 
    }
    else {
      log_prior_fixed <- - colSums(working_v[2:p_z,]^2)/(2*var_fixed_covariate) - (p_z-1)*log(2*pi)/2 - (p_z-1)*log(var_fixed_covariate)/2 
    }
    log_prior <- log_prior + log_prior_fixed
  }
  
  # beta likelihood
  eta <- J_gamma %*% working_v
  
  llh_v <- colSums(d*(log(k) + k*eta + (k-1)*log(t)) - (t*exp(eta))^k)
  
  posterior_v <- llh_v + log_prior
  
  # don't need to devide by n_particles since this is a common factor
  # ALL$llh <- log(sum(exp(posterior_v - density_v))/n_particles)
  
  # aviod numerical outflow
  diff_v <- posterior_v - density_v
  min_value <- - max(diff_v)
  ALL$llh <- log(sum(exp(diff_v + min_value))) - min_value
  ALL$lmp <- hyper_par$log_m_prior(p_gam, hyper_par$h, hyper_par$p)
  
  ALL$eta <- optimals$eta
  
  
  return(ALL)
  
}
