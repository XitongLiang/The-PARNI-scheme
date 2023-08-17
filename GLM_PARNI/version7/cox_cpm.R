# correlated pseudo-marginal for logistic regression
cox_cpm <- function(ALL, hyper_par, ...) {
  
  curr <- ALL$curr
  p_gam <-ALL$p_gam
  v_indx <- ALL$v_indx
  v <- ALL$v
  u <- ALL$u
  g <- ALL$g
  
  
  t <- hyper_par$t
  d <- hyper_par$d
  
  Z <- hyper_par$Z
  p_z <- hyper_par$p_z
  inv_variances <- hyper_par$inv_variances
  # var_intercept <- hyper_par$var_intercept
  # var_fixed_covariate <- hyper_par$var_fixed_covariate
  
  
  n_particles <- hyper_par$n_particles
  rho <- hyper_par$rho
  
  new_v_indx <- which(curr == 1)
  
  p_theta <- p_z + p_gam
  
  if (p_theta > 0) {
    
    
    initn <- is.null(v) & is.null(u)
    
    if (initn) {
      if (p_gam == 0) {
        v <- NULL
      }
      else {
        v <- matrix(rnorm(p_gam*n_particles), ncol = n_particles)
      }
      
      if (p_z > 0) {
        u <- matrix(rnorm(p_z*n_particles), ncol = n_particles)
      }
      
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
      if (p_z > 0){
        u <- rho * u + sqrt(1-rho^2) * matrix(rnorm(p_z*n_particles), ncol = n_particles)
      }
      
    }
    
    
    # update particles
    ALL$v <- v
    ALL$u <- u
    ALL$v_indx <- new_v_indx
    
    X_gamma <- as.matrix(hyper_par$X[, new_v_indx], ncol = p_gam)
    
    
    optimals <- cox_NR_coxph(t, Z, X_gamma, d, p_gam, p_z, g, inv_variances)
    
    optimal_betahat <- optimals$betahat
    optimal_inv_hessian <- optimals$inv_hessian 
    
    working_vs <- mvnorm_trans(rbind(u,v), optimal_betahat, optimal_inv_hessian)
    
    working_v <- working_vs$v
    density_v <- working_vs$density
    # print(density_v)
    # numerators of improtance sampling
    # beta prior
    
    if (p_z > 0) {
      
      if (p_z > 1) {
        
        if (n_particles == 1) {
          log_prior_fixed <- - sum(working_v[(1:p_z),]^2*inv_variances)/2 - p_z*log(2*pi)/2 + sum(log(inv_variances))/2
        }
        else {
          log_prior_fixed <- - colSums(working_v[(1:p_z),]^2*inv_variances)/2 - p_z*log(2*pi)/2 + sum(log(inv_variances))/2
        }
        
      }
      else {
        
        log_prior_fixed <- - (working_v[1,]^2*inv_variances)/2 - log(2*pi)/2 + log(inv_variances)/2
        
      }
      
      
      
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
      
      log_prior <- log_prior + log_prior_fixed
      
    }
    else {
      if (p_gam > 1) {
        if (n_particles == 1) {
          log_prior <- - sum(working_v^2)/(2*g) - p_gam*log(2*pi)/2 - p_gam*log(g)/2
        }
        else {
          log_prior <- - colSums(working_v^2)/(2*g) - p_gam*log(2*pi)/2 - p_gam*log(g)/2
        }
        
        # log_prior <- - colSums(working_v[-(1:p_z),]^2)/(2*g) - (p_gam+1)*log(2*pi)/2 - (p_gam+1)*log(g)/2
      }
      else if (p_gam == 0) {
        log_prior <- 0
      }
      else {
        log_prior <- - working_v^2/(2*g) - log(2*pi)/2 - log(g)/2
      }
    }
    
    
    
    # beta likelihood
    eta <- optimals$J %*% working_v
    
    t2_v <- sapply(1:n_particles, 
                   function(i, exp_eta, t) {
                     return(cox_compute_t2(exp_eta[,i], t))
                   }, 
                   exp_eta = exp(eta), t = t)
    
    llh_v <- colSums(d*(eta - log(t2_v)))
    
    
    posterior_v <- llh_v + log_prior
    
    # don't need to devide by n_particles since this is a common factor
    # ALL$llh <- log(sum(exp(posterior_v - density_v))/n_particles)
    
    # aviod numerical outflow
    diff_v <- posterior_v - density_v
    min_value <- - max(diff_v)
    ALL$llh <- log(sum(exp(diff_v + min_value))) - min_value
    
    ALL$eta <- optimals$eta
    
    
    
  }
  else {
    
    ALL$llh <- - cox_log_post_beta(NULL, t, NULL, d)
    ALL$eta <- rep(0, hyper_par$n)
    
  }
  
  ALL$lmp <- hyper_par$log_m_prior(p_gam, hyper_par$h, hyper_par$p)
  
  
  return(ALL)
  
}

