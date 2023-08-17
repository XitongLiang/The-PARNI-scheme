# ALA about 0
cox_ALA <- function(ALL, hyper_par, ...) {
  
  curr <- ALL$curr
  p_gam <-ALL$p_gam
  g <- ALL$g
  
  
  t <- hyper_par$t
  d <- hyper_par$d
  
  Z <- hyper_par$Z
  p_z <- hyper_par$p_z
  inv_variances <- hyper_par$inv_variances
  
  W <- hyper_par$W
  
  p_theta <- p_gam + p_z
  
  if (p_theta >= 1) {
    
    if (p_z > 0) {
      
      if (p_gam > 0) {
        
        X <- as.matrix(hyper_par$X[, which(curr==1)], ncol = p_gam)
        J <- cbind(Z, X)
        full_inv_variances <- c(rep(inv_variances, p_z), rep(1/g, p_gam))
        cox_formula <- formula(Surv(t,d) ~ ridge(Z, theta = inv_variances, scale = FALSE) + 
                                 ridge(X, theta = 1/g, scale = FALSE))
        
      }
      else {
        J <- Z
        full_inv_variances <- inv_variances
        cox_formula <- formula(Surv(t,d) ~ ridge(Z, theta = inv_variances, scale = FALSE))
      }
    }
    else {
      X <- as.matrix(hyper_par$X[, which(curr==1)], ncol = p_gam)
      J <- X
      full_inv_variances <- 1/g
      cox_formula <- formula(Surv(t,d) ~ ridge(X, theta = 1/g, scale = FALSE))
    }
    
    
    
    # initial point
    tilde_beta <- rep(0, p_theta) 
    
    
    # new beta
    
    ALL$tilde_betahat <- tilde_beta
    
    eta <- as.vector(J %*% tilde_beta) 
    
    cox_fit <- coxph(cox_formula, control = coxph.control(iter.max = 0),
                     init = tilde_beta)
    
    inv_B <- cox_fit$var
    
    L_invB <- chol(inv_B)
    
    lambda <- exp(eta)            
    tilde_y <- - cox_pseudores(lambda, t, d)
    
    
    
    
    if (p_theta == 1) {
      log_det <- log(L_invB)
    }
    else {
      log_det <- sum(log(diag(L_invB)))
    }
    
    # print(log_det)
    
    # print(sum((L_invB %*% t(J) %*% (y-mu))^2)/2)
    
    grad_beta <- t(J) %*% tilde_y + tilde_beta * full_inv_variances
    
    
    ALL$approx_llh <- log_det + p_theta*log(2*pi)/2 - 
      cox_log_post_beta(tilde_beta, t, J, d, g, inv_variances, p_z, p_gam)  +
      sum((L_invB %*% grad_beta)^2)/2
    
    
  }
  else {
    
    ALL$approx_llh <- - cox_log_post_beta(NULL, t, NULL, d)
    
  }
  
 
  
  ALL$lmp <- hyper_par$log_m_prior(p_gam, hyper_par$h, hyper_par$p)
  
  return(ALL)
  
}


