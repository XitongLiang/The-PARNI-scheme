log_llh_g_L <- function(gamma, hyper_par){
  
  
  p_gam <- sum(gamma)
  
  n <- hyper_par$n
  p <- hyper_par$p
  h <- hyper_par$h
  g <- hyper_par$g
  yty <- hyper_par$yty
  
  includes <- which(gamma == 1)
  
  if(p_gam == 0){
    A <- yty
  }
  else{
    
    ytX <- hyper_par$ytX
    
    Xg <- hyper_par$X[, includes]
    L_Xg <- chol(t(Xg) %*% Xg)
    
    ytXgFXgty <- sum((forwardsolve(t(L_Xg), diag(p_gam)) %*% 
                        as.matrix(ytX[includes]))^2)
    
    A <-  yty - g/(1+g) * ytXgFXgty
    
  }
  
  log_llh <-  - p_gam/2 * log(1+g) - n * log(A)/2
  
  # log_m <- hyper_par$log_m_prior(p_gam, h, p)
  
  # log_post <- log_llh + log_m
  
  return(log_llh)
  
}



log_llh_ind_L <- function(gamma, hyper_par){
  
  p_gam <- sum(gamma)
  
  n <- hyper_par$n
  p <- hyper_par$p
  h <- hyper_par$h
  # h_alpha <- hyper_par$h_alpha
  # h_beta <- hyper_par$h_beta
  g <- hyper_par$g
  yty <- hyper_par$yty
  
  
  includes <- which(gamma == 1)
  
  if(p_gam == 0){
    A <- yty
    L_Vg <- NULL
    sqrt_det_Vg <- 0
  }
  else{
    
    X <- hyper_par$X
    # V <- hyper_par$V
    ytX <- hyper_par$ytX
    
    # Vg <- V[includes, includes]
    Xg <- X[, includes]
    Vg <- t(Xg) %*% Xg
    diag(Vg) <- diag(Vg) + 1/g
    L_Vg <- chol(Vg)
    
    ytXgFXgty <- sum((forwardsolve(t(L_Vg), diag(p_gam)) %*% 
                        as.matrix(ytX[includes]))^2)
    
    A <-  yty - ytXgFXgty
    
    if (p_gam == 1){
      sqrt_det_Vg <- log(L_Vg)
    }
    else {
      sqrt_det_Vg <- sum(log(diag(L_Vg)))
    }
    
  }
  
  log_llh <- - p_gam/2 * log(g) - sqrt_det_Vg - n * log(A)/2
  
  # fixed h
  # log_m <- p_gam * (log(h) - log(1-h))
  
  # binomial-beta prior
  # log_m <- sum(log(1:p_gam)) - sum(log(p - p_gam + h_beta - 0:p_gam))
  # log_m <- hyper_par$log_m_prior(p_gam, h, p) # lbeta(p_gam+h_alpha, p-p_gam+h_beta)
  
  # log_post <- log_llh + log_m
  
  return(log_llh)
}


