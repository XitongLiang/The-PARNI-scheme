# bayes factor

# g-prior
bf_g_L <- function(gamma, hyper_par){
  
  includes <- which(gamma == 1)
  
  n <- hyper_par$n
  p <- hyper_par$p
  yty <- hyper_par$yty
  # XtX <- hyper_par$XtX
  ytX <- hyper_par$ytX
  g <- hyper_par$g
  X <- hyper_par$X
  diag_V <- hyper_par$diag_V
  # y <- hyper_par$y
  # h <- hyper_par$h
  
  g_ratio <- g/(g+1)
  inv_sqrtg1 <- 1/sqrt(g+1)
  
  n_power <- n/2
  
  BF <- rep(NA, p)
  p_gam <- sum(gamma)
  
  # empty model
  if(p_gam == 0){
    A <- yty
    for(j in 1:p){
      tilda_A <- yty - g_ratio * ytX[j]^2/diag_V[j]
      
      A_ratio <- A/tilda_A
      
      BF[j] <- (A_ratio)^n_power * inv_sqrtg1
    }
  }
  else{
    
    zj <- 0
    
    Xg <- X[, includes]
    XgtX <- t(Xg) %*% X
    XgtXg <- XgtX[, includes]  # XtX[includes, includes]
    L_Xg <- chol(XgtXg)
    
    ytXg <- matrix(ytX[includes],nrow = 1,ncol = p_gam)
    
    F <- solve(XgtXg)
    L_F <- chol(F)
    
    ytXgFXgty <- sum((solve(t( L_Xg )) %*% t(ytXg))^2)
    
    A <- yty - ytXgFXgty * g_ratio
    
    # XgtX <- matrix(XtX[includes,], nrow = p_gam, ncol = p)
    # XgtX <- XtX[includes,]
    
    d_vec <- 1 / (diag_V - rowSums((t(XgtX) %*% t(L_F))^2))
    
    ytXFXtxj_vec <- ytXg %*% F %*% XgtX
    tilda_A_vec <- A - d_vec * (ytXFXtxj_vec - ytX)^2 * g_ratio
    
    # print(( d_vec * (ytXFXtxj_vec - ytX)^2 * g_ratio)[1])
    
    BF <- (A / tilda_A_vec)^n_power * inv_sqrtg1
    
    for (j in includes) {
      
      if(p_gam == 1){
        
        tilde_A <- yty
        A_ratio <- tilde_A/A
        
      }
      else{
        
        zj <- zj+1
        
        # dj <- 1/F[zj,zj]
        # tilde_A <- A + (ytXg %*% F[,zj] )^2 * g_ratio * dj
        # A_ratio <- tilde_A/A
        
        A_ratio <- 1 + (ytXg %*% F[,zj] )^2 * g_ratio / (A * F[zj,zj])
        
      }
      
      
      BF[j] <- A_ratio^n_power * inv_sqrtg1
      
    }
  }
  
  return(BF)
}


# independent prior
bf_ind_L <- function(gamma, hyper_par, L_gamma = NULL){
  
  includes <- which(gamma == 1)
  
  n <- hyper_par$n
  p <- hyper_par$p
  yty <- hyper_par$yty
  XtX <- hyper_par$XtX
  ytX <- hyper_par$ytX
  g <- hyper_par$g
  # y <- hyper_par$y
  # h <- hyper_par$h
  X <- hyper_par$X
  # diag_XtX <- hyper_par$diag_XtX
  diag_V <- hyper_par$diag_V # + inv_g
  inv_g <- 1/g
  inv_sqrt_g <- sqrt(inv_g)
  n_power <- n/2
  
  BF <- rep(NA, p)
  p_gam <- sum(gamma)

  # empty model
  if(p_gam == 0){
    A <- yty
    for(j in 1:p){
      
      dj <- 1 / diag_V[j]
      
      tilda_A <- yty - ytX[j]^2 * dj
      
      A_ratio <- A/tilda_A
      
      BF[j] <- (A_ratio)^n_power * sqrt(dj) * inv_sqrt_g
    }
  }
  else{
    
    zj <- 0
    
    # Vg <- V[includes, includes]
    
    if (is.null(L_gamma)) {
      Xg <- X[, includes]
      XgtX <- t(Xg) %*% X
      Vg <- XgtX[, includes]
      
      if (p_gam == 1){
        Vg <- Vg + inv_g
      }
      else {
        diag(Vg) <- diag(Vg) + inv_g
      }
      L_Vg <- chol(Vg)
    }
    else {
      L_Vg <- L_gamma
    }
    
    ytXg <- matrix(ytX[includes],nrow = 1,ncol = p_gam)
    
    F <- solve(Vg)
    L_F <- chol(F)
    
    ytXgFXgty <- sum((solve(t( L_Vg )) %*% t(ytXg))^2)
    
    A <- yty - ytXgFXgty
    
    # XgtX <- matrix(XtX[includes,], nrow = p_gam, ncol = p)
    # XgtX <- XtX[includes,]
    
    
    inv_d_vec <- 1 / (diag_V - rowSums((t(XgtX) %*% t(L_F))^2))
    # inv_d_vec <- 1 / (diag(XtX) + inv_g - diag(t(XgtX) %*% F %*%  XgtX))
    inv_d_vec[includes] <- 0
    
    #if (any(inv_d_vec < 0)){
    #  print(inv_d_vec)
    #}
    
    ytXFXtxj_vec <- ytXg %*% F %*% XgtX
    tilda_A_vec <- A - inv_d_vec * (ytXFXtxj_vec - ytX)^2 
    
    
    BF <- sqrt(inv_d_vec) * inv_sqrt_g * (A / tilda_A_vec)^n_power 
    
    
    for (j in includes) {
      
      if(p_gam == 1){
        
        tilde_A <- yty
        A_ratio <- tilde_A/A
        
        BF[j] <-  sqrt(F) * inv_sqrt_g * A_ratio^n_power
        
      }
      else{
        
        zj <- zj+1
        
        inv_dj <- F[zj,zj]
        
        tilde_A <- A + (ytXg %*% F[,zj] )^2 / inv_dj
        A_ratio <- tilde_A/A
        
        # A_ratio <- 1 + (ytXg %*% F[,zj] )^2 / (A * F[zj,zj])
        
        BF[j] <- sqrt(inv_dj) * inv_sqrt_g * A_ratio^n_power
        
      }
      
      
    }
  }
  
  return(BF)
}

