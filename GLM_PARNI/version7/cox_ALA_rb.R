cox_ALA_rb <- function(ALL, hyper_par) {
  
  p <- hyper_par$p
  X <- hyper_par$X
  d <- hyper_par$d
  t <- hyper_par$t
  h <- hyper_par$h_exp
  
  Z <- hyper_par$Z
  p_z <- hyper_par$p_z
  inv_variances <- hyper_par$inv_variances
  
  curr <- ALL$curr
  p_gam <- ALL$p_gam
  g <- ALL$g
  
  
  W <- hyper_par$W
  tilde_y <- hyper_par$tilde_y
  
  p_theta <- p_z + p_gam
  
  if (p_theta >= 1) {
    
    includes <- which(curr == 1)
    
    if (p_z > 0) {
      if (p_gam > 0) {
        J_gamma <- cbind(Z, as.matrix(hyper_par$X[, includes], ncol = p_gam))
      }
      else {
        J_gamma <- Z
      }
    }
    else {
      J_gamma <- as.matrix(hyper_par$X[, includes], ncol = p_gam)
    }
    
    
    W_J_gamma <- W %*% J_gamma
    WX <- W %*% X
    
    # B <- t(W_J_gamma) %*% W_J_gamma
    B <- t(W_J_gamma) %*% J_gamma
    # print(inv_variances)
    diag(B) <- diag(B) + c(rep(inv_variances, p_z), rep(1/g, p_gam))
    inv_B <- solve(B)
    L_invB <- chol(inv_B)
    
    
    tyT_Jg_invB <- t(tilde_y) %*% J_gamma %*% inv_B
    tXT_tJg_LinvB <- t(X) %*% W_J_gamma %*% t(L_invB)
    tyT_Jg_LinvB <- as.numeric(t(tilde_y) %*% J_gamma %*% t(L_invB))
    
    tyT_X <- as.numeric(t(tilde_y) %*% X)
    
    tilde_diag_XtX <- colSums( WX * X )
    
    # Bayes_factor <- rep(NA, p)
    
    d_vec <- tilde_diag_XtX + 1/(g) - rowSums(tXT_tJg_LinvB^2)
    
    d_vec[includes] <- 1
    Bayes_factor <- d_vec^{-1/2}*g^{-1/2}*exp((tyT_X - colSums(tyT_Jg_LinvB * t(tXT_tJg_LinvB)))^2/(2*d_vec))
    
    # print(tilde_diag_XtX[1573])
    # print(rowSums(tXT_tJg_LinvB^2)[1573])
    
    q_j <- 0
    
    for (j in includes) {
      q_j <- q_j + 1
      d_j <- 1/inv_B[(p_z+q_j),(p_z+q_j)]
      Bayes_factor[j] <- d_j^{-1/2}*g^{-1/2}*exp(d_j*tyT_Jg_invB[(p_z+q_j)]^2/2)
    }
    
    
  }
  else {
    
    WX <- W %*% X
    
    d_vec <- colSums(WX * X) + 1/g
    Bayes_factor <- d_vec^{-1/2}*g^{-1/2}*exp(as.vector(t(tilde_y) %*% X)^2/(2*d_vec))
    
  }
  
  
  
  PIPs <- h*Bayes_factor/(1-h+h*Bayes_factor)
  
  # print(Bayes_factor[1])
  
  return(PIPs)
  
}








