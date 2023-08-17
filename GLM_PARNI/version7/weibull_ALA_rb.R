weibull_ALA_rb <- function(ALL, hyper_par) {
  
  p <- hyper_par$p
  p_z <- hyper_par$p_z
  
  X <- hyper_par$X
  inv_variances <- hyper_par$inv_variances
  d <- hyper_par$d
  t <- hyper_par$t
  h <- hyper_par$h_exp
  
  curr <- ALL$curr
  p_gam <- ALL$p_gam
  k <- ALL$k
  g <- ALL$g
  
  tilda_y <- k*(t^k - d)
  W <- k^2*t^k
  
  includes <- which(curr == 1)
  J_gamma <- cbind(hyper_par$Z, X[, includes])
  
  tilda_J_gamma <- J_gamma * sqrt(W)
  tilda_X <- X * sqrt(W)
  
  B <- t(tilda_J_gamma) %*% tilda_J_gamma
  diag(B) <- diag(B) + c(inv_variances, rep(1/g, p_gam))
  inv_B <- solve(B)
  L_invB <- chol(inv_B)
  
  
  tyT_Jg_invB <- t(tilda_y) %*% J_gamma %*% inv_B
  tXT_tJg_LinvB <- t(tilda_X) %*% tilda_J_gamma %*% t(L_invB)
  tyT_Jg_LinvB <- as.numeric(t(tilda_y) %*% J_gamma %*% t(L_invB))
  
  tyT_X <- as.numeric(t(tilda_y) %*% X)
  
  tilda_diag_XtX <- colSums(tilda_X^2)
  
  # Bayes_factor <- rep(NA, p)
  
  d_vec <- tilda_diag_XtX + 1/(g) - rowSums(tXT_tJg_LinvB^2)
  
  d_vec[includes] <- 1
  Bayes_factor <- d_vec^{-1/2}*g^{-1/2}*exp((tyT_X - colSums(tyT_Jg_LinvB * t(tXT_tJg_LinvB)))^2/(2*d_vec))
  
  q_j <- 0
  
  for (j in includes) {
    q_j <- q_j + 1
    d_j <- 1/inv_B[(p_z+q_j),(p_z+q_j)]
    Bayes_factor[j] <- d_j^{-1/2}*g^{-1/2}*exp(d_j*tyT_Jg_invB[(p_z+q_j)]^2/2)
  }
  
  PIPs <- h*Bayes_factor/(1-h+h*Bayes_factor)
  
  return(PIPs)
  
}








